#include <stdio.h>

#include "consts.h"
#include "params.h"
#include "poly.h"
#include "polyvec.h"
#include "randombytes.h"

#include "aes256ctr.h"

typedef struct {
    uint8_t *sx;
    uint8_t neg_start;
    uint8_t cnt;
} sppoly; // sparse poly

void samplevec_hwt(int8_t *s, uint8_t *seed, uint16_t nonce, uint16_t hmwt);
void smaug_matrix_vector_mul(polyvec *r, const polyvec a[KEM_K],
                             const int8_t *svec);
void smg_poly_add(uint16_t res[2 * KEM_N], const poly *op1, const uint8_t deg);
void smg_poly_sub(uint16_t res[2 * KEM_N], const poly *op1, const uint8_t deg);
void smg_poly_reduce(poly *res, const uint16_t temp[2 * KEM_N]);
void smg_poly_mult_add(poly *res, const poly *op1, const sppoly *op2);
uint8_t convToIdx(uint8_t *res, const uint8_t res_length, const uint8_t *op,
                  const size_t op_length);

static void poly_naivemul(poly *c, const poly *a, const poly *b) {
    unsigned int i, j;
    int16_t t[2 * KEM_N] = {0};

    for (i = 0; i < KEM_N; ++i)
        for (j = 0; j < KEM_N; ++j)
            t[i + j] =
                (t[i + j] + (int32_t)a->coeffs[i] * b->coeffs[j]) % KEM_Q;

    for (i = KEM_N; i < 2 * KEM_N; i++)
        c->coeffs[i - KEM_N] = (t[i - KEM_N] - t[i]) % KEM_Q;
}

void samplevec_hwt(int8_t *s, uint8_t *seed, uint16_t nonce, uint16_t hmwt) {
    uint32_t i = 0, pos = 0;
    __attribute__((aligned(32))) uint32_t buf[1088 / 4] = {0};

    aes256ctr_prf((uint8_t *)buf, sizeof(buf), seed, nonce);

    for (i = KEM_N * KEM_K - hmwt; i < KEM_N * KEM_K; ++i) {
        uint64_t deg_tmp = 0;
        deg_tmp = (uint64_t)buf[pos] * (i + 1);
        uint32_t deg = (uint32_t)(deg_tmp >> 32);

        s[i] = s[deg];
        s[deg] = ((buf[(hmwt + (pos >> 4))] >> (pos & 0x0f)) & 0x02) - 1;
        pos++;
    }
}

void smg_poly_add(uint16_t res[2 * KEM_N], const poly *op1, const uint8_t deg) {
    for (size_t i = 0; i < KEM_N; ++i) {
        res[deg + i] += op1->coeffs[i];
    }
}

void smg_poly_sub(uint16_t res[2 * KEM_N], const poly *op1, const uint8_t deg) {
    for (size_t i = 0; i < KEM_N; ++i) {
        res[deg + i] -= op1->coeffs[i];
    }
}

void smg_poly_reduce(poly *res, const uint16_t temp[2 * KEM_N]) {
    for (size_t j = 0; j < KEM_N; ++j) {
        res->coeffs[j] += temp[j] - temp[j + KEM_N];
        res->coeffs[j] %= KEM_Q;
    }
}

void smg_poly_mult_add(poly *res, const poly *op1, const sppoly *op2) {
    uint16_t temp[KEM_N * 2] = {0};

    for (size_t j = 0; j < op2->neg_start; ++j) {
        smg_poly_add(temp, op1, (op2->sx)[j]);
    }

    for (size_t j = op2->neg_start; j < op2->cnt; ++j) {
        smg_poly_sub(temp, op1, (op2->sx)[j]);
    }
    smg_poly_reduce(res, temp);
}

uint8_t convToIdx(uint8_t *res, const uint8_t res_length, const uint8_t *op,
                  const size_t op_length) {
    uint8_t index = 0, b = 0;
    uint8_t index_arr[2] = {0, res_length - 1}; // 0 for positive, 1 for
                                                // negative

    for (size_t i = 0; i < op_length; ++i) {
        index = ((op[i] & 0x80) >> 7) & 0x01;
        b = (-(uint64_t)op[i]) >> 63;
        res[index_arr[index]] ^= (-b) & (res[index_arr[index]] ^ i);
        index_arr[index] += op[i];
    }

    return index_arr[0];
}

void smaug_matrix_vector_mul(polyvec *r, const polyvec a[KEM_K],
                             const int8_t *svec) {
    //convToIdx
    size_t cnt_arr_idx = 0;
    uint8_t cnt_arr[KEM_K] = {0};

    for (size_t i = 0; i < KEM_N * KEM_K; ++i) {
        cnt_arr_idx = ((i & 0x700) >> 8) & (-(svec[i] & 0x01));
        cnt_arr[cnt_arr_idx] += (svec[i] & 0x01);
    }

    sppoly sp_vec[KEM_K];
    for (size_t i = 0; i < KEM_K; ++i) {
        sp_vec[i].cnt = cnt_arr[i];
        sp_vec[i].sx = (uint8_t *)calloc(cnt_arr[i], sizeof(uint8_t));
        sp_vec[i].neg_start =
            convToIdx(sp_vec[i].sx, sp_vec[i].cnt,
                      (const uint8_t *)svec + (i * KEM_N), KEM_N);
    }

    for (size_t i = 0; i < KEM_K; ++i)
        printf("%hu ", sp_vec[i].neg_start);
    printf("\n");

    // smaug mult
    for (int i = 0; i < KEM_K; i++) {
        for (int j = 0; j < KEM_K; j++) {
            smg_poly_mult_add(&(r->vec[i]), &(a[i].vec[j]), &(sp_vec[j]));
        }
    }
}

int main(void) {
    int i, j, err, hmwt = 0;
    uint16_t nonce = 0;
    uint8_t seed[POLYMUL_SYMBYTES];
    polyvec a[KEM_K];
    polyvec s, t, u;
    int8_t svec[KEM_N * KEM_K] = {0};
    poly tmp;
    err = 0;

#ifdef LIGHTSABER
    hmwt = 140;
#elif defined(SABER)
    hmwt = 198;
#elif defined(FIRESABER)
    hmwt = 176;
#endif

    randombytes(seed, POLYMUL_SYMBYTES);
    for (i = 0; i < POLYMUL_SYMBYTES; ++i) {
        seed[i] = i;
    }
    for (i = 0; i < KEM_K; i++)
        polyvec_uniform(&a[i], seed, nonce++);
    // polyvec_noise(&s, seed, nonce++); // replace with hwt sampling
    samplevec_hwt(svec, seed, nonce++, hmwt);
    for (i = 0; i < KEM_K; ++i) {
        for (j = 0; j < KEM_N; ++j) {
            s.vec[i].coeffs[j] = svec[KEM_N * i + j];
        }
    }

    for (i = 0; i < KEM_K; i++) {
        poly_naivemul(&t.vec[i], &a[i].vec[0], &s.vec[0]);
        for (j = 1; j < KEM_K; j++) {
            poly_naivemul(&tmp, &a[i].vec[j], &s.vec[j]);
            poly_add(&t.vec[i], &t.vec[i], &tmp);
        }
    }

    // replace with smaug mult including convToIdx
    // saber_matrix_vector_mul(&u, a, &s);
    smaug_matrix_vector_mul(&u, a, svec);
    for (i = 0; i < KEM_K; i++)
        for (j = 0; j < KEM_N; j++)
            if ((u.vec[i].coeffs[j] - t.vec[i].coeffs[j]) % KEM_Q) {
                printf("ERROR1: %d, %d, %d, %d\n", i, j, t.vec[i].coeffs[j],
                       u.vec[i].coeffs[j]);
                err = 1;
            }

    polyvec_matrix_vector_mul(&u, a, &s, 0); // not need to touch
    for (i = 0; i < KEM_K; i++)
        for (j = 0; j < KEM_N; j++)
            if ((u.vec[i].coeffs[j] - t.vec[i].coeffs[j]) % KEM_Q) {
                printf("ERROR2: %d, %d, %d, %d\n", i, j, t.vec[i].coeffs[j],
                       u.vec[i].coeffs[j]);
                err = 1;
            }
    if (!err)
        printf("ALL GOOD.\n");
    return err;
}
