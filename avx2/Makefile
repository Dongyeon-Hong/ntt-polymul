CC ?= /usr/bin/cc
CFLAGS += -Wall -Wextra -Wpedantic -Wmissing-prototypes -Wredundant-decls \
  -Wshadow -Wpointer-arith -mavx2 -mbmi2 -mpopcnt -maes \
  -march=native -mtune=native -O3 -fomit-frame-pointer -flto
SABERFLAGS = -march=native -mtune=native -O3 -fomit-frame-pointer
RM = /bin/rm

HEADERS = params.h polyvec.h poly.h cbd.h aes256ctr.h randombytes.h consts256.h consts512.h consts1024.h fq.inc shuffle.inc
SOURCES = poly.c cbd.c aes256ctr.c randombytes.c
SOURCES256N = $(SOURCES) ntt256n.S invntt256n.S basemul256x1.S
OBJECTS = aes256ctr.o
OBJECTS256NX2 = $(OBJECTS) cbd256.o ntt256n.o invntt256n.o basemul256x1x2.o
OBJECTS256NX3 = $(OBJECTS) cbd256.o ntt256n.o invntt256n.o basemul256x1x3.o
OBJECTS256NX4 = $(OBJECTS) cbd256.o ntt256n.o invntt256n.o basemul256x1x4.o

.PHONY: all clean

all: \
  test_ntt256n \
  test_nttrange256n \
  test_basemul256n \
  test_speed256nx2 \
  test_speed256nx3 \
  test_speed256nx4 \
  test_range256n \
  test_lightsabermul \
  test_sabermul \
  test_firesabermul \

lib: \
  liblightsabermul.a \
  libsabermul.a \
  libfiresabermul.a \

%.o: %.c $(HEADERS)
	$(CC) $(CFLAGS) -c -o $@ $<

%.o: %.S $(HEADERS)
	$(CC) $(CFLAGS) -c -o $@ $<

cbd256.o: cbd.c $(HEADERS)
	$(CC) $(CFLAGS) -DKEM_N=256 -c -o $@ $<

poly_lightsaber.o: poly.c $(HEADERS)
	$(CC) $(CFLAGS) -DLIGHTSABER -c -o $@ $<

poly_saber.o: poly.c $(HEADERS)
	$(CC) $(CFLAGS) -DSABER -c -o $@ $<

poly_firesaber.o: poly.c $(HEADERS)
	$(CC) $(CFLAGS) -DFIRESABER -c -o $@ $<

polyvec_lightsaber.o: polyvec.c $(HEADERS)
	$(CC) $(CFLAGS) -DLIGHTSABER -c -o $@ $<

polyvec_saber.o: polyvec.c $(HEADERS)
	$(CC) $(CFLAGS) -DSABER -c -o $@ $<

polyvec_firesaber.o: polyvec.c $(HEADERS)
	$(CC) $(CFLAGS) -DFIRESABER -c -o $@ $<

sabermul/lightsabermul.o: sabermul/consts.h sabermul/matrix.c sabermul/scm_avx.c sabermul/sabermul.c sabermul/toom-cook_4way.c $(HEADERS)
	$(CC) $(SABERFLAGS) -DLIGHTSABER sabermul/sabermul.c -c -o $@

sabermul/sabermul.o: sabermul/consts.h sabermul/matrix.c sabermul/scm_avx.c sabermul/sabermul.c sabermul/toom-cook_4way.c $(HEADERS)
	$(CC) $(SABERFLAGS) -DSABER sabermul/sabermul.c -c -o $@

sabermul/firesabermul.o: sabermul/consts.h sabermul/matrix.c sabermul/scm_avx.c sabermul/sabermul.c sabermul/toom-cook_4way.c $(HEADERS)
	$(CC) $(SABERFLAGS) -DFIRESABER sabermul/sabermul.c -c -o $@

basemul256x1x2.o: basemul256x1.S $(HEADERS)
	$(CC) -c $(CFLAGS) -DKEM_K=2 -c -o $@ $<

basemul256x1x3.o: basemul256x1.S $(HEADERS)
	$(CC) -c $(CFLAGS) -DKEM_K=3 -c -o $@ $<

basemul256x1x4.o: basemul256x1.S $(HEADERS)
	$(CC) -c $(CFLAGS) -DKEM_K=4 -c -o $@ $<

liblightsabermul.a: poly_lightsaber.o polyvec_lightsaber.o $(OBJECTS256NX2)
	$(AR) -r $@ $^

libsabermul.a: poly_saber.o polyvec_saber.o $(OBJECTS256NX3)
	$(AR) -r $@ $^

libfiresabermul.a: poly_firesaber.o polyvec_firesaber.o $(OBJECTS256NX4)
	$(AR) -r $@ $^

test_ntt256n: $(SOURCES256N) $(HEADERS) test_nttn.c cpucycles.c cpucycles.h
	$(CC) $(CFLAGS) $(SOURCES256N) -DLIGHTSABER test_nttn.c cpucycles.c -o $@

test_nttrange256n: $(SOURCES) $(HEADERS) test_nttn.c nttrange.c nttrange256n.c consts256n7681.c consts256n10753.c nttrange.h
	$(CC) $(CFLAGS) $(SOURCES) -DFIRESABER -DSCALAR_NTT test_nttn.c nttrange.c nttrange256n.c consts256n7681.c consts256n10753.c -o $@ -ldb

test_basemul256n: $(SOURCES256N) $(HEADERS) test_basemuln.c
	$(CC) $(CFLAGS) $(SOURCES256N) -DLIGHTSABER test_basemuln.c -o $@

test_lightsabermul: $(SOURCES256N) $(HEADERS) test_polyvecmul.c polyvec.c sabermul/lightsabermul.o
	$(CC) $(CFLAGS) $(SOURCES256N) -DLIGHTSABER test_polyvecmul.c polyvec.c sabermul/lightsabermul.o -o $@

test_sabermul: $(SOURCES256N) $(HEADERS) test_polyvecmul.c polyvec.c sabermul/sabermul.o
	$(CC) $(CFLAGS) $(SOURCES256N) -DSABER test_polyvecmul.c polyvec.c sabermul/sabermul.o -o $@

test_firesabermul: $(SOURCES256N) $(HEADERS) test_polyvecmul.c polyvec.c sabermul/firesabermul.o
	$(CC) $(CFLAGS) $(SOURCES256N) -DFIRESABER test_polyvecmul.c polyvec.c sabermul/firesabermul.o -o $@

test_speed256nx2: $(SOURCES256N) $(HEADERS) test_speed.c polyvec.c sabermul/lightsabermul.o cpucycles.c cpucycles.h speed_print.c speed_print.h
	$(CC) $(CFLAGS) $(SOURCES256N) -DLIGHTSABER test_speed.c polyvec.c sabermul/lightsabermul.o cpucycles.c speed_print.c -o $@

test_speed256nx3: $(SOURCES256N) $(HEADERS) test_speed.c polyvec.c sabermul/sabermul.o cpucycles.c cpucycles.h speed_print.c speed_print.h
	$(CC) $(CFLAGS) $(SOURCES256N) -DSABER test_speed.c polyvec.c sabermul/sabermul.o cpucycles.c speed_print.c -o $@

test_speed256nx4: $(SOURCES256N) $(HEADERS) test_speed.c polyvec.c sabermul/firesabermul.o cpucycles.c cpucycles.h speed_print.c speed_print.h
	$(CC) $(CFLAGS) $(SOURCES256N) -DFIRESABER test_speed.c polyvec.c sabermul/firesabermul.o cpucycles.c speed_print.c -o $@

test_range256n: test_range256n.c nttrange.c nttrange256n.c consts256n7681.c consts256n10753.c nttrange.h consts256.h
	$(CC) $(CFLAGS) -DFIRESABER test_range256n.c nttrange.c nttrange256n.c consts256n7681.c consts256n10753.c -o $@ -ldb

clean:
	-$(RM) -rf *.o sabermul/*.o *.a *.so
	-$(RM) -f access.db
	-$(RM) -f test_ntt256n
	-$(RM) -f test_nttrange256n
	-$(RM) -f test_basemul256n
	-$(RM) -f test_speed256nx2
	-$(RM) -f test_speed256nx3
	-$(RM) -f test_speed256nx4
	-$(RM) -f test_range256n
	-$(RM) -f test_lightsabermul
	-$(RM) -f test_sabermul
	-$(RM) -f test_firesabermul
