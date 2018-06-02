CC=/usr/bin/gcc
LDFLAGS=-lcrypto

Test: test_NTT.c fips202.c pre_ntt.c ntt.c poly.c precomp.c reduce.c rng.c  pre_ntt.h fips202.h ntt.h poly.h reduce.h rng.h
	$(CC) test_NTT.c fips202.c pre_ntt.c ntt.c poly.c  reduce.c rng.c -o test $(LDFLAGS)

speed: speed.c cpucycles.c pre_ntt.c fips202.c ntt.c poly.c precomp.c reduce.c rng.c cpucycles.h pre_ntt.h fips202.h ntt.h poly.h reduce.h rng.h
	$(CC) speed.c cpucycles.c pre_ntt.c fips202.c ntt.c precomp.c poly.c  reduce.c rng.c -o speed $(LDFLAGS)
.PHONY: clean

clean:
	-rm PQCgenKAT_kem
