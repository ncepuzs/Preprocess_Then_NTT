CC=/usr/bin/gcc
LDFLAGS=-lcrypto
CFLAGS    =  -O3 -fomit-frame-pointer -march=native 

Test: test_NTT.c fips202.c pt_ntt.c ntt.c poly.c precomp.c reduce.c rng.c  pt_ntt.h fips202.h ntt.h poly.h reduce.h rng.h
	$(CC) test_NTT.c fips202.c pt_ntt.c ntt.c poly.c  reduce.c rng.c -w -o test $(LDFLAGS)

Speed: speed.c cpucycles.c pt_ntt.c fips202.c ntt.c poly.c precomp.c reduce.c rng.c cpucycles.h pt_ntt.h fips202.h ntt.h poly.h reduce.h rng.h
	$(CC) $(CFLAGS) speed.c cpucycles.c pt_ntt.c fips202.c ntt.c precomp.c poly.c  reduce.c rng.c -w -o speed $(LDFLAGS)
.PHONY: clean

clean:
	-rm PQCgenKAT_kem
