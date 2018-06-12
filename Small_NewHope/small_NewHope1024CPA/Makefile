CC=/usr/bin/gcc
LDFLAGS=-lcrypto
CFLAGS    = -O3 -fomit-frame-pointer -march=native

SOURCES= cpapke.c fips202.c kem.c ntt.c poly.c precomp.c pt_ntt.c reduce.c rng.c verify.c
HEADERS= api.h cpapke.h fips202.h ntt.h params.h poly.h pt_ntt.h reduce.h rng.h verify.h 

PQCgenKAT_kem: PQCgenKAT_kem.c $(HEADERS) $(SOURCES)
	$(CC) -w -o $@ PQCgenKAT_kem.c $(SOURCES) $(LDFLAGS)
Speed: speed.c cpucycles.c $(HEADERS) $(SOURCES)
	$(CC) speed.c cpucycles.c $(SOURCES) -w -o speed $(LDFLAGS) $(CFLAGS)
.PHONY: clean

clean:
	-rm PQCgenKAT_kem
