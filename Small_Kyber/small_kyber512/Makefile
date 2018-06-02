CC=/usr/bin/gcc
LDFLAGS=-lcrypto
CFLAGS    = -O3 -fomit-frame-pointer -march=native

SOURCES= cbd.c fips202.c indcpa.c kem.c ntt.c poly.c polyvec.c precomp.c pt_ntt.c reduce.c rng.c verify.c
HEADERS= api.h cbd.h fips202.h indcpa.h ntt.h params.h poly.h polyvec.h pt_ntt.h reduce.h rng.h verify.h

PQCgenKAT_kem: $(HEADERS) $(SOURCES) PQCgenKAT_kem.c
	$(CC) -o $@ $(SOURCES) PQCgenKAT_kem.c $(LDFLAGS)
Speed512: $(SOURCES) $(HEADERS) cpucycles.h cpucycles.c speed.c
	$(CC) $(CFLAGS) -DKYBER_K=2 $(SOURCES)  cpucycles.c speed.c -o  speed $(LDFLAGS)
	
.PHONY: clean

clean:
	-rm PQCgenKAT_kem

