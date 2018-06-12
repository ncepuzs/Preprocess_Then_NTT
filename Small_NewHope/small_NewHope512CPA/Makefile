CC=/usr/bin/gcc
LDFLAGS=-lcrypto
CFLAGS    = -O3 -fomit-frame-pointer -march=native


SOURCES= cpapke.c fips202.c kem.c ntt.c poly.c precomp.c reduce.c rng.c verify.c
HEADERS= api.h cpapke.h fips202.h ntt.h params.h poly.h reduce.h rng.h verify.h 

PQCgenKAT_kem: PQCgenKAT_kem.c $(HEADERS) $(SOURCES) pt_ntt.c pt_ntt.h
	$(CC) -w -o $@ PQCgenKAT_kem.c $(SOURCES) pt_ntt.c $(LDFLAGS)

Speed: speed.c pt_ntt.c cpucycles.c pt_ntt.h $(HEADERS) $(SOURCES)
	$(CC) speed.c pt_ntt.c cpucycles.c $(SOURCES) -w -o speed $(LDFLAGS) $(CFLAGS)
.PHONY: clean

clean:
	-rm PQCgenKAT_kem test speed
