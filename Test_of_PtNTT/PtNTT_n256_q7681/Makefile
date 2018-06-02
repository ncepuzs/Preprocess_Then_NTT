CC 		    = /usr/bin/cc
CFLAGS    =  -O3 -fomit-frame-pointer -march=native
RM 		    = /bin/rm


Speed:  poly.c fips202.c reduce.c cbd.c precomp.c ntt.c pt_ntt.c randombytes.c randombytes.h pt_ntt.h cpucycles.h cpucycles.c speed.c 
	$(CC) $(CFLAGS) poly.c fips202.c reduce.c pt_ntt.c cbd.c randombytes.c precomp.c ntt.c cpucycles.c speed.c -w -o speed

Test: test_ntt.c poly.c fips202.c reduce.c cbd.c pt_ntt.c precomp.c ntt.c randombytes.c randombytes.h params.h poly.h pt_ntt.h reduce.h cbd.h ntt.h
	$(CC)  test_ntt.c poly.c fips202.c reduce.c cbd.c pt_ntt.c randombytes.c precomp.c ntt.c -w -o test -lcrypto


.PHONY: clean test

clean:
	-$(RM) *.o
	-$(RM) -r test
	-$(RM) -r speed

