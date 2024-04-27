CC=/usr/bin/g++
CFLAGS += -march=core-avx2 -fomit-frame-pointer
#LDFLAGS=-lcrypto

SOURCES= Math.cpp readwrite.cpp ntt.cpp SIMD.cpp main.cpp 
HEADERS= Math.h params.h readwrite.h ntt.h SIMD.h
main: $(HEADERS) $(SOURCES)
	$(CC) $(CFLAGS) -o $@ $(SOURCES) $(LDFLAGS)

.PHONY: clean

clean:
	-rm main

