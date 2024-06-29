CC=/usr/bin/g++
CFLAGS += -march=core-avx2 -fomit-frame-pointer -std=c++17 -O2
#LDFLAGS=-lcrypto

SOURCES= Math.cpp SIMD.cpp readwrite.cpp ntt.cpp  main.cpp 
HEADERS= ../params.h Math.h SIMD.h readwrite.h ntt.h
main: $(HEADERS) $(SOURCES)
	$(CC) $(CFLAGS) -o $@ $(SOURCES) $(LDFLAGS)

.PHONY: clean

clean:
	-rm main

