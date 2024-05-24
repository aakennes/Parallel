CC=/usr/bin/g++
CFLAGS += -march=core-avx2 -fomit-frame-pointer -std=c++17 -lpthread 
#LDFLAGS=-lcrypto

SOURCES= SIMD.cpp 
HEADERS= SIMD.h
main: $(HEADERS) $(SOURCES)
	$(CC) $(CFLAGS) -o $@ $(SOURCES) $(LDFLAGS)

.PHONY: clean

clean:
	-rm main

