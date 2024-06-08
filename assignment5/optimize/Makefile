CC=mpic++
CFLAGS += -march=core-avx2 -fomit-frame-pointer -std=c++17 -lpthread -O2
#LDFLAGS=-lcrypto

SOURCES= SIMD.cpp 
HEADERS= SIMD.h ../params.h
main: $(HEADERS) $(SOURCES)
	$(CC) $(CFLAGS) -o $@ $(SOURCES) $(LDFLAGS)

.PHONY: clean

clean:
	-rm main

