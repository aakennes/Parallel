CC=mpic++
CFLAGS += -fomit-frame-pointer -std=c++17 -fopenmp -O2
#LDFLAGS=-lcrypto

SOURCES= Math.cpp readwrite.cpp ntt.cpp main.cpp 
HEADERS= Math.h ../params.h readwrite.h ntt.h
main: $(HEADERS) $(SOURCES)
	$(CC) $(CFLAGS) -o $@ $(SOURCES) $(LDFLAGS)

.PHONY: clean

clean:
	-rm main

