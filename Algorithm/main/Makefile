CC=/usr/bin/g++
CFLAGS += -fomit-frame-pointer -std=c++17
#LDFLAGS=-lcrypto

SOURCES= Math.cpp readwrite.cpp ntt.cpp  main.cpp 
HEADERS= ../params.h Math.h readwrite.h ntt.h
main: $(HEADERS) $(SOURCES)
	$(CC) $(CFLAGS) -o $@ $(SOURCES) $(LDFLAGS)

.PHONY: clean

clean:
	-rm main

