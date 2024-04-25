CC=/usr/bin/g++
CFLAGS += -march=native -fomit-frame-pointer
#LDFLAGS=-lcrypto

SOURCES= Math.cpp readwrite.cpp main.cpp
HEADERS= Math.h params.h readwrite.h
main: $(HEADERS) $(SOURCES)
	$(CC) $(CFLAGS) -o $@ $(SOURCES) $(LDFLAGS)

.PHONY: clean

clean:
	-rm main

