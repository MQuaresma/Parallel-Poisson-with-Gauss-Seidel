CC=gcc
SRC=src
BIN=bin
CFLAGS=-O0

default: all

all:
	$(CC) $(CFLAGS) $(SRC)/main.c -o $(BIN)/main

clean:
	rm -rf bin/
