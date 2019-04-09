CC=gcc
SRC=src
BIN=bin
CFLAGS=-O0

default: all

all: seq pseudo_par parallel
	$(CC) $(CFLAGS) $(SRC)/main.c -o $(BIN)/main

seq:
	@echo "Not implemented!"

pseudo_par:
	@echo "Not implemented!"
	
parallel:
	@echo "Not implemented!"

clean:
	rm -rf bin/
