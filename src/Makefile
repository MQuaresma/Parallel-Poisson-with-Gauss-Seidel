CC=gcc
MPICC=mpicc
SEQ=SEQ
RB-SEQ=RB-SEQ
RB-PAR=RB-PAR
UTILS=UTILS
BIN=bin
CFLAGS=-O0 -lm
IFLAGS=-I $(UTILS)

default: all

all: seq pseudo_par parallel
	$(CC) $(CFLAGS) main.c -o $(BIN)/main

seq: $(SEQ)/seq.c $(UTILS)/*.[hc]
	$(CC) $(CFLAGS) $(IFLAGS) $(SEQ)/seq.c $(UTILS)/utils.c -o $(BIN)/seq

pseudo_par: $(RB-SEQ)/rb_seq.c $(UTILS)/*.[hc]
	$(CC) $(CFLAGS) $(IFLAGS) $(RB-SEQ)/rb_seq.c $(UTILS)/utils.c -o $(BIN)/rb_seq
	
parallel: $(RB-PAR)/rb_par.c $(UTILS)/*.[hc]
	$(MPICC) $(CFLAGS) $(IFLAGS) $(RB-PAR)/rb_par.c $(UTILS)/utils.c -o $(BIN)/rb_par

clean:
	rm -rf bin/*
