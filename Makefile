CC     = gcc
CFLAGS = -g -D_GNU_SOURCE
LFLAGS =
OFILES = collect.o process.o util.o
HEADERS = clink.h

FILES = COPYRIGHT Makefile clink.doc clink.h collect.c process.c util.c

%.o: %.c $(HEADERS)
	$(CC) $(CFLAGS) -c $<

clink: $(OFILES)
	$(CC) $(LFLAGS) -o clink $(OFILES) -lm

tar:
	rm -f clink.1.0/*
	cp $(FILES) clink.1.0
	tar -cvf clink.1.0.tar clink.1.0
	tar -czf clink.1.0.tar.gz clink.1.0