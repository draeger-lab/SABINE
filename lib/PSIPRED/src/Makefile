######################################################################
#	Makefile for PSIpred2                                        #
######################################################################

#### CC		Command to execute C compiler
#### CFLAGS	Flags to pass to C compiler.

INCDIR		= .
CC		= gcc

LIBS		= -lm

all:		psipred psipass2 seq2mtx pfilt

install:
		/bin/cp psipred psipass2 seq2mtx pfilt ../bin

clean:
		/bin/rm -f psipred psipass2 seq2mtx pfilt

psipred:	sspred_avpred.c ssdefs.h sspred_net.h
		$(CC) sspred_avpred.c $(LIBS) -o psipred

psipass2:	sspred_hmulti.c ssdefs.h sspred_net2.h
		$(CC) sspred_hmulti.c $(LIBS) -o psipass2

seq2mtx:	seq2mtx.c
		$(CC) seq2mtx.c $(LIBS) -o seq2mtx

pfilt:		pfilt.c
		$(CC) pfilt.c $(LIBS) -o pfilt
