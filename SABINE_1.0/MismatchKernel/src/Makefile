
# MITRA MakeFile


ifeq ($(MODE),debug)
	#Debug flags
	COMPFLAGS =  -c -g -pg -Wall
	CC = cc
	CLIBS = -lm -g -pg
else 
ifeq ($(MODE),memwatch)
	#Memwatch flags
	COMPFLAGS = -DMEMWATCH -DMW_STDIO -c -g -Wall -pg
	CC = cc
	CLIBS = -DMEMWATCH -DMW_STDIO -lm  -g -Wall -pg
else
	#Efficient Flags
	COMPFLAGS = -c -O2 
	CC = gcc
	CLIBS = -lm  -O2
endif
endif


INCLUDES = HashTable.h mitra.h MiscUtil.h SymbolTable.h SymbolTable.h Globals.h  Input.h



string-kernel 	: string-kernel.o HashTable.o MiscUtil.o SymbolTable.o Globals.o  Input.o memwatch.o
	${CC} ${CLIBS} -o  string-kernel string-kernel.o HashTable.o MiscUtil.o SymbolTable.o Globals.o  Input.o memwatch.o

%.o	: %.c
	${CC} ${COMPFLAGS} $< -o $@




TAGS  :
	etags *.c *.h

clean :
	rm string-kernel *.o


