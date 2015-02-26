#============================================================================#
#           GeCo & GeDe 2014-2015 IEETA/DETI, UNIVERSITY OF AVEIRO           #
#============================================================================#
BIN    = .
CC     = gcc
CPLP   = -g
#-----------------------------------------------------------------------------
CFLAGS = -O3 -Wall $(CPLP)
#-----------------------------------------------------------------------------
LIBS   = -lm
DEPS   = defs.h
PROGS  = $(BIN)/GeCo $(BIN)/GeDe
OBJS   = mem.o common.o buffer.o context.o bitio.o arith.o \
         arith_aux.o 
#-----------------------------------------------------------------------------
all:
	$(MAKE) progs
progs: $(PROGS)
$(BIN)/GeCo: geco.c $(DEPS) $(OBJS)
	$(CC) $(CFLAGS) -o $(BIN)/GeCo geco.c $(OBJS) $(LIBS)
$(BIN)/GeDe: gede.c $(DEPS) $(OBJS)
	$(CC) $(CFLAGS) -o $(BIN)/GeDe gede.c $(OBJS) $(LIBS)
mem.o: mem.c mem.h $(DEPS)
	$(CC) -c $(CFLAGS) mem.c
buffer.o: buffer.c buffer.h $(DEPS)
	$(CC) -c $(CFLAGS) buffer.c
common.o: common.c common.h $(DEPS)
	$(CC) -c $(CFLAGS) common.c
context.o: context.c context.h $(DEPS)
	$(CC) -c $(CFLAGS) context.c
#alpha.o: alpha.c alpha.h $(DEPS)
#	$(CC) -c $(CFLAGS) alpha.c
#gfcm.o: gfcm.c gfcm.h $(DEPS)
#	$(CC) -c $(CFLAGS) gfcm.c
bitio.o: bitio.c bitio.h
	$(CC) -c $(CFLAGS) bitio.c
arith.o: arith.c arith.h
	$(CC) -c $(CFLAGS) arith.c
arith_aux.o: arith_aux.c arith_aux.h
	$(CC) -c $(CFLAGS) arith_aux.c

clean:
	/bin/rm -f *.o
cleanall:
	/bin/rm -f *.o $(PROGS)
#-----------------------------------------------------------------------------
