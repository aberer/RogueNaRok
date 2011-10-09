CC = gcc 

CFLAGS =  -D_GNU_SOURCE -Wall -DNDEBUG #    -DNEW_SUBSET  #   #
LFLAGS =   -lm

#  -funroll-loops 
ifeq ($(mode), debug)
 CFLAGS += -g
else
 CFLAGS += -march=native -O3
ifeq ($(mode), profile)
 CFLAGS += -pg -g
endif
endif
ifeq ($(mode), parallel)
CFLAGS += -DPARALLEL -DPORTABLE_PTHREADS 
LFLAGS += -pthread 
endif
ifeq ($(mode), parallelDebug)
CFLAGS += -DPARALLEL -g
LFLAGS += -pthread 
endif

RM = rm -fr

TARGETS = RogueNaRok rnr-prune rnr-lsi  rnr-tii  rnr-mast 
# TESTS=rnr-test

all :  $(TARGETS)

rnr-objs = common.o RogueNaRok.o  Tree.o BitVector.o HashTable.o List.o Array.o  Dropset.o ProfileElem.o legacy.o newFunctions.o parallel.o
lsi-objs = rnr-lsi.o common.o Tree.o BitVector.o   HashTable.o legacy.o newFunctions.o List.o
tii-objs = rnr-tii.o common.o BitVector.o Tree.o HashTable.o List.o legacy.o newFunctions.o 
mast-objs = rnr-mast.o common.o List.o Tree.o BitVector.o HashTable.o legacy.o newFunctions.o
prune-objs = rnr-prune.o common.o Tree.o BitVector.o HashTable.o  legacy.o newFunctions.o List.o
# rnr-test-objs = Dropset.o ProfileElem.o  rnr-test.o List.o HashTable.o common.o  BitVector.o

rnr-lsi: $(lsi-objs)
	$(CC) $(LFLAGS) -o $@   $^ $(CFLAGS)
rnr-tii: $(tii-objs)
	$(CC) $(LFLAGS) -o $@   $^ $(CFLAGS)
rnr-mast: $(mast-objs)
	$(CC) $(LFLAGS) -o $@   $^ $(CFLAGS)
rnr-prune: $(prune-objs)
	$(CC) $(LFLAGS) -o $@   $^ $(CFLAGS)
RogueNaRok: $(rnr-objs)
	$(CC) $(LFLAGS) -o $@   $^ $(CFLAGS)

# rnr-test : $(rnr-test-objs)
# 	$(CC) $(LFLAGS) -o $@   $^ $(CFLAGS) -lcunit


%.o : %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

clean : 
	$(RM) $(rnr-objs) $(lsi-objs) $(tii-objs) $(mast-objs) $(prune-objs) $(TARGETS)  $(TESTS) $(rnr-test-objs)

# test : rnr-test
