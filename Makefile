CC=gcc
LD=gcc

OPT_LEVEL=-O0 -g
#OPT_LEVEL=-flto -mavx2 -Ofast -fomit-frame-pointer

CC_OPTS=$(OPT_LEVEL) -fPIC -Wall -Wno-unused-variable -Werror -Wfatal-errors
LD_OPTS=$(OPT_LEVEL) -lm -lc

HEADERS:=$(wildcard *.h) Makefile

.SUFFIXES:

TARGETS=spmv_test
all : $(TARGETS)

%.o: %.c $(HEADERS)
	$(CC) -c $< -o $@ $(CC_OPTS)

%: %.o spmv.o mmul.o rand.o
	$(LD) $^ -o $@ $(LD_OPTS)

clean:
	rm -rf $(TARGETS) *.o
