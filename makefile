CC=gcc
DLDSRC= src/*.c
CFLAGS=
DEBUG= -g

all:
ifeq ($(CFLAGS),-DUSE_OPENMP)
	$(CC) -O3 -fopenmp $(CFLAGS) $(DEBUG) -Wall $(DLDSRC) -o dld -lm -llapack
else ifeq ($(CFLAGS),-DMINGW32)
	$(CC) -O3 -fopenmp $(CFLAGS) $(DEBUG) -Wall $(DLDSRC) win/*.dll -o win/dld -lm
else
	$(CC) -O3 $(DEBUG) -Wall $(DLDSRC) -o dld -lm -llapack
endif

