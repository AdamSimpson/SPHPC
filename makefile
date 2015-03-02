CC=mpicc
CFLAGS=-g -O3 -Wall
LDLIBS=-lm

SRC_DIR=../src
BIN_DIR=../bin
SRC_FILES=geometry.c neighbors.c fileio.c communication.c fluid.c simulation.c
OBJ_FILES=$(subst .c,.o, $(SRC_FILES))

all: sph

%.o: $(SRC_DIR)/%.c
	$(CC) -c $(CFLAGS) -o $@ $^

sph: $(OBJ_FILES)
	mkdir -p $(BIN_DIR)
	$(CC) $(CFLAGS) $(LDLIBS) -o $(BIN_DIR)/$@ $^
