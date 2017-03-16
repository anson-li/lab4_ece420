###
# Designating default compiler and compiler values
###
CC = gcc
CFLAGS = -g -Wall
LIBS = -lpthread -lm

###
# Setting default make
###
default: datatrim serialtester main
all: datatrim serialtester main

###
# Support
###
datatrim:
	$(CC) -o datatrim datatrim.c

serialtester:
	$(CC) -o serialtester serialtester.c Lab4_IO.c

main:
	mpicc -g -Wall -o main main.c LAB4_IO.c 

###
# Clean process
###
clean:
	rm -f datatrim serialtester main
