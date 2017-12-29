CC=gcc
CFLAGS=-O2 -gdwarf-2

read-close-gaps : read-close-gaps.c
	echo "Making read-close-gaps..."
	$(CC) $(CFLAGS) -o read-close-gaps read-close-gaps.c

describe-assembly : describe-assembly.c
	echo "Making describe-assembly..."
	$(CC) $(CFLAGS) -o describe-assembly describe-assembly.c -lm 
