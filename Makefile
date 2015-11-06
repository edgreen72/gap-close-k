CC=gcc
CFLAGS=-g -gdwarf-2

read-close-gaps : read-close-gaps.c
	echo "Making read-close-gaps..."
	$(CC) $(CFLAGS) -o read-close-gaps read-close-gaps.c

