CFLAGS = -Wall -pedantic -O3
LIBS = -lpthread -lgmp -lm
SOURCE = ep2.c
EXEC = ep2
CC = gcc

ep2: $(SOURCE)
	$(CC) -o $(EXEC) $(CFLAGS) $(SOURCE) $(LIBS) 

.PHONY: clean
clean:
	rm -rf $(EXEC) *~
