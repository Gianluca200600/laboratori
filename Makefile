CC = gcc
CFLAGS = -Wall -Wextra

LIBRARY_SRCS = library/blas.c library/odesolver.c library/rkm.c library/solver.c library/examples.c
LIBRARY_OBJS = $(LIBRARY_SRCS:.c=.o)
LIBRARY_A = library/liblibrary.a

# All targets
TARGETS = l_2_1 l_2_2 l_2_3 l_3_1 l_3_2 l_3_3 l_4_1

# Default: build all
all: $(TARGETS)

# Create static library
$(LIBRARY_A): $(LIBRARY_OBJS)
	ar rcs $@ $^

# Compile library sources to objects
library/%.o: library/%.c
	$(CC) $(CFLAGS) -c $< -o $@

# Lab 2
l_2_1: l02/l_2_1.c $(LIBRARY_A)
	$(CC) $(CFLAGS) $< -L library/ -llibrary -o ./l02/l_2_1

l_2_2: l02/l_2_2.c $(LIBRARY_A)
	$(CC) $(CFLAGS) $< -L library/ -llibrary -o ./l02/l_2_2

l_2_3: l02/l_2_3.c $(LIBRARY_A)
	$(CC) $(CFLAGS) $< -L library/ -llibrary -o ./l02/l_2_3

# Lab 3
l_3_1: l03/l_3_1.c $(LIBRARY_A)
	$(CC) $(CFLAGS) $< -L library/ -llibrary -o ./l03/l_3_1

l_3_2: l03/l_3_2.c $(LIBRARY_A)
	$(CC) $(CFLAGS) $< -L library/ -llibrary -o ./l03/l_3_2

l_3_3: l03/l_3_3.c $(LIBRARY_A)
	$(CC) $(CFLAGS) $< -L library/ -llibrary -o ./l03/l_3_3

# Lab 4
l_4_1: l04/l_4_1.c $(LIBRARY_A)
	$(CC) $(CFLAGS) $< -L library/ -llibrary -o ./l04/l_4_1

l_4_2: l04/l_4_2.c $(LIBRARY_A)
	$(CC) $(CFLAGS) $< -L library/ -llibrary -o ./l04/l_4_2


clean:
	rm -f library/*.o library/*.a $(TARGETS) $(TARGETS:=.exe)

.PHONY: all clean
