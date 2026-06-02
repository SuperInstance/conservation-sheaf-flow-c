CC = gcc
CFLAGS = -Wall -Wextra -std=c11 -O2 -D_POSIX_C_SOURCE=200809L
LDFLAGS = -lm

SRCDIR = src
INCDIR = include
TESTDIR = tests

SRCS = $(SRCDIR)/conservation_sheaf.c $(SRCDIR)/flow_dynamics.c $(SRCDIR)/transport.c $(SRCDIR)/csf_theorem.c
OBJS = $(SRCS:.c=.o)

TARGET = libcsf.a
TEST_BIN = test_csf

.PHONY: all test clean

all: $(TARGET)

$(TARGET): $(OBJS)
	ar rcs $@ $^

$(SRCDIR)/%.o: $(SRCDIR)/%.c
	$(CC) $(CFLAGS) -I$(INCDIR) -c $< -o $@

test: $(TESTDIR)/test_csf.c $(TARGET)
	$(CC) $(CFLAGS) -I$(INCDIR) $< $(TARGET) $(LDFLAGS) -o $(TEST_BIN)
	./$(TEST_BIN)

clean:
	rm -f $(SRCDIR)/*.o $(TARGET) $(TEST_BIN)
