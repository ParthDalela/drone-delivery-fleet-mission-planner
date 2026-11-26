CC = gcc
CFLAGS = -Wall -Wextra -std=c17 -O3 -Iinclude

SRC = src/algorithms.c src/simulator.c src/main.c
OBJ = $(SRC:.c=.o)
TARGET = drone_planner

all: $(TARGET)

$(TARGET): $(OBJ)
	$(CC) $(CFLAGS) -o $@ $^ -lm

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -f src/*.o $(TARGET)

run: all
	./$(TARGET)

benchmark: all
	./$(TARGET) --benchmark
.PHONY: all clean run benchmark