#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct fifo {
	void *buffer;
  void *buffer_end;

  size_t capacity;
  size_t count;
  size_t sz;

  void *head;
  void *tail;
} fifo_t;

fifo_t *fifo_init(size_t capacity, size_t sz);
int	fifo_push(fifo_t *fifo, const void *item);
int	fifo_pop(fifo_t *fifo, void *item);