#include "fifo.h"

fifo_t *fifo_init(size_t capacity, size_t sz)
{
	fifo_t *fifo = (fifo_t *) malloc(sizeof(fifo_t));

	fifo->buffer = (void *) malloc(capacity * sz);
	fifo->buffer_end = fifo->buffer + capacity * sz;

	fifo->capacity = capacity;
	fifo->sz = sz;

	fifo->head = fifo->buffer;
	fifo->tail = fifo->buffer;

	return fifo;
}

void fifo_delete(fifo_t *fifo)
{
	free(fifo->buffer);
	fifo->buffer = NULL;
	free(fifo);
}

int fifo_push(fifo_t *fifo, const void *item)
{
	if(fifo->count == fifo->capacity) {
		return -1;
	}

	memcpy(fifo->head, item, fifo->sz);

	fifo->head = (char *) fifo->head + fifo->sz;

	if(fifo->head == fifo->buffer_end)
		fifo->head = fifo->buffer;

	fifo->count++;

	return 0;
}

int fifo_pop(fifo_t *fifo, void *item)
{
	if(fifo->count == 0)
		return -1;

	memcpy(item, fifo->tail, fifo->sz);

	fifo->tail = (char *)fifo->tail + fifo->sz;

	if(fifo->tail == fifo->buffer_end)
		fifo->tail = fifo->buffer;

	fifo->count--;

	return 0;
}