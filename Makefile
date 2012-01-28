CC=gcc
CFLAGS=-Wall
LDFLAGS=-lsndfile -lfftw3
OBJS=fifo.o watermark.o
TARGET=watermark

$(TARGET): $(OBJS)
	$(CC) -o $(TARGET) $(OBJS) $(LDFLAGS)

test: fifo.o test.o
	$(CC) -o test test.o fifo.o $(LDFLAGS)


clean:
	rm -f $(OBJS) $(TARGET)
