CFLAGS=-O
LDFLAGS=-O
LIBS=-lsndfile -lsamplerate -lm

noaa-apt:	noaa-apt.o
		$(CC) $(LDFLAGS) -o $@ $< $(LIBS)
