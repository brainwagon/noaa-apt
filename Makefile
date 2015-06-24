CFLAGS=-O -I/opt/local/include
LDFLAGS=-O
LIBS=-L/opt/local/lib -lsndfile -lsamplerate -lm

noaa-apt:	noaa-apt.o
		$(CC) $(LDFLAGS) -o $@ $< $(LIBS)
