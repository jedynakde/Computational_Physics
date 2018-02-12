CFLAGS = -g -std=gnu99  -I $(HOME)/include 
LDFLAGS = -g -L$(HOME)/lib -L/usr/X11R6/lib 
#CFLAGS = -O3 -std=gnu99  -I $(HOME)/include 
#LDFLAGS = -O3 -L$(HOME)/lib -L/usr/X11R6/lib 
LIBS = -lgraph -lm -lX11

% : %.c

% : %.o 
	$(CC) $(LDFLAGS) $<  $(LIBS)  -o $@

%.o : %.c
	$(CC) $(CFLAGS) -c $<  -o $@
