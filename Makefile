all: MolGen skimmer 206 283 apertures orifice transcooling
CFLAGS = -Wall -O3 -I/usr/include/libxml2 -I/usr/include 
LIBS = -L/usr/lib/x86_64-linux-gnu -lgsl -lgslcblas -lm -lpthread -lconfig 

MolGen : MolGen.c
	gcc -o MolGen MolGen.c $(CFLAGS) $(LIBS)

skimmer : skimmer.c
	gcc -o skimmer skimmer.c $(LIBS)

206 : 206.c
	gcc -o 206 206.c $(LIBS)

283 : 283.c
	gcc -o 283 283.c $(LIBS)

apertures : apertures.c
	gcc -o apertures apertures.c $(LIBS)

orifice : orifice.c
	gcc -o orifice orifice.c $(LIBS)

transcooling: transcooling.c
	gcc -o transcooling transcooling.c $(LIBS)
