all: MolGen skimmer 206 283 apertures orifice
CFLAGS = -Wall -O3 -I/usr/include/libxml2 -I/usr/include
LIBS = -lgsl -lgslcblas -lm -lpthread -L/usr/lib/x86_64-linux-gnu -lconfig

clean:
	rm *.o trajectoriesswitched trajectoriestime

MolGen : MolGen.c
	g++ -o MolGen MolGen.c $(CFLAGS) $(LIBS)

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
