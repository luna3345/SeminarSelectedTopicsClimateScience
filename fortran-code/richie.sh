#compile the main program richie
f77 -c richie.f

#compile the subroutines. 
f77 -c getdat.f tablot.f impfil.f hfilt1.f ricdat.f dolph.f

#load the programs
f77 -o richie.x                  richie.o            \
                                 getdat.o            \
                                 tablot.o            \
                                 impfil.o            \
                                 hfilt1.o            \
                                 dolph.o             \
                                 ricdat.o

# run it
#richie.x > richie.lpt
