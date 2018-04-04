#Original codes are in /Volumes/Pogoshare/Shared/src/SlideAndOverlay
FC=gfortran
TARGETS=slide-and-match2
BINDIR=../../bin
all: $(TARGETS)
install: all
	install $(TARGETS) $(BINDIR)
%.o: %.F90
	$(FC) -c $(FFLAGS) $< -o $@
%.o: %.f90
	$(FC) -c $(FFLAGS) $< -o $@
slide-and-match: slide-and-match.o mt95.o
	$(FC) $(FFLAGS) $^ -o $@ $(LDFLAGS)
slide-and-match2: slide-and-match2.o mt95.o
	$(FC) $(FFLAGS) $^ -o $@ $(LDFLAGS)
clean:
	-rm $(TARGETS)
	-rm *.o *~ *.mod
