FC     = gfortran
FCFLAGS = -Dstandalone -fdefault-real-8 -fdefault-double-8

SRC=machine.F90 ../module_sf_mynnsfc.F90 module_sf_mynnsfc_tests.F90 driver.F90
OBJ=machine.o module_sf_mynnsfc.F90.o module_sf_mynnsfc_tests.o driver.o

%.o: %.F90
	$(FC) $(FCFLAGS) -c $(SRC)

driver: $(OBJ)
	$(FC) $(FFLAGS) -o $@ $(OBJ)

clean:
	@rm -f *.mod *.o driver
