EXEC = run

FC = gfortran

FFLAGS = -Wall 

SRC = constantes.f90 Legendre.f90 functions.f90 matrix.f90 main.f90
OBJ = $(SRC:.f90=.o)

all: $(EXEC)

$(EXEC): $(OBJ)
	$(FC) $(FFLAGS) -o $@ $(OBJ)

%.o: %.f90
	$(FC) $(FFLAGS) -c $<

clean:
	rm -f *.o *.mod $(EXEC)
