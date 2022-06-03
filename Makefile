FC = gfortran
EXE = run

OP_LINK =
OP_COMP =

SRC = main.f90 \
constant.f90

OBJ = main.o \
constant.o

all: $(EXE)

compile: $(OBJ)

$(EXE): $(OBJ)
	$(FC) -o $(EXE) $(OBJ) $(OP_LINK)

$(OBJ):
	$(FC) -c $(OP_COMP) $(SRC_DIR)$(@:.o=.f90) -o $@

.PHONY = clean

clean:
	rm $(OBJ) *.mod

cleaner:
	$(clean)
	rm $(EXE)

constant.o: constant.f90

main.o: main.f90 constant.o
