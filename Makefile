FC = gfortran
EXE = run

OP_LINK =
OP_COMP =

SRC = main.f90 \
constant.f90 \
declaration.f90 \
read_input.f90 \
boundary.f90

OBJ = main.o \
constant.o \
declaration.o \
read_input.o \
boundary.o

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

declaration.o: declaration.f90 constant.o

read_input.o: read_input.f90 constant.o declaration.o

boundary.o: boundary.f90 constant.o declaration.o

main.o: main.f90 constant.o declaration.o read_input.o boundary.o
