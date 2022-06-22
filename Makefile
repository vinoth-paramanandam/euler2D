FC = gfortran
EXE = run

OP_LINK = -O3
OP_COMP = -O3

SRC = main.f90 \
constant.f90 \
declaration.f90 \
read_input.f90 \
grid.f90 \
boundary.f90 \
initial.f90 \
misc.f90 \
odesolver.f90 \
flux.f90

OBJ = main.o \
constant.o \
declaration.o \
read_input.o \
grid.o \
boundary.o \
initial.o \
misc.o \
odesolver.o \
flux.o

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

grid.o: grid.f90 constant.o declaration.o

boundary.o: boundary.f90 constant.o declaration.o

initial.o: initial.f90 constant.o declaration.o misc.o

misc.o: misc.f90 constant.o declaration.o

odesolver.o: odesolver.f90 constant.o declaration.o

flux.o: flux.f90 constant.o declaration.o

main.o: main.f90 constant.o declaration.o read_input.o grid.o \
	initial.o misc.o boundary.o odesolver.o flux.o
