ifort -qopenmp -o test LapackRoutines.f90 Hamiltonian.f90 plotBands.f90 -qmkl -lpthread -lm
