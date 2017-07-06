CXXFLAGS=-O3 -fopenmp -march=native -mtune=native -lgmpxx -lgmp -g

.PHONY: all
all: find-possible-multiplicities gen-krylov minpoly krylov-to-eigenspaces trace-a2 trace-a3-a4
