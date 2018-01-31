#!/bin/bash

#gfortran -fcheck=bounds -cpp -ffree-line-length-200 -O0 -g -c data_global.f90
#gfortran -fcheck=bounds -cpp -ffree-line-length-200 -O0 -g  -c mod_grid.f90 
#gfortran -fcheck=bounds -cpp -ffree-line-length-200 -O3 -g  mod_grid.o data_global.o mg.f90

gfortran -cpp -ffree-line-length-200 -O3 -c data_global.f90
gfortran -cpp -ffree-line-length-200 -O3 -c mod_grid.f90 
gfortran -cpp -ffree-line-length-200 -O3  mod_grid.o data_global.o mg.f90
