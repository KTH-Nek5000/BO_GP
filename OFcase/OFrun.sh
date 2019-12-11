#!/bin/bash

mpirun -np $1 simpleFoam -parallel > log &
wait $!
