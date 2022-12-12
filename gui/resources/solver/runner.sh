#!/usr/bin/env bash

#./solver input.txt output.txt explicit
mpiexec -n 4 solver input.txt output.txt explicit
