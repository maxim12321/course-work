#! /usr/bin/env sh

mkdir -p build
cd ./build
cmake .. -B ./
make
cp solver ../solver
cp solver ../../gui/resources/solver
