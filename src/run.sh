#!/bin/bash

RESULT_DIR=$PWD/Results
rm -rf $RESULT_DIR
mkdir -p $RESULT_DIR

if [ $1 == "cluster" ]; then
    module load gcc/5.3.0
    module load gnu/openmpi_eth/1.8.2
    cd src
fi


N=10
while [ $N -lt 10000 ]; do
    RESULT_FILE=$RESULT_DIR/n_$N.txt

    make clean
    sed -i '' 's/#define N [0-9]*/#define N '"$N"'/g' UTILS/utils.h
    make

    ./bin/seq >> $RESULT_FILE
    echo " " >> $RESULT_FILE

    ./bin/rb_seq >> $RESULT_FILE
    echo " " >> $RESULT_FILE

    mpirun -np 8 bin/rb_par >> $RESULT_FILE

    let N=N*10
done
