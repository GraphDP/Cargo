#!/bin/bash -x

#code            data             eps ItrNum 
./cargo /data/${1} 0.5 30
./cargo /data/${1} 1 30
./cargo /data/${1} 1.5 30
./cargo /data/${1} 2 30
./cargo /data/${1} 2.5 30
./cargo /data/${1} 3 30
cd ../
