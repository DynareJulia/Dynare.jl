#!/bin/bash

export OMP_NUM_THREADS=1

julia -t 20 --project=../../../work-1.10 bench_irbc.jl
#julia -t 19 --project=../../../work-1.10 bench_irbc.jl

#julia -t 10 --project=../../../work-1.10 bench_irbc.jl
