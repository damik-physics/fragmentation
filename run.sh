#!/bin/bash

for i in {1..20}
do 
	make -f makefile_fragmentation SITES=$1 K=1 && ulimit -s unlimited && ./frag_L$1_k1
	mv ./ti_energies* ti_states* energies* states* *entropy* times* ipr* lvl_stat_* ./output
done 

