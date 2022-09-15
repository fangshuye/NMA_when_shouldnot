#!/bin/bash

for r in {1..4}
do
	sbatch TSA_V5.sh $r
done

