#!/usr/bin/env bash

for i in $(seq 1 32);
do
    printf "\nCORES = $i" >> time.txt
	\time -o time.txt -a make test_flag t=$i
done
