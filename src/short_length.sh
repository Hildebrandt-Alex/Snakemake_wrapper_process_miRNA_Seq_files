#!/bin/sh

paste - - - - |
awk -F'[\t ]' '!x[$3]++ {print length($3)}' |
awk '{dups[$0]++} END{for (num in dups) {print num,dups[num]}}' -  |
sort -k1 -n - ;
