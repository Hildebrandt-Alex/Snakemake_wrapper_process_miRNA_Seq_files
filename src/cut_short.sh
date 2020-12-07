#!/bin/sh

grep -A 14 -P  '\t0\t' - | 
sed 's/\t//';