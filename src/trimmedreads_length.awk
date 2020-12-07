#!/usr/bin/awk -f
BEGIN{}
f && !NF{exit} /length/ {f=1} f
END{}