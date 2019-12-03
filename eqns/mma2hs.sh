#!/bin/bash

# Convert a Mathematica equation to a Haskell equations

# Start with "InputForm" Mathematica output

perl -p -e 's/\\\[Gamma\]/g/g' |
perl -p -e 's/\[\]//g' |
perl -p -e 's/\[u, v\]//g' |
perl -p -e 's/Derivative\[1, 0\]\[([a-z]*)\]/\1u/g' |
perl -p -e 's/Derivative\[0, 1\]\[([a-z]*)\]/\1v/g' |
perl -p -e 's/Derivative\[1, 1\]\[([a-z]*)\]/\1uv/g' |
perl -p -e 's/Derivative\[2, 0\]\[([a-z]*)\]/\1uu/g' |
perl -p -e 's/Derivative\[0, 2\]\[([a-z]*)\]/\1vv/g' |
tr '\n' ' '
echo
