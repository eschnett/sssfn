#!/bin/bash

# Convert a Mathematica equation to a Haskell equations

# Start with "InputForm" Mathematica output

tr '\n' ' ' |
perl -p -e 's/\\\[Gamma\]/g/g' |
perl -p -e 's/\[\]//g' |
perl -p -e 's/\[u,\s*v\]//g' |
perl -p -e 's/Derivative\[1, 0\]\[([a-z]*)\]/\1u/g' |
perl -p -e 's/Derivative\[0, 1\]\[([a-z]*)\]/\1v/g' |
perl -p -e 's/Derivative\[1, 1\]\[([a-z]*)\]/\1uv/g' |
perl -p -e 's/Derivative\[2, 0\]\[([a-z]*)\]/\1uu/g' |
perl -p -e 's/Derivative\[0, 2\]\[([a-z]*)\]/\1vv/g' |
perl -p -e 's/r1du\[([a-z]*)\]/r1\1u/g' |
perl -p -e 's/r1dv\[([a-z]*)\]/r1\1v/g' |
perl -p -e 's/r1dr\[([a-z]*)\]/r1\1r/g' |
perl -p -e 's/r2\[([a-z]*)\]/r2\1/g' |
perl -p -e 's/z\[([a-z]*)\]/z\1/g'
echo
