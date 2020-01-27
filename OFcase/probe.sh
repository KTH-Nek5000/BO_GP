#! /bin/sh

#deleting all parenthesis everywhere in the file, using a temporary "tmp" file
cat postProcessing/probes/0/U | tr -d "()" > postProcessing/probes/0/tmp

# call plot file
gnuplot probe.plt
