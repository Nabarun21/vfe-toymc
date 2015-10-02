#!/bin/bash

for f in all_plots/*;
  do
     [ -d $f ] && cd "$f" && echo Entering into $f and removing all pdf files
     rm *.pdf
     cd ../..
  done;