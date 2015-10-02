#!/bin/bash

for f in all_plots/*;
do
    [ -d $f ] && cd "$f" && echo Entering into $f and copying
    cp -f *.png /afs/cern.ch/user/n/ndev/www/logbooks/vfe-logbook
    cd ../..
done;