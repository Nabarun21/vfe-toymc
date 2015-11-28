#!/bin/bash

usage='Usage: -pt <pulse type>  -pw <pulse width> -ns <num_samples> -nf <sampling frequency> -ph <phase position> (-cpt <old pulse type> -cns <old no. of samples> -cnf<old sampling freq>) '

args=`getopt rdlp: -- "$@"`
if test $? != 0
     then
         echo $usage
         exit 1
fi

eval set -- "$args"


for i
 do
    case "$i" in
      -pt) shift; pulse_type=$2;shift;;
      -ns) shift; num_samples=$2;shift;;
      -nf) shift; sampl_freq=$2;shift;;
      -ph) shift; phase_pos=$2;shift;;
      -rp) shift; risepos=$2;shift;;
      -cpt) shift; old_pulse_type=$2;shift;;
      -cns) shift; old_pulse_nsample=$2;shift;;
      -cnf) shift; old_pulse_nfreq=$2;shift;;
    esac
done


echo "replacing pulse shape definitions infile"

if [ "X"${old_pulse_type} != "X" ]
then
    add2=$( echo $(echo $old_pulse_type) | awk '{print toupper($0)}')
    repl2=$( echo $(echo $pulse_type) | awk '{print toupper($0)}')
    sed  -i.bak -e "s/$add2/$repl2/g" extractsigma.C 
    
fi

echo "check"

if [ "X"${old_pulse_nsample} != "X" ]
then
    add=$old_pulse_nsample","$old_pulse_nfreq
    repl=$num_samples","$sampl_freq
    sed  -i.bak -e "s/$add/$repl/g" extractsigma.C
fi


echo "recompiling"

g++ -o Pu_MyScript_new pu_myscript.cc PulseChiSqSNNLS.cc -std=c++11 `root-config --cflags --glibs`


echo "calling simulator"

./Pu_MyScript_new $pulse_type $num_samples $sampl_freq $phase_pos $risepos