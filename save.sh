#!/bin/bash

if [ $# -ne 1 ]; then
   echo "usage: ./save.sh saveDirectoryName"
   exit 1
fi

target="../data/$1"
if [ ! -d $target ]; then
   mkdir $target
else
   rm $target/*
fi

echo "save to " $target

cp figs $target -r
cp OFpost/beta.dat $target
cp gpOptim/workDir/gpList.dat $target
cp log $target

