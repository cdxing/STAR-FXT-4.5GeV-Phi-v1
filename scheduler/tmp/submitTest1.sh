#!/bin/sh
#nJob1=$1
nJob1=0
maxJob1=$1
#nJob2=$3
# nJob2=0
# maxJob2=$2
DEBUG=true
while [ $nJob1 -le $maxJob1 ]
do
    inputPar1=$nJob1

    echo $inputPar1
    if ${DEBUG}; then
    star-submit-template -template submitPicoDstJobs2.xml -entities myFileName=vertexZcut_v2_$nJob1,cutTestPar1=$inputPar1

    fi
    ((nJob1++))
done
#echo $nJob1
