#!/bin/sh
#nJob1=$1
nJob1=0
maxJob1=$1
#nJob2=$3
nJob2=0
maxJob2=$2

nJob3=$3
nNumber=0

DEBUG=true
while [ $nJob1 -le $maxJob1 ]
do
    echo $nJob1
    nJob2=0
    while [ $nJob2 -le $maxJob2 ]
    do
        echo $nJob2
        if ${DEBUG}; then
        echo "star-submit-template -template submitPicoDstJobs3.xml -entities cutID=$nJob1,variationID=$nJob2,versionID=$nJob3"
        star-submit-template -template submitPicoDstJobs3.xml -entities cutID=$nJob1,variationID=$nJob2,versionID=$nJob3
        ((nNumber++))
        echo "total submits: $nNumber"
        fi
        if(($nJob1 == 0 && $nJob2 ==0)); then
        break
        fi
        if(($nJob1 == 2 && $nJob2 ==4)); then
        break
        fi
        if(($nJob1 == 3 && $nJob2 ==4)); then
        break
        fi
        if(($nJob1 == 20 && $nJob2 ==4)); then
        break
        fi
        if(($nJob1 == 21 && $nJob2 ==4)); then
        break
        fi
        ((nJob2++))
    done

    ((nJob1++))
done
#echo $nJob1
