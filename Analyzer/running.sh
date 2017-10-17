#!/bin/bash

echo "cmsRun = "$[$[`ps aux | grep cmsRun | grep anal | wc -l`]/2]

dir=$1
sum1=0
sum2=0
for file in `ls $dir/log*`
do
  n1=`grep NEVENTS $file | tail -n1 | awk '{print $4}'`
  n2=`grep NEVENTS $file | tail -n1 | awk '{print $8}'`
  sum1=$[$n1+$sum1]
  sum2=$[$n2+$sum2]
  echo n=$n1 K [${n2}]
done
echo sum=$sum1 K [$sum2]
