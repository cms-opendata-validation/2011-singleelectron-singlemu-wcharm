#!/bin/bash

dir=$1
sum1=0
sum2=0
for file in `ls $dir/log*`
do
  n1=`grep Processed $file | tail -n1 | awk '{print $2}'`
  n2=`grep Processed $file | tail -n1 | awk '{print $5}'`
  sum1=$[$n1+$sum1]
  sum2=$[$n2+$sum2]
  echo n=$n1 [${n2}]
done
echo sum=$sum1 [$sum2]
