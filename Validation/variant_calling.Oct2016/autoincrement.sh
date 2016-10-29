#!/bin/bash    

#This files increments the variable n automatically after each run. It would be handy to identify consequent runs of the workflow.

n=9;#the variable that I want to be incremented
next_n=$[$n+1]
sed -i "/#the variable that I want to be incremented$/s/=.*#/=$next_n;#/" ${0}
echo $n




