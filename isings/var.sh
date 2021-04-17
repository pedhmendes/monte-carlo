#!/bin/bash

echo "" >> $1
echo "" >> $1

#variancia
awk 'BEGIN{sum=0; summ=0}{arr[NR]=$1; sum+=$1}END{media=sum/NR; for(i=1; i<=NR; i++){summ+=(arr[i]-media)^2}; summ=sqrt(summ/NR); print media "\t" summ}' $1 >> $1
