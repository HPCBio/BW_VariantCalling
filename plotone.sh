#!/bin/bash
umask 0027

USAGE="Usage : $0 <datafile>"

if [ "$#" -eq "0" ] 
then
   echo $USAGE
   exit 1
fi

input=$1
data=${input}.dat
tmp=${input}.tmp
plot=${input}.gp
pngfig=${input}.png

tail -n +7 $input | tr "," " " > $tmp
sed ':a;N;$!ba;s/buffers\n//g' $tmp > $tmp.new
sed ':a;N;$!ba;s/cached\n//g' $tmp.new > $tmp
#mv $tmp.new $tmp
rm $tmp.new 
size=`wc $tmp | awk '{ print $1 }'`
start_t=`head -n 1 $tmp | awk '{ print $1 }'`
end_t=`tail -n 1 $tmp | awk '{ print $1 }'`
title=`grep cmd $input | awk '{ print $2 }'`
title=`echo "$input : $title"`
cat $tmp | awk -v start_t=$start_t '{ print ($1-start_t) " " $2 " " $3 " " $4 " " $5 " " $6 " " $7 }' > $data

xr_max=`tail -n 1 $data | awk '{ print $1 }'`


echo "#gplot file" > $plot
echo "set term png size 1024,680 enhanced " >> $plot
echo "set output \"$pngfig\"" >> $plot
#echo "set ytics font \", 8\"" >> $plot

echo "set multiplot layout 3,1 title \"${title}\" " >> $plot
echo "set tmargin 0" >> $plot
echo "unset xtics" >> $plot
echo "set xrange [0:$xr_max]" >> $plot
echo "set lmarg 13" >> $plot
echo "set rmarg 10" >> $plot

#echo "set title \"Memory Usage\"" >> $plot
echo "set key center right nobox" >> $plot
echo "set ytics" >> $plot
echo "set y2tics" >> $plot
echo "set ylabel \"Memory usage (kB)\"" >> $plot
#echo "set y2label \"Memory usage (kB)\"" >> $plot
echo "plot '$data' using 1:4 with lp title \"VmSize (kB)\", '$data' using 1:5 with lp title \"VmRSS (kB)\"" >> $plot

#echo "set title \"I/O\"" >> $plot
echo "set key center right nobox" >> $plot
echo "set ytics" >> $plot
echo "set y2tics" >> $plot
echo "set ylabel \"I/O (bytes/sec)\"" >> $plot
#echo "set y2label \"I/O (bytes/sec)\"" >> $plot
echo "plot '$data' using 1:6 with lp title \"rchar(bytes/sec)\", '$data' using 1:7 with lp title \"wchar(bytes/sec)\"" >> $plot

#echo "set title \"threads\"" >> $plot
echo "set key center right nobox" >> $plot
echo "set ytics nomirror" >> $plot
echo "set y2tics" >> $plot
echo "set yrange [0:100]" >> $plot
echo "set ylabel \"# of threads\"" >> $plot
#echo "set autoscale  y" >> $plot
echo "set autoscale y2" >> $plot
echo "set y2label \"CPU load (%)\"" >> $plot
echo "set xtics nomirror rotate by -45 " >> $plot
echo "set xlabel \"time (s)\"" >> $plot
echo "plot '$data' using 1:2 with lp title \"# threads\" axes x1y1, '$data' using 1:3 with lp title \"% CPU\" axes x1y2" >> $plot

echo "unset multiplot" >> $plot

gnuplot $plot
