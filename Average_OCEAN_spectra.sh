#!/bin/bash
## To use: Average_OCEAN_spectra.sh  $folder $absorbing_atom $starting_NO $ending_NO $n,l quantum numbers
## e.g. Average_OCEAN_spectra.sh ./Fe2O3 Fe 1 6 1s
## 3 photon files are included (3 incidentor momentum transfer directions)
here=$PWD
cd $1
ctr_atm=$2  
[ ${#ctr_atm} -eq 1 ] && ctr_atm=$2\_
gaga=`printf "%02d" $3`
cat absspct_$ctr_atm.00$gaga\_$5\_01 | awk '{print $1}' > Eng1.dat
#touch Combo_$ctr_atm_$$ctr_atmabs1.dat  Combo_$ctr_atm_$$ctr_atmabs2.dat
for ii in `seq -w $3 $4`
do
  echo $ii
  i=`printf "%02d" ${ii#0}`
  [ -f buffer$i.dat ] && rm buffer$i.dat
  touch buffer$i.dat
  [ -f Combo_$2\_$i.dat ] && rm Combo_$2\_$i.dat
  # Averaging over photon files 1 2 3
  for j in 01 02 03
  do
    cat absspct_$ctr_atm.00$i\_$5\_$j | awk '{print $2}' > tmp1.dat
    paste buffer$i.dat tmp1.dat > Combo_$2\_$i.dat
    cp Combo_$2\_$i.dat buffer$i.dat
  done
  cat Combo_$2\_$i.dat | awk '{print ($1+$2+$3)*0.3333333}' > ave_$2\_$i.dat
  paste Eng1.dat ave_$2\_$i.dat > gaga
  mv gaga ave_$2\_$i.dat
done

[ -f buffer.dat ] && rm buffer.dat
touch buffer.dat
for jj in `seq -w $3 $4`
do
  j=`printf "%02d" ${jj#0}`
  cat ave_$2\_$j.dat | awk '{print $2}' > tmp1.dat
  paste buffer.dat tmp1.dat > Combo_$2\_all.dat
  cp Combo_$2\_all.dat buffer.dat
done

$here/Calc_Ocean_Ave.exe $2
paste Eng1.dat Ave_$2\_all.dat > gaga
mv gaga Ave_$2\_all.dat
rm buffer*.dat

paste Eng1.dat Combo_$2\_all.dat > gaga
mv gaga Combo_$2\_all.dat
