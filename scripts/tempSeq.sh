#!/bin/sh
# A script to do kMC simulation at temperature sequence.

for i in $(seq 150 +10 700)
do
mkdir ${i}
done

for i in $(seq 150 +10 700)
do
cp -r  sub.kmc Thermo config.txt ./${i}
done

for i in $(seq 150 +10 700)
do
cd ${i}
echo "7c Temperature = ${i}"
sed -i  "7c Temperature = ${i}" config.txt
cd ..
done

for i in $(seq 210 +10 390)
do
cd ${i}
qsub sub.kmc
cd ..
done