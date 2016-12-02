#!/bin/sh
# A script to do kMC simulation at temperature sequence.

for i in $(seq 150 +10 700)
do
mkdir ${i}
done

for i in $(seq 500 +10 700)
do
cp -r  sub.kmc ./${i}
done

for i in $(seq 150 +10 700)
do
cd ${i}
echo "5c Temperature = ${i}"
sed -i  "5c Temperature = ${i}" config.txt
cd ..
done

for i in $(seq 150 +10 700)
do
cd ${i}
qsub sub.kmc
cd ..
done