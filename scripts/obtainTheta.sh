#!/bin/sh

for i in $(seq 150 +10 700)
do
cd ${i}
tail -5 kmc.log | head -1 | awk '{printf "%.3f \n", $7 }'
cd ..
done
