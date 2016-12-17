#!/bin/sh

for i in $(seq 230 +10 700)
do
cd ${i}
tail -11 kmc.log | head -1 | awk '{printf "%.3f \n", $9 }'
cd ..
done
