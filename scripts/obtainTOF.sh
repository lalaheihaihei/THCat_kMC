#!/bin/sh

for i in $(seq 350 +10 700)
do
cd ${i}
tail -10 kmc.log | head -1 | awk '{printf "%.3f \n", $3 }'
cd ..
done
