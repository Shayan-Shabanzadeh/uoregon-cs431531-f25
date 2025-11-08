#!/bin/bash

echo "Threads,Serial,P1_Time,P2_Time" > prefix_results.csv

for f in output/prefix_*threads.out; do
    threads=$(echo $f | sed -E 's/.*prefix_([0-9]+)threads.*/\1/')

    serial=$(grep "O(N-1)" "$f" | awk -F':' '{print $2}' | tr -d ' ()s')
    p1=$(grep "O(NlogN)" "$f" | awk -F':' '{print $2}' | tr -d ' ()s')
    p2=$(grep "2(N-1)" "$f" | awk -F':' '{print $2}' | tr -d ' ()s')

    echo "$threads,$serial,$p1,$p2" >> prefix_results.csv
done

echo "prefix_results.csv generated."