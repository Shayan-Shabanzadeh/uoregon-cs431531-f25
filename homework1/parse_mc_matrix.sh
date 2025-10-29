#!/bin/bash
out="results_mc_matrix.csv"
echo "threads,guesses,mean_time,mean_abs_err" > $out

for f in output/mc_matrix_*.out; do
    threads=$(grep "BEGIN MATRIX" "$f" \
              | sed 's/.*threads=\([0-9]*\).*/\1/')

    grep -n "BEGIN BLOCK" "$f" | while read -r line; do
        lineno=$(echo "$line" | cut -d: -f1)
        guesses=$(echo "$line" | sed 's/.*guesses=\([0-9]*\).*/\1/')

        block=$(sed -n "$((lineno+1)),$ p" "$f" | sed '/BEGIN /Q')

        times=$(echo "$block" | grep "Time to calculate Pi in //" | grep "guesses" | awk '{print $NF}')
        [ -z "$times" ] && continue
        mean_time=$(echo "$times" | awk '{s+=$1} END{printf "%.6f", s/NR}')

        pi_vals=$(echo "$block" | grep "^Pi is" | awk '{print $NF}')
        mean_err=$(echo "$pi_vals" | awk -v pi=3.141592653589793 '{d=$1-pi; if(d<0)d=-d; s+=d} END{printf "%.10f",s/NR}')

        echo "$threads,$guesses,$mean_time,$mean_err" >> $out
    done
done

echo "Wrote $out"