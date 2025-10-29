#!/bin/bash

out="results_integ_critical.csv"
echo "method,variant,threads,steps,mean_time,mean_abs_err" > $out

for f in output/pi_critical_matrix_*.out; do
    read m v t <<<$(grep "BEGIN MATRIX" "$f" | awk '{print $3,$4,$5}')
    method=${m#method=}
    variant=${v#variant=}
    threads=${t#threads=}

    [ "$variant" != "critical" ] && continue

    grep -n "BEGIN BLOCK" "$f" | while read -r line; do
        lineno=$(echo "$line" | cut -d: -f1)
        steps=$(echo "$line" | sed 's/.*steps=\([0-9]*\).*/\1/')
        block=$(sed -n "$((lineno+1)),$ p" "$f" | sed '/BEGIN /Q')

        times=$(echo "$block" | grep "Time to calculate Pi in //" | awk '{print $NF}')
        [ -z "$times" ] && continue
        mean_time=$(echo "$times" | awk '{s+=$1} END{printf "%.6f", s/NR}')

        pi_vals=$(echo "$block" | grep "^Pi is" | awk '{print $NF}')
        mean_err=$(echo "$pi_vals" | awk -v pi=3.141592653589793 \
                      '{d=$1-pi; if(d<0)d=-d; s+=d} END{printf "%.10f", s/NR}')

        echo "$method,$variant,$threads,$steps,$mean_time,$mean_err" >> $out
    done
done

echo "Wrote $out"