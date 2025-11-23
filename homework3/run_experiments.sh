#!/bin/bash

MATRIX="cant/cant.mtx"
VECTOR="cant/b.mtx"
OUT="myans.mtx"

THREADS=("1" "2" "4" "8" "16" "32" "64" "128")

CSV="results.csv"
echo "threads,load,convert,lock_init,coo_par,csr_par,coo_ser,csr_ser,store" > $CSV

for t in "${THREADS[@]}"; do
    export OMP_NUM_THREADS=$t
    echo "Running with $t threads..."

    OUTPUT=$(./spmv $MATRIX $VECTOR $OUT)

    LOAD=$(echo "$OUTPUT"      | grep -E "^Load[[:space:]]"       | awk '{print $2}')
    CONVERT=$(echo "$OUTPUT"   | grep -E "^Convert[[:space:]]"    | awk '{print $2}')
    LOCK=$(echo "$OUTPUT"      | grep -E "^Lock Init[[:space:]]"  | awk '{print $3}')

    COO_PAR=$(echo "$OUTPUT" | grep "^COO SpMV" | grep -v "Serial" | awk '{print $NF}')
    CSR_PAR=$(echo "$OUTPUT" | grep "^CSR SpMV" | grep -v "Serial" | awk '{print $NF}')
    COO_SER=$(echo "$OUTPUT" | grep "^COO SpMV Serial" | awk '{print $NF}')
    CSR_SER=$(echo "$OUTPUT" | grep "^CSR SpMV Serial" | awk '{print $NF}')

    STORE=$(echo "$OUTPUT"     | grep -E "^Store[[:space:]]"      | awk '{print $2}')

    echo "$t,$LOAD,$CONVERT,$LOCK,$COO_PAR,$CSR_PAR,$COO_SER,$CSR_SER,$STORE" >> $CSV
done

echo "finished!"