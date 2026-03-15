for n in 10 20 50 100 200 500 1000 2000 5000 10000; do
    echo "Running input_cross_${n}nm.toml"
    OMP_NUM_THREADS=36 ./build/ntmc "calculation-cross/input_cross_${n}nm.toml"
done