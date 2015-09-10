#!/usr/bin/env 

for f in `ls *.vtk`; do 
    if [ -f "test_results/$f" ]; then
        diff $f test_results/$f
        if [ $? != 0 ]; then
            echo "ERROR: $f and test_results/$f are different"
            exit 1
        fi
    fi
done

for f in `ls *.ply`; do
    if [ -f "test_results/$f" ]; then
        diff $f test_results/$f
        if [ $? != 0 ]; then
            echo "ERROR: $f and test_results/$f are different"
            exit 1
        fi
    fi
done

