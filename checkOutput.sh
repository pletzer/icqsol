#!/usr/bin/env 

good="yes"
for f in `ls *.vtk`; do 
    if [ -f "test_results/$f" ]; then
        diff $f test_results/$f > /dev/null
        if [ $? != 0 ]; then
            echo "ERROR: $f and test_results/$f are different"
            good="no"
        fi
    fi
done

for f in `ls *.ply`; do
    if [ -f "test_results/$f" ]; then
        diff $f test_results/$f > /dev/null
        if [ $? != 0 ]; then
            echo "ERROR: $f and test_results/$f are different"
            good="no"
        fi
    fi
done

if [ $good == "no" ]; then
    exit 1
fi