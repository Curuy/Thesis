#!/bin/bash

echo ""
echo "========================================="
echo "          COMPILING TESTS                "
echo "========================================="


make clean
make tests

TEST_DIR="test"

echo ""
echo "========================================="
echo "          RUNNING TESTS                  "
echo "========================================="
SECONDS=0
for file in "$TEST_DIR"/test_*; do
    if [ -f "$file" ] && [ -x "$file" ]; then
        echo "-> Running $(basename "$file")"
        "$file"
        echo "-----------------------------------------"
    fi
done

duration=$SECONDS
mins=$((duration / 60))
secs=$((duration % 60))


echo "Total Execution Time: ${mins}m ${secs}s"
echo "========================================="


echo ""
echo "========================================="
echo "          CLEANING                       "
echo "========================================="

make clean

