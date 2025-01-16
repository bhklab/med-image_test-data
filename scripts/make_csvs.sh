#!/bin/bash

# Create the test_data directory if it doesn't exist
mkdir -p test_data

# Generate 10 CSV files
for i in {1..10}; do
  # Create the CSV file with 10 rows and 10 columns
  echo "Col_1,Col_2,Col_3,Col_4,Col_5,Col_6,Col_7,Col_8,Col_9,Col_10" > test_data/test_file_$i.csv
  for j in {1..10}; do
    # Generate a row with random numbers between 1 and 100
    echo "$(shuf -i 1-100 -n 10 | tr '\n' ',' | sed 's/,$//')" >> test_data/test_file_$i.csv
  done
done

echo "10 CSV files have been created in the 'test_data' directory."