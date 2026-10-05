#!/bin/bash

# Run max_jobs copies of skeletor at the same time in order to process many runs.
# $1 = list of run numbers

#Set these variables:
# don't need a final slash in your paths below:
# where to find your stg1 files:
INDIR=/data
#where to write your skeletor output root files
OUTDIR=/data/skeletor
#Choices: GV, KS, or blank
TOFtype=GV
# Max number of running copies at the same time
# NB: Full isochrones v9 take 6GB/run, so set the number of runs to execute 
# simultaneously accordingly so that you don't run out of memory.
max_jobs=35

if [ ! -f "$1" ]; then
    echo "Error: File $1 not found."
    exit 1
fi

# Find length of list of runs:
run_count=$(wc -l < $1)
echo "Processing $run_count runs..."

# Loop through input files from 1 to 40 (or change pattern as needed)
# Loop through each line of the file
while IFS= read -r line; 
  do
      echo "Processing: $line"
      fname=run_$line.stg1.root
      infile=$INDIR/$fname
      outfile=$OUTDIR/skeletor.$line.root
      echo Processing $infile into $outfile...
      ./bin/skeletor.exe --ckf -c -o $outfile $infile &> $OUTDIR/$line.out &

      # Check how many jobs are running right now
      while [ $(jobs -r | wc -l) -ge $max_jobs ]; do
          # Wait a short moment before checking again
          sleep 60
      done

  done < $1

# Wait for all background jobs to finish before the script exits
wait

echo "All tasks are done!"
