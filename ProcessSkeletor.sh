#!/bin/bash

# pass the txt file of run numbers for processing by skeletor as $1

#Set these variables:
# don't need a final slash in your paths below:
# where to find your stg1 files:
INDIR=/data
#where to write your skeletor output root files
OUTDIR=/data/skeletor
#Choices: GV, KS, or blank
#TOFtype=GV
TOFtype=

# Check if the file actually exists before reading
if [ ! -f "$1" ]; then
    echo "Error: File $1 not found."
    exit 1
fi

# Loop through each line of the file
while IFS= read -r line; do
    # Process the line (e.g., print it out)
    echo "Processing: $line"
    fname=run_$line.stg1.root
    infile=$INDIR/$fname
    outfile=$OUTDIR/skeletor.$line.root
    echo Processing $infile into $outfile...
#Step 1: Make TOF tree. Input: stg1 file, output skeletor.runnum.root
    ./bin/skeletor.exe --ckf -c -o $outfile $infile &> $OUTDIR/$line.out
#Step 2: Make histograms and profiles. Output: skeletor.runnum.AnalTOFSkel.root
    echo Processing $outfile with AnalTOFSkel.C...
    root -l -q -b 'AnalTOFSkel.C("'"${outfile}"'")'
#Step 3: Fit profiles, make pdfs of plots. Input skeletor.runnum.AnalTOFSkel.root, output: MakeSkeletorPlotsBatch.runnum.pdf and several txt files with fit params
    analfile=$OUTDIR/skeletor.$line.AnalTOFSkel.root
    echo Processing $analfile with MakeSkeletorPlotsBatch.C...
    root -l -q -b 'MakeSkeletorPlotsBatch.C("'"${analfile}"'","'"${line}"'","'"${TOFtype}"'")'
done < $1
