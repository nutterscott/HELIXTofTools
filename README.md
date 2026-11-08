These tools run in the HELIX environment.

# Files
## Analysis scripts for root:
- **AnalTOFSkel.C**: Takes skeletor output root file and creates histograms
- **PSUThicknessTools.C**: Routines for returning thickness of scintillator as function of position.
- **PSUTofTools12.C**: Various routines to assist in analysis; reads in v13 of tofcal.root.
- **MakeSkeletorPlotsBatch.C**: Takes AnalTOFSkel.C output and makes pdf of histograms. Also does some fits (DCTY vs TimeDiff and DCTY vs ChargeDivFrac) and writes fit results to text files. NB: Full isochrones v9 take 6GB/run, so set the number of runs to execute simultaneously accordingly so that you don't run out of memory.
- **WriteRunNumbers.C**: Use with files in Analysis/Runlists/ to generate your own runlists.
- Plus more C scripts that I am too lazy to document here, hoping that you can figure out what they do based on their names.

## Lists of runs
- **AllRuns.txt**: Duh. All runs, even not good ones.
- **2024_Kiruna_Float_Boff.txt**: Official list
- **2024_Kiruna_Float_Bon.txt**: Official list
- **AllGoodRuns.txt**: Combination of 2024 Bon and Boff above


## Bash scripts to run analysis. Takes list of runs as input.
- **ProcessSkeletor.sh**: Start to finish processing of a run from run list, including skeletor, AnalTOFSkel.C, and MakeSkeletorPlotsBatch.C
- **MultiSkeletor.ss**: Run multiple copies of skeletor simultaneously from input run list.

# Processing:
For each run, do this in order:
- skeletor
- AnalTOFSkel.C
- MakeSkeletorPlotsBatch.C or other C scripts

Check codes for input parameters, or example in ProcessSkeletor.sh.

# Owner:
Scott Nutter
nutters@nku.edu
