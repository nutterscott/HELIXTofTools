These tools run in the HELIX environment.

# Files
## Analysis scripts for root:
- **AnalTOFSkel.C**: Takes skeletor output root file and creates histograms
- **PSUThicknessTools.C**: Routines for returning thickness of scintillator as function of position.
- **PSUTofTools12.C**: Various routines to assist in analysis; reads in v12 of tofcal.root.
- **MakeSkeletorPlotsBatch.C**: Takes AnalTOFSkel.C output and makes pdf of histograms. Also does some fits (DCTY vs TimeDiff and DCTY vs ChargeDivFrac) and writes fit results to text files.
- **WriteRunNumbers.C**: Use with files in Analysis/Runlists/ to generate your own runlists.

## Lists of runs
- **AllRuns.txt**: Duh. All runs, even not good ones.
- **2024_Kiruna_FLoat_Boff.txt**: Official list
- **2024_Kiruna_FLoat_Bon.txt**: Official list


## Bash scripts to run analysis. Takes list of runs as input.
- **ProcessSkeletor.sh**: Start to finish processing of a run from run list, including skeletor, AnalTOFSkel.C, and MakeSkeletorPlotsBatch.C
- **MultiSkeletor.ss**: Run multiple copies of skeletor simultaneously from input run list.

Owner:
Scott Nutter
nutters@nku.edu
