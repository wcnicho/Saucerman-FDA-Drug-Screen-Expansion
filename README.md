# Saucerman-FDA-Drug-Screen-Expansion
- Expands upon the Drug Screen created by [Saucerman et al.](https://doi.org/10.1002/psp4.12599) to include more drug targets
- Uses and updated Sex-differentiated model created by Kelsey Watts et al.
 - Previously expanded by [Jesse Rogers et al.](https://doi.org/10.7554/eLife.62856)
 
 # Files
- tblParse.m
  - combines together drug targets found in the DrugBank database
- Drug Actions.xlsx
  - list of drug name and their action (Competitive/Noncompetitive) as found on pubmed
- SimulatedDrugsUpdate.xlsx
  - final updated drug target sheet to be used in the Saucerman drug code
- drugbank_resultsComb.xlsx
 - combined spreadsheet of drug targets from 'drugbank_results (1).csv', 'drugbank_results (2).csv', and 'drugbank_results.csv'
- makeSimulatedSheet.m
  - code creating the formatted 'SimulatedDrugsUpdate.xlsx' spreadsheet
- targets.txt
  - list of the expanded model's gene targets
