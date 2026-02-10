
## 
This code implements the Fast Multipole Method (FMM) for evaluating 
secular equations arising in eigenvalue problems.

## Requirements
- MATLAB R2022b or later

## Usage
Run the main MATLAB script:

FMM_run.m file
Load "FMM_test_data" matlab data file from the test_data folder in to FMM_run.m and run the script.

## Note 1
Vary number of expansion coefficients pp in FMM_run.m file for obtaining different results.
## Note 2
Multipole-to-local translations are avoided in derivative evaluation
to reduce round-off errors.

## Assumptions and Limitations
- The current implementation assumes a uniform distribution of
  sources, leading to the use of a uniform binary tree.
- This assumption is appropriate for the test problems considered in
  the manuscript.
- Performance and accuracy for highly nonuniform point distributions
  are not addressed here and are left for future work.

License
For academic use only.








