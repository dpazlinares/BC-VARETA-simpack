# BC-VARETA
Includes the routines of the Brain Connectivity Variable Resolution Tomographic Analysis (BC-VARETA), a Simulation Package for MEEG  and other well established Methods for comparison purpose. BC-VARETA and the other Methods within the package extract the Source Activity and Connectivity given a single frequency component in the Fourier Transform Domain of an Individual MEEG Data.

***First download the complementary file 'data.rar' from the link below, extract the contained subfolder and copy it along with the routines into a common folder***

https://lstneuro-my.sharepoint.com/:u:/g/personal/deirel_paz_neuroinformatics-collaboratory_org/EVCqmmZ9d9dLrOTb37bpfrwB8yU3xAIk6AXJB7vzzkiz2g?e=Oa2Eqc

- Main (**execute this routine for demosntration**): generates simulation of 4 cortical connected points and reproduces the results of 
  BC-VARETA, sLORETA and LCMV.
- surfpatch_v1: creates parches around the cortical points for visualization of the connectivity in a reduced space  
- mkpinknoise: generates pink noise in source and sensor spaces
- xspectrum: computes the spectra of the simulated scalp activity 
- scalp_topography: crates plots with the topographic maps of activity and connectivity
- bcvareta: executes BC-VARETA method
- bcvareta_initial_values: computes 'bcvareta' initialization
- screening_ssbl: extracts the posibly active generators as part of 'bcvareta_initial_values', using the Elastic Net Structured Sparse
  Bayesian Learning
- screening: applies a smoothing to the outputs of 'screening_ssbl'
- mkfilt_eloreta: computes eLORETA method
- mkfilt_lcmv: computes LCMV method
