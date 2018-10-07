# BC-VARETA-simpack
Includes the routines of the Brain Connectivity Variable Resolution Tomographic Analysis (BC-VARETA), a Simulation Package for MEEG  and other well established Methods for comparison purpose. BC-VARETA and the other Methods within the package extract the Source Activity and Connectivity given a single frequency component in the Fourier Transform Domain of an Individual MEEG Data. See the pdf file "Brief of Theory and Results" for an insight to this methodology.
References:
Paz-Linares, D., Gonzalez-Moreira, E., Martinez-Montes, E. and Valdes-Sosa, P.A., 2018. Note on the Estimation of Embedded Hermitian Gaussian Graphical Models for MEEG Source Activity and Connectivity Analysis in the Frequency Domain. Part I: Single Frequency Component and Subject. arXiv preprint arXiv:1810.01174.
https://arxiv.org/abs/1810.01174
Paz-Linares, D., Gonzalez-Moreira, E., Martinez-Montes, E., Valdes-Hernandez, P.A., Bosch-Bayard, J., Bringas-Vega, M.L. and Valdes-Sosa, P.A., 2018. Caulking the Leakage Effect in MEEG Source Connectivity Analysis. arXiv preprint arXiv:1810.00786.
https://arxiv.org/abs/1810.00786

***First download the complementary file 'data.rar' from the link below, extract the contained subfolder and copy it along with the routines into a common folder***

https://lstneuro-my.sharepoint.com/:u:/g/personal/deirel_paz_neuroinformatics-collaboratory_org/EVCqmmZ9d9dLrOTb37bpfrwB8yU3xAIk6AXJB7vzzkiz2g?e=Oa2Eqc

BC-VARETA-master:
- Main (**execute this routine for demosntration**): generates simulation of 4 cortical connected points and reproduces the results of 
  BC-VARETA, sLORETA and LCMV.
- surfpatch_v1: creates parches around the cortical points for visualization of the connectivity in a reduced space  
- bcvareta: executes BC-VARETA method
- bcvareta_initial_values: computes 'bcvareta' initialization
- screening_ssbl: extracts the posibly active generators as part of 'bcvareta_initial_values', using the Elastic Net Structured Sparse
  Bayesian Learning
- trascendent_term: nonlinear function for regularization parameters estimation within the function 'screening_ssbl'     
- screening: applies a smoothing to the outputs of 'screening_ssbl'
- mkfilt_eloreta: computes eLORETA method
- mkfilt_lcmv: computes LCMV method
- data: folder containing leadfield, surfaces, colormaps, etc 
