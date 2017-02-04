The explanation of the code.


This code is based on  manuscript "Inferring interactions in complex microbial communities from
nucleotide sequence data and environmental parameters" that was submitted to PLOS one.


This folder contains two R-script files, four R-data files.  The "Abundscale.Rdata" and "Parascale.Rdata" are the original data sets used in our manuscript. They will take a long time to compute. The "Abundscale_short.Rdata" and "Parascale_short.Rdata" are two smaller data sets which take much less computation time for reader to test. 

The interaction analysis workflow is explained in file “test_git.R”. The functions which are used in “test_git.R” are defined in the file “InteractionAnalysis_git.R”.


The code “test_git.R” contains three parts, the first part is used for our original data set and takes a long time for computation; the second part is used for the smaller dataset and repeats the same calculation as the first part; the third part contains the additional tools which are not presented in our manuscript, but, might be useful in the future analysis.  


The basic roadmap of the workflow can be summarized as follows:

(1)    Load the data frame  of species abundance and the environmental parameter. Scale the data if it is necessary. 

(2)    Filter the abundance and environmental parameter data frame out the “NA” and “0” values. Use the function: FilterNA(), and Filterrepeat()

(3)    Do the test of singularity on both of the abundance and environmental parameter data frame.

(4)    Define the function parameters which are needed for the numerical calculation 

(5)    Calculate the rate of change of species abundance with respect to environmental parameter (p^k_i in our paper). use the function:  rate_change(), or rate_change02(),

(6)    Calculate the interaction values for each species in the sense of environmental parameter(\beta_{ij\alpha}^{k} in our paper). use the function: interinf(), or interinf02()

(7)    use the summarizing strategy to calculate the global interaction matrix from (\beta_{ij\alpha}^{k}  . Use the function: interMatrix02()

(8)    Extract the six examples for robust test:  strong positive, median positive, low positive, strong negative, median negative, low negative. Use the function: six_test_robust()

(9)     Do the robust test 


The detailed description of these function can be found in the file “InteractionAnalysis_git.R”.


ps: additional tools (it is not presented in the manuscript):   Figure out the dominant parameters. Use the function: main_para()












