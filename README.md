The explanation of the code.

This folder contains two R-script files, two R-data files which are the species abundance information, and a txt files which is the environmental parameter table. 
The interaction analysis workflow is explained in file “test_git.R”. The functions which are used in “test_git.R” are defined in the file “InteractionAnalysis_git.R”.
The basic roadmap of the workflow can be summarized as follows:
(1)	Load the dataframe  of species abundance and environmental parameter. Scale the data if it is necessary. 
(2)	Filter the abundance and environmental parameter dataframe out the “NA” and “0” values. Use the function: FilterNA(), and Filterrepeat()
(3)	Do the test of singularity on both of the abundance and environmental parameter dataframe.
(4)	Define the function parameters which are needed in the numerical calculation 
(5)	Calculate the rate of change of species abundance with respect to environmental parameter (p^k_i in our paper). use the function:  rate_change(), or rate_change02(),
(6)	Calculate the interaction values for each species in the sense of environmental parameter(\beta_{ij\alpha}^{k} in our paper). use the function: interinf(), or interinf02()
(7)	use the summarizing strategy to calculate the global interaction matrix from (\beta_{ij\alpha}^{k}  . Use the function: interMatrix02()
(8)	Figure out the dominant parameters. Use the function: main_para()
(9)	Extract the six examples for robust test:  strong positive, median positive, low positive, strong negative, median negative, low negative. Use the function: six_test_robust()
(10)	 Do the robust test 
The detailed description of these function can be found in the file “InteractionAnalysis_git.R”.








