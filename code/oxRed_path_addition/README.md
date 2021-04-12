#Description
Scripts for incorporation of the oxido-reductive pathway for lactose utilization into `CintGEM`. Additionally, gene IDs are mapped to those used in the available RNAseq datasets for **C. intermedia** to enable straightforward data incorporation into our model.

#Prerequisites
Run the scripts in the folder `code/RNAseq_analysis` prior to any script in this folder.

#Instructions
Run the scripts in the following order.
 
1. `oxidopathwayaddition.m`: Integration of the pathway into the model and evaluation of flux capabilities. 
2. `mapGeneIDs.m`: Mapping of gene IDs between model and RNAseq data (using results from orthofinder runs).
3. `lxr4_exploration.m`: explore presence of lxr4 in the model and data. Result: lxr4 was already present in our model but assigned to other reaction(s).

#Results
The resulting functional model can be found in `models/candida_intermedia/cintGEM_oxido.mat`