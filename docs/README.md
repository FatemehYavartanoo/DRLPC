## Documentation Folder

This folder contains additional resources and documentation for the DRLPC project. Below is a description of the files:

- **DRLPC_UserGuide.pdf**: A comprehensive user guide that explains how to use the DRLPC function, including:
   - Preparing input data.
   - Running the DRLPC algorithm.
   - Interpreting the outputs.

- **DRLPC_algorithm.png**: A figure illustrating the DRLPC algorithm.

- **DRLPC_final_results.rds**: The results are processed after running the DRLPC algorithm. This file contains the reduced dimensions and results for the analyzed genes.

To view the detailed instructions, refer to [DRLPC_UserGuide.pdf](DRLPC_UserGuide.pdf).

To view the DRLPC algorithm figure, click on [DRLPC-algorithm.png](DRLPC-algorithm.png).


### DRLPC Final Results  
The complete results file, **`DRLPC_final_results.rds`**, can be found in the `docs` folder.  
To load the file into R for further analysis, use the following command:  
```R
DRLPC_results <- readRDS("docs/DRLPC_final_results.rds")


