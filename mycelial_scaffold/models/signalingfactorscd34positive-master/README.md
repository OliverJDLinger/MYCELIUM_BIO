# **NicheHotSpotter**

## Scripts to analyse CD34 positive and negative muscle stem cells to identify niche induced signaling factors that maintain these stem cell phenotypes.

### Instruction to use in a local computer


1.  Download the files in a local directory.
2.  From an RStudio terminal, install required packages by typing:
    ```R
    install.packages( c("topGO", "plyr", "igraph", "Matrix", "reshape2", "RSpectra", "data.table", "dplyr", "snow", "DT",  "openxlsx", "markdown", "visNetwork") )
    ```
3. Run the the markdown file musclestemcells.Rmd to produce the results
4. Network visualizations are made using [Cytoscape](https://cytoscape.org/) from the generated network output files 
 
#### Authors
NicheHotSpotter was developed in the [Computational Biology Group](https://wwwfr.uni.lu/lcsb/research/computational_biology) by
- [Srikanth Ravichandran](https://wwwen.uni.lu/lcsb/people/srikanth_ravichandran)
- [Antonio del Sol](https://wwwfr.uni.lu/lcsb/people/antonio_del_sol_mesa)