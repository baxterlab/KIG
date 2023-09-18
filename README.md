# Known Ionomic Gene List pipeline (v1.0) 

##### Version 2: the updated list with OrthoFinder orthologs can be found [here](https://github.com/danforthcenter/KIG_v2/)

RScript that takes the primary genes from the Known Ionomic Genes list and finds orthologs from other species and appends them to the list as inferred genes. Each species will have its own csv output file that contains all of its primary genes and genes that have been inferred as orthologs to other species' primary genes, along with each gene's corresponding orthologs in every other species. Current primary species, dependent upon submissions to the KnownIonomeList, include Arabidopsis thaliana, Oryza sativa, Zea mays, and Medicago truncatula. Current inferred species include all the organisms listed in the primary list, as well as Sorghum bicolor, Glycine max, Setaria viridis, and Setaria italica.

R version 3.5.0 (2018-04-23)

[Whitt et al. 2020](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7576880/)

#### R Packages: The KIG List script runs with these additional packages:
    1.1   plyr and dplyr (Wickham 2011 & Wickham et al. 2018)
    1.2   data.table (Dowle & Srinivasan 2018)
    1.3   doParallel (Microsoft & Weston 2017)
    1.4   readr (Wickham et al. 2017)

#### Input list
The pipeline takes in the known ionome list of submissions from collaborators as csv input to start the pipeline. The file should have one gene per line, with no duplicates, and include information about the gene ID, the species, elements the gene was associated with, the tissue type of the analysis and a citation of the study characterizing this gene/element interaction. Gene name, closest ortholog species and additional comments are optional, but if included will be copied along with the gene to its respective species table in the output. 

#### Function
The beginning of the pipeline reads in the primary list, and sorts it by species and gene ID. White spaces are removed from the elements file to prepare for string manipulations later in the pipeline. The geneFilePrep function reads in Phytozome InParanoid files to determine orthologs for each of the genes in our primary list of genes. After the an ortholog search has been executed for each gene in the list, the orthologs of all the genes are sorted into data frames by species. The function noDups is run to remove multiple entries of an ortholog. This happens when multiple primary genes infer the same gene as an ortholog. noDups removes the multiple entries and creates a string of all the genes that inferred this ortholog and writes it to the "Inferred from" column of the final output table. noDups also condenses all of the unique elements from the multiple primary genes into a string that is written to the column "Elements" of the final output table. geneFilePrep is run again to find the orthologs of these inferred genes. A list of all the orthologs for a gene is condensed into lists by species and appended to their species column in the list. For example, if you wanted to see all the G. max orthologs for AT1G01580, you would go to the A.thaliana_knownIonomicsGenesWOrthologs table, look at the row for AT1G01580 and the column "G. max orthologs". All of the species in the pipeline will have their primary and inferred genes combined into the same data table and saved to a csv with the suffix "_knownIonomicsGenesWOrthologs".

#### Output files
Each output file (species.name_knownIonomicsGenesWOrthologs.csv), will retain the same format as the primary list, except it will only include genes for its respective species, and it will have additional columns from the pipeline. After the original columns you will find the generated columns "Primary/Inferred", "Inferred from", and ortholog columns for each species.

#### Downloading new phytozome tables
For the purposes of finding orthologs and building a known gene list, only the default columns are needed from the phytozome ortholog comparison tool. If there are other things you want included in the list output (ie. GO annotations, protein similarity scores, etc.) they can be added in and merged on phytozome using the edit-query tool - though it can be a bit tricky if the servers are busy or faulty. If you're having trouble using their web browser to download comparison files, try using their biomart interface (https://phytozome.jgi.doe.gov/biomart/martview). When downloading files to add species, make sure to get files comparing each of the current species (or those being analyzed in the input list) and the new species (ex: A.thaliana_new.species.csv, G.max_new.species.csv, etc.).

#### References

Dowle, M. & Srinivasan, A. 2018. data.table: Extension of `data.frame`. R package version 1.11.4. https://CRAN.R-project.org/package=data.table

Microsoft Corporation & Weston, S. 2017. doParallel: Foreach Parallel Adaptor for the 'parallel' Package. Rpackage version 1.0.11. https://CRAN.R-project.org/package=doParallel

Wickham, H. 2011. The Split-Apply-Combine Strategy for Data Analysis. Journal of Statistical Software, 40(1), 1-29. URL http://www.jstatsoft.org/v40/i01/.

Wickham, H., François, R., Henry, L. and Müller, K. (2018). dplyr: A Grammar of Data Manipulation. R package version 0.7.5. https://CRAN.R-project.org/package=dplyr

Wickham, H., Hester, J. & Francois, R. 2017. readr: Read Rectangular Text Data. R package version 1.1.1. https://CRAN.R-project.org/package=readr
