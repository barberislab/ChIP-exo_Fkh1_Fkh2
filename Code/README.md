# Python notebooks

The 4 Jupyter notebooks in this folder contain all Python ode to perform all analyses discussed in the text with the exception of the binding motif analysis. 


# Targets mapped on KEGG pathways through R Pathview

We used the R library Pathview (http://pathview.r-forge.r-project.org/) to map our annotated ChIP-seq dataset onto KEGG pathway maps. 

## Original Pathview reference
Luo W, Brouwer C. Pathview: an R/Biocondutor package for pathway-based data integration and visualization. Bioinformatics, 2013, 29(14):1830-1831, doi: 10.1093/bioinformatics/btt285

## Organization
2 folders:
Mondeel2017 considers only those targets suggested by our study. 
Consensus considers a set of targets derived from our study and 3 previous studies (Venters et al, Hu et al. MacIsaac et al. and Ostrow et al.)

Each folder consist of subfolders for various KEGG pathway maps. 

## Joint images for Fkh1 and Fkh2
Each enzyme or gene in the KEGG maps will have two associated states by virtue of dividing the node in half: the left half reflects results for Fkh1, the right half reflects results for Fkh2. 

## Colors
Mondeel2018: targets identified by our study are coloured in a green hue.

## Consensus:
Green - targets verified in this study and at least 1 ChIP-chip study
Red - targets verified in this study but none of the ChIP-chip studies
Yellow - Targets unverified in this study but found in at least 1 ChIP-chip study

## Note of warning:
Some KEGG pathway maps contain nodes related to several genes/protein/isoenzymes. When Pathview projects data on such maps it maps all isoenzymes to that particular node and by default sums their values. We opted to change this so that the colour shown is based on the maximum value over all genes mapped to that particular node. For instance, in the cell cycle map you will notice in the uncoloured map that there is a node labeled: Clb1/2. In the images including our data this will be labeled ‘Club’, but in fact the data will relate to both Clb1 and Clb2. If, hypothetically, Clb1 would have a signal-to-noise ratio of 1 and Clb2 a value of 3. The image will colour the node according to the maximum, i.e. 3. 
