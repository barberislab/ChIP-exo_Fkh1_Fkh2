[![Binder](https://mybinder.org/badge.svg)](https://mybinder.org/v2/gh/barberislab/ChIP-exo_Fkh1_Fkh2/master)

# ChIP-exo analysis highlights Fkh1 and Fkh2 transcription factors as hubs that integrate multi-scale networks in budding yeast
This repository contains all experimental and computational data and the subsequent outputs for the publication titled "ChIP-exo analysis highlights Fkh1 and Fkh2 transcription factors as hubs that integrate multi-scale networks in budding yeast" by Thierry D.G.A. Mondeel, Petter Holland, Jens Nielsen and Matteo Barberis. 

The article is currently in revision. A link to the full article will appear here once it is accepted and publicly available.

# Repository organization

The "Data" folder contains the raw experimental data, the processed experimental data (GEM, MACE and maxPeak), binding motif analysis (MEME, DREME and CentriMo), a version of the SQL database from [GEMMER](http://gemmer.barberislab.com) and data from the literature we have used in our analysis.

The "Code" folder contains all Python code in Jupyter notebooks and R scripts to generate to figures and tables presented in the manuscript starting from the "Data" folder. 

The output of the scripts and notebooks in the "Code" folder is saved in the "Figures" and "Tables" folders. Only these outputs were used for the figures and tables presented in the manuscript and the Supplementary material. 

See the READMEs in the individual folders for further details.

# maxPeak peak detection methodology 

The "maxPeak_PDM" folder contains a fully functional R script illustrating the use of the maxPeak peak detection method which we introduced in the publication. The resulting .bed files are the ones we utilised for all subsequent analyses in the text. See the README in the "maxPeak_PDM" folder and the main text for further details.

# Reproduce our analysis offline and in the cloud

Most of the analysis covered in the publication is reproducible using the jupyter notebooks of Python code in this repository through Binder or by downloading the files and using your local Python installation. To start, click the "Binder" button above. From there navigate to the "Code" folder and you will be able to execute the various Jupyter notebooks. 

To reproduce the rest of the analysis simply use the .py and .r scripts in the "Code" and "Data" folders on your personal computer. 
