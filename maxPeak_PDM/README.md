# maxPeak PDM

This folder contains a fully functional example of the maxPeak peak detection method. 

Input to maxPeak comprises processed .wig files. To get these, the number of reads on both + and – DNA strands was summed up genome wide for each nucleotide position. At this stage, biological duplicates were averaged.

Running the maxPeak method is performed by simply running the maxPeak.r script. This will generate .bed files for each input .wig file. The .bed files will list the gene name, chromosome, a 10 bp window around the peak bp within the 1000 bp window upstream of the ORF of the gene and a signal to noise ratio (SNR). The SNR is calculated as the number of reads in the peak location normalized by the 65th percentile of the peaks across all genes.

To apply the maxPeak method to your own data simply change the .wig file names defined within the R script.

See the publication for further details on the method and its results.