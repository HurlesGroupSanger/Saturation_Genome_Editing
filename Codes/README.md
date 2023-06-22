Codes used in DDX3X SGE analysis.

This code is for reproducing the analysis in "Saturation genome editing of DDX3X clarifies pathogenicity of germline and
somatic variation". 

Preprint:medRxiv (2022), [doi:10.1101/2022.06.10.22276179](https://www.medrxiv.org/content/10.1101/2022.06.10.22276179v1)

We first generate the oligos library with the code in [SGE oligo generation](https://github.com/HurlesGroupSanger/Saturation_Genome_Editing/tree/main/Codes/SGE%20oligo%20generation).

After we have the sequencing result, we perform the read counting with the script in [SGE Fastq-to-Count](https://github.com/HurlesGroupSanger/Saturation_Genome_Editing/tree/main/Codes/SGE%20Fastq-to-Count).

Then, we calculate the LFC and LFC-trend with the R code in [SGE LFC and LFC-trend calculation](https://github.com/HurlesGroupSanger/Saturation_Genome_Editing/tree/main/Codes/SGE%20LFC%20and%20LFC-trend%20calculation/LFC-trend).

For reproducing the figures, the Python codes are in the [jupyternotebooks for main figures in DDX3X paper](https://github.com/HurlesGroupSanger/Saturation_Genome_Editing/tree/main/Codes/jupyternotebooks%20for%20main%20figures%20in%20DDX3X%20paper).

For variant effect visualization (Figure S3 and Figure S6), the R codes are in the [sge variant effect visualization](https://github.com/HurlesGroupSanger/Saturation_Genome_Editing/tree/main/Codes/sge%20variant%20effect%20visualization).

