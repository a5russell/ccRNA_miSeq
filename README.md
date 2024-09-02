# ccRNA_miSeq

Repository for analyzing a multiplexed PCR with template-switching, followed by Gibson circularization and miSeq, to identify ccRNA.

Original effort included SNP-containing templates, but their abundance was too low to reliably use to correct for incorrect template-switching.

We instead assess chimera rate, whcih was ~<0.05%, and so unlikely to explain our observations.

## Directories

The following directories, and their purpose, exist within this repository.

- <b>Database</b>       Repository for influenza genomic sequences. A/WSN/1933 BLAST database and STAR indicies provided. also contains TSO adapter sequence and matches for     5' and 3' ends of vRNA segments.
- <b>Results</b>        Final datafiles for analyses after processing.
- <b>Scripts</b>        Short scripts written for this analysis. Seperated from jupyter notebooks for readability and portability.
- <b>Sequencing</b>     Folder that would contain NGS samples.

  ## Dependencies

  - <b>Python</b>      run with version 3.7. Available from https://www.python.org/downloads/release/python-370/
- <b>Trimmomatic</b> run with version 0.39. Available from http://www.usadellab.org/cms/?page=trimmomatic
- <b>STAR</b>        run with version 2.7.1.a. Available from https://github.com/alexdobin/STAR
- <b>Samtools</b>    run with version 1.9. Available from http://www.htslib.org/
- <b>FastQC</b>      run with version 0.11.8. Available from https://www.bioinformatics.babraham.ac.uk/projects/fastqc/

The following python packages and versions were used. All were installed using mamba. (https://mamba.readthedocs.io/en/latest/installation/mamba-installation.html)
- <b>numpy</b>       run with version 1.26.4. (https://numpy.org/)
- <b>matplotlib</b>  run with version 3.8.4. (https://matplotlib.org/)
- <b>seaborn</b>     run with version 0.13.2. (https://seaborn.pydata.org/)
- <b>pandas</b>      run with version 2.2.2.(https://pandas.pydata.org/)
- <b>scipy</b>       run with version 1.13.0. (https://www.scipy.org/)
- <b>statsmodels</b> run with version 0.14.1. (https://www.statsmodels.org/stable/index.html)

