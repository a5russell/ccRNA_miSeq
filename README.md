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

  

