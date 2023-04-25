## Containers for BENGAL nextflow pipeline

Yuyao Song <ysong@ebi.ac.uk>

We provide a docker container for running SCCAF assessment. 

For cluster execuition, we convert the docker container into a singularity container. 

Pull SCCAF containers as follows:

`singularity pull sccaf.sif docker://yysong123/intgpy:sccaf`

This container is built for linux/amd64 machines.
