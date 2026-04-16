[← Back to main page](index.md)

# FAQ

## Pipeline initialisation

### Workflow defines configfile config_chip.yml but it is not present or accessible.

This error occurs when the pipeline is run without a config file present in the working directory. Ensure that seqnado-config has been run before starting the pipeline and that you are in the new directory created by seqnado-config.

Follow the [Configuration Guide](configuration.md) instructions to create a config file.


## Singularity configuration

### Workflow Error

Failed to pull singularity image from library://asmith151/seqnado/seqnado_pipeline:latest:  
FATAL: Unable to get library client configuration:  
remote has no library client (see https://apptainer.org/docs/user/latest/endpoint.html#no-default-remote)

Fix:

re-run seqnado init: See the [Initialisation Guide](initialisation.md)

or

```bash
apptainer remote add --no-login SylabsCloud cloud.sylabs.io  
apptainer remote use SylabsCloud  
```

## Optional configuration

### Can I merge multiple samples into a single sample?

Yes, you can merge multiple samples into a single sample to generate merged bigWig files and consensus peaks. To do this, you need to create a design file that specifies the samples to be merged. The design file should have a column named "merge" that specifies the samples to be merged e.g.:


| sample | r1 | r2 | merge |
|--------|----|----|--------|
| sample-control-rep1 | /path/to/sample-control-rep1_R1.fastq.gz | /path/to/sample-control-rep1_R2.fastq.gz | control | 
| sample-control-rep2 | /path/to/sample-control-rep2_R1.fastq.gz | /path/to/sample-control-rep2_R2.fastq.gz | control | 
| sample-control-rep3 | /path/to/sample-control-rep3_R1.fastq.gz | /path/to/sample-control-rep3_R2.fastq.gz | control | 
| sample-treated-rep1 | /path/to/sample-treated-rep1_R1.fastq.gz | /path/to/sample-treated-rep1_R2.fastq.gz | treated |
| sample-treated-rep2 | /path/to/sample-treated-rep2_R1.fastq.gz | /path/to/sample-treated-rep2_R2.fastq.gz | treated |
| sample-treated-rep3 | /path/to/sample-treated-rep3_R1.fastq.gz | /path/to/sample-treated-rep3_R2.fastq.gz | treated |
