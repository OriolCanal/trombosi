# trombosi

## Pipeline description

To run this pipeline we need to have the plugins installed, and the cache file of the vep program also installed

### Command line arguments

Mandatory:

-d, --input_directory: input directory where fastq or fastq.gz files are stored

Optional:

-r, --reference: path where reference genome file is stored along with its bwa index files


### Algorithm steps:

* First of all the pipeline decompress the fastq.gz files that are in the input directory. Once the files are decompressed, it lists all the fastq files and match them according to the RB code (creating a list of lists containing paired end reads fastq files, and a list of single fastq files).
* The paired-end reads are analysed using fastq and give the output parameters to a quality dictionary.


```
sudo docker pull broadinstitute/picard
sudo docker pull ensemblorg/ensembl-vep
sudo docker pull broadinstitute/gatk
sudo docker pull quay.io/biocontainers/mosdepth
```

docker images:

|REPOSITORY |  TAG |  IMAGE ID | CREATED | SIZE|
| --- | --- | --- | --- | --- |
|broadinstitute/picard |latest | 7d402db81b4b | 3 weeks ago | 1.29GB |
|broadinstitute/picard  | <none> |  d6f85b73a21a | 2 months ago | 1.29GB |
|perl | latest | 0c14309df6d3 | 2 months ago | 894MB |
|ensemblorg/ensembl-vep | latest | 60d876219751 | 2 months ago | 704MB |
|broadinstitute/gatk | latest | d98b4a3aacfc | 3 months ago | 4.51GB |
|hello-world | latest | feb5d9fea6a5 | 16 months ago | 13.3kB |
|broadinstitute/gatk | 4.1.3.0 | 0ef3d00baa09 | 3 years ago | 3.72GB |
|quay.io/biocontainers/mosdepth | 0.2.4--he527e40_0 |912bd2aff907 | 4  years ago |  81.9MB|
