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
* The paired-end reads are analysed using fastqc and give the output parameters to a quality dictionary.
* Then, fastq files are trimmed using [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic), we also drop the reads that are not paired.
* The paired reads, are then aligned to the reference grch37.13 genome using [bwa mem algorithm](https://github.com/lh3/bwa). Obtaining a sam file.
* The sam file is then converted into a sorted bam file along with its index file. It is also calculate the % enrichment parameter, (which is the number of regions aligned into the bed file divided by the number of primary reads aligned), the call rate (percentage of exons that are covered at a 1X, 10X, 20X, 30X, 100X, 200X depth coverage) and the exon lost (total number of exons that have not been covered at a given coverage) along with many other BAM stats. These parameters are stored in the quality dictionary.
* The bam files are processed using gatk [AddOrReplaceReadGroups](https://gatk.broadinstitute.org/hc/en-us/articles/360037226472-AddOrReplaceReadGroups-Picard-) (that assigns read groups in the SAM record), [MarkDuplicates](https://gatk.broadinstitute.org/hc/en-us/articles/360037052812-MarkDuplicates-Picard-) (to mark the duplicate reads) and [BaseRecalibrator](https://gatk.broadinstitute.org/hc/en-us/articles/360036898312-BaseRecalibrator) (to create a recalibration table) and [ApplyBQSR](https://gatk.broadinstitute.org/hc/en-us/articles/360036856671-ApplyBQSR)  (to recalibrate the quality scores using the obtained recalibration table) obtaining a ready-bam file along with its index file.
* Then the variant calling is executed using [GATKhaplotypeCaller](https://gatk.broadinstitute.org/hc/en-us/articles/360037225632-HaplotypeCaller) to call the variants of the sample obtaining a vcf file as a output.
* The obtained variants are annotated using the [Ensembl Variant Effect Predictor (vep)](https://www.ensembl.org/info/docs/tools/vep/index.html).
* Finally, the rare variants (population frequency < 0,01%) along with some annotations are converted into an excel file.
* To check the quality of the reads, it is also creted an excel with different parameters values that have been obtained during the pipeline.

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
