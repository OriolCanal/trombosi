import pysam
import subprocess
from subprocess import PIPE
import re
from trombosis import *
from errors import *
import pandas as pd
import io
from bs4 import BeautifulSoup
import json
import vcf



#common variables to define the reference genome path and file
reference_genome_path = "/home/ocanal/vep_data/fasta_file/grch37/"
abs_reference_genome_file = f"{reference_genome_path}GRCh37.p13.genome.fa"
reference_genome_file = "GRCh37.p13.genome.fa"


class FASTQ_file:
    """
    Class with different methods to work with BAM files
    """



    def __init__(self, fastQ_file):
        self.fastq_file = fastQ_file
        self.fastq_absPath = f"{directory}{self.fastq_file}"

    def get_ID(self):
        """
        get the RB of the file
        """
        ID_regex_exp = "[A-Z]{2}[0-9]{5}"
        ID = re.search(ID_regex_exp, self.fastq_file).group(0)
        
        #Raising error if the ID nomenclature is not valid
        if (ID):
            return (ID)
        else:
            errmsg = f"No typical nomenclature with RBXXXXX where X are digits has been found for file {self.bam_file}"
            raise NotTypicalNomenclature(self.bam_file, errmsg)

    def get_lane_number(self):
        """
        obtaine the lane number based on fastq naming convention
        """
        lane_regex_exp = "[L]{2}[0-9]{3}"
        lane = re.search(lane_regex_exp, self.fastq_file).group(0)

        if (lane):
            return(lane)
        else:
            logging.critical(f"The format of the file {self.fastq_file} is not correct, the lane LXXX where X are digits have not been found.")
            
    def quality_check_fastq(self):
        """
        Perform a quality check of the fastq data and generate HTML report
        """

        #define ouptut directory and file
        output_directory = f"{directory}Quality_report"
        output_file = f"{output_directory}/{self.get_ID()}"

        logging.info(f"Running fastQ to perform the quality check. HTML report will be stored in {output_file}")


        #create directory where output will be stored if it doesn't exist
        dirExists = os.path.exists(output_directory)
        if dirExists:
            pass
        else:
            os.makedirs(output_directory)
            print(f"A new directory to store html from fastqc has been created in: {output_directory}")

        # run fastQC
        subprocess.run(["fastqc", "-o", output_file, self.fastq_file])


    def trim_fastq(self):
        """ 
        Trimmomatic will remove low-quality base calls from the beginning and end of the reads, apply a sliding window trimming
        algorithm to remove low-quality bases and discard reads that are shorter that 50 bases after trimming
        """

        #define ouptut file and folder
        output_folder = f"{directory}FastQ_Trimmed"
        trimmed_fastq = f"{output_folder}/{self.get_ID}_trim.fastq"

        #Set the Trimmomatic parameters
        trimmomatic_params = "LEADING:20", "TRAILING:20", "SLIDINGWINDOW:4:20", "MINLEN:50"

        #create folder if it doesn't exist
        dirExists = os.path.exists(output_folder)
        if dirExists:
            pass
        else:
            os.makedirs(output_folder)
            print(f"A new directory to store trimmed fastq has been created in: {output_folder}")

        #run trimmomatic using subprocess.run
        subprocess.run(["trimmomatic", "SE", "-phred33", self.fastq_file, trimmed_fastq, trimmomatic_params])

        return (trimmed_fastq)


    def align_to_reference_hg38 (self, trimmed_fastq):
        """
        Align the fastq file to a reference genome (hg38)
        """
        
        #define ouptut file and folder
        output_folder = f"{directory}SAM"
        sam_file = f"{output_folder}/{self.get_ID()}.sam"

        #create folder if it doesn't exist
        dirExists = os.path.exists(output_folder)
        if dirExists:
            pass
        else:
            os.makedirs(output_folder)
            logging.info(f"A new directory to store sam files has been created in: {output_folder}")
        logging.info("aligning to reference genome")      

        #the docker container of the bwa: biocontainer/bwa can't be installed so I did it using subprocess.run
        subprocess.run(["bwa", "mem", "-M", "-t", "4", abs_reference_genome_file, trimmed_fastq, "-o", sam_file])
        logging.info("alignment finished!")
        return(sam_file)

    def sam_to_sorted_bam (self, sam):
        """
        giving a sam file as input, it creates a sorted bam file with its index
        """
        #define ouptut file and folder
        bam_folder = f"{directory}BAM"
        bam_file = f"{bam_folder}/{self.get_ID()}.bam"
        sorted_bam_file = f"{bam_folder}/{self.get_ID()}_sorted.bam"

        #create folder if it doesn't exist
        dirExists = os.path.exists(bam_folder)
        if dirExists:
            pass
        else:
            os.makedirs(bam_folder)
            print(f"A new directory to store bam files has been created in: {bam_folder}")

        #I can't download biocontainer/samtools so I will run it on command line with subprocess run 

        #create the bam file from the sam file
        cmd1 = f"samtools view -S -b {sam} > {bam_file}"
        logging.info(f"Compressing the SAM file into a BAM file : \n{cmd1}")
        subprocess.run(cmd1, shell=True)

        #sort the bam file
        cmd2 = f"samtools sort {bam_file} -o {sorted_bam_file}"
        logging.info(f"Sorting the BAM file: \n {cmd2}")
        subprocess.run(cmd2, shell=True)

        #create the index bam file
        cmd3 = f"samtools index {sorted_bam_file}"
        logging.info(f"Creating the index file: \n {cmd3}")
        subprocess.run(cmd3, shell=True)

        return(sorted_bam_file)

class Paired_FASTQ_file(FASTQ_file):
    """
    class for fastq paired end reads files
    """

    def __init__(self, fastq_file, paired_fastq_file):
        self.fastq_file = fastq_file
        self.fastq_absPath = f"{directory}{self.fastq_file}"
        self.paired_fastq_file = paired_fastq_file
        self.paired_fastq_absPath = f"{directory}{self.paired_fastq_file}"


    def quality_check_fastq(self):
        """
        Perform a quality check of the fastq data and generate HTML report
        """
        logging.info(f"Performing the quality check of the FastQ files {self.fastq_file, self.paired_fastq_file} ")
        #define ouptut directory and file
        output_directory = f"{directory}/Quality_report"


        #create directory where output will be stored if it doesn't exist
        dirExists = os.path.exists(output_directory)
        if dirExists:
            pass
        else:
            os.makedirs(output_directory)
            print(f"A new directory to store html from fastqc has been created in: {output_directory}")
        print (f"fastqc", "-o", f"{output_directory}/", self.fastq_absPath, self.paired_fastq_absPath)
        # run fastQC
        subprocess.run(["fastqc", "-o", f"{output_directory}/", self.fastq_absPath, self.paired_fastq_absPath])

        return (0)


    def trim_fastq(self):
        """ 
        Trimmomatic will remove low-quality base calls from the beginning and end of the reads, apply a sliding window trimming
        algorithm to remove low-quality bases and discard reads that are shorter that 50 bases after trimming
        """

        #define ouptut file and folder
        output_folder = f"{directory}FastQ_Trimmed"
        forward_paired_trimmed_fastq = f"{output_folder}/{self.get_ID()}_paired_F.fastq"
        reverse_paired_trimmed_fastq = f"{output_folder}/{self.get_ID()}_paired_R.fastq"
        forward_unpaired_trimmed_fastq = f"{output_folder}/{self.get_ID()}_unpaired_F.fastq"
        reverse_unpaired_trimmed_fastq = f"{output_folder}/{self.get_ID()}_unpaired_R.fastq"
        
        #Set the Trimmomatic parameters
        trimmomatic_params = "LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:50"

        #create folder if it doesn't exist
        dirExists = os.path.exists(output_folder)
        if dirExists:
            pass
        else:
            os.makedirs(output_folder)
            print(f"A new directory to store trimmed fastq has been created in: {output_folder}")
        
        trim_cmd = f"trimmomatic PE -phred33 \
        {self.fastq_absPath} \
        {self.paired_fastq_absPath} \
        {forward_paired_trimmed_fastq} \
        {forward_unpaired_trimmed_fastq} \
        {reverse_paired_trimmed_fastq} \
        {reverse_unpaired_trimmed_fastq} \
        {trimmomatic_params}"

        #run trimmomatic using subprocess.run  
        subprocess.run(trim_cmd, shell= True)

        return (forward_paired_trimmed_fastq, forward_unpaired_trimmed_fastq, reverse_paired_trimmed_fastq, reverse_unpaired_trimmed_fastq)


    def align_to_reference_hg38(self, forward_paired_trimmed_fastq, reverse_paired_trimmed_fastq):
        """
        align the paired reads to the reference genome to create a bam file
        """

        #define ouptut file and folder
        sam_folder = f"{directory}SAM"
        sam = f"{sam_folder}/{self.get_ID()}.sam"

        #create folder if it doesn't exist
        dirExists = os.path.exists(sam_folder)

        if dirExists:
            pass
        else:
            os.makedirs(sam_folder)
            print(f"A new directory to store html from fastqc has been created in: {sam_folder}")

        #defining the output files of the indexing necessary to run bwa-mem
        output_index_files = [f"{reference_genome_path}{reference_genome_file}.amb",
        f"{reference_genome_path}{reference_genome_file}.ann",
        f"{reference_genome_path}{reference_genome_file}.pac"]

        #checking if the reference genome have been indexed to run bwa-mem
        for output_index_f in output_index_files:
            if os.path.exists(output_index_f):
                pass
            else:
                logging.critical(f"The index file for the reference genome in order to run bwa-mem has not been found. Please execute the following command: \
                    bwa index {abs_reference_genome_file}")


        #command to run bwa mem for the alignment of the reads to the reference genome hg38
        cmd = ["bwa", "mem",
            "-t 4", #use 3 threads
            "-M", #mark secondary alignments
            "-S", #output in sam format
            abs_reference_genome_file,
            forward_paired_trimmed_fastq,
            reverse_paired_trimmed_fastq,
            #"-o",
            #sam
            ]
        
        print(cmd)

        result = subprocess.run(cmd, cwd=sam_folder, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

        if result.returncode == 0:
            logging.info("Alignment was successful!")
            #redirect the output to the output aligned file
            with open(sam, "wb") as output_file:
                output_file.write(result.stdout)

        
        return (sam)

class Bam_file:
    """
    Class with different methods to work with BAM files
    """

    def __init__(self,bam_file):
        self.bam_absPath = bam_file
        self.bam_file = self.bam_absPath.replace(directory,"")
        
    def get_bamfile (self):
        return (self.bam_file)

    def get_grouped_bamfile (self):
        grouped_dir_list = self.bam_file.split("/")[-2:]
        grouped_dir = "/".join(str(path_item) for path_item in grouped_dir_list)
        
        return (grouped_dir)

    def get_ID(self):
        """
        get the ID of the file. Ex: 
            file: RB22106_9999999.rmdup.bam
            return: RB22106
        """

        regex_exp = "[A-Z]{2}[0-9]{5}"
        ID = re.search(regex_exp, self.bam_file).group(0)
        
        #Raising error if the ID nomenclature is not valid
        if (ID):
            return (ID)
        else:
            errmsg = f"No typical nomenclature with RBXXXXX where X are digits has been found for file {self.bam_file}"
            raise NotTypicalNomenclature(self.bam_file, errmsg)


    def get_full_ID(self):
        """	
        get the full ID of the file. Ex: 
            file: RB22106_9999999.rmdup.bam
            return: RB22106_9999999
        """
        regex_exp ="[A-Z]{2}[0-9]{5}"
        print (self.bam_file)
        full_ID = re.search(regex_exp, self.bam_file).group(0)

        if (full_ID):
            return (full_ID)
        else:
            errmsg = f"No typical nomenclature with RBXXXXX_XXXXXXX where X are digits has been found for file {self.bam_file}"
            raise NotTypicalNomenclature(self.bam_file, errmsg)

    def parse_bam (self):
        bamfile = pysam.AlignmentFile(self.bam_absPath,"rb")
        return(bamfile)
    
    def check_integrity(self):
        """
        testing if files are not corrupted. Data in the middle of the file is not read (matter of time cost) but it
        is useful for testing that files are not truncated
        """

        integrity_cmd = f"samtools quickcheck -v {self.bam_absPath}"
        cmd_output = subprocess.run(integrity_cmd, shell = True)
        cmd_stderr= cmd_output.stderr.strip()
        if cmd_stderr:
            logging.critical (f"Truncated BAM file : {self.bam_file}") 
        else: 
            return (True)

    def get_coverage (self):
        """
        run mosdepth to obtain depth calculations of the BAM file. 
        """
        abs_output_dir = f"{directory}/mosdepth/"
        rel_output_dir = "mosdepth/"
        dir_exists = os.path.exists(abs_output_dir)
        if dir_exists:
            pass 
        else:
            os.mkdir(abs_output_dir)

        mosdepth = "quay.io/biocontainers/mosdepth:0.2.4--he527e40_0"
        mosdepth_outname = f"{self.get_ID()}_mosdepth"
        mosdepth_rel_f = f"{rel_output_dir}{mosdepth_outname}"


        abs_bed = "/home/ocanal/vep_data/tromb_bed/"
        bed_file = "THROMBOSIS.V1.61.CDS.bed"

        cmd = f"docker run -v {directory}:/work_data -v {abs_bed}:/bed {mosdepth} mosdepth --fast-mode --thresholds 1,10,20,30,100,200 --by /bed/{bed_file}  /work_data/{mosdepth_rel_f} /work_data/{self.get_bamfile()}"
        print(f"RUNNING MOSDEPTH TO OBTAIN THE COVERAGE:\n {cmd}")
        subprocess.run(cmd, shell=True)

        #converting the output to a pandas dataframe
        thresholds_output = f"{abs_output_dir}{mosdepth_outname}.thresholds.bed.gz"
        
        colnames = ["chrom", "start", "end", "region", "1X", "10X", "20X", "30X", "100X", "200X"]
        mosdepth_pddf = pd.read_csv(thresholds_output, compression="gzip", names = colnames, sep="\t", comment="#")
        print(mosdepth_pddf)
        return(mosdepth_pddf)

    def reads_stats(self):
        """
        read different stats from the bam file
        """
        cmd_samtools = f"samtools flagstat {self.bam_absPath}"
        samtools_output = subprocess.run(cmd_samtools, shell=True, stdout = PIPE, stderr = PIPE)
        samtools_stderr = samtools_output.stderr
        samtools_stdout = samtools_output.stdout
        
        #print the output of the command to a stats_txt file
        with open(f"{directory}/BAM/{self.get_ID()}_BAM_stats.txt", "a") as stats_f:
            print(samtools_stdout, file=stats_f)

        print (samtools_stdout)
        #taking the different output values and assigning it to variables
        output_lines = []
        for line in samtools_stdout.splitlines():
            output_lines.append(str(line))

            word_list = line.split(" ")
            pattern = re.compile(r"^\d+")
            
            if "total" in line:
                QC_passed_reads = word_list[0]
                QC_failed_reads = word_list[2]
                print (f"QC_passed_reads : {QC_passed_reads}")
                total_reads = QC_failed_reads + QC_passed_reads
                #print (f"passed: {QC_passed_reads}\n failed: {QC_failed_reads}")

            elif "secondary" in line:
                # look for digits in the line
                secondary= word_list[0]
                secondary_failed = word_list[2]

            
            elif "supplementary" in line:
                supplementary = word_list[0]
                supplementary_failed = word_list[2]

            elif "duplicates" in line:
                duplicates = word_list[0]
                duplicates_failed = word_list[2]

            elif "mapped (" in line:
                mapped = word_list[0]
                mapped_failed = word_list[2]
                percentage = word_list[4].replace("(","")

            elif "paired in sequencing" in line:
                paired_in_seq = word_list[0]
                paired_in_seq_failed = word_list[2]

            elif "read1" in line:
                read1 = word_list[0]
                read1_failed = word_list[2]

            elif "read2" in line:
                read2 = word_list[0]
                read2_failed = word_list[2]

            elif "properly paired" in line:
                properly_paired = word_list[0]
                properly_paired_failed = word_list[2]
                properly_paired_percentage = word_list[4].replace("(","")

            elif "with itself and mate mapped" in line:
                both_reads_mapped = word_list[0]
                both_reads_mapped_failed = word_list[2]

            elif "singletons" in line:
                singletons = word_list[0]
                singletons_failed = word_list[2]
                singletons_percentage = word_list[4].replace("(","")

            elif "with mate mapped to a different chr (mapQ>=5)" in line:
                pair_id_dif_chr_mapQ = word_list[0]
                pair_id_dif_chr_mapQ_failed = word_list[2]

            elif "with mate mapped to a different chr" in line:
                pair_in_dif_chr =  word_list[0]
                pair_in_dif_chr_failed =  word_list[2]

            #The format will be type: passed QCfilters + not passed QCfilters
        metrics_dict = {
        "RB" : self.get_ID(),
        "Total": QC_passed_reads + " + " + QC_failed_reads,
        "Secondary": secondary + " + " + secondary_failed,
        "Supplementary": supplementary + " + " + supplementary_failed,
        "Duplicates": duplicates + " + " + duplicates_failed,
        "Mapped" : mapped + " + " + mapped_failed + " + " + percentage + "%",
        "Paired_in_seq" : paired_in_seq + " + " + paired_in_seq_failed,
        "Read1" : read1 + " + " + read1_failed,
        "Read2" : read2 + " + " + read2_failed,
        "Properly_paired" : properly_paired + " + " + properly_paired_failed + " + " + properly_paired_percentage + "%",
        "Both_reads_mapped" : both_reads_mapped + " + " + both_reads_mapped_failed,
        "Singletons" : singletons + " + " + singletons_failed + " + " + singletons_percentage + "%",
        "Pair_in_dif_chr" : pair_in_dif_chr + " + " + pair_in_dif_chr_failed,
        }

        if samtools_stderr:
            logging.warning(f" when analysing file {self.bam_file} STDERR is the following:\n {samtools_stderr}")
        
        print(metrics_dict)
        return (metrics_dict)


    def group_reads (self):
        """
        execute gatk AddOrReplaceReadGroups from docker image
        """

        output_directory = f"{directory}/grouped_reads"
        gatk_version = "broadinstitute/gatk:4.1.3.0"
        output_file = f"grouped_reads/{self.get_ID()}.gr"

        dir_exists = os.path.exists(output_directory)
        if dir_exists:
            pass 
        else:
            os.mkdir(output_directory)
#            -LB Tromb                         \
#          -ID {self.get_ID()}                   \
#          -PL ILLUMINA                          \
#          -PU unit1                             \
#          -SM {self.get_ID()}"

        #replace gropus: RGPL: Platform, RGSM: sample name
        cmd = f"docker run -v {directory}:/work_data -it {gatk_version} gatk AddOrReplaceReadGroups \
          -I /work_data/{self.get_bamfile()}    \
          -O /work_data/{output_file}           \
		  -RGID 4 \
		  -RGLB {self.get_ID()} \
		  -RGPL ILLUMINA \
		  -RGPU unit1 \
		  -RGSM 20"

        subprocess.run(cmd, shell=True)
        path_output_file=f"{output_directory}/{self.get_ID()}.gr"
        return(path_output_file)

    def mark_duplicates (self):
        """
        already done in the files so we won't apply this filter, the reads that are marked as duplicated will be ignored by the gatk haplotype caller.
        """

        picard_path = "/home/ocanal/Desktop/INSTALLATION/picard"
        duplicates_folder = f"{directory}/marked_duplicates_BAM"
        gr_marked_bam_file = f"{duplicates_folder}/{self.get_ID()}markDuplicates.bam" 

        #we will use the grouped bam file as the input file (output of addorreplacegroups)
        gr_bam = f"{directory}grouped_reads/{self.get_ID()}.gr"

        #Create the output folder if it doesn't existreference_ge
        dirExists = os.path.exists(duplicates_folder)
        if dirExists:
            pass
        else:
            os.makedirs(duplicates_folder)
            print(f"A new directory to store html from fastqc has been created in: {duplicates_folder}")


        metrics_file = f"{directory}/BAM/{self.get_ID()}_metrics.txt"

        duplicates_cmd = f"java -jar {picard_path}/picard.jar MarkDuplicates \
            -I {gr_bam}\
            -O {gr_marked_bam_file}\
            -M {metrics_file}"

        subprocess.run(duplicates_cmd, shell= True)

        return (gr_marked_bam_file)


    def quality_recalibration (self):

        #input folder where marked duplicates files is stored
        duplciates_folder = f"{directory}marked_duplicates_BAM"
        #abs path of the input marked dupolicated bam file
        gr_marked_bam_file = f"{duplciates_folder}/{self.get_ID()}markDuplicates.bam" 
        #relative path of marked duplicates file 
        rel_duplicates_file = f"marked_duplicates_BAM/{self.get_ID()}markDuplicates.bam"
        #version of gatk that will be used by docker
        gatk_version = "broadinstitute/gatk:4.1.3.0"
        #output folder
        quality_recab_folder = f"{directory}ready_BAM"
        #abs path of the metrics file output
        rel_metrics_f = f"ready_BAM/{self.get_ID()}_metrics.table"
        #abs path of the final outut quality recalibrated file
        rel_ready_bam = f"ready_BAM/{self.get_ID()}_recab.bam"
        abs_ready_bam = f"{directory}{rel_ready_bam}"
        #reference dbsnp common variant file
        gatk_common_variants_file = "/home/ocanal/vep_data/gatk_resources/hg19/v150/"
        
        dir_exists = os.path.exists(quality_recab_folder)
        if dir_exists:
            pass 

        else:
            os.mkdir(quality_recab_folder)   

        #create metrics file for base recalibrator (If I don't open the file an error appears: the metrics file doesn't exist) 

        recal_cmd = f"docker run \
        -v {directory}:/work_data \
        -v $HOME/vep_data:/opt/vep/.vep:Z \
        -v {gatk_common_variants_file}:/dbSNP \
        -it {gatk_version} gatk BaseRecalibrator \
        -I /work_data/{rel_duplicates_file} \
        -R /opt/vep/.vep/fasta_file/grch37/GRCh37.p13.genome.fa \
        --known-sites /dbSNP/00-All.vcf.gz\
        -O /work_data/{rel_metrics_f}"
        print(recal_cmd)
        subprocess.run(recal_cmd, shell=True)

        #apply base recalibration based on metrix file
        apply_recal_cmd = f"docker run \
        -v {directory}:/work_data \
        -v $HOME/vep_data:/opt/vep/.vep:Z \
        -it {gatk_version} gatk ApplyBQSR \
        -R /opt/vep/.vep/fasta_file/grch37/GRCh37.p13.genome.fa \
        -I /work_data/{rel_duplicates_file} \
        --bqsr-recal-file /work_data/{rel_metrics_f} \
        -O /work_data/{rel_ready_bam}"

        subprocess.run(apply_recal_cmd, shell = True)

        return(abs_ready_bam)

    def index_bam (self):
        """
        create a index bam in the grouped folder
        """
        rel_ready_bam = f"ready_BAM/{self.get_ID()}_recab.bam"
        abs_ready_bam = f"{directory}{rel_ready_bam}"
        abs_ready_bam_ind = f"{abs_ready_bam}.bai"

        cmd = f"samtools index {abs_ready_bam} {abs_ready_bam_ind}"
        subprocess.run(cmd, shell=True)

    def haplotype_caller (self):
        """
        Executing haplotype caller (standard) with hg.38 as reference. Indel realignment is already performed by the haplotype caller since version 4
        """

        output_directory = f"{directory}/vcf"
        dir_exists = os.path.exists(output_directory)
        rel_output_file = f"vcf/{self.get_full_ID()}.vcf.gz"
        rel_ready_bam = f"ready_BAM/{self.get_ID()}_recab.bam"
        
        if dir_exists:
            pass 

        else:
            os.mkdir(output_directory)

        reference_genome = "GRCh37.p13.genome.fa" 
        gatk_version = "broadinstitute/gatk:4.1.3.0"
        #bed file is where the regions of interest are defined
        #add these 2 lines to define the regions of interest of the genome (not used because are defined in grch19 and we are using grch38)
        #bed_file = "THROMBOSIS.V1.61.CDS.bed"
        #          -L /reference/{bed_file} \

        cmd = f"docker run \
        -v {directory}:/work_data \
        -v {reference_genome_path}:/reference \
        -it {gatk_version} gatk HaplotypeCaller \
        -R /reference/{reference_genome} \
        -I /work_data/{rel_ready_bam} \
        -O /work_data/{rel_output_file}"

        print (f"Running Haplotype caller: \n {cmd}")
        subprocess.run(cmd, shell=True )


class Vcf_Class:
    """
    Class with different methods to work with VCF files
    """

    def __init__(self,vcf_file):
        self.vcf_file = vcf_file
        self.vcf_absPath = directory + "vcf/" + self.vcf_file
        self.abs_directory = directory 
        self.rel_path = f"vcf/{self.vcf_file}"
        self.abs_output_dir = f"{directory}annotations"
        self.abs_vep_file = f"{self.abs_output_dir}/{self.get_full_ID()}"
        self.excel_name = f"{self.get_full_ID()}.xlsx"
        self.excel_path = f"{self.abs_directory}/excels/"
        self.excel_file_path = f"{self.excel_path}{self.excel_name}"

    def get_full_ID(self):
        regex_exp ="[A-Z]{2}[0-9]{5}"
        full_ID = re.search(regex_exp, self.vcf_file).group(0)
        return (full_ID)

    def variant_annotation (self):
        """
        we need the cache saved in /opt/vep/.vep which can be obtained with the commands: 

            cd $HOME/.vep
            curl -O https://ftp.ensembl.org/pub/release-108/variation/indexed_vep_cache/homo_sapiens_vep_108_GRCh38.tar.gz
            tar xzf homo_sapiens_vep_108_GRCh38.tar.gz
        """ 

        #INSTALL DOCKER VEP PLUGINS
        #check if the folder already exists:
        
        dir_exists = os.path.exists(self.abs_output_dir)
        os.system(f"chmod 777 {self.abs_output_dir}")
        if dir_exists:
            pass  

        else:
            os.mkdir(self.abs_output_dir)
            os.system(f"chmod 777 {self.abs_output_dir}")

        vep = "ensemblorg/ensembl-vep"
        rel_ouptut_file = f"annotations/{self.get_full_ID()}"
 
        # REVEL plugin adds the REVEL score for missense variants to VEP output
        path_revel_data = "/opt/vep/.vep/REVEL/new_tabbed_revel.tsv.gz"
        #gnomAD includes the frequencies of the alternative allele in different populations
        path_gnomAD_data = "/opt/vep/.vep/gnomAD/gnomad.genomes.tabbed.tsv.gz"
        #maxEntScan get splice sites predictions
        maxEntScan_data = "/opt/vep/.vep/maxEntScan/maxent_data/"
        #SplicaAI data
        path_SpliceAI_snv = "/opt/vep/.vep/SpliceAI/ghrc38/spliceai_scores.raw.snv.hg38.vcf.gz"
        path_SpliceAI_indel = "/opt/vep/.vep/SpliceAI/ghrc38/spliceai_scores.raw.indel.hg38.vcf.gz"

        #CADD path
        path_snvCADD = "/opt/vep/.vep/CADD/gnomad.genomes.r2.1.1.snv.tsv.gz"
        path_indelsCADD = "/opt/vep/.vep/CADD/gnomad.genomes.r2.1.1.indel.tsv.gz"

        #The :Z in the docker volume is to give root privilegies to the mounted volume (to read/write files)

        #              --plugin gnomADc,{path_gnomAD_data}\
        #          --plugin REVEL,{path_revel_data}\
        #          --plugin SpliceAI,snv={path_SpliceAI_snv},indel={path_SpliceAI_indel}\
        #          --plugin REVEL,{path_revel_data}\
        #          --plugin CADD,{path_snvCADD},{path_indelsCADD}\
      

        cmd = f"docker run \
        -v {self.abs_directory}:/work_data \
        -v {reference_genome_path}:/reference \
        -v $HOME/vep_data:/opt/vep/.vep:Z \
        -it {vep} vep \
        --cache --offline \
        --fasta /reference/{reference_genome_file} \
        --format vcf \
        --tab \
        --assembly GRCh37 \
        --hgvsg --everything --force_overwrite \
        -i /work_data/{self.rel_path} \
        --plugin MaxEntScan,{maxEntScan_data},SWA \
        -o /work_data/{rel_ouptut_file}"

        print(f"RUNNING VEP:\n {cmd}")
        run = subprocess.run(cmd, shell=True)


    def vep_parser (self):
        """ 
        Parse the information obtained by the vep annotated file and returns a pandas dataframe with the annotated variants
        where each row is a variant and each columns a property of the variant (name, poblational frequency etc...) (108 columns)
        obtained from the annotated file (/annotations/ID)
        """
        #Open the ouptut vep file
        with open(self.abs_vep_file) as f:
            # it contains the variable position information (excludint the header from the file)
            lines = [l for l in f if not l.startswith("##")]
        #From the annotation lines return a pandas dataframe rows:variable_sites, columns:different_annotation.
        return (pd.read_csv(
            io.StringIO(''.join(lines)), sep='\t').rename(columns={'#Uploaded_variation': 'Uploaded_variation'}))
        
        #Insert a column named RB at the 0th column of the dataframe with the RB identifier of the sample
        #final_pd = pd_dataframe.insert(0, 'RB' ,self.get_full_ID())
        #print(final_pd)
        #return (final_pd)


    def add_RB_column (self, pd_df):
        

        pd_df["RB"] = self.get_full_ID()
        print(pd_df)
        return(pd_df)

    def filtering_df_by_AlleFreq (self, pd_df, freq):
        """ 
        filters the annotated variants by the poblational frequency of the variant
        """

        #Extracting variants that have a not null value on maximal poblational frequency
        pd_df_FreqNotNull = pd_df[pd_df.MAX_AF != '-']
        
        #converting MAX_AF values into floats to perform the comparasion with freq
        pd_df_FreqNotNull['MAX_AF'] = pd_df_FreqNotNull['MAX_AF'].astype(float)
    
        return(pd_df_FreqNotNull.loc[pd_df_FreqNotNull['MAX_AF'] <= freq])


    def vcf_to_pd (self):
        """
        obtain the Quality of a variant from the location of the variant
        """

        #column names to give to the pandas dataframe (using the same as the vcf file)
        colnames = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "20"]

        return(pd.read_csv(self.vcf_absPath, compression="gzip", names=colnames, sep="\t", comment="#"))


    def list_uploaded_variation (self, pd_df):
        """From a pandas dataframe, it obtains a list of the Uploaded_variation values"""
        chr_position_list = []
        
        for row in pd_df.index:
            chr_position_list.append(pd_df["Uploaded_variation"][row])
        
        return(chr_position_list)
        
    def split_uploaded_variation (self, chr_position_list):
        """From a list of uploaded_variation values, we create a list of list containing: [chr1, position, alleles] for each variant"""
        chr_position = []
        chr_position_listoflist = []
        for var in chr_position_list:
            items = var.split("_")
            chr_position = [items[0], int(items[1])]
            print (f"chr position  {chr_position}")
            chr_position_listoflist.append(chr_position)
        
        return(chr_position_listoflist)

    def extract_vcf_column (self, pd_df, property = str, chrom = str, pos = int):
        """ 
        extract a column of the vcf file variant givent the position of the variant in the following format:
        chrom = chr1
        pos = 1554362
        We should give to the function one of the colnames to extract the information
        """

        #qual = pd_df["QUAL"].where((pd_df["POS"] == "1554362") & (pd_df["CHROM"] == "chr1"))
        #print(f"pos = {pos}")
        return(pd_df.loc[(pd_df["CHROM"] ==chrom) & (pd_df["POS"] == pos), property])
    
    def extract_variant_coverage (self, mosdepth_pddf, chrom=str, pos= int):
        x1 = (mosdepth_pddf.loc[(mosdepth_pddf["chrom"] == chrom) & (mosdepth_pddf["start"]<= pos) & (mosdepth_pddf["end"]>pos), "1X"])
        x10 = (mosdepth_pddf.loc[(mosdepth_pddf["chrom"] == chrom) & (mosdepth_pddf["start"]<= pos) & (mosdepth_pddf["end"]>pos), "10X"])
        x20 = (mosdepth_pddf.loc[(mosdepth_pddf["chrom"] == chrom) & (mosdepth_pddf["start"]<= pos) & (mosdepth_pddf["end"]>pos), "20X"])
        x30 = (mosdepth_pddf.loc[(mosdepth_pddf["chrom"] == chrom) & (mosdepth_pddf["start"]<= pos) & (mosdepth_pddf["end"]>pos), "30X"])
        x100 = (mosdepth_pddf.loc[(mosdepth_pddf["chrom"] == chrom) & (mosdepth_pddf["start"]<= pos) & (mosdepth_pddf["end"]>pos), "100X"])
        x200 = (mosdepth_pddf.loc[(mosdepth_pddf["chrom"] == chrom) & (mosdepth_pddf["start"]<= pos) & (mosdepth_pddf["end"]>pos), "200X"])
        start = int((mosdepth_pddf.loc[(mosdepth_pddf["chrom"] == chrom) & (mosdepth_pddf["start"]<= pos) & (mosdepth_pddf["end"]>pos), "start"]))
        end = int((mosdepth_pddf.loc[(mosdepth_pddf["chrom"] == chrom) & (mosdepth_pddf["start"]<= pos) & (mosdepth_pddf["end"]>pos), "end"]))
        bases = end - start
        print(x1, x10, x20, x30, x100, x200, start, end, bases)
        if x200 >= bases:
            return ("200x")
        elif x100 >= bases:
            return ("100x")
        elif x30 >= bases:
            return ("30x")
        elif x20 >= bases:
            return("20x")
        elif x10 >= bases:
            return("10x")
        elif x1 >= bases:
            return("1x")
        else:
            return ("not_covered")

    def add_coverage_df (pddf, coverage_list):
        pddf["coverage"]= coverage_list
        return(pddf)

    def extract_ref_alt_reads ( self, values_list):
        """
        Get the total amount of reads of a variant position,
        the total amount of alternative reads supporting a variant,
        and the % of reads supporting a variant"""

        alt_reads_list = []
        total_reads_list = []
        percent_alt_list = []
        for values in values_list:
            params = values[0].split(":")

            #AD where it is stored ref_reads,alt_reads is situated in params[1]
            params_list=params[1].split(",")
            ref_reads = int(params_list[0])
            alt_reads = int(params_list[1])

            total_reads = alt_reads + ref_reads
            percent_alt = (alt_reads/(total_reads))*100

            alt_reads_list.append(alt_reads)
            total_reads_list.append(total_reads)
            percent_alt_list.append(percent_alt)

        return(alt_reads_list, total_reads_list, percent_alt_list)


    def add_df_qual( self, pd_df, qual):
        """
        Adds the quality of the reads extracted from the vcf file to the pandas dataframe
        """
        pd_df["QUAL"] = qual
        return(pd_df)


    def reads_counts_to_df (self, pd_df, alt_reads_list, total_reads_list, percent_alt_list):
        
        pd_df["count_alternative_reads"] = alt_reads_list
        pd_df["count_total_reads"] = total_reads_list
        pd_df["%_alt_reads"] = percent_alt_list
            
        return(pd_df)
        


    def extracting_columns (self, pandas_df):
        """
        From the pandas dataframe create an excel with the columns required
        """
        # the ones that I have to add once I have applied the pluggins:
        # "SpliceAI_pred", "CADD_PHRED", "CADD_RAW", "REVEL",
        #selecting the columns of interest to extract from the pandas dataframe 
        columns = ["RB", "Uploaded_variation", "coverage", "QUAL", "%_alt_reads", "count_total_reads", "count_alternative_reads", "Consequence", "Gene", "SWISSPROT","cDNA_position", "CDS_position", "Protein_position",
        "HGVSc", "HGVSp", "HGVSg", "MANE_SELECT", "MANE_PLUS_CLINICAL", "CLIN_SIG", "Existing_variation" ,"UNIPROT_ISOFORM", "MAX_AF", "PHENO",
        "PUBMED", "Feature_type", "MES-SWA_acceptor_alt", "MES-SWA_acceptor_diff", "MES-SWA_acceptor_ref", "MES-SWA_acceptor_ref_comp",
        "MES-SWA_donor_alt", "MES-SWA_donor_diff", "MES-SWA_donor_ref", "MES-SWA_donor_ref_comp"]
        #creating a pd dataframe with the columns of interest
        print(pandas_df[columns])
        return(pandas_df[columns])

    def df_to_excel (self, final_pd):
        """Create an excel in /directory/excels/ with the data extracted from vep with the selecting columns from extracting_columns function"""

        #create the excel directory if it does not exist
        dir_exists = os.path.exists(self.excel_path)
        if dir_exists:
            pass  

        else:
            os.mkdir(self.excel_path)

        #extracting the pandas dataframe to an excel

        final_pd.to_excel(self.excel_file_path)
        logging.info(f"Results are written to excel {self.excel_file_path} successfully!")








