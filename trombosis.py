import os
import subprocess
from classes import *
from errors import *
from multiprocessing import Pool
import argparse
import logging


parser = argparse.ArgumentParser(description="Pipeline to automatize the process from fastq files to vcf annotations")

parser.add_argument("-d", "--input_directory", action="store", help="input directory where bam files are stored")
parser.add_argument("-r", "--reference", required=False, default="/home/ocanal/vep_data/fasta_file/grch37/GRCh37.p13.genome.fa", action = "store", help="input directory where the fq reference files with its index and dict files are stored")
parser.add_argument("-x", "--remove_tmp_output", required=False, action="store_true", help="remove temporary outputs (Sam, bam, ready_bam, vcf ...) when the program has finished")
parser.add_argument("-b", "--bed_file", required=False, default="/home/ocanal/vep_data/tromb_bed/THROMBOSIS.V1.61.CDS.bed", help="path to your bed file")
parser.add_argument("-np", "--no_pluggins", required=False, action="store_false", help="if -np is specified, no pluggins will be used for the annotation process using vep, otherwise, CADD, SpliceAI, REVEL, maxEntScan and gnomAD will be used.")

args = parser.parse_args()
abs_bed_file = args.bed_file
bed_file = os.path.basename(abs_bed_file)
bed_folder = os.path.dirname(abs_bed_file)
directory = args.input_directory
if not directory.endswith("/"):
    directory += "/"

reference = args.reference
remove_tmp_output = args.remove_tmp_output
pluggins_activated = args.no_pluggins




root_logger = logging.getLogger()
root_logger.setLevel(logging.DEBUG)
handler = logging.FileHandler(f"{directory}/logging.log", "a", "utf-8")
# logging.getLogger().addHandler(logging.StreamHandler())

# checking existance of reference genome and bed file
input_files = [abs_bed_file, reference]

for input_file in input_files:
    if not os.path.exists(input_file):
        raise FilesNotFound(f"The reference genome or bed file did not exists {input_file}, \
                            you can use the -b or -r flags to indicate where they are located")

root_logger.addHandler(handler)

run_name = os.path.basename(os.path.normpath(directory))



def uncompress_files(directory):
    """ From a given directory it decompress all .tgz files"""
    cmpr_dirs = []

    subdirs = os.listdir(directory)
    for subdir in subdirs:
        subdir_full_path = f"{directory}{subdir}"
        if subdir_full_path.endswith("fastq.gz"):
            cmpr_dirs.append(subdir_full_path)
            uncmpr_cmd = f"gzip -d {subdir_full_path}"
            logging.info(f"uncompressing file {subdir}")
            subprocess.run(uncmpr_cmd, shell=True)
    return (cmpr_dirs)

def compress_files(directory):

    dirs = os.listdir(directory)
    for direct in dirs:
        direct = f"{directory}{direct}"
        if direct.endswith(".fastq"):
            cmpr_file_cmd = f"gzip {direct}"
            logging.info(f"compressing file: {direct}")
            subprocess.run(cmpr_file_cmd, shell=True)


def list_fastq_files(directory):
    """
    From a given directory it lists you all the fastq files
    """
    fastq_files = []
    files = [f for f in os.listdir(directory) if os.path.isfile(os.path.join(directory,f))]
    if not files:
        raise (FilesNotFound("The input directory contains no files on it"))
    for file in files:
        if file.endswith(".fastq"):
            fastq_files.append(file)
    if not fastq_files:
        logging.critical = f"No fastq files found in the input directory: {directory}"        

    return (fastq_files)


def match_paired_fastq(fastq_files):
    """
    From a list of fastq files, it mathes the ones with the same RB creating a dictionary:
    paired_fq[RB]= 1 / 2

    2 if there are paired reads, otherwise 1.

    If there are more than 2 fastq files with the same ID, it returns an error
    """
    regex_exp = "[A-Z]{2}[0-9]{5}"

    fq_dictionary = {}

    for fq_file in fastq_files:

        if re.match(regex_exp, fq_file): 

            # obtain the RB of the file
            ID = re.search(regex_exp, fq_file).group(0)

            # if the RB have not already been added to the dictionary
            if ID not in fq_dictionary:

                # create a dictionary with RB = 1
                fq_dictionary[ID] = 1

            # If the RB already exists in the dictionary
            else:

                # add 1 to the current dictionary value
                fq_dictionary[ID] += 1
                if fq_dictionary[ID] > 2:
                    logging.critical(f"Be aware, you have more than 2 files \
                                        with the following RB {ID}")
        else:
            logging.warning(f"The file {fq_file} have not the typical \
                            nomenclature and won't be analysed")
    paired_reads_IDs = [key for key, value in fq_dictionary.items() if value == 2]
    single_reads_IDs = [key for key, value in fq_dictionary.items() if value == 1]

    return (single_reads_IDs, paired_reads_IDs)


def single_reads_IDs_to_fqfile(fq_file_list, ID):

    match = [fq_file for fq_file in fq_file_list if re.search(ID, fq_file)]
    se_fq_file = str(match[0])
    return (se_fq_file)


def paired_reads_IDs_to_fqfile(fq_file_list, ID):

    match = [fq_file for fq_file in fq_file_list if re.search(ID, fq_file)]
    pe_fq_files = list(match)

    return (pe_fq_files)


def list_readybam_files(directory):
    """
    From a given directory it lists you all the bam files
    """
    ready_bam_files = []
    ready_bam_dir = f"{directory}/grouped_reads"
    files = [f for f in os.listdir(ready_bam_dir) if os.path.isfile(os.path.join(ready_bam_dir,f))]
    for file in files:
        if file.endswith(".bam"):
            ready_bam_files.append(file)

    return (ready_bam_files)

#UNCOMMENT TO ANALYSE SINGLE READS
# def single_reads_fastQ_to_bam(fq_file):
#     """
#     Different steps that are run by the pipeline.
#     """

#     logging.info(f"Starting to analyse the file {fq_file}")

#     #Create a FastQ class 
#     FastQ_class = FASTQ_file(fq_file)

#     #Perform a quality check and generate a HTML report

#     FastQ_class.quality_check_fastq()

#     #Use Trimmomatic to remove low-quality base calls and discard short reads (<30)
#     trimmed_fqfile = FastQ_class.trim_fastq()

#     #Align the trimmed file to the reference genome Hg38
#     sam_file = FastQ_class.align_to_reference_hg37(trimmed_fqfile)

#     #Calculate %enrichment = number of raeds aligned in bed region / total number of reads aligned
#     FastQ_class.enrichment(sam_file)

#     #Convert the sam file to a sorted bam file along with its index
#     sorted_bam = FastQ_class.sam_to_sorted_bam(sam_file)

#     logging.info(f"Sorted Bam file have been created for the input file {fq_file}")

#     return (sorted_bam)


def paired_reads_fastq_to_bam(fq_file_1, fq_file_2, qual_dict):
    """
    giving a paired reads fastq it will create a bam
    """

    logging.info(f"Starting to analyse the paired reads files {fq_file_1, fq_file_2}")

    # create a paired FASTQ class
    Paired_FastQ_class = Paired_FASTQ_file(fq_file_1, fq_file_2)

    # Perform a quality check and generate a HTML report
    Paired_FastQ_class.quality_check_fastq(qual_dict)

    # Use Trimmomatic to remove low quality base calls and discard short reads (<30).
    # output is trimmed_fq that is a tuple where:
    # trimmed_fq[0] = forward paired trimmed fastQ file
    # trimmed_fq[1] = forward unpaired trimmed fastQ file
    # trimmed_fq[2] = reverse paired trimmed fastQ file
    # trimmed_fq[3] = reverse unpaired trimmed fastQ file
    # trimmed_fq[4] = quality_dict with trimmomatic output quality parameters
    trimmed_fq = Paired_FastQ_class.trim_fastq(qual_dict)
    qual_dict = trimmed_fq[4]

    for key, value in qual_dict.items():
        print(key, value)

    # align the paired reads to the reference genome
    sam_file = Paired_FastQ_class.align_to_reference_hg37(trimmed_fq[0], trimmed_fq[2])

    # Convert the sam file to a sorted bam file along with its index
    sorted_bam = Paired_FastQ_class.sam_to_sorted_bam(sam_file)

    # Calculate %enrichment = number of raeds aligned in bed region / total number of reads aligned
    qual_dict = Paired_FastQ_class.enrichment(sam_file, qual_dict)

    logging.info(f"Sorted Bam file have been created for the input file {fq_file_1, fq_file_2}")

    return (sorted_bam, qual_dict)


def bam_to_final_excel(sorted_bam, qual_dict):
    """
    Function that applies the different functions created in classes.py to obtain an excel with the 
    rare variants from a sorted bam file
    """

    logging.info(f"Variant annotation have started for file {sorted_bam}")

    # Bam class is created
    Bam_class = Bam_file(sorted_bam)
    # print(Bam_class.bam_file)

    # The integrity of the BAM file is checked
    # Bam_class.check_integrity()

    # dictionary of the reads stats of the BAM file is created
    qual_dict = Bam_class.reads_stats(qual_dict)

    mosdepth_pddf = Bam_class.get_coverage()
    # AddOrReplaceReadGroups is executed and output files are stored in /grouped_reads/
    # The Bam file is ready to be proceeded by haplotype caller
    Bam_gr = Bam_class.group_reads()

    # mark duplicates
    marked_bam = Bam_class.mark_duplicates()

    ready_bam= Bam_class.quality_recalibration()

    # We take the ready_bam as the working file from now on and create a class for this file  
    ready_Bam_class = Bam_file(ready_bam)

    # Index the bam file (.bai) to run haplotype caller 
    ready_Bam_class.index_bam()

    # Perform the variant calling
    ready_Bam_class.haplotype_caller()

    # Create a vcf class
    vcf_file = f"{ready_Bam_class.get_full_ID()}.vcf.gz"

    Vcf_file_Class = Vcf_Class(vcf_file)

    # run variant annotations 
    Vcf_file_Class.variant_annotation()

    # Create a pandas dataframe (variants in rows and parameters as columns)
    vep_pddf = Vcf_file_Class.vep_parser()

    # We add an identifier column where it will be stored the RB of the sample
    vep_pddf_RB = Vcf_file_Class.add_RB_column(vep_pddf)

    # Filter the rare variants (population frequency < 0.1%) We take a frequency of 0.001 considering that vep gives us the 
    vep_filtered = Vcf_file_Class.filtering_df_by_AlleFreq(vep_pddf_RB, 0.001)

    # creting a dataframe to store all the information from the vcf file
    vcf_pd = Vcf_file_Class.vcf_to_pd()

    # From a pandas dataframe, it obtains a list of the Uploaded_variation values
    vep_variation = Vcf_file_Class.list_uploaded_variation(vep_filtered)

    # From a list of uploaded_variation values, we create a list of list containing: [chr1, position, alleles] for each variant
    vep_variation_list = Vcf_file_Class.split_uploaded_variation(vep_variation)

    # calculating call rate and exon lost
    qual_dict = Vcf_file_Class.call_rate_perc(mosdepth_pddf, qual_dict)

    # comparing the vep vs vcf position of variants as in deletions may be different
    vep_adj_variation_list = Vcf_file_Class.comparing_vcf_vep(vcf_pd, vep_variation_list)

    # get mean coverage
    Vcf_file_Class.get_mean_coverage(qual_dict)

    # List where the quality of rare variants are stored
    rare_variant_qual = []
    rare_variant_info = []
    coverage_variant = []

    # For every rare variant, we obtain the quality of the variant and append to rare_variant_list   
    for vep_variation in vep_adj_variation_list:
        print(f"uploaded_variation[0]: \n {vep_variation[0]} \n uploadedvariation[1]: {vep_variation[1]}")
        # print(f"vcf pd: \n {vcf_pd}")
        rare_variant_qual.append(list(Vcf_file_Class.extract_vcf_column(vcf_pd, "QUAL" , chrom = vep_variation[0], pos = vep_variation[1])))
        print("rare_variant_qual")
        rare_variant_info.append(list(Vcf_file_Class.extract_vcf_column(vcf_pd, "THROMBOSI" , chrom = vep_variation[0], pos = vep_variation[1])))
        coverage_variant.append(str(Vcf_file_Class.extract_variant_coverage(mosdepth_pddf, chrom = vep_variation[0], pos = vep_variation[1])))

    print(f"rare variants info: {rare_variant_info}")

    # getting list of reads supporting alternative allele, total reads covering a varaint and % reads supporting a variant
    alt_reads_list, total_reads_list, percent_alt_list = Vcf_file_Class.extract_ref_alt_reads(rare_variant_info)

    # Add the qual parameter to the vep dataframe
    qual_df = Vcf_file_Class.add_df_qual(vep_filtered, rare_variant_qual)

    # Adding the total count of allele reads, the reads supporting the variant and the % of reads supporting the variant
    reads_df =Vcf_file_Class.reads_counts_to_df(qual_df, alt_reads_list, total_reads_list, percent_alt_list)

    coverage_df = Vcf_file_Class.add_coverage_df(reads_df,coverage_variant)

    mane_pddf = Vcf_file_Class.create_mane_pddf("/home/ocanal/vep_data/mane_data/MANE.GRCh38.v1.0.summary.txt")
    feature_list = Vcf_file_Class.extract_feature(coverage_df)
    mane_list = []
    for feature in feature_list:
        mane = Vcf_file_Class.extract_mane(mane_pddf, feature)
        mane_list.append(mane)

    complete_df = Vcf_file_Class.add_mane_to_pddf(mane_list, coverage_df)
    # Extracting the columns that we want in the final excel
    reduced_pandas = Vcf_file_Class.extracting_columns(complete_df)
    # print(coverage_df)

    # Creating the resulting excel
    Vcf_file_Class.df_to_excel(reduced_pandas)
    print(f"qual_dict: \n {qual_dict}\n end qual dict")
    logging.info(f"The annotations have been performed correctly and the rare variants have been annotated to the excel\n\n")
    return (qual_dict)


def append_qual_pddf(qual_pddf, qual_dict):

    qual_pddf = qual_pddf.append(qual_dict, ignore_index = True)
    return (qual_pddf)


def join_excels(parent_dir):

    # Initialize an empty list to store the data
    all_data = []

    # Iterate over the subdirectories
    for subdir in os.listdir(parent_dir):
        subdir_path = os.path.join(parent_dir, subdir)

        # Check if the current item is a directory
        if os.path.isdir(subdir_path):
            # Check if the subdirectory contains an "Excel" directory
            excel_dir_path = os.path.join(subdir_path, 'excels')
            if os.path.isdir(excel_dir_path):
                # Iterate over the Excel files in the excels directory
                for file in os.listdir(excel_dir_path):
                    # Check if the file is an Excel file
                    if file.endswith('.xlsx'):
                        file_path = os.path.join(excel_dir_path, file)
                        # Read the Excel file into a pandas DataFrame
                        df = pd.read_excel(file_path)
                        # Append the data to the list
                        all_data.append(df)

    # create a joined_excels folder, if it haven't already exists

    output_directory = f"{parent_dir}/joined_excels"
    dir_exists = os.path.exists(output_directory)
    if dir_exists:
        pass 

    else:
        os.mkdir(output_directory)

    # Concatenate all the data into a single DataFrame
    all_data_df = pd.concat(all_data, ignore_index=True)
    all_data_sorted = all_data_df.sort_values('Uploaded_variation')
    # Save the DataFrame to an Excel file
    all_data_df.to_excel(f'{output_directory}/variants_joined.xlsx', index=False)
    all_data_sorted.to_excel(f'{output_directory}/variants_joined_sorted.xlsx', index = False)
    logging.info("Joined excel have been created successfully")


def rearrange_df_cols(qual_df):
    for col in qual_df.columns:
        print(col)
    cols = ["Run",
            "RB",
            "fastqc_version",
            "encoding",
            "numb_exons",
            "seq_len_R1",
            "seq_len_R2",
            "GC_R1_perc",
            "GC_R2_perc",
            "total_sequences_R1",
            "total_sequences_R2",
            "total_read_pairs",
            "surviving_reads_pairs",
            "surviving_reads_pairs_perc",
            "dropped_trim",
            "dropped_perc_trim",
            "only_forward_drop",
            "only_forward_drop_perc",
            "only_reverse_drop",
            "only_reverse_drop_perc",
            "dropped_trim",
            "dropped_perc_trim",
            "GC_cont_dist_R1",
            "GC_cont_dist_R2",
            "adapter_content_R1",
            "adapter_content_R2",
            "base_qual_cont_R1",
            "base_qual_cont_R2",
            "overrepresented_seq_R1",
            "overrepresented_seq_R2",
            "per_base_N_cont_R1",
            "per_base_N_cont_R2",
            "per_base_seq_qual_R1",
            "per_base_seq_qual_R2",
            "per_tile_seq_qual_R1",
            "per_tile_seq_qual_R2",
            "seq_duplic_level_R1",
            "seq_duplic_level_R2",
            "seq_len_dist_R1",
            "seq_len_dist_R2",
            "seq_qual_score_R1",
            "seq_qual_score_R2",
            "seq_poor_quality_R1",
            "seq_poor_quality_R2",
            "reads_aligned_in_bed",
            "Total_mapped_reads",
            "enrichment_factor_%",
            "mean_coverage",
            "x1_call_percent",
            "x10_call_percent",
            "x20_call_percent",
            "x30_call_percent",
            "x100_call_percent",
            "x200_call_percent",
            "x1_exon_lost",
            "x10_exon_lost",
            "x20_exon_lost",
            "x30_exon_lost",
            "x100_exon_lost",
            "x200_exon_lost"]
    
    qual_df = qual_df[cols]
    return (qual_df)


def remove_temporary_files (directory):
    subdirs_to_remove = ["BAM", "FastQ_Trimmed", "grouped_reads", "marked_duplicates_BAM", "mosdepth", "Quality_report", "SAM"]
    for dir in subdirs_to_remove:
        final_dir = os.path.join(directory,dir)
        try:
            os.rmdir(final_dir)
        except OSError as e:
            print(f"Error: {final_dir} {e.strerror}")


if __name__ == "__main__":

    # Step 1: decompress files and csorted_bamreate a list of the existing BAM files
    uncompress_files(directory)

    # Step 2 detect the fastq files 
    fq_files = list_fastq_files(directory)

    # and seperate in 2 different lists the single end reads and the paired end reads
    ID_fq_single_reads_list, ID_fq_paired_reads_list = match_paired_fastq(fq_files)

    #UNCOMMENT TO ANALYSE SINGLE READS
    #if there are some single end reads fastq files
    # if ID_fq_single_reads_list:
        
    #     #analyse each fastq file
    #     for ID_single_reads in ID_fq_single_reads_list:
    #         logging.info (f"The fastQ file: {fq_file} will be analysed by the pipeline")
            
    #         #match the single reads IDs with its files
    #         fq_file = single_reads_IDs_to_fqfile(fq_file, ID_fq_single_reads_list)
            
    #         #perform the analysis
    #         sorted_bam = single_reads_fastQ_to_bam(fq_file)
    #         bam_to_final_excel(sorted_bam)

    # if there are some paired end reads fastq files:
    if ID_fq_paired_reads_list:
        
        # define the pandas dataframe where qualityparams will be stored
        qual_pddf = pd.DataFrame()
        # analyse each paired end reads
        for ID_paired_reads in ID_fq_paired_reads_list:
            qual_dict = {}
            # match the paired reads files associated with the ID
            paired_list = paired_reads_IDs_to_fqfile(fq_files, ID_paired_reads)
            print(paired_list)
            sorted_bam, qual_dict = paired_reads_fastq_to_bam(paired_list[0], paired_list[1],qual_dict)
            qual_dict2 = bam_to_final_excel(sorted_bam, qual_dict)
            qual_dict.update(qual_dict2)
            # striping new line characters from dict values
            for key, value in qual_dict.items():
                print(f"key {key} \n value: {value}")
                if type(value) == str:
                    qual_dict[key] = value.rstrip()
            
            qual_pddf = append_qual_pddf(qual_pddf, qual_dict)
            print(qual_pddf)
            qual_pddf_rearranged = rearrange_df_cols(qual_pddf)
            qual_pddf_rearranged.to_excel(f"{directory}{run_name}_QC.xlsx", index = False)
    # remove temporary intermediate (BAM, SAM, vcf ...) directories if the -x command line argument is set
    if remove_tmp_output != False:
        remove_temporary_files(directory)


    
    
    
    #WHEN ANALYSISNG FROM BAM FILES 
    #bams_to_analyse =[]
    #Checking if bam files have been already analysed. If so don't analyse again. CHANGE TO CHECK THE LAST STEP OF THE PIPELINE!

    # for bam_file in bam_files:
    #     print(f"HEEY! {bam_file}")
    #     Bam_class = Bam_file(bam_file)
    #     output_directory = f"{directory}vcf"
    #     abs_output_file = f"{output_directory}/{Bam_class.get_full_ID()}.vcf.gz"
    #     print (abs_output_file)

    #     if os.path.exists(abs_output_file):
    #         #bam_files.remove(bam_file)
    #         print(f"file {bam_file} won't be analyzed")

    #     else:
    #         bams_to_analyse.append(bam_file) 
    #         print(f"FILE {bam_file} will be analyzed by the pipeline")

    #     pipeline_process(bam_file)
    # print(f"The following bam files will be analysed:\n {bams_to_analyse}")
    # parerallize the pipeline process to process one file in each CPU
    #the process pool is closed automatically when it finishes
    #with Pool() as pool:
    #	pool.map(pipeline_process, (bm_file for bm_file in bam_files))
    #CHANGE BAM_FILES FOR BAMS_TO_ANALYSE WHEN WORKING


        


        