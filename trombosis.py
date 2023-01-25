import os
import subprocess
from classes import *
from multiprocessing import Pool
import argparse
import logging

parser = argparse.ArgumentParser(description = "Pipeline to automatize the process from fastq files to vcf annotations")
group = parser.add_mutually_exclusive_group(required=True)
group.add_argument ("-d", "--input_directory", action = "store" ,help="input directory where bam files are stored")
group.add_argument ("-r", "--reference", required = False, default = "/home/ocanal/Desktop/reference_genomes/trombosi/hg38/Homo_sapiens_assembly38.fasta", action = "store", help="input directory where the fq reference files with its index and dict files are stored")
args = parser.parse_args()
directory = args.input_directory
reference = args.reference

root_logger = logging.getLogger()
root_logger.setLevel(logging.DEBUG)
handler = logging.FileHandler(f"{directory}/log.log", "w", "utf-8")
root_logger.addHandler(handler)

if not directory.endswith("/"):
    directory += "/"

def uncompress_files (directory):
    """ From a given directory it decompress all .tgz files"""
    cmpr_dirs = []
    dirs = os.listdir(directory)
    for direct in dirs:
        if direct.endswith(".tgz"):
            cmpr_dirs.append(direct)
            uncmpr_cmd = f"tar -zxvf {direct}"
            subprocess.run(uncmpr_cmd, shell=True)



def list_fastq_files (directory):
    """From a given directory it lists you all the bam files"""
    fastq_files=[]
    files = [f for f in os.listdir(directory) if os.path.isfile(os.path.join(directory,f))]
    for file in files:
        if file.endswith(".fastq"):
            fastq_files.append(file)
    if not fastq_files:
        logging.critical = f"No fastq files found in the input directory: {directory}"
        

    return(fastq_files)

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

        #obtain the RB of the file
        ID = re.search(regex_exp, fq_file).group(0)

        # if the RB have not already been added to the dictionary
        if ID not in fq_dictionary:

            #create a dictionary with RB = 1
            fq_dictionary[ID] = 1

        #If the RB already exists in the dictionary
        else:

            #add 1 to the current dictionary value
            fq_dictionary[ID] += 1
            if fq_dictionary[ID] > 2:
                logging.critical (f"Be aware, you have more than 2 files with the following RB {ID}")

    paired_reads_IDs = [key for key, value in fq_dictionary.items() if value == 2]
    single_reads_IDs = [key for key, value in fq_dictionary.items() if value == 1]



    return (single_reads_IDs, paired_reads_IDs)
    
def single_reads_IDs_to_fqfile (fq_file_list, ID):

    match = [fq_file for fq_file in fq_file_list if re.search(ID, fq_file)]
    se_fq_file = str(match[0])
    return(se_fq_file)

def paired_reads_IDs_to_fqfile(fq_file_list, ID):

    match = [fq_file for fq_file in fq_file_list if re.search(ID, fq_file)]
    pe_fq_files = list(match)


    return(pe_fq_files)



    
def list_readybam_files (directory):
    """From a given directory it lists you all the bam files"""
    ready_bam_files=[]
    ready_bam_dir = f"{directory}/grouped_reads"
    files = [f for f in os.listdir(ready_bam_dir) if os.path.isfile(os.path.join(ready_bam_dir,f))]
    for file in files:
        if file.endswith(".bam"):
            ready_bam_files.append(file)

    return(ready_bam_files)

def single_reads_fastQ_to_bam(fq_file):
    """
    Different steps that are run by the pipeline.
    """

    logging.info(f"Starting to analyse the file {fq_file}")

    #Create a FastQ class 
    FastQ_class = FASTQ_file(fq_file)

    #Perform a quality check and generate a HTML report

    FastQ_class.quality_check_fastq()

    #Use Trimmomatic to remove low-quality base calls and discard short reads (<30)
    trimmed_fqfile = FastQ_class.trim_fastq()

    #Align the trimmed file to the reference genome Hg38
    sam_file = FastQ_class.align_to_reference_hg38(trimmed_fqfile)

    #Convert the sam file to a sorted bam file along with its index
    sorted_bam = FastQ_class.sam_to_sorted_bam(sam_file)
    
    logging.info(f"Sorted Bam file have been created for the input file {fq_file}")

    return (sorted_bam)

def paired_reads_fastq_to_bam (fq_file_1, fq_file_2):
    """giving a paired reads fastq it will create a bam"""
    
    logging.info(f"Starting to analyse the paired reads files {fq_file_1, fq_file_2}")

    #create a paired FASTQ class
    Paired_FastQ_class = Paired_FASTQ_file(fq_file_1, fq_file_2)

    #Perform a quality check and generate a HTML report
    Paired_FastQ_class.quality_check_fastq()

    #Use Trimmomatic to remove low quality base calls and discard short reads (<30).
    #output is trimmed_fq that is a tuple where:
    #trimmed_fq[0] = forward paired trimmed fastQ file
    #trimmed_fq[1] = forward unpaired trimmed fastQ file
    #trimmed_fq[2] = reverse paired trimmed fastQ file
    #trimmed_fq[3] = reverse unpaired trimmed fastQ file
    trimmed_fq = Paired_FastQ_class.trim_fastq()

    #align the paired reads to the reference genome
    sam_file = Paired_FastQ_class.align_to_reference_hg38(trimmed_fq[0], trimmed_fq[2])
    
    #Convert the sam file to a sorted bam file along with its index
    sorted_bam = Paired_FastQ_class.sam_to_sorted_bam(sam_file)
    
    logging.info(f"Sorted Bam file have been created for the input file {fq_file_1, fq_file_2}")

    return (sorted_bam)

def bam_to_final_excel (sorted_bam):
    """
    Function that applies the different functions created in classes.py to obtain an excel with the 
    rare variants from a sorted bam file
    """

    logging.info(f"Variant annotation have started for file {sorted_bam}")

    #Bam class is created
    Bam_class = Bam_file(sorted_bam)
    #print(Bam_class.bam_file)

    #The integrity of the BAM file is checked
    #Bam_class.check_integrity()

    #dictionary of the reads stats of the BAM file is created
    #Bam_class.reads_stats()

    mosdepth_pddf=Bam_class.get_coverage()
    #AddOrReplaceReadGroups is executed and output files are stored in /grouped_reads/
    #The Bam file is ready to be proceeded by haplotype caller
    Bam_gr = Bam_class.group_reads()

    #mark duplicates
    marked_bam = Bam_class.mark_duplicates()

    #UNCOMMENT THE FOLLOWING LINE!
    ready_bam= Bam_class.quality_recalibration()


    #We take the ready_bam as the working file from now on and create a class for this file  
    ready_Bam_class = Bam_file(ready_bam)
    print(ready_Bam_class.get_ID())
    print(re)
    #Index the bam file (.bai) to run haplotype caller 
    ready_Bam_class.index_bam()

    #UNCOMMENT THE FOLLOWING LINE
    #ready_Bam_class.haplotype_caller()

    #Create a vcf class
    vcf_file = f"{ready_Bam_class.get_full_ID()}.vcf.gz"
    
    Vcf_file_Class = Vcf_Class(vcf_file)

    #run variant annotations 

    Vcf_file_Class.variant_annotation()

    #Create a pandas dataframe (variants in rows and parameters as columns)
    pandasdt = Vcf_file_Class.vep_parser()

    #We add an identifier column where it will be stored the RB of the sample
    pd_dt_RB = Vcf_file_Class.add_RB_column(pandasdt)

    #Filter the rare variants (population frequency < 0.1%) We take a frequency of 0.001 considering that vep gives us the 
    pandas_filtered = Vcf_file_Class.filtering_df_by_AlleFreq(pd_dt_RB, 0.001)

    #creting a dataframe to store all the information from the vcf file
    vcf_pd=Vcf_file_Class.vcf_to_pd()

    #From a pandas dataframe, it obtains a list of the Uploaded_variation values
    uploaded_variation = Vcf_file_Class.list_uploaded_variation(pandas_filtered)
    #From a list of uploaded_variation values, we create a list of list containing: [chr1, position, alleles] for each variant
    uploaded_variation_list = Vcf_file_Class.split_uploaded_variation(uploaded_variation)
    #List where the quality of rare variants are stored
    rare_variant_qual = []
    rare_variant_info = []
    coverage_variant = []
    #For every rare variant, we obtain the quality of the variant and append to rare_variant_list 
    
    for uploaded_variation in uploaded_variation_list:
        rare_variant_qual.append(float(Vcf_file_Class.extract_vcf_column(vcf_pd, "QUAL" , chrom = uploaded_variation[0], pos = uploaded_variation[1])))
        rare_variant_info.append(list(Vcf_file_Class.extract_vcf_column(vcf_pd, "20" , chrom = uploaded_variation[0], pos = uploaded_variation[1])))
        #coverage_variant.append(list(Vcf_file_Class.extract_variant_coverage(mosdepth_pddf, chrom=uploaded_variation[0], pos=uploaded_variation[1])))
    alt_reads_list, total_reads_list, percent_alt_list = Vcf_file_Class.extract_ref_alt_reads(rare_variant_info)
    
    #Add the qual parameter to the vep dataframe
    qual_df = Vcf_file_Class.add_df_qual(pandas_filtered, rare_variant_qual)
    reads_df =Vcf_file_Class.reads_counts_to_df(qual_df, alt_reads_list, total_reads_list, percent_alt_list)
    #coverage_df = Vcf_file_Class.add_coverage_df(reads_df,coverage_variant)
    #Extracting the columns that we want in the final excel
    reduced_pandas = Vcf_file_Class.extracting_columns(coverage_df)
    #print(coverage_df)

    #Creating the resulting excel
    Vcf_file_Class.df_to_excel(reduced_pandas)

    logging.info(f"The annotations have been performed correctly and the rare variants have been annotated to the excel")

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

    #create a joined_excels folder, if it haven't already exists

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
    logging.info(f"Joined excel have been created successfully")


if __name__ =="__main__":

    #Step 1: decompress files and create a list of the existing BAM files
    uncompress_files(directory)

    #Step 2 detect the fastq files 
    fq_files = list_fastq_files(directory)

    #and seperate in 2 different lists the single end reads and the paired end reads
    ID_fq_single_reads_list, ID_fq_paired_reads_list = match_paired_fastq(fq_files)

    #if there are some single end reads fastq files
    if ID_fq_single_reads_list:
        
        #analyse each fastq file
        for ID_single_reads in ID_fq_single_reads_list:
            logging.info (f"The fastQ file: {fq_file} will be analysed by the pipeline")
            
            #match the single reads IDs with its files
            fq_file = single_reads_IDs_to_fqfile(fq_file, ID_fq_single_reads_list)
            
            #perform the analysis
            sorted_bam = single_reads_fastQ_to_bam(fq_file)
            bam_to_final_excel(sorted_bam)

    
    #if there are some paired end reads fastq files:
    if ID_fq_paired_reads_list:

        #analyse each paired end reads
        for ID_paired_reads in ID_fq_paired_reads_list:

            #match the paired reads files associated with the ID
            paired_list = paired_reads_IDs_to_fqfile(fq_files, ID_paired_reads)
            print(paired_list)
            sorted_bam = paired_reads_fastq_to_bam(paired_list[0], paired_list[1])
            bam_to_final_excel(sorted_bam)
    
    join_excels(directory)
    
    
    
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


        


        