import vcf
import os
import pandas as pd
import os
import subprocess
from multiprocessing import Pool
import argparse
import logging
import pysam
import subprocess
from subprocess import PIPE
import re
import pandas as pd
import io
from bs4 import BeautifulSoup
import json
import vcf


directory = "/home/ocanal/Desktop/trombosi_fastq/"

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
        #print(pd_df)
        return(pd_df)


    def extracting_columns (self, pandas_df):
        """
        From the pandas dataframe create an excel with the columns required
        """
        # the ones that I have to add once I have applied the pluggins:
        # "SpliceAI_pred", "CADD_PHRED", "CADD_RAW", "REVEL",
        # selecting the columns of interest to extract from the pandas dataframe 
        columns = ["RB", "Uploaded_variation", "QUAL", "%_alt_reads", "count_total_reads","count_alternative_reads", "Consequence", "Gene", "SWISSPROT","cDNA_position", "CDS_position", "Protein_position",
                    "HGVSc", "HGVSp", "HGVSg", "MANE_SELECT", "MANE_PLUS_CLINICAL", "CLIN_SIG", "Existing_variation" ,"UNIPROT_ISOFORM", "MAX_AF", "PHENO",
                    "PUBMED", "Feature_type", "MES-SWA_acceptor_alt", "MES-SWA_acceptor_diff", "MES-SWA_acceptor_ref", "MES-SWA_acceptor_ref_comp",
                    "MES-SWA_donor_alt", "MES-SWA_donor_diff", "MES-SWA_donor_ref", "MES-SWA_donor_ref_comp"]
        # creating a pd dataframe with the columns of interest
        # print(pandas_df[columns])
        return(pandas_df[columns])

    def df_to_excel (self, final_pd):
        """Create an excel in /directory/excels/ with the data extracted from vep with the selecting columns from extracting_columns function"""

        # create the excel directory if it does not exist
        dir_exists = os.path.exists(self.excel_path)
        if dir_exists:
            pass  

        else:
            os.mkdir(self.excel_path)

        # Extracting the pandas dataframe to an excel

        final_pd.to_excel(self.excel_file_path)
        logging.info(f"Results are written to excel {self.excel_file_path} successfully!")

    def get_coverage(self):

        abs_output_dir = f"{directory}/mosdepth/"        
        
        colnames = ["chrom", "start", "end", "region", "1X", "10X", "20X", "30X", "100X", "200X"]
        mosdepth_pddf = pd.read_csv(f"{abs_output_dir}RB31799_mosdepth.thresholds.bed.gz", compression="gzip", names = colnames, sep="\t", comment="#")
        print(mosdepth_pddf)
        return(mosdepth_pddf)



        
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
        vcf_abs_path = "/home/ocanal/Desktop/trombosi_fastq/vcf/RB31799.vcf.gz"

        #column names to give to the pandas dataframe (using the same as the vcf file)
        colnames = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "20"]

        return(pd.read_csv(vcf_abs_path, compression="gzip", names=colnames, sep="\t", comment="#"))



    def list_uploaded_variation (self, pd_df):
        """From a pandas dataframe, it obtains a list of the Uploaded_variation values"""
        chr_position_list = []
        
        for ind in pd_df.index:
            chr_position_list.append(pd_df["Uploaded_variation"][ind])
        
        return(chr_position_list)

    def split_uploaded_variation (self, chr_position_list):
        """From a list of uploaded_variation values, we create a list of list containing: [chr1, position, alleles] for each variant"""
        chr_position = []
        chr_position_listoflist = []
        for var in chr_position_list:
            items = var.split("_")
            chr_position = [items[0], int(items[1])]
            
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
    
    def add_df_qual( self, pd_df, qual):

        pd_df["vcf_QUAL"] = qual
        return(pd_df)

    def add_df_coverage(self, pddf, coverage_variant):
        """
        add the coverage of the SNP to the final dataframe"""
        pddf["coverage"] = coverage_variant
        return(pddf)

    def call_rate_perc(self, mosdepth_pddf, qual_dict):
        """determines the % of regions specified in the trombosi bed file are covered 
        by a certain number of reads (1 10 20 30 100 200)"""
        
        mosdepth_pddf = mosdepth_pddf.reset_index() #make sure indexes pair with number of rows

        #setting all the values of coverage to 0
        coverage_dic = {}
        x1_region_call = x10_region_call = x20_region_call = x30_region_call = x100_region_call = x200_region_call = 0
        x1_exon_lost = x10_exon_lost = x20_exon_lost = x30_exon_lost = x100_exon_lost = x200_exon_lost = 0
    
       
        for index, row in mosdepth_pddf.iterrows():

            start = int(row["start"])
            end = int(row["end"])
            length = end - start
            x1_bases = int(row["1X"])
            x10_bases = int(row["10X"])
            x20_bases = int(row["20X"])
            x30_bases = int(row["30X"])
            x100_bases = int(row["100X"])
            x200_bases = int(row["200X"])

            if x200_bases >= length:
                x200_region_call += 1

            elif x100_bases >= length:
                x100_region_call += 1

            elif x30_bases >= length:
                x30_region_call += 1
            
            elif x20_bases >= length:
                x20_region_call += 1
            
            elif x10_bases >= length:
                x10_region_call += 1

            elif x1_bases >= length:
                x1_region_call += 1

            else:
                #logging.critical(f"Be aware! the region {str(row["region"])} has not been covered for the sample {self.getID}")
                pass

            if x1_bases < length:
                x1_exon_lost += 1

            elif x10_bases < length:
                x10_exon_lost += 1

            elif x20_bases < length:
                x20_exon_lost += 1

            elif x30_bases < length:
                x30_exon_lost += 1

            elif x100_bases < length:
                x100_exon_lost += 1

            elif x200_bases < length:
                x200_exon_lost += 1
            

            
        #number of regions is equal to the number of dataframe rows
        number_exons = mosdepth_pddf.shape[0]

        #a lower call rate will also be covered if a higher call rate has been called
        #e.g. the number of regions covered by 200x will also be covered by the others call rates
        x100_region_call += x200_region_call
        x30_region_call += x100_region_call
        x20_region_call += x30_region_call
        x10_region_call += x20_region_call
        x1_region_call += x10_region_call 

        #Otherwise, in exon lost is the contrary, if a exon is not covered at 1x,
        #it means that it won't be covered at higher coverages
        x10_exon_lost += x1_exon_lost
        x20_exon_lost += x10_exon_lost
        x30_exon_lost += x20_exon_lost
        x100_exon_lost += x30_exon_lost
        x200_exon_lost += x100_exon_lost

        qual_dict["x1_exon_lost"] = x1_exon_lost
        qual_dict["x10_exon_lost"] = x10_exon_lost
        qual_dict["x20_exon_lost"] = x20_exon_lost
        qual_dict["x30_exon_lost"] = x30_exon_lost
        qual_dict["x100_exon_lost"] = x100_exon_lost
        qual_dict["x200_exon_lost"] = x200_exon_lost
        

        #converting call rates to % and adding to qual dict
        qual_dict["x200_call_percent"] = (x200_region_call / number_exons) * 100
        qual_dict["x100_call_percent"] = (x100_region_call / number_exons) * 100
        qual_dict["x30_call_percent"] = (x30_region_call / number_exons) * 100
        qual_dict["x20_call_percent"] = (x20_region_call / number_exons) * 100
        qual_dict["x10_call_percent"] = (x10_region_call / number_exons) * 100
        qual_dict["x1_call_percent"] = (x1_region_call / number_exons) * 100

        qual_dict["numb_exons"] = number_exons

        qual_dict["RB"] = self.get_full_ID()

        return(qual_dict)
    

    def get_mean_coverage (self, qual_dict):
        full_path = f"{directory}mosdepth/"
        file = f"{self.get_full_ID()}_mosdepth.regions.bed.gz"
        full_file = full_path + file

        colnames = ["chr", "start", "end", "region", "coverage"]

        mosdepth_pddf = pd.read_csv(f"{full_file}", compression="gzip", names = colnames, sep="\t", comment="#")
        
        #number of exons is equal to the number of rows of the dataframe:
        number_exons = mosdepth_pddf.shape[0]

        #sum all the coverages
        sum_coverage = mosdepth_pddf["coverage"].sum()

        mean_coverage = sum_coverage/number_exons

        qual_dict["mean_coverage"] = mean_coverage

        for key,value in qual_dict.items():
            print (key, value)

        return(qual_dict)

    def extract_variant_coverage (self, mosdepth_pddf, chrom=str, pos= int):
        """
        It takes the rare variant from the mosdepth output (comparing that the SNP position is
        between the position interval of mosdepth), it takes the 
        minimum x (x10, x20, x30... meaning that all bases are paired at least by this number of reads)
        and returns the maximum xnumber that all bases are covered by x number of reads.
        """
        print(chrom, pos)
        try:
            x1 = (mosdepth_pddf.loc[(mosdepth_pddf["chrom"] == chrom) & (mosdepth_pddf["start"]<= pos) & (mosdepth_pddf["end"]>pos), "1X"]).astype(int).item()
            x10 = (mosdepth_pddf.loc[(mosdepth_pddf["chrom"] == chrom) & (mosdepth_pddf["start"]<= pos) & (mosdepth_pddf["end"]>pos), "10X"]).astype(int).item()
            x20 = (mosdepth_pddf.loc[(mosdepth_pddf["chrom"] == chrom) & (mosdepth_pddf["start"]<= pos) & (mosdepth_pddf["end"]>pos), "20X"]).astype(int).item()
            x30 = (mosdepth_pddf.loc[(mosdepth_pddf["chrom"] == chrom) & (mosdepth_pddf["start"]<= pos) & (mosdepth_pddf["end"]>pos), "30X"]).astype(int).item()
            x100 = (mosdepth_pddf.loc[(mosdepth_pddf["chrom"] == chrom) & (mosdepth_pddf["start"]<= pos) & (mosdepth_pddf["end"]>pos), "100X"]).astype(int).item()
            x200 = (mosdepth_pddf.loc[(mosdepth_pddf["chrom"] == chrom) & (mosdepth_pddf["start"]<= pos) & (mosdepth_pddf["end"]>pos), "200X"]).astype(int).item()
            start = (mosdepth_pddf.loc[(mosdepth_pddf["chrom"] == chrom) & (mosdepth_pddf["start"]<= pos) & (mosdepth_pddf["end"]>pos), "start"]).astype(int).item()
            end = (mosdepth_pddf.loc[(mosdepth_pddf["chrom"] == chrom) & (mosdepth_pddf["start"]<= pos) & (mosdepth_pddf["end"]>pos), "end"]).astype(int).item()
            bases = end - start
        except:
            x1 = x10 = x20 = x30 = x100 = x200 = 0
            start = 0
            end = 1
            bases = 1
            
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
            return ("not_into_bed_interval")

    def qual_excel (qual_dict):

        excel_filename= f"{directory}Quality_excel"
        qual_pddf = pd.DataFrame(qual_dict, index = [0])
        qual_pddf.to_excel(excel_filename, index = False, header=False,mode= 'a')

        return(0)

    def extract_ref_alt_reads ( self, values_list):
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

    def reads_counts_to_df (self, pd_df, alt_reads_list, total_reads_list, percent_alt_list):
        
        pd_df["count_alternative_reads"] = alt_reads_list
        pd_df["count_total_reads"] = total_reads_list
        pd_df["%_alt_reads"] = percent_alt_list
            
        return(pd_df)
        














#print(pd_df)

#qual = pd_df["QUAL"].where((pd_df["POS"] == "1554362") & (pd_df["CHROM"] == "chr1"))

#Create a vcf class
vcf_file = f"/home/ocanal/Desktop/trombosi_fastq/RB31799.vcf.gz"

Vcf_file_Class = Vcf_Class(vcf_file)

pandasdt = Vcf_file_Class.vep_parser()



#We add an identifier column where it will be stored the RB of the sample
pd_dt_RB = Vcf_file_Class.add_RB_column(pandasdt)
#Filter the rare variants (population frequency < 0.1%) We take a frequency of 0.001 considering that vep gives us the 
pandas_filtered = Vcf_file_Class.filtering_df_by_AlleFreq(pd_dt_RB, 0.001)
#print(pandas_filtered)
vcf_pd=Vcf_file_Class.vcf_to_pd()
uploaded_variation = Vcf_file_Class.list_uploaded_variation(pandas_filtered)
uploaded_variation_list = Vcf_file_Class.split_uploaded_variation(uploaded_variation)
mosdepth_pddf = Vcf_file_Class.get_coverage()
qual_dict = {}
qual_dict = Vcf_file_Class.call_rate_perc(mosdepth_pddf, qual_dict)
Vcf_file_Class.get_mean_coverage(qual_dict)
#print(uploaded_variation)
rare_variant_qual = []
rare_variant_info = []
coverage_variant= []
#print(uploaded_variation_list)
for uploaded_variation in uploaded_variation_list:
    rare_variant_qual.append(float(Vcf_file_Class.extract_vcf_column(vcf_pd, "QUAL" , chrom = uploaded_variation[0], pos = uploaded_variation[1])))
    rare_variant_info.append(list(Vcf_file_Class.extract_vcf_column(vcf_pd, "20" , chrom = uploaded_variation[0], pos = uploaded_variation[1])))
    coverage_variant.append(str(Vcf_file_Class.extract_variant_coverage(mosdepth_pddf, chrom=uploaded_variation[0], pos=uploaded_variation[1])))
print (f"coverage variant: \n {coverage_variant}\n final coverage variant ")
#print (f"info list = {rare_variant_info}")
#print(f"dfa {str(rare_variant_info[0][0])}")
#print (f"qual list = {rare_variant_qual}")
alt_reads_list, total_reads_list, percent_alt_list = Vcf_file_Class.extract_ref_alt_reads(rare_variant_info)

qual_df = Vcf_file_Class.add_df_qual(pandas_filtered, rare_variant_qual)



reads_df =Vcf_file_Class.reads_counts_to_df(qual_df, alt_reads_list, total_reads_list, percent_alt_list)

#Extracting the columns that we want in the final excel
reduced_pandas = Vcf_file_Class.extracting_columns(pandas_filtered)
Vcf_file_Class.qual_excel(qual_dict)
#Creating the resulting excel
Vcf_file_Class.df_to_excel(reduced_pandas)

