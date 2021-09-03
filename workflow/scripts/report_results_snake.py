#!/usr/bin/env python3
# -*- coding: utf-8 -*-
""" report_results.py

    Script to create a report of the results for every sample,
    from the three programs used in this pipeline: Aldy, Optitype 
    and GATK. It also states if every necessary position to properly
    genotype has sufficient coverage.
    
    
    ****************************************************************
    *   Julia Corell (@jucosie)                                    *                                                   
    *   Version 1.0                                                *                                      
    ****************************************************************
"""
import os
import sys
import glob
import pandas as pd
import vcf
import numpy as np
import re

def get_hla_genotype(path_to_file, sample, DF_results):
    """
    Extracts the HLA results from a given file and adds them to the DF
    of results.

    Parameters
    ----------
    path_to_file : String
        Absolute path to HLA results file.
    sample : String
        The name of the sample.
    DF_results: Pandas Dataframe
        Dataframe where is stored every result of a given sample.

    Returns
    -------
    None.
    """
    try:
        hla_df = pd.read_csv(path_to_file, sep="\t")
        for letter in ["A", "B", "C"]:
            hla_res = list(hla_df[letter + "1"]) + list(hla_df[letter + "2"])
            DF_results.loc["HLA-" + letter, "Result_1"] = hla_res[0]
            DF_results.loc["HLA-" + letter,"Result_2"] = hla_res[1]
            print("The results of HLA-{} have been extracted from {}".format(letter,\
                (path_to_file)))
    except OSError:
        print ("HLA file of sample {} not found.".format(sample))

def assign_genotype(genotype, gene_sym, rs_id, DF_results):
    """
    Translates the diplotype results of a given VCF instance and adds
    them to the DF of results.

    Parameters
    ----------
    genotype : String
        Encoded allele values.
    gene_sym : String
        The gene symbol of the gene.
    rs_id : String
        If available, the rsID of the instance.    
    DF_results: Pandas Dataframe
        Dataframe where is stored every result of a given sample.

    Returns
    -------
    None.
    """
    if genotype in ["0/0", "0|0"]:
        DF_results.loc[gene_sym, "Result_1"] = "REF"
        DF_results.loc[gene_sym, "Result_2"] = "REF"
    if genotype in ["0/1", "1/0", "0|1", "1|0"]:
        DF_results.loc[gene_sym, "Result_1"] = rs_id
        DF_results.loc[gene_sym, "Result_2"] = "REF"
    if genotype in ["1/1", "1|1"]:
        DF_results.loc[gene_sym, "Result_1"] = rs_id
        DF_results.loc[gene_sym, "Result_2"] = rs_id
    if genotype in ["./.", ".|."]:
        DF_results.loc[gene_sym, "Result_1"] = "NOT_CALLED"
        DF_results.loc[gene_sym, "Result_2"] = "NOT_CALLED"

def assign_genotype_2(genotype, rs_id, DF_results_mult):
    """
    Translates the diplotype results of a given VCF instance and adds
    them to the DF of multiple SNP results.

    Parameters
    ----------
    genotype : String
        Encoded allele values.
    rs_id : String
        If available, the rsID of the instance.    
    DF_results_mult: Pandas Dataframe
        Dataframe where is stored every result of those genes
        which can have more than one variant, of a given sample.

    Returns
    -------
    None.
    """
    i = DF_results_mult[DF_results_mult["rsid"] == rs_id].index

    if genotype in ["0/0", "0|0"]:
        DF_results_mult.loc[i, "Result_1"] = "REF"
        DF_results_mult.loc[i, "Result_2"] = "REF"
    if genotype in ["0/1", "1/0", "0|1", "1|0"]:
        DF_results_mult.loc[i, "Result_1"] = rs_id
        DF_results_mult.loc[i, "Result_2"] = "REF"
    if genotype in ["1/1", "1|1"]:
        DF_results_mult.loc[i, "Result_1"] = rs_id
        DF_results_mult.loc[i, "Result_2"] = rs_id
    if genotype in ["./.", ".|."]:
        DF_results_mult.loc[i, "Result_1"] = "NOT_CALLED"
        DF_results_mult.loc[i, "Result_2"] = "NOT_CALLED"

def extract_SNP_from_merged(path_to_vcf , sample, variants_table, DF_results):
    """
    Extracts the results of the Variant Calling from the VCF of the filtered 
    SNVs and INDELs, and adds them to the DF of results.

    Parameters
    ----------
    path_to_vcf: String
        Absolute path to the filtered VCF file.
    sample: String
        The name of the sample.
    variants_table: Pandas Dataframe
        Auxiliary table with information about every position queried
        in the VC.
    DF_results: Pandas Dataframe
        Dataframe where is stored every result of a given sample.

    Returns
    -------
    None.
    """
    #First we want to get every variant which has been called and has passed the filters
    try:
        merged_vcf = vcf.Reader(open(path_to_vcf, "r"))
        for record in merged_vcf:
            if record.FILTER == []:
                if variants_table["end"].astype(str).str.contains(str(record.POS)).any(): #For the SNPs
                    gene = variants_table.loc[(variants_table.end == record.POS) & \
                        (variants_table.chr == record.CHROM) , "gene"].item()
                    rsid = variants_table.loc[(variants_table.end == record.POS) & \
                        (variants_table.chr == record.CHROM), "rsid"].item()
                    if gene in DF_results.index:
                        print("The results of {} have been extracted from {}".format(gene,\
                            (path_to_vcf)))
                        assign_genotype(record.genotype(sample)["GT"], gene, rsid, DF_results)
                else:
                    for n in range(len(variants_table.index)): ##For those positions defined by an interval, e.g. rs11322783
                            if variants_table.loc[n,"start"] <= record.POS < variants_table.loc[n, "end"]:
                                gene = variants_table.loc[n, "gene"]
                                rsid = variants_table.loc[n, "rsid"]
                                if gene in DF_results.index:
                                    print("The results of {} have been extracted from {}".format(gene,\
                                        (path_to_vcf)))
                                    assign_genotype(record.genotype(sample)["GT"], gene, rsid, DF_results)
            else:
                print("{} was detected in sample {} but didn't pass the filters".format(record.ID, sample))
                gene = variants_table.loc[(variants_table.rsid == record.ID) , "gene"].item()
                DF_results.loc[gene, "Result_1"] = "FAILED_QA"
                DF_results.loc[gene, "Result_2"] = "FAILED_QA"
    except OSError:
        print ("File merged.vcf not accesible.")


def check_missing_genes (DF_results, variants_table):
    """
    Checks which genes didn't appear in the filtered VCF (because there
    wasn't any call). 

    Parameters
    ----------
    variants_table: Pandas Dataframe
        Auxiliary table with information about every position queried
        in the VC.
    DF_results: Pandas Dataframe
        Dataframe where is stored every result of a given sample.

    Returns
    -------
    list_of_genes: List
        List of genes which have to be checked in the unfiltered VCF
    """
    genes_missing = DF_results[DF_results["Result_1"].isna()]
    list_of_genes = []
    for g in genes_missing.index:
        if g in list(variants_table["gene"]):
            list_of_genes.append(g)
    return list_of_genes

def extract_SNP_from_output(path_to_vcf, sample, variants_table, missing_genes, DF_results):
    """
    Extracts the results of the Variant Calling from the VCF which stores
    every position queried, and adds them to the DF of results.

    Parameters
    ----------
    path_to_vcf: String
        Absolute path to the unfiltered VCF file.
    sample: String
        The name of the sample.
    variants_table: Pandas Dataframe
        Auxiliary table with information about every position queried
        in the VC.
    missing_genes: List
        Output of the function check_missing_genes
    DF_results: Pandas Dataframe
        Dataframe where is stored every result of a given sample.

    Returns
    -------
    None.
    """
    try:
        output_vcf = vcf.Reader(open(path_to_vcf, "r"))
        for rec in output_vcf:
            if variants_table["end"].astype(str).str.contains(str(rec.POS)).any():
                gene_symbol = variants_table.loc[(variants_table.end == rec.POS) & \
                    (variants_table.chr == rec.CHROM), "gene"].item() 
                rsID = variants_table.loc[(variants_table.end == rec.POS) & \
                    (variants_table.chr == rec.CHROM), "rsid"].item() #Get the rsID of the  position from the variants table
                if gene_symbol in missing_genes:
                    print("The results of {} variant {}, have been extracted from {}".format(gene_symbol, \
                                rsID, (path_to_vcf)))
                    assign_genotype(rec.genotype(sample)["GT"],gene_symbol, rsID, DF_results)
            else:
                    for n in range(len(variants_table.index)): ##For those positions defined by an interval, e.g. rs11322783
                            if variants_table.loc[n,"start"] <= rec.POS < variants_table.loc[n, "end"]:
                                gene = variants_table.loc[n, "gene"]
                                rsID = variants_table.loc[n, "rsid"]
                                if gene in DF_results.index:
                                    print("The results of {} have been extracted from {}".format(gene_symbol,\
                                        (path_to_vcf)))
                                    assign_genotype(rec.genotype(sample)["GT"], gene, rsID, DF_results)
    except OSError:
        print ("File output.vcf is not accesible")
        

def check_filter(path_to_merged, record_ID):
    """
    Checks the status of the FILTER flag of a given VCF record.

    Parameters
    ----------
    path_to_merged: String
        Absolute path to the filtered VCF file.
    DF_results: Pandas Dataframe
        Dataframe where is stored every result of a given sample.

    Returns
    -------
    list_of_genes: List
        List of genes which have to be checked in the unfiltered VCF
    """
    try:
        merged_vcf = vcf.Reader(open(path_to_merged, "r"))
        for instance in merged_vcf:
            if instance.ID == record_ID:
                if instance.FILTER == []:
                    return True
                else:
                    return False
    except OSError:
        print ("File merged.vcf not accesible.")


def extract_multipleSNP_from_output(path_to_vcf, path_to_merged, sample, variants_bed_sub, DF_results_mult):
    """
    Extracts the results of the Variant Calling from the VCF which stores
    every position queried, and adds them to the DF of multiple SNP results.
    If the possition has a variant, it checks its FILTER status in the filtered
    VCF.

    Parameters
    ----------
    path_to_vcf: String
        Absolute path to the unfiltered VCF file.
    path_to_merged: String
        Absolute path to the filtered VCF file.
    sample: String
        The name of the sample.
    variants_bed_sub: Pandas Dataframe
        Subset of the variants_table just with the genes for which we have
        called more than one position
    DF_results_mult: Pandas Dataframe
        Dataframe where is stored every result of those genes
        which can have more than one variant, of a given sample.

    Returns
    -------
    None.
    """
    try:
        output_vcf = vcf.Reader(open(path_to_vcf, "r"))
        for rec in output_vcf:
            if rec.ID:
                if variants_bed_sub["rsid"].str.contains(rec.ID).any():
                    if check_filter(path_to_merged, rec.ID):
                        gene_symbol = variants_bed_sub.loc[(variants_bed_sub.rsid == rec.ID), \
                            "gene"].item()
                        print("The results of {} variant {}, have been extracted from VC_calls/merged.vcf".format(gene_symbol, \
                            rec.ID, ))
                        assign_genotype_2(rec.genotype(sample)["GT"],rec.ID,DF_results_mult)
                    else:
                        print("{} was detected but didn't pass the filters".format(rec.ID))
                        ind = DF_results_mult[DF_results_mult["rsid"] == rec.ID].index
                        DF_results_mult.loc[ind, "Result_1"] = "FAILED_QA"
                        DF_results_mult.loc[ind, "Result_2"] = "FAILED_QA"
            elif variants_bed_sub["end"].astype(str).str.contains(str(rec.POS)).any():
                gene_symbol = variants_bed_sub.loc[(variants_bed_sub.end == rec.POS) & \
                    (variants_bed_sub.chr == rec.CHROM), "gene"].item() 
                rsID = variants_bed_sub.loc[(variants_bed_sub.end == rec.POS) & \
                    (variants_bed_sub.chr == rec.CHROM), "rsid"].item() #Get the rsID of the  position from the variants table
                print("The results of {} variant {}, have been extracted from {}".format(gene_symbol, \
                            rsID, (path_to_vcf)))
                assign_genotype_2(rec.genotype(sample)["GT"],rsID,DF_results_mult)
    except OSError:
        print ("File output.vcf not accesible.")


def extract_aldy_solution(file):
    """
    Function to extract the genotype called by Aldy of a given .aldy file.

    Parameters
    ----------
    file: String
        Absolute path to the aldy file.

    Returns
    -------
    major_1, major_2: String
        The diplotype assigned by Aldy.
    """
    try:
        if os.stat(file).st_size != 0:
            with open(file) as f:
                for line in f:
                    if line.startswith("#Solution 1"):
                        sol = line	                
            chunks = sol.split(":")
            solution_1 = chunks[1].split(",")[0].strip()
            solution_2 = chunks[1].split(",")[1].strip()
            major_1 = find_major_allele(solution_1)
            major_2 = find_major_allele(solution_2)
            return major_1, major_2
        else:
            print ("File {} is not accesible".format(file))
            major_1 = "NOT_CALLED"
            major_2 = "NOT_CALLED"
            return major_1, major_2
    except OSError:
        print ("File {} is not accesible".format(file))
        major_1, major_2 = "NOT_CALLED"
        return major_1, major_2

def find_major_allele(string):
    """
    Extract the major star allele from a string of the Aldy output.

    Parameters
    ----------
    string: String
        Aldy result.

    Returns
    -------
    found: String
        The major star allele.
    """
    try:
        found = re.search(r'^(\*[0-9]+|\*[a-zA-Z]+[0-9]*)', string).group(1)
    except AttributeError:
        found = np.NaN
    return found

def assign_aldy_solutions(list_aldy_files, DF_results, sample):
    """
    Given a list of aldy files, extract the results for each gene and 
    add them to the DF of results. 

    Parameters
    ----------
    list_aldy_files: List
        List of Aldy files of a given sample.
    sample: String
        The name of the sample.
    DF_results: Pandas Dataframe
        Dataframe where is stored every result of a given sample.
    Returns
    -------
    """
    if list_aldy_files:
        for aldy in list_aldy_files:
            file_name = aldy.split("/")[-1]
            gene_name = file_name.split(".")[1] #Be careful if you change the 
                                                #naming pattern of aldy files
            sol_1,sol_2 = extract_aldy_solution(aldy)
            DF_results.loc[gene_name, "Result_1"] = sol_1
            DF_results.loc[gene_name,"Result_2"] = sol_2
            print("The results of {} have been extracted from {}".format(gene_name,\
                (aldy)))
    else:
        print("We cannot show the results of Aldy for {}".format(sample))

def check_coverage(cov_df, gene, threshold = 20):
    """
    Checks the coverage of every position of interest of a given gene. 
    If the coverage is greater than the threshold it returns a "PASS"
    flag, otherwise returns "FAILED". For the HLAs always returns
    "CHECK".

    Parameters
    ----------
    cov_df: Pandas Dataframe
        Dataframe which stores the coverage every position queried.
    gene: String
        Gene symbol.
    threshold: Integer
        Threshold to assign the flag.
    Returns
    -------
    found: String
        The major star allele.
    """
    cov_df_subset = cov_df[cov_df["Gene"] == gene]
    if gene in ["HLA-A","HLA-B", "HLA-C", "HLA-DRB1"]:
    	flag = "CHECK"
    elif all(cov_df_subset.loc[:, "Coverage"] > threshold):
        flag = "PASS"
    else:
        flag = "FAILED"
    return flag

def assign_coverage_df_results(path_to_cov_df, sample, DF_results):
    """
    Adds the coverage flag to every gene in DF results.

    Parameters
    ----------
    path_to_cov_df: String
        Absolute path to the coverage file.
    sample: String
        The name of the sample.
    DF_results: Pandas Dataframe
        Dataframe where is stored every result of a given sample.
    
    Returns
    -------
    None.
    """
    try:
        cov_df = pd.read_csv(path_to_cov_df, sep = "\t", header = None)
        cov_df.columns = ["Chromosome", "Start", "End", "Gene", "Coverage"]
        for gene in DF_results.index:
            flag_cov = check_coverage(cov_df, gene,20) #Change threshold if needed
            DF_results.loc[gene, "Coverage"] = flag_cov
    except OSError:
        print("Mosdepth file for sample {} not accesible".format(sample))

def assign_coverage_df_multipleSNP(path_to_cov_df, variants_bed_sub, sample, DF_results_mult):
    """
    Adds the coverage flag to every gene in DF of multiple results.

    Parameters
    ----------
    path_to_cov_df: String
        Absolute path to the coverage file.
    sample: String
        The name of the sample.
    variants_bed_sub: Pandas Dataframe
        Subset of the variants_table just with the genes for which we have
        called more than one position.
    DF_results_mult: Pandas Dataframe
        Dataframe where is stored every result of those genes
        which can have more than one variant, of a given sample.
    
    Returns
    -------
    None.
    """
    try:
        cov_df = pd.read_csv(path_to_cov_df, sep = "\t", header = None)
        cov_df.columns = ["Chromosome", "Start", "End", "Gene", "Coverage"]
        for id in variants_bed_sub["rsid"]:
            flag_cov = check_coverage_single_pos(cov_df,variants_bed_sub,id,20)
            i = DF_results_mult[DF_results_mult["rsid"] == id].index
            DF_results_mult.loc[i, "Coverage"] = flag_cov
    except OSError:
        print ("Mosdepth file for sample {} not accesible".format(sample))

def check_coverage_single_pos(cov_df, variants_bed_sub, rsid, threshold = 20):
    """
    Checks the coverage of a single position. If the coverage is greater 
    than the threshold it returns a "PASS" flag, otherwise returns "FAILED". 
    For the HLAs always returns "CHECK".

    Parameters
    ----------
    cov_df: Pandas Dataframe
        Dataframe which stores the coverage every position queried.
    variants_bed_sub: Pandas Dataframe
        Subset of the variants_table just with the genes for which we have
        called more than one position.
    rsid : String
        The rsID of the instance.
    threshold: Integer
        Threshold to assign the flag.

    Returns
    -------
    found: String
        The major star allele.
    """
    position = variants_bed_sub.loc[variants_bed_sub.rsid == rsid, "end"].item()
    if position > threshold:
        flag = "PASS"
    else:
        flag = "FAILED"
    return flag

def create_df(gene_sym):
    """
    Creates a Pandas dataframe in which the results of every gene will 
    be stored. The row names will be the list of gene symbols (extracted
    from the coverage file).

    Parameters
    ----------
    gene_sym: List 
        List of gene symbols.

    Returns
    -------
    DF_results: Pandas Dataframe
        Dataframe where is stored every result of a given sample.
    """
    DF_results = pd.DataFrame(data=np.NaN, index = gene_sym, \
        columns=["Result_1","Result_2","Coverage"])
    #Delete those genes which have more than one instance (will be reported 
    # in a separated file)
    DF_results = DF_results.drop(["RYR1","CACNA1S","ITPA"], axis=0)
    return DF_results

def create_df_mult(variants_table):
    """
    Creates the dataframe to store the results of the genes which have
    been queried multiple times. We will have as many rows as positions 
    we are querying per gene. There will be an extra column with the rsid 
    we are querying.

    Parameters
    ----------
    variants_table: Pandas Dataframe
        Auxiliary table with information about every position queried
        in the VC.
        
    Returns
    -------
    DF_results_mult: Pandas Dataframe
        Dataframe where is stored every result of those genes
        which can have more than one variant, of a given sample.
    """

    DF_results_mult = variants_table[(variants_table["gene"] == "RYR1") | \
        (variants_table["gene"] == "CACNA1S") | (variants_table["gene"] == "ITPA")]
    DF_results_mult = DF_results_mult[["gene","rsid"]]
    DF_results_mult = DF_results_mult.reset_index(drop=True)
    DF_results_mult["Result_1"] = np.NaN
    DF_results_mult["Result_2"] = np.NaN
    DF_results_mult["Coverage"] = np.NaN
    return DF_results_mult

def main():
    print(__doc__)
    #Get sample name
    sample = snakemake.params.sample
    
    sys.stderr = open(snakemake.log[0], "w")
    
    optitype_file = snakemake.input[0]
    output_vcf= snakemake.input[1]
    merged_vcf = snakemake.input[2]
    mosdepth_file = snakemake.input[3]
    list_files = snakemake.input[4:]
    
    coverage_df = pd.read_csv(mosdepth_file, sep = "\t", header = None)
    coverage_df.columns = ["Chromosome", "Start", "End", "Gene", "Coverage"]
    gene_symbols = list(coverage_df["Gene"].unique())

    #Read variants file
    variants_bed = pd.read_csv(snakemake.params.variants, sep="\t")

    #We also need a subset of the variants_file just with the genes for which we have
    #called more than one position
    variants_bed_subset = variants_bed[(variants_bed["gene"] == "RYR1") | \
        (variants_bed["gene"] == "CACNA1S") | (variants_bed["gene"] == "ITPA")]
    variants_bed_subset = variants_bed_subset.reset_index(drop=True)

    df_results = create_df(gene_symbols)
    df_results_multi = create_df_mult(variants_bed)
    get_hla_genotype(optitype_file,sample,df_results)
    extract_SNP_from_merged(merged_vcf, sample,variants_bed, df_results)
    miss_genes = check_missing_genes(df_results, variants_bed)
    extract_SNP_from_output(output_vcf, sample, variants_bed, miss_genes, df_results)
    extract_multipleSNP_from_output(output_vcf, merged_vcf, sample, variants_bed_subset, df_results_multi)
    assign_aldy_solutions(list_files, df_results, sample)
    assign_coverage_df_results(mosdepth_file, sample, df_results)
    assign_coverage_df_multipleSNP(mosdepth_file, variants_bed_subset, sample, df_results_multi)
    df_results.to_csv(snakemake.output[0], sep ="\t", index_label = "Gene")
    df_results_multi.to_csv(snakemake.output[1], sep ="\t")
    
    print("\n FINISHED!")
if __name__ == '__main__':
    main()
