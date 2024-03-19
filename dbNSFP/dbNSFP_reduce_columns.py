import pandas as pd

columns = ["#chr", "pos(1-based)", "ref", "alt", "rs_dbSNP", "genename", "HGVSp_ANNOVAR", "HGVSp_snpEff", "Ensembl_geneid", "Ensembl_transcriptid", "Ensembl_proteinid", "SIFT4G_score", "Polyphen2_HDIV_score", "MutationTaster_score", "FATHMM_score", "REVEL_score", "AlphaMissense_score", "CADD_phred", "GERP++_RS", "phyloP100way_vertebrate", "phastCons100way_vertebrate", "1000Gp3_AF", "ExAC_AF", "gnomAD_exomes_AF", "gnomAD_exomes_NFE_AF", "gnomAD_genomes_AF", "gnomAD_genomes_NFE_AF", "clinvar_id", "clinvar_clnsig", "clinvar_trait", "clinvar_OMIM_id", "GTEx_V8_eQTL_gene", "GTEx_V8_eQTL_tissue"]

chromosomes = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y", "M"]

for chromosome in chromosomes:
    df = pd.read_csv("dbNSFP4.7a_variant.chr" + chromosome + ".gz", compression = "gzip", header = 0, sep = "\t", usecols=columns, low_memory=False)
    
    if chromosome == "M":
        chromosome = "MT"
        df = df.replace({"#chr": {"M" : "MT"}})
    
    df = df.rename(columns={"#chr" : "chr", "pos(1-based)": "pos"})
    df.to_csv("dbNSFP4.7a_variant.chr" + chromosome + "_small.txt", index = False, sep = "\t")
    print("chromosome " + chromosome + " done")
