import argparse
import numpy as np
import pandas as pd
import toml


parser = argparse.ArgumentParser(
    prog="GenerateJob",
    description="Generate SLURM job for samples based on a TOML project config file.",
)

parser.add_argument("--genome", type=str, required=True, help="Reference genome.")
parser.add_argument("--output", type=str, required=True, help="Output path.")
parser.add_argument("--sample", type=str, required=True, help="Sample name.")

args = parser.parse_args()

sample = args.sample
genome = args.genome
path = args.output

if genome == "grch37" or genome == "grch38":
    # dbNSFP
    chromosomes = [
        "1",
        "2",
        "3",
        "4",
        "5",
        "6",
        "7",
        "8",
        "9",
        "10",
        "11",
        "12",
        "13",
        "14",
        "15",
        "16",
        "17",
        "18",
        "19",
        "20",
        "21",
        "22",
        "X",
        "Y",
        "MT",
    ]

    var = pd.read_csv(
        path + "/" + sample + "_annotated.txt", header=0, sep="\t", low_memory=False
    )
    var = var.rename(columns={"CHROM": "CHROM_" + genome, "POS": "POS_" + genome})

    var = var.astype({"CHROM_" + genome: str, "POS_" + genome: str})

    appended_data = []
    ref = "/lustre09/project/6019267/shared/tools/main_pipelines/short-read/dbNSFP"

    for chromosome in chromosomes:
        print(chromosome)
        db = pd.read_csv(
            ref + "/dbNSFP4.7a_variant.chr" + chromosome + "_small.txt",
            header=0,
            sep="\t",
            low_memory=False,
        )

        db = db.astype({"CHROM_" + genome: str, "POS_" + genome: str})

        var_chr = var[var["CHROM_" + genome] == chromosome]

        m = pd.merge(
            var_chr,
            db,
            how="left",
            on=["CHROM_" + genome, "POS_" + genome, "REF", "ALT"],
        )
        appended_data.append(m)

    final = pd.concat(appended_data)
    final.to_csv(path + "/" + sample + "_annotated_dbNSFP.txt", sep="\t", index=False)

    n_rows = len(final.index)
    rows_ann = final["rsID"].count()
    percentage = round(rows_ann / n_rows * 100, 2)

    print(
        "\nFinal\nannotated variants: ",
        rows_ann,
        "\ntotal variants:",
        n_rows,
        "\npercentage dbNSFP annotation: ",
        percentage,
        "%",
    )

    # Formating
    final = final.rename(
        columns={
            "REF": "Ref",
            "ALT": "Alt",
            "QUAL": "Quality",
            # "TYPE": "Variation",
            "HOM": "Zygosity",
            "DP": "Read_depth",
            "ANN[*].GENE": "Gene_name",
            "ANN[*].GENEID": "Gene_id",
            "ANN[*].FEATUREID": "Transcript",
            "ANN[*].EFFECT": "Effect",
            "ANN[*].IMPACT": "Impact",
            "ANN[*].BIOTYPE": "Biotype",
            "ANN[*].HGVS_C": "Codon_change",
            "ANN[*].HGVS_P": "Protein_change",
        }
    )

    final["Position_" + genome] = (
        final["CHROM_" + genome].astype(str) + ":" + final["POS_" + genome].astype(str)
    )
    final["Quality"] = round(final["Quality"], 2)
    final = final.replace({"Zygosity": {True: "Hom", False: "Het"}})

    columns = [
        "Gene_name",
        "Gene_id",
        "Transcript",
        "Effect",
        "Impact",
        "Biotype",
        "Codon_change",
        "Protein_change",
    ]

    for index, row in final.iterrows():
        infos = []

        for c in columns:
            cell_value = str(row[c])
            str_split = [sub.replace('"', "") for sub in cell_value.split("|")]

            # Append to infos, or take the first split value if there's only one
            infos.append(str_split if len(str_split) > 1 else [cell_value])

            # Update the final DataFrame with the first part of split value
            final.at[index, c] = str_split[0]

            # Transpose and concatenate corresponding elements from each list
            if len(infos) > 1:
                info_concat = [
                    "|".join([info[i] for info in infos]) for i in range(len(infos[0]))
                ]
                final.at[index, "Infos"] = "; ".join(info_concat)
            else:
                final.at[index, "Infos"] = "; ".join(infos[0])

    final = final[
        [
            "Position_" + genome,
            "Ref",
            "Alt",
            "Quality",
            "Read_depth",
            "Zygosity",
            "rsID",
            # "Variation",
            "Gene_name",
            "Gene_id",
            "Transcript",
            "Effect",
            "Impact",
            "Biotype",
            "Codon_change",
            "Protein_change",
            "Infos",
            "genename",
            "Ensembl_geneid",
            "Ensembl_transcriptid",
            "Ensembl_proteinid",
            "codon_change",
            "protein_change",
            "SIFT_score",
            "PolyPhen2_HDIV_score",
            "GERP_score",
            "phyloP_score",
            "phastCons_score",
            "MutationTaster_score",
            "FATHMM_score",
            "REVEL_score",
            "AlphaMissense_score",
            "CADD_phred_score",
            "1000G_AF",
            "ExAC_AF",
            "gnomAD_exomes_AF",
            "gnomAD_exomes_NFE_AF",
            "gnomAD_genomes_AF",
            "gnomAD_genomes_NFE_AF",
            "clinvar_id",
            "clinvar_trait",
            "clinvar_clnsig",
            "OMIM",
        ]
    ]

    # count "Het" for each gene
    final["het_count"] = final.groupby("Gene_name")["Zygosity"].transform(
        lambda x: (x == "Het").sum()
    )

    # Remplace "Het" for "Multiple-het" if het_count > 1
    final.loc[(final["het_count"] > 1) & (final["Zygosity"] == "Het"), "Zygosity"] = (
        "Multiple-het"
    )

    # Remove het_count
    final = final.drop(columns=["het_count"])

    final = final.replace(".", np.nan)
    final.to_csv(path + "/" + sample + "_variants_all.txt", sep="\t", index=False)

    ## Filtering
    # https://www.htslib.org/workflow/filter.html
    # https://jp.support.illumina.com/content/dam/illumina-support/help/Illumina_DRAGEN_Bio_IT_Platform_v3_7_1000000141465/Content/SW/Informatics/Dragen/QUAL_QD_GQ_Formulation_fDG.htm

    df_filtered = final[(final["Quality"] > 20) & (final["Read_depth"] > 5)]

    df_filtered.to_csv(
        path + "/" + sample + "_variants_filtered.txt", sep="\t", index=False
    )

    print(">>> Filters:")
    print("\t>>> Quality Phred > 20")
    print("\t>>> Number of total reads > 5")

    all_var = len(final.index)
    filtered_var = len(df_filtered.index)
    percentage = round(filtered_var / all_var * 100, 2)

    print(
        "\nNumber of variants in all: ",
        all_var,
        "\nNumber of variants in filtered: ",
        filtered_var,
        "\npercentage filtered variants: ",
        percentage,
        "%",
    )

else:
    final = pd.read_csv(
        path + "/" + sample + "_annotated.txt", header=0, sep="\t", low_memory=False
    )

    final = final.rename(
        columns={
            "REF": "Ref",
            "ALT": "Alt",
            "QUAL": "Quality",
            "TYPE": "Variation",
            "HOM": "Zygosity",
            "DP": "Read_depth",
            "ANN[*].GENE": "Gene_name",
            "ANN[*].GENEID": "Gene_id",
            "ANN[*].FEATUREID": "Transcript",
            "ANN[*].EFFECT": "Effect",
            "ANN[*].IMPACT": "Impact",
            "ANN[*].BIOTYPE": "Biotype",
            "ANN[*].HGVS_C": "Codon_change",
            "ANN[*].HGVS_P": "Protein_change",
        }
    )

    final["Position"] = final["CHROM"].astype(str) + ":" + final["POS"].astype(str)
    final["Quality"] = round(final["Quality"], 2)
    final = final.replace({"Zygosity": {True: "Hom", False: "Het"}})

    columns = [
        "Gene_name",
        "Gene_id",
        "Transcript",
        "Effect",
        "Impact",
        "Biotype",
        "Codon_change",
        "Protein_change",
    ]

    for index, row in final.iterrows():
        infos = []

        for c in columns:
            cell_value = str(row[c])
            str_split = [sub.replace('"', "") for sub in cell_value.split("|")]

            # Append to infos, or take the first split value if there's only one
            infos.append(str_split if len(str_split) > 1 else [cell_value])

            # Update the final DataFrame with the first part of split value
            final.at[index, c] = str_split[0]

            # Transpose and concatenate corresponding elements from each list
            if len(infos) > 1:
                info_concat = [
                    "|".join([info[i] for info in infos]) for i in range(len(infos[0]))
                ]
                final.at[index, "Infos"] = "; ".join(info_concat)
            else:
                final.at[index, "Infos"] = "; ".join(infos[0])

    final = final[
        [
            "Position",
            "Ref",
            "Alt",
            "Quality",
            "Read_depth",
            "Zygosity",
            "Variation",
            "Gene_name",
            "Gene_id",
            "Transcript",
            "Effect",
            "Impact",
            "Biotype",
            "Codon_change",
            "Protein_change",
            "Infos",
        ]
    ]

    # count "Het" for each gene
    final["het_count"] = final.groupby("Gene_name")["Zygosity"].transform(
        lambda x: (x == "Het").sum()
    )

    # Remplace "Het" for "Multiple-het" if het_count > 1
    final.loc[(final["het_count"] > 1) & (final["Zygosity"] == "Het"), "Zygosity"] = (
        "Multiple-het"
    )

    # Remove het_count
    final = final.drop(columns=["het_count"])

    final = final.replace(".", np.nan)
    print(final)

    final.to_csv(path + "/" + sample + "_variants_all.txt", sep="\t", index=False)

    df_filtered = final[(final["Quality"] > 20) & (final["Read_depth"] > 5)]
    df_filtered.to_csv(
        path + "/" + sample + "_variants_filtered.txt", sep="\t", index=False
    )

    print(">>> Filters:")
    print("\t>>> Quality Phred > 20")
    print("\t>>> Number of total reads > 5")

    all_var = len(final.index)
    filtered_var = len(df_filtered.index)
    percentage = round(filtered_var / all_var * 100, 2)

    print(
        "\nNumber of variants in all: ",
        all_var,
        "\nNumber of variants in filtered: ",
        filtered_var,
        "\nPercentage filtered variants: ",
        percentage,
        "%",
    )
