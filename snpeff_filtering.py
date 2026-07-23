import argparse
import pandas as pd
import numpy as np

print(">>> Filtering snpEff output\n")

parser = argparse.ArgumentParser(
    prog="Filter snpEff output.",
    description="Filter snpEff output.",
)

parser.add_argument("--path", type=str, required=True, help="Output path.")
parser.add_argument("--sample", type=str, required=True, help="Sample name.")

args = parser.parse_args()

path = args.path
sample = args.sample

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
