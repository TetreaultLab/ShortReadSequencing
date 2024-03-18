import pandas as pd

name = "Mus_musculus.GRCm39.105"
genome = "grcm39"
lines = []
with open(
    "/lustre03/project/6019267/shared/tools/PIPELINES/References/" + name + ".gtf", "r"
) as f:
    for line in f:
        line = line.strip()
        if not line.startswith("#"):
            line_strip = line.split("\t")
            matchers = ["gene_id"]
            matching = [s for s in line_strip if any(xs in s for xs in matchers)]
            line_gene = matching[0].split(";")
            myList = [i.strip().split(" ") for i in line_gene]
            matchers2 = ["gene_id", "gene_name", "gene_biotype"]
            final = [s for s in myList if any(xs in s for xs in matchers2)]
            final = [i[1].replace('"', "") for i in final]
            # print(final)
            if len(final) == 3:
                lines.append(final)

df = pd.DataFrame(lines, columns=["gene_id", "gene_name", "gene_biotype"])
df = df.drop_duplicates()
df = df.sort_values(by="gene_id", ascending=True)
df.to_csv(genome + "_gene_association.txt", sep="\t", index=False)
