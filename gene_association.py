import pandas as pd
import subprocess

# ["Homo_sapiens.GRCh38.105", "grch38"]
# ["Homo_sapiens.GRCh37.87", "grch37"]
# ["Mus_musculus.GRCm39.105", "grcm39"]
# ["Danio_rerio.GRCz11.105", "grcz11"]
# ["Caenorhabditis_elegans.WBcel235.105", "wbcel1235"]

name = "Caenorhabditis_elegans.WBcel235.105"
genome = "wbcel1235"

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
            # print(myList)
            matchers2 = ["gene_id", "gene_name"]
            final = []
            for m in myList:
                # print(m)
                if "gene_id" in m:
                    final.append(m)
                if "gene_name" in m:
                    final.append(m)

            # final = [s for s in myList if any(xs in s for xs in matchers2)]
            final = [i[1].replace('"', "") for i in final]
            # print(final)
            if len(final) == 2:
                lines.append(final)
            if len(final) == 1:
                final = [i for i in final for r in range(2)]
                lines.append(final)


df = pd.DataFrame(lines, columns=["gene_id", "gene_name"])
df = df.drop_duplicates()
# df = df.sort_values(by="gene_id", ascending=True)
df.to_csv(genome + "_gene_association_not_sorted.txt", sep="\t", index=False)

with open(genome + "_gene_association.txt", "w") as outfile:
    subprocess.run(
        ["sort", genome + "_gene_association_not_sorted.txt"],
        stdout=outfile,
    )

subprocess.run(["rm", genome + "_gene_association_not_sorted.txt"])
