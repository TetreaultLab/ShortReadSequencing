# Make directories
mkdir -p /lustre03/project/6019267/shared/tools/{REFERENCE_GENOMES,ANNOTATED_GTF}/{human,mouse,worm,zebrafish}
mkdir -p /lustre03/project/6019267/shared/tools/{REFERENCE_GENOMES,ANNOTATED_GTF}/human/{hg19_GRCh37,hg38_GRCh38}/{index_star,index_bowtie2,index_bwa,index_salmon}
mkdir -p /lustre03/project/6019267/shared/tools/{REFERENCE_GENOMES,ANNOTATED_GTF}/mouse/mm39_GRCm39/{index_star,index_bowtie2,index_bwa,index_salmon}
mkdir -p /lustre03/project/6019267/shared/tools/{REFERENCE_GENOMES,ANNOTATED_GTF}/worm/ce11_WBcel235/{index_star,index_bowtie2,index_bwa,index_salmon}
mkdir -p /lustre03/project/6019267/shared/tools/{REFERENCE_GENOMES,ANNOTATED_GTF}/zebrafish/danRer11_GRCz11/{index_star,index_bowtie2,index_bwa,index_salmon}

# Download ENSEMBL fasta/annotation files
# human - GRCh37
rsync -avzzP rsync://ftp.ensembl.org/ensembl/pub/grch37/current/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz REFERENCE_GENOMES/human/hg19_GRCh37 # GCA_000001405.14
rsync -avzzP rsync://ftp.ensembl.org/ensembl/pub/grch37/current/gtf/homo_sapiens/Homo_sapiens.GRCh37.87.gtf.gz ANNOTATED_GTF/human/hg19_GRCh37
# human - GRCh38
rsync -avzzP rsync://ftp.ensembl.org/ensembl/pub/release-110/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz REFERENCE_GENOMES/human/hg38_GRCh38 # GCA_000001405.29
rsync -avzzP rsync://ftp.ensembl.org/ensembl/pub/release-110/gtf/homo_sapiens/Homo_sapiens.GRCh38.110.gtf.gz ANNOTATED_GTF/human/hg38_GRCh38
rsync -avzzP rsync://ftp.ensembl.org/ensembl/pub/release-110/regulation/homo_sapiens/homo_sapiens.GRCh38.Regulatory_Build.regulatory_features.20221007.gff.gz ANNOTATED_GTF/human/hg38_GRCh38
**rsync -avzzP rsync://ftp.ensembl.org/ensembl/pub/release-110/variation/gvf/homo_sapiens/1000GENOMES-phase_3.gvf.gz ANNOTATED_GTF/human/hg38_GRCh38
**https://ftp.ensembl.org/pub/release-110/variation/gvf/homo_sapiens/

# mouse
rsync -avzzP rsync://ftp.ensembl.org/ensembl/pub/release-110/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz REFERENCE_GENOMES/mouse/mm39_GRCm39 # GCA_000001635.9
rsync -avzzP rsync://ftp.ensembl.org/ensembl/pub/release-110/gtf/mus_musculus/Mus_musculus.GRCm39.110.gtf.gz ANNOTATED_GTF/mouse/mm39_GRCm39
rsync -avzzP rsync://ftp.ensembl.org/ensembl/pub/release-110/regulation/mus_musculus/mus_musculus.GRCm39.Regulatory_Build.regulatory_features.20221007.gff.gz ANNOTATED_GTF/mouse/mm39_GRCm39
rsync -avzzP rsync://ftp.ensembl.org/ensembl/pub/release-110/variation/gvf/mus_musculus/mus_musculus.gvf.gz ANNOTATED_GTF/mouse/mm39_GRCm39

# worm
rsync -avzzP rsync://ftp.ensembl.org/ensembl/pub/release-110/fasta/caenorhabditis_elegans/dna/Caenorhabditis_elegans.WBcel235.dna.toplevel.fa.gz REFERENCE_GENOMES/worm/ce11_WBcel235 # GCA_000002985.3
rsync -avzzP rsync://ftp.ensembl.org/ensembl/pub/release-110/gtf/caenorhabditis_elegans/Caenorhabditis_elegans.WBcel235.110.gtf.gz ANNOTATED_GTF/worm/ce11_WBcel235

# zebrafish
rsync -avzzP rsync://ftp.ensembl.org/ensembl/pub/release-110/fasta/danio_rerio/dna/Danio_rerio.GRCz11.dna.primary_assembly.fa.gz REFERENCE_GENOMES/zebrafish/danRer11_GRCz11 # GCA_000002035.4
rsync -avzzP rsync://ftp.ensembl.org/ensembl/pub/release-110/gtf/danio_rerio/Danio_rerio.GRCz11.110.gtf.gz ANNOTATED_GTF/zebrafish/danRer11_GRCz11
rsync -avzzP rsync://ftp.ensembl.org/ensembl/pub/release-110/variation/gvf/danio_rerio/danio_rerio.gvf.gz ANNOTATED_GTF/zebrafish/danRer11_GRCz11

# transcriptome (cDNA)
rsync -avzzP rsync://ftp.ensembl.org/ensembl/pub/grch37/current/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh37.cdna.all.fa.gz REFERENCE_GENOMES/human/hg19_GRCh37
rsync -avzzP rsync://ftp.ensembl.org/ensembl/pub/release-110/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz REFERENCE_GENOMES/human/hg38_GRCh38
rsync -avzzP rsync://ftp.ensembl.org/ensembl/pub/release-110/fasta/mus_musculus/cdna/Mus_musculus.GRCm39.cdna.all.fa.gz REFERENCE_GENOMES/mouse/mm39_GRCm39
rsync -avzzP rsync://ftp.ensembl.org/ensembl/pub/release-110/fasta/danio_rerio/cdna/Danio_rerio.GRCz11.cdna.all.fa.gz REFERENCE_GENOMES/zebrafish/danRer11_GRCz11
rsync -avzzP rsync://ftp.ensembl.org/ensembl/pub/release-110/fasta/caenorhabditis_elegans/cdna/Caenorhabditis_elegans.WBcel235.cdna.all.fa.gz REFERENCE_GENOMES/worm/ce11_WBcel235



#!/bin/bash
#
## Load environment
#
module load StdEnv/2020 star/2.7.9a bowtie2/2.5.1 bwa/0.7.17 gcc/9.3.0 openmpi/4.0.3 salmon/1.4.0
#
## Genome Generate
#
for genome in hg19 hg38 mm39 ce11 danRer11
do
    ref_genome=""
    fasta=""
    annotation=""
    case "$genome" in
        "hg19")
        ref_genome="human/hg19_GRCh37"
        fasta="Homo_sapiens.GRCh37.dna.primary_assembly.fa"
        annotation="Homo_sapiens.GRCh37.87.gtf"
        transcriptome="Homo_sapiens.GRCh37.cdna.all.fa"
        SA="14"
        ;;
        "hg38")
        ref_genome="human/hg38_GRCh38"
        fasta="Homo_sapiens.GRCh38.dna.primary_assembly.fa"
        annotation="Homo_sapiens.GRCh38.110.gtf"
        transcriptome="Homo_sapiens.GRCh38.cdna.all.fa"
        SA="14"
        ;;
        "mm39")
        ref_genome="mouse/mm39_GRCm39"
        fasta="Mus_musculus.GRCm39.dna.primary_assembly.fa"
        annotation="Mus_musculus.GRCm39.110.gtf"
        transcriptome="Mus_musculus.GRCm39.cdna.all.fa"
        SA="14"
        ;;
        "ce11")
        ref_genome="worm/ce11_WBcel235"
        fasta="Caenorhabditis_elegans.WBcel235.dna.toplevel.fa"
        annotation="Caenorhabditis_elegans.WBcel235.110.gtf"
        transcriptome="Caenorhabditis_elegans.WBcel235.cdna.all.fa"
        SA="12"
        ;;
        "danRer11")
        ref_genome="zebrafish/danRer11_GRCz11"
        fasta="Danio_rerio.GRCz11.dna.primary_assembly.fa"
        annotation="Danio_rerio.GRCz11.110.gtf"
        transcriptome="Danio_rerio.GRCz11.cdna.all.fa"
        SA="14"
        ;;
    esac
    
    

    gunzip -c -k /lustre03/project/6019267/shared/tools/REFERENCE_GENOMES/${ref_genome}/${fasta}.gz > /lustre03/project/6019267/shared/tools/REFERENCE_GENOMES/${ref_genome}/${fasta}
    gunzip -c -k /lustre03/project/6019267/shared/tools/ANNOTATED_GTF/${ref_genome}/${annotation}.gz > /lustre03/project/6019267/shared/tools/ANNOTATED_GTF/${ref_genome}/${annotation}
    #
    STAR --runThreadN 8 \
        --runMode genomeGenerate \
        --genomeDir /lustre03/project/6019267/shared/tools/REFERENCE_GENOMES/${ref_genome}/index_star/ \
        --genomeFastaFiles /lustre03/project/6019267/shared/tools/REFERENCE_GENOMES/${ref_genome}/${fasta} \
        --sjdbGTFfile /lustre03/project/6019267/shared/tools/ANNOTATED_GTF/${ref_genome}/${annotation} \
        --sjdbOverhang 99
    #
    # Bowtie2
    #
    bowtie2-build --threads 8 \
                /lustre03/project/6019267/shared/tools/REFERENCE_GENOMES/${ref_genome}/${fasta} \
                /lustre03/project/6019267/shared/tools/REFERENCE_GENOMES/${ref_genome}/index_bowtie2/${genome}
    #
    # BWA
    #
    bwa index -p /lustre03/project/6019267/shared/tools/REFERENCE_GENOMES/${ref_genome}/index_bwa/${genome} \
            -a bwtsw \
            /lustre03/project/6019267/shared/tools/REFERENCE_GENOMES/${ref_genome}/${fasta}
    #
    # Salmon
    # https://combine-lab.github.io/alevin-tutorial/2019/selective-alignment/
    #
    grep "^>" /lustre03/project/6019267/shared/tools/REFERENCE_GENOMES/${ref_genome}/${fasta} | cut -d " " -f 1 > /lustre03/project/6019267/shared/tools/REFERENCE_GENOMES/${ref_genome}/decoys.txt
    sed -i.bak -e 's/>//g' /lustre03/project/6019267/shared/tools/REFERENCE_GENOMES/${ref_genome}/decoys.txt
    cat /lustre03/project/6019267/shared/tools/REFERENCE_GENOMES/${ref_genome}/${transcriptome}.gz /lustre03/project/6019267/shared/tools/REFERENCE_GENOMES/${ref_genome}/${fasta}.gz > /lustre03/project/6019267/shared/tools/REFERENCE_GENOMES/${ref_genome}/gentrome.fa.gz
    
    salmon index -p 8 \
                 -t /lustre03/project/6019267/shared/tools/REFERENCE_GENOMES/${ref_genome}/gentrome.fa.gz \
                 -i /lustre03/project/6019267/shared/tools/REFERENCE_GENOMES/${ref_genome}/index_salmon/ \
                 -d /lustre03/project/6019267/shared/tools/REFERENCE_GENOMES/${ref_genome}/decoys.txt
    #
    ## Remove none GZ files
    #
    rm /lustre03/project/6019267/shared/tools/REFERENCE_GENOMES/${ref_genome}/${fasta}
    rm /lustre03/project/6019267/shared/tools/ANNOTATED_GTF/${ref_genome}/${annotation}
    #
    # DONE
    #
    echo -n "GENERATE INDEXES DONE"
done