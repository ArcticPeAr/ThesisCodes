 
projectTfDict = {"BRCA-US": ["FULL"]}


PROJECT = []
TFACTOR = []
for key, value in projectTfDict.items():
    for v in value:
        PROJECT.append(key)
        TFACTOR.append(v)

print(PROJECT)
print(TFACTOR)

############################################################################
rule all:
    input:
        expand("/media/veracrypt10/Biotin/CistopicDir/Output/{project}/MethylTable_{project}_{tfactor}.csv.xz", project=PROJECT, tfactor=TFACTOR), \
        expand("/media/veracrypt10/Biotin/CistopicDir/Output/{project}/{project}_{tfactor}_universe.bed", project=PROJECT, tfactor=TFACTOR), \
        expand("/media/veracrypt10/Biotin/CistopicDir/Output/{project}/Complete_{project}_{tfactor}.pdf", project=PROJECT, tfactor=TFACTOR)

############################################################################

rule FromRawToCis:
    input:
        CpG = "/media/veracrypt10/Biotin/CistopicDir/{project}.RDS",
        coords = "/media/veracrypt10/Biotin/CistopicDir/{project}_{tfactor}.bed.hg19.wgEncodeHaibMethyl450CpgIslandDetails_emap.probes.bed",
        premeta = "/media/veracrypt10/Biotin/CistopicDir/Meta/sampleinfo_TCGA_RNA_seq_cluster.txt"
    output:
        pdfpdf = "/media/veracrypt10/Biotin/CistopicDir/Output/{project}/{project}_{tfactor}.pdf",
        methTable2GZ = "/media/veracrypt10/Biotin/CistopicDir/Output/{project}/MethylTable_{project}_{tfactor}.csv.xz",
        topicAssigToPatientOut = "/media/veracrypt10/Biotin/CistopicDir/Output/{project}/topicAssigToPatient_{project}_{tfactor}.csv",
        RegScrPrtopicOut = "/media/veracrypt10/Biotin/CistopicDir/Output/{project}/RegScrPrtopic_{project}_{tfactor}.csv",
        RegAssigUnormalOut = "/media/veracrypt10/Biotin/CistopicDir/Output/{project}/RegAssigUnormal_{project}_{tfactor}.csv",
        ctoOut = "/media/veracrypt10/Biotin/CistopicDir/Output/{project}/CTO_{project}_{tfactor}.rds",
        binaOut = "/media/veracrypt10/Biotin/CistopicDir/Output/{project}/bina_{project}_{tfactor}.csv",
        metaOut = "/media/veracrypt10/Biotin/CistopicDir/Output/{project}/meta_{project}_{tfactor}.csv",
        hmatOut = "/media/veracrypt10/Biotin/CistopicDir/Output/{project}/hmat_{project}_{tfactor}.csv",
        ClusterDFOut = "/media/veracrypt10/Biotin/CistopicDir/Output/{project}/ClusterDF_{project}_{tfactor}.csv"
    priority:
        100
    params:
        string = "{project}"
    threads: workflow.cores
    shell:
        "Rscript full.R --threads {threads} {params.string_param} {input.CpG} {input.coords} {input.premeta} {output.pdfpdf} {output.methTable2GZ} {output.topicAssigToPatientOut} {output.RegScrPrtopicOut} {output.RegAssigUnormalOut} {output.ctoOut} {output.binaOut} {output.metaOut} {output.hmatOut} {output.ClusterDFOut}"

# Not sure if --threads is needed in above shell command

rule fICA:
    input:
        ctoIn = "/media/veracrypt10/Biotin/CistopicDir/Output/{project}/CTO_{project}_{tfactor}.rds"
    output:
        fICApdf = "/media/veracrypt10/Biotin/CistopicDir/Output/{project}/fICAplots_{project}_{tfactor}.pdf"
    priority:
        100
    shell:
        "Rscript fICA.R {input.ctoIn} {output.fICAPDF}"

rule UMAP: 
    input:
        ctoIn = "/media/veracrypt10/Biotin/CistopicDir/Output/{project}/CTO_{project}_{tfactor}.rds"
    output:
        UMAPpdf = "/media/veracrypt10/Biotin/CistopicDir/Output/{project}/UMAPplots_{project}_{tfactor}.pdf"
    priority:
        99
    shell:
        "Rscript UMAP.R {input.ctoIn} {output.UMAPPDF}"
    
rule heatmapAna:
    input: 
        ctoIn = "/media/veracrypt10/Biotin/CistopicDir/Output/{project}/CTO_{project}_{tfactor}.rds"
    output:
        