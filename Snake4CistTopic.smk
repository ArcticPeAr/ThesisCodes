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
        CpG = "/media/veracrypt10/Biotin/CistopicDir/{project}-small.RDS",
        coords = "/media/veracrypt10/Biotin/CistopicDir/{project}_{tfactor}.bed.hg19.wgEncodeHaibMethyl450CpgIslandDetails_emap.probes.bed",
        premeta = "/media/veracrypt10/Biotin/CistopicDir/Meta/sampleinfo_TCGA_RNA_seq_cluster.txt"
    output:
        methTable2GZ = "/media/veracrypt10/Biotin/CistopicDir/Output/{project}/MethylTable_{project}_{tfactor}.csv.xz",
        pdfpdf = "/media/veracrypt10/Biotin/CistopicDir/Output/{project}/{project}_{tfactor}.pdf",
        ctoOut = "/media/veracrypt10/Biotin/CistopicDir/Output/{project}/CTO_{project}_{tfactor}.rds",
        topicAssigToPatientOut = "/media/veracrypt10/Biotin/CistopicDir/Output/{project}/topicAssigToPatient_{project}_{tfactor}.csv",
        RegScrPrtopicOut = "/media/veracrypt10/Biotin/CistopicDir/Output/{project}/RegScrPrtopic_{project}_{tfactor}.csv",
        RegAssigUnormalOut = "/media/veracrypt10/Biotin/CistopicDir/Output/{project}/RegAssigUnormal_{project}_{tfactor}.csv",
        binaOut = "/media/veracrypt10/Biotin/CistopicDir/Output/{project}/bina_{project}_{tfactor}.csv",
    params:
        string = "{project}"
    threads: workflow.cores
    shell:
        "Rscript CisTopic.R {threads} {params.string} {input.CpG} {input.coords} {input.premeta} {output.methTable2GZ} {output.pdfpdf} {output.topicAssigToPatientOut} {output.RegScrPrtopicOut} {output.RegAssigUnormalOut} {output.ctoOut} {output.binaOut}"

#Use Indepentent Component Analysis for decomposing the data into independent components. This can reveal hidden structures or patterns that are not apparent in original data
rule fICA:
    input:
        ctoIn = "/media/veracrypt10/Biotin/CistopicDir/Output/{project}/CTO_{project}_{tfactor}.rds"
    output:
        fICApdf = "/media/veracrypt10/Biotin/CistopicDir/Output/{project}/fICAplots_{project}_{tfactor}.pdf"
    priority:
        100
    shell:
        "Rscript fICA.R {input.ctoIn} {output.fICAPDF}"

#Dimensionality reduction with UMAP
rule UMAP: 
    input:
        ctoIn = "/media/veracrypt10/Biotin/CistopicDir/Output/{project}/CTO_{project}_{tfactor}.rds"
    output:
        UMAPpdf = "/media/veracrypt10/Biotin/CistopicDir/Output/{project}/UMAPplots_{project}_{tfactor}.pdf"
    priority:
        99
    shell:
        "Rscript UMAP.R {input.ctoIn} {output.UMAPPDF}"

#Supposed to find clusters for each topic and make a pdf with the clusters.
rule clusters:
    input: 
        ctoIn = "/media/veracrypt10/Biotin/CistopicDir/Output/{project}/CTO_{project}_{tfactor}.rds"
    output:
        clusterPDF = "/media/veracrypt10/Biotin/CistopicDir/Output/{project}/Clusterplots_{project}_{tfactor}.pdf"
    priority:
        98
    shell:
        "Rscript Clusters.R {input.ctoIn} {output.clusterPDF}"



#Make a bed file with the CpG coordinates and the topic scores for each CpG
#Python find highest contributing patient and CpG to topics and make violinplots with seaborn.
rule seabornViolin:
    input:
        methTable2GZ = "/media/veracrypt10/Biotin/CistopicDir/Output/{project}/MethylTable_{project}_{tfactor}.csv.xz",
        RegScrPrtopicOut = "/media/veracrypt10/Biotin/CistopicDir/Output/{project}/RegScrPrtopic_{project}_{tfactor}.csv",
        RegAssigUnormalOut = "/media/veracrypt10/Biotin/CistopicDir/Output/{project}/RegAssigUnormal_{project}_{tfactor}.csv",
        topicAssigToPatientOut = "/media/veracrypt10/Biotin/CistopicDir/Output/{project}/topicAssigToPatient_{project}_{tfactor}.csv",
        binaOut = "/media/veracrypt10/Biotin/CistopicDir/Output/{project}/bina_{project}_{tfactor}.csv",
    output:
        ProbeTopicScoreCSV = "/media/veracrypt10/Biotin/CistopicDir/Output/{project}/ProbeTopicScore_{project}_{tfactor}.csv",    
        ViolinPDF = "/media/veracrypt10/Biotin/CistopicDir/Output/{project}/Violin_{project}_{tfactor}.pdf"
    shell:
        "python3 seabornViolin.py {input.methTable2GZ} {input.RegScrPrtopicOut} {input.RegAssigUnormalOut} {input.topicAssigToPatientOut} {input.binaOut} {output.ViolinPDF} {output.ProbeTopicScoreCSV}"

#Make heatmap in ggplot2 from top patients and CpGs found in 'rule seabornViolin'
rule probeHeatmap:
    input:
        "/media/veracrypt10/Biotin/CistopicDir/Output/{project}/ProbeTopicScore_{project}_{tfactor}.csv"
    output:
        "/media/veracrypt10/Biotin/CistopicDir/Output/{project}/PrbTpcScrHeatmap_{project}_{tfactor}.pdf"
    shell:
        "Rscript Heatmap.R {input} {output}" 


#CSVs can be big so this rule compresses the CSV to XZ to save space.
rule MethTabl_CSVtoXZ:
    input:
        CSV = "/media/veracrypt10/Biotin/CistopicDir/Output/{project}/MethylTable_{project}_{tfactor}.csv"
    output:
        XZ = "/media/veracrypt10/Biotin/CistopicDir/Output/{project}/MethylTable_{project}_{tfactor}.csv.xz"
    shell:
        "xz --threads=0 -9v {input.CSV}"


#Make outputfolder and add variable number of bedfiles there:
checkpoint get_bed:
    input:
        "/media/veracrypt10/Biotin/CistopicDir/Output/{project}/CTO_{project}_{tfactor}.rds"
    output:
        directory("/media/veracrypt10/Biotin/CistopicDir/Output/{project}/bedfiles_{project}_{tfactor}")
    priority:
        99
    shell:
        "Rscript CisBed.R /media/veracrypt10/Biotin/CistopicDir/Output/{wildcards.project}/CTO_{wildcards.project}_{wildcards.tfactor}.rds /media/veracrypt10/Biotin/CistopicDir/Output/{wildcards.project}/bedfiles_{wildcards.project}_{wildcards.tfactor}"


#Liftover bed from HG19 to HG38 and then delete the HG19 (To easier use snakemake with unibind)
rule liftover:
    input:
        HG19Bed = "/media/veracrypt10/Biotin/CistopicDir/Output/{project}/bedfiles_{project}_{tfactor}/Topic_{bednumber}.bed"
    output:
        HG38Bed = "/media/veracrypt10/Biotin/CistopicDir/Output/{project}/bedfiles_{project}_{tfactor}/Topic_{bednumber}_HG38.bed"
    priority:
        98
    shell:
        "/home/petear/Documents/liftOverDir/liftOver {input.HG19Bed} /media/veracrypt10/Biotin/CistopicDir/hg19ToHg38.over.chain.gz {output.HG38Bed} unMapped"





#make function to continue to work on sepearate bedfiles
def cont_work_bed(wildcards):
    checkpoint_output = checkpoints.get_bed.get(**wildcards).output[0]  #Collect each output (from output[0]) from checkpoint
    file_names = expand("/media/veracrypt10/Biotin/CistopicDir/Output/{project}/bedfiles_{project}_{tfactor}/Topic_{bednumber}_HG38.bed",
                        project = wildcards.project,
                        tfactor = wildcards.tfactor,
                        bednumber = glob_wildcards(os.path.join(checkpoint_output, "Topic_{bednumber}.bed")).bednumber)
    return file_names


#Create universe (Background)
rule universe:
    input:
        cont_work_bed
    output:
        Universe = "/media/veracrypt10/Biotin/CistopicDir/Output/{project}/{project}_{tfactor}_universe.bed"
    shell:
        "cat /media/veracrypt10/Biotin/CistopicDir/Output/{wildcards.project}/bedfiles_{wildcards.project}_{wildcards.tfactor}/*_HG38.bed | bedtools sort | bedtools merge > {output}"


#Run unibind enrichment analysis on each bed file from checkpoint get_bed. (Should this also be a checkpoint?) .
rule unibind:
    input:
        topicBed = "/media/veracrypt10/Biotin/CistopicDir/Output/{project}/bedfiles_{project}_{tfactor}/Topic_{bednumber}_HG38.bed",
        Universe = "/media/veracrypt10/Biotin/CistopicDir/Output/{project}/{project}_{tfactor}_universe.bed"
    output:
        "/media/veracrypt10/Biotin/CistopicDir/Output/{project}/bedfiles_{project}_{tfactor}/UB_output_Topic{bednumber}/allEnrichments_swarm.pdf"
    shell:
        "bash /media/veracrypt10/Biotin/CistopicDir/unibind_enrichment/bin/UniBind_enrich.sh oneSetBg /media/veracrypt10/Biotin/CistopicDir/UnibindFolder/hg38_robust_UniBind_LOLA.RDS {input.topicBed} {input.Universe} /media/veracrypt10/Biotin/CistopicDir/Output/{wildcards.project}/bedfiles_{wildcards.project}_{wildcards.tfactor}/UB_output_Topic{wildcards.bednumber}"


#make function to continue to work on sepearate bedfiles
def unibindFunc(wildcards):
    checkpoint_output = checkpoints.get_bed.get(**wildcards).output[0]  #Collect each output (from output[0]) from checkpoint
    file_names = expand("/media/veracrypt10/Biotin/CistopicDir/Output/{project}/bedfiles_{project}_{tfactor}/UB_output_Topic{bednumber}/allEnrichments_swarm.pdf",
                        project = wildcards.project,
                        tfactor = wildcards.tfactor,
                        bednumber = glob_wildcards(os.path.join(checkpoint_output, "Topic_{bednumber}.bed")).bednumber)
    return file_names

rule mergeUBPDF:
    input:
        unibindFunc
    output:
        "/media/veracrypt10/Biotin/CistopicDir/Output/{project}/UBallEnrichmentsSwarmplot_{project}_{tfactor}.pdf"
    shell:
        "pdfunite {input} {output}"

rule mergeCisTPDF:
    input:
        cistopicPDF = "/media/veracrypt10/Biotin/CistopicDir/Output/{project}/{project}_{tfactor}.pdf",
        UBpdf = "/media/veracrypt10/Biotin/CistopicDir/Output/{project}/UBallEnrichmentsSwarmplot_{project}_{tfactor}.pdf"
    output:
        "/media/veracrypt10/Biotin/CistopicDir/Output/{project}/merged_{project}_{tfactor}.pdf"
    shell:
        "pdfunite {input.cistopicPDF} {input.UBpdf} {output}"

rule testFunc:
    input:
        unibindFunc
    output:
        "test_{project}_{tfactor}.txt"
    shell:
        "echo 'done' > {output}"

# #&& rm Topic_{wildcards.bednumber}.bed

    ########################################
rule rGREAT:
    input:
        HG38Bed = "/media/veracrypt10/Biotin/CistopicDir/Output/{project}/bedfiles_{project}_{tfactor}/Topic_{bednumber}_HG38.bed"
    output:
        GreatPDF = "/media/veracrypt10/Biotin/CistopicDir/Output/{project}/bedfiles_{project}_{tfactor}/Topic_{bednumber}_HG38.bed_GREAT.pdf"
    shell:
        "Rscript GREAT.R /media/veracrypt10/Biotin/CistopicDir/Output/{wildcards.project}/bedfiles_{wildcards.project}_{wildcards.tfactor}/Topic_{wildcards.bednumber}_HG38.bed /media/veracrypt10/Biotin/CistopicDir/Output/{wildcards.project}/bedfiles_{wildcards.project}_{wildcards.tfactor}/{wildcards.bednumber}_HG38.bed_GREAT.pdf"

def rGREATFunc(wildcards):
    checkpoint_output = checkpoints.get_bed.get(**wildcards).output[0]  #Collect each output (from output[0]) from checkpoint
    file_names = expand("/media/veracrypt10/Biotin/CistopicDir/Output/{project}/bedfiles_{project}_{tfactor}/Topic_{bednumber}_HG38.bed_GREAT.pdf",
                        project = wildcards.project,
                        tfactor = wildcards.tfactor,
                        bednumber = glob_wildcards(os.path.join(checkpoint_output, "Topic_{bednumber}.bed")).bednumber)
    return file_names

#Merge PDFs from rGREAT
rule MergeGREATPDFs:
    input:
        rGREATFunc
    output:
        "/media/veracrypt10/Biotin/CistopicDir/Output/{project}/GREAT_merged_{project}_{tfactor}.pdf"
    shell:
        "pdfunite {input} {output}"

#Merge all PDFs to make a final PDF with full report and all plots
rule completePDF:
    input:
        cisTPDF =   "/media/veracrypt10/Biotin/CistopicDir/Output/{project}/merged_{project}_{tfactor}.pdf",
        GreatPDF = "/media/veracrypt10/Biotin/CistopicDir/Output/{project}/GREAT_merged_{project}_{tfactor}.pdf"
    output:
        "/media/veracrypt10/Biotin/CistopicDir/Output/{project}/Complete_{project}_{tfactor}.pdf"
    shell:
        "pdfunite {input.cisTPDF} {input.GreatPDF} {output}"
 
