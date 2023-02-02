projectTfDict = {"BRCA-US": ["FULL"]}

###  DEFs   ################################################################
# inputs = [expand('{project}_{tfactor}',
    #                  project=key, tfactor=value)
#           for key, value in projectTfDict.items()]

############################################################################


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
        expand("/home/petear/CistopicDir/Output/{project}/MethylTable_{project}_{tfactor}.csv.xz", project=PROJECT, tfactor=TFACTOR), \
        expand("/home/petear/CistopicDir/Output/{project}/{project}_{tfactor}_universe.bed", project=PROJECT, tfactor=TFACTOR), \
        expand("/home/petear/CistopicDir/Output/{project}/Complete_{project}_{tfactor}.pdf", project=PROJECT, tfactor=TFACTOR)

        ############################################################################



# #Rename probes to fit current project
# if TFACTOR == "FULL":
#     rule renameFullProbeslist:
#         input:
#             "FULL.bed.hg19.wgEncodeHaibMethyl450CpgIslandDetails_emap.probes.bed"
#         output:
#             "{project}_FULL.bed.hg19.wgEncodeHaibMethyl450CpgIslandDetails_emap.probes.bed"
#         shell:
#             "cp *.bed.hg19.wgEncodeHaibMethyl450CpgIslandDetails_emap.probes.bed {wildcards.project}_FULL.bed.hg19.wgEncodeHaibMethyl450CpgIslandDetails_emap.probes.bed"

#Remove rows(probes) with more than 50% NAs and impute the rest (Change FROM AMELIA to METHYLIMP when finished):
#metadata with PAM50 and ER.Status is in same location with this name:: sampleinfo_TCGA_RNA_seq_cluster.txt
#metadata with icgc_donor_id etc is at same location with this name::  {project}_sample_Info_260620.tsv
rule FromRawToCis:
    input:
        CpG = "/home/petear/CistopicDir/{project}.RDS",
        coords = "/home/petear/CistopicDir/{project}_{tfactor}.bed.hg19.wgEncodeHaibMethyl450CpgIslandDetails_emap.probes.bed",
        premeta = "/home/petear/CistopicDir/Meta/sampleinfo_TCGA_RNA_seq_cluster.txt"
    output:
        pdfpdf = "/home/petear/CistopicDir/Output/{project}/{project}_{tfactor}.pdf",
        methTable2GZ = "/home/petear/CistopicDir/Output/{project}/MethylTable_{project}_{tfactor}.csv.xz",
        topicAssigToPatientOut = "/home/petear/CistopicDir/Output/{project}/topicAssigToPatient_{project}_{tfactor}.csv",
        RegScrPrtopicOut = "/home/petear/CistopicDir/Output/{project}/RegScrPrtopic_{project}_{tfactor}.csv",
        RegAssigUnormalOut = "/home/petear/CistopicDir/Output/{project}/RegAssigUnormal_{project}_{tfactor}.csv",
        ctoOut = "/home/petear/CistopicDir/Output/{project}/CTO_{project}_{tfactor}.rds",
        binaOut = "/home/petear/CistopicDir/Output/{project}/bina_{project}_{tfactor}.csv",
        metaOut = "/home/petear/CistopicDir/Output/{project}/meta_{project}_{tfactor}.csv",
        hmatOut = "/home/petear/CistopicDir/Output/{project}/hmat_{project}_{tfactor}.csv",
        ClusterDFOut = "/home/petear/CistopicDir/Output/{project}/ClusterDF_{project}_{tfactor}.csv"
    priority:
        100
    shell:
        "Rscript full.R {input.CpG} {input.coords} {input.premeta} {output.pdfpdf} {output.methTable2GZ} {output.topicAssigToPatientOut} {output.RegScrPrtopicOut} {output.RegAssigUnormalOut} {output.ctoOut} {output.binaOut} {output.metaOut} {output.hmatOut} {output.ClusterDFOut}"

#Python find highest contributing patient and CpG to topics and make violinplots with seaborn.
rule topPatientCpG:
    input:
        methTab = "/home/petear/CistopicDir/Output/{project}/MethylTable_{project}_{tfactor}.csv",
        regScrNorm = "/home/petear/CistopicDir/Output/{project}/RegScrPrtopic_{project}_{tfactor}.csv",
        regScrUnrm = "/home/petear/CistopicDir/Output/{project}/RegAssigUnormal_{project}_{tfactor}.csv",
        topicAssig = "/home/petear/CistopicDir/Output/{project}/topicAssigToPatient_{project}_{tfactor}.csv",
        binFrame = "/home/petear/CistopicDir/Output/{project}/binarized_{project}_{tfactor}.csv"
    output:
        "/home/petear/CistopicDir/Output/{project}/Violin_{project}_{tfactor}.pdf"
        "/home/petear/CistopicDir/Output/{project}/ProbeTopicScore_{project}_{tfactor}.csv"
    script:
        "seabornViolin.py"


#Make heatmap in ggplot2 from top patients and CpGs found in 'rule seabornViolin'
rule probeHeatmap:
    input:
        "/home/petear/CistopicDir/Output/{project}/ProbeTopicScore_{project}_{tfactor}.csv"
    output:
        "/home/petear/CistopicDir/Output/{project}/PrbTpcScrHeatmap_{project}_{tfactor}.pdf"
    shell:
        "Rscript Heatmap.R {input} {output}"

#CSVs can be big so this rule compresses the CSV to XZ to save space.
rule MethTabl_CSVtoXZ:
    input:
        CSV = "/home/petear/CistopicDir/Output/{project}/MethylTable_{project}_{tfactor}.csv"
    output:
        XZ = "/home/petear/CistopicDir/Output/{project}/MethylTable_{project}_{tfactor}.csv.xz"
    shell:
        "xz --threads=0 -9v {input.CSV}"


#Make outputfolder and add variable number of bedfiles there:
checkpoint get_bed:
    input:
        "/home/petear/CistopicDir/Output/{project}/CTO_{project}_{tfactor}.rds"
    output:
        directory("/home/petear/CistopicDir/Output/{project}/bedfiles_{project}_{tfactor}")
    priority:
        99
    shell:
        "Rscript CisBed.R /home/petear/CistopicDir/Output/{wildcards.project}/CTO_{wildcards.project}_{wildcards.tfactor}.rds /home/petear/CistopicDir/Output/{wildcards.project}/bedfiles_{wildcards.project}_{wildcards.tfactor}"


#Liftover bed from HG19 to HG38 and then delete the HG19 (To easier use snakemake with unibind)
rule liftover:
    input:
        HG19Bed = "/home/petear/CistopicDir/Output/{project}/bedfiles_{project}_{tfactor}/Topic_{bednumber}.bed"
    output:
        HG38Bed = "/home/petear/CistopicDir/Output/{project}/bedfiles_{project}_{tfactor}/Topic_{bednumber}_HG38.bed"
    priority:
        98
    shell:
        "liftOver {input.HG19Bed} /home/petear/CistopicDir/hg19ToHg38.over.chain.gz {output.HG38Bed} unMapped"





#make function to continue to work on sepearate bedfiles
def cont_work_bed(wildcards):
    checkpoint_output = checkpoints.get_bed.get(**wildcards).output[0]  #Collect each output (from output[0]) from checkpoint
    file_names = expand("/home/petear/CistopicDir/Output/{project}/bedfiles_{project}_{tfactor}/Topic_{bednumber}_HG38.bed",
                        project = wildcards.project,
                        tfactor = wildcards.tfactor,
                        bednumber = glob_wildcards(os.path.join(checkpoint_output, "Topic_{bednumber}.bed")).bednumber)
    return file_names


#Create universe (Background)
rule universe:
    input:
        cont_work_bed
    output:
        Universe = "/home/petear/CistopicDir/Output/{project}/{project}_{tfactor}_universe.bed"
    shell:
        "cat /home/petear/CistopicDir/Output/{wildcards.project}/bedfiles_{wildcards.project}_{wildcards.tfactor}/*_HG38.bed | bedtools sort | bedtools merge > {output}"


#Run unibind enrichment analysis on each bed file from checkpoint get_bed. (Should this also be a checkpoint?) .
rule unibind:
    input:
        topicBed = "/home/petear/CistopicDir/Output/{project}/bedfiles_{project}_{tfactor}/Topic_{bednumber}_HG38.bed",
        Universe = "/home/petear/CistopicDir/Output/{project}/{project}_{tfactor}_universe.bed"
    output:
        "/home/petear/CistopicDir/Output/{project}/bedfiles_{project}_{tfactor}/UB_output_Topic{bednumber}/allEnrichments_swarm.pdf"

    shell:
        "bash /home/petear/CistopicDir/UnibindFolder/UnibindMaster/bin/UniBind_enrich.sh oneSetBg /home/petear/CistopicDir/UnibindFolder/LolaUpdated/UniBind_LOLA.RDS {input.topicBed} {input.Universe} /home/petear/CistopicDir/Output/{wildcards.project}/bedfiles_{wildcards.project}_{wildcards.tfactor}/UB_output_Topic{wildcards.bednumber}"


#make function to continue to work on sepearate bedfiles
def unibindFunc(wildcards):
    checkpoint_output = checkpoints.get_bed.get(**wildcards).output[0]  #Collect each output (from output[0]) from checkpoint
    file_names = expand("/home/petear/CistopicDir/Output/{project}/bedfiles_{project}_{tfactor}/UB_output_Topic{bednumber}/allEnrichments_swarm.pdf",
                        project = wildcards.project,
                        tfactor = wildcards.tfactor,
                        bednumber = glob_wildcards(os.path.join(checkpoint_output, "Topic_{bednumber}.bed")).bednumber)
    return file_names

rule mergeUBPDF:
    input:
        unibindFunc
    output:
        "/home/petear/CistopicDir/Output/{project}/UBallEnrichmentsSwarmplot_{project}_{tfactor}.pdf"
    shell:
        "pdfunite {input} {output}"

rule mergeCisTPDF:
    input:
        cistopicPDF = "/home/petear/CistopicDir/Output/{project}/{project}_{tfactor}.pdf",
        UBpdf = "/home/petear/CistopicDir/Output/{project}/UBallEnrichmentsSwarmplot_{project}_{tfactor}.pdf"
    output:
        "/home/petear/CistopicDir/Output/{project}/merged_{project}_{tfactor}.pdf"
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
        HG38Bed = "/home/petear/CistopicDir/Output/{project}/bedfiles_{project}_{tfactor}/Topic_{bednumber}_HG38.bed"
    output:
        GreatPDF = "/home/petear/CistopicDir/Output/{project}/bedfiles_{project}_{tfactor}/Topic_{bednumber}_HG38.bed_GREAT.pdf"
    shell:
        "Rscript GREAT.R /home/petear/CistopicDir/Output/{wildcards.project}/bedfiles_{wildcards.project}_{wildcards.tfactor}/Topic_{wildcards.bednumber}_HG38.bed /home/petear/CistopicDir/Output/{wildcards.project}/bedfiles_{wildcards.project}_{wildcards.tfactor}/{wildcards.bednumber}_HG38.bed_GREAT.pdf"

def rGREATFunc(wildcards):
    checkpoint_output = checkpoints.get_bed.get(**wildcards).output[0]  #Collect each output (from output[0]) from checkpoint
    file_names = expand("/home/petear/CistopicDir/Output/{project}/bedfiles_{project}_{tfactor}/Topic_{bednumber}_HG38.bed_GREAT.pdf",
                        project = wildcards.project,
                        tfactor = wildcards.tfactor,
                        bednumber = glob_wildcards(os.path.join(checkpoint_output, "Topic_{bednumber}.bed")).bednumber)
    return file_names


rule MergeGREATPDFs:
    input:
        rGREATFunc
    output:
        "/home/petear/CistopicDir/Output/{project}/GREAT_merged_{project}_{tfactor}.pdf"
    shell:
        "pdfunite {input} {output}"


rule completePDF:
    input:
        cisTPDF =   "/home/petear/CistopicDir/Output/{project}/merged_{project}_{tfactor}.pdf",
        GreatPDF = "/home/petear/CistopicDir/Output/{project}/GREAT_merged_{project}_{tfactor}.pdf"
    output:
        "/home/petear/CistopicDir/Output/{project}/Complete_{project}_{tfactor}.pdf"
    shell:
        "pdfunite {input.cisTPDF} {input.GreatPDF} {output}"
