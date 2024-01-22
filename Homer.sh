bedFile=$1

# Determine hg version based on file name
if [[ $bedFile == *"HG38.bed" ]]; then
    hg="hg38"
else
    hg="hg19"
fi

annotatedFile=${bedFile%.*}_annotatedPeaks.txt
genesFile=${bedFile%.*}_genes.txt
motifsDir=${bedFile%.*}_Motifs

# Annotate peaks with Homer
echo "Annotating peaks with Homer"
annotatePeaks.pl $bedFile $hg > $annotatedFile

# Extract genes from annotated file
echo "Extracting genes from annotated file"
cut -f4 $annotatedFile | sort | uniq > $genesFile

# Find motifs in bed file
echo "Finding motifs in bed file"
findMotifsGenome.pl $bedFile $hg $motifsDir -size 200 -mask

# Find motifs in genes
echo "Finding motifs in genes"
findMotifs.pl $genesFile $hg $motifsDir -size 200 -mask
