#!/bin/bash

if [ $# -ne 6 ]
then
    echo -e "\nusage: `basename $0` <genome> <fasta.in> <min_length> <flank> <end_left> <end_right>\n"
    echo -e "DESCRIPTION: This script takes a fasta sequence <fasta.in>, blasts it to the genome, recovers locations with alingment length > <min_length>, prints them as bed file, extends bed coordinates <flank> bases in each direction, and makes fasta from that BED.\n"

    echo -e "INPUT:     <genome>      location of genome in fasta format"
    echo -e "           <fasta.in>    query sequence if fasta format"
    echo -e "           <min_length>  min length of the blast hit. If set to 0, min length = (length of query)/2"
    echo -e "           <flank>       number of bases to extend the genome coordinates of the matched locus\n"
    
    echo -e "OUTPUT:    produces a <fasta.in.bed> which is the blast results BED file file"
    echo -e "           produces a <fasta.in.blast.flank.bed> which is the extended BED"
    echo -e "           produces a <fasta.in.blast.flank.bed.fa> bedtools getfasta \n"
    
    echo -e "REQUIRES:  Installation of BEDTOOLS, SAMTOOLS and BLAST+ accessible from PATH\n"

    
     
    exit
fi

genome=$1
fasta_in=$2
out=`basename $2`
# out=`echo "$fasta_in" | cut -d'.' -f1`
min_length=$3
flank_l=$4
flank_r=$4
end_left=$5
end_right=$6


# look for blast database, if not found, make it
FILE=$genome.nin
if [ ! -f "$FILE" ]; then
    echo -e "\n\nBlast database doesn't exist. Running makeblastdb, this can take some time\n\n"
    makeblastdb -in $genome -dbtype nucl
    echo -e "\n\n"
fi

# look for genome.fasta.length file, if not found, make it
FILE=$genome.length
if [ ! -f "$FILE" ]; then
    echo -e "\n\nFile with genome lengths not found. Making it now, this can take some time"
    samtools faidx $genome
    awk '{OFS="\t"; print $1,$2}' < $genome.fai > $genome.length
	echo -e "\n\n"
fi



# if the value 0 is entered as min_length by the user, redifine min_length as half the length of  fasta.in 
if [ $3 == 0 ]
then 
	min_length=`grep -v ">" $2 | wc | awk '{print $3/2}'`
else
	min_length=$3
fi


echo "query sequence" $out
echo "minimum length of blast hit = " $min_length
echo "the hit locus will be extended" = $flank "bases in each direction"

# run blast # evalue cutoff is 1e-20 - this can be changed to a lower value ( 1e-10 for example) to capture less similar hits 
echo "#qseqid sseqid pident length mismatch qstart qend sstart send sstrand" > $out.blast.o # this step writes the column names in the out put file
blastn -query $fasta_in -db $genome -outfmt "6 qseqid sseqid pident length mismatch qstart qend sstart send sstrand" -evalue 1e-20 | awk -v "ml=$min_length" '{OFS="\t"; if ($4 > ml) {print $0}}' >> $out.blast.o

# parse blast result into a bed file. 
awk '{OFS="\t"; if ($1~/^\#/) {} else { if ($10~/plus/) {print $2, $8, $9, $1, $3, "+"} else {print $2, $9, $8, $1, $3, "-"}}}' < $out.blast.o > $out.blast.bed

# extend boundaries "flank" bases up and down of the blast hit locations
# Adding 300 on completed ends to allow the main script to find the edges again.

if [ "$end_left" = TRUE ] ; then
	flank_l=300
fi	
if [ "$end_right" = TRUE ] ; then
	flank_r=300	
fi

bedtools slop -s  -i $out.blast.bed  -g $genome.length -l $flank_l -r $flank_r > $out.blast.flank.bed


# extract fasta sequence from the reference genome

bedtools getfasta -fi  $genome -fo $out.blast.bed.fa  -bed $out.blast.flank.bed -s 

# how many sequences ended up in the multi fasta 
fasta_count=`grep -c ">" $out.blast.bed.fa`

echo "the fasta has "$fasta_count " sequences"

# make alignment. If multiple processors are available, `--threads n` can be increases by changing c to a higher number. For a standard laptop computer this could be 2 or 3. 

mafft --reorder  --thread -1 $out.blast.bed.fa > $out.maf.fa


# remove redundant files
rm $out.blast.o $out.blast.bed *.blast.flank.bed
