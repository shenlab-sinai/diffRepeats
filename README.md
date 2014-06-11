diffRepeats
===========

Quantify repeat element enrichment with next-generation sequencing data

Li Shen
Mount Sinai School of Medicine
New York, NY
June 2012


INTRODUCTION

This is diffRepeats, the little brother of diffReps. It is a program used to 
quantify repeat elements. It takes fastq files as inputs and uses BWA to perform 
short read alignment. It then parses the alignment results and summarize a count
table for each repeat element as a row and each fastq file as a column. Any 
differential analysis program can be used on this count table to find the 
significant elements.


PREREQUISITES

You must install BWA short read aligner and make sure "bwa" is in your PATH. 
diffRepeats also uses Samtools to read BAM files so you should install that too.
In order to enable parallel processing, diffRepeats requires Parallel::ForkManager 
which can be installed from CPAN. 


INSTALLATION

Installation is rather easy. Just copy diffRepeats.pl to one of your searchable 
directory such as ~/bin. 


REPEAT ELEMENTS

Repeat elements can be downloaded from Repbase (http://www.girinst.org/repbase/). 
Download them as fasta format for the species desired, such as: Mus_musculus.fa.
Then issue a command at console:

bwa index Mus_musculus.fa

to build an index on the fasta file. It will generate a few files which facilitate
BWA to find an alignment for each shrot read.


DISCUSSION

Differential analysis for repeat elements is unique and has to be taken care of 
specially. This is because each repeat has multiple copies on the reference genome.
Even though these copies may have evolved into slightly different versions with
variable regions, they are still highly similar to each other. It is therefore 
extremely difficult for a short read aligner to tell them apart. With the current 
read length of around 100bp, it is nearly impossible to quanfity one specific copy 
of a repeat element. We therefore take a less ambitious goal to quantify each 
repeat element without specifying where they exactly locate. This is achieved by 
aligning the short reads to a database of repeat sequences. These repeat sequences
constitute a so called "repetitive genome" which is used in sequence alignment just 
like a reference genome. 










