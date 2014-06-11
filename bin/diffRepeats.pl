#!/usr/bin/perl -w
#
#	diffRepeats: Generate count table for repeat elements, which
#	can be used for differential analysis.
#
#	- Li Shen, MSSM
#	Oct. 2011
#
#	Program arguments:
#	Input: RepBase fasta, fastq file names.
#	Output: BAM for each fastq file, count table.
#


use strict;
use Getopt::Long;
use Parallel::ForkManager;

# Program arguments.
my($repbase, @fq, $tbl);
my $mqual = 20;	# map quality requirement.
my $qillu = 0;	# quality encoding.
my $nproc = 1;	# number of threads.
##
if(@ARGV < 1 or !GetOptions("--repbase=s" => \$repbase,
							"--fq=s{1,}" => \@fq,
							"--tbl=s" => \$tbl,
							"--mqual:i" => \$mqual,
							"--qillu" => \$qillu,
							"--nproc:i" => \$nproc)){
	print "#### diffRepeats: Generate count table for repeat elements. ####\n";
	print "Usage: $0 --repbase repbase.fa --fq fq1 [fq2]...[fqN] --tbl output.txt\n";
	print "--repbase  Fasta represents DNA sequences of repeats downloaded from RepBase.\n";
	print "--fq       Short read libraries in fastq format.\n";
	print "--tbl      Output read count table.\n";
	print "Options:\n";
	print "--mqual    Map quality requirement (Default=20).\n";    
	print "--qillu    Use illumina [1.3, 1.8) quality encoding (Default=Sanger).\n";
	print "--nproc    Number of threads (Default=1).\n";
	print "################################################################\n";
	exit 0;
}

# Read fasta downloaded from RepBase and extract the repeat elements.
# Info to extract: name, type and origin.
open REPBASE, "<$repbase" or die "Read repbase fasta error: $!\n";
my %rep_tbl;
while(<REPBASE>){
	chomp;
	next unless /^\>/;	# ignore DNA sequences.
	my($name,$type,$origin) = split /\t/;
	if($name =~ /^\>(.+)$/){
		$name = $1;
	}else{
		warn "Misformatted repeat name for line: $0\n";
		next;
	}
	$rep_tbl{$name} = {"type" => $type, "origin" => $origin};
}
close REPBASE;

# Find eligible sequence files.
my @fqe;
my $sux;
foreach my $q(@fq){
	if($q =~ /(.+)\.(txt|fq|fastq)/){
		push @fqe, $1;
		$sux = $2;
	}else{
		warn "Unrecognized sequence file name: $q. Ignored.\n";
	}
}

# Align fastq files to the repeat fasta and generate BAM files.
# bwa aln [-I] repbase.fa reads.fq|bwa samse repbase.fa - reads.fq|samtools view -b -S -o file -q $mqual
my $pm = new Parallel::ForkManager($nproc);
foreach my $q(@fqe){
	# Determine the basename.
	my $infq = $q . '.' . $sux;
	my $outbam = $q . '.bam';
	my $qstr = $qillu ? '-I' : '';
	my $cmd = "bwa aln $qstr $repbase $infq|bwa samse $repbase - $infq|samtools view -b -S -o $outbam -q $mqual -";
	my $pid = $pm->start and next;
	system($cmd);
	$pm->finish;
}
$pm->wait_all_children;

# Read BAM and fill in read counts for the repeat elements.
# Setup counters for all entries in the repeat table.
while(my($name,$entry) = each %rep_tbl){
	my @z;
	for(my $i=0; $i < @fqe; $i++){
		push @z, 0;
	}
	$rep_tbl{$name}->{'count'} = \@z;
}
my $i = 0;	# iterator for eligible sequence file.
foreach my $q(@fqe){
	# Determine the basename.
	my $outbam = $q . '.bam';
	open BAM, "samtools view $outbam |" or die "Open bam for streaming error:$!\n";
	while(<BAM>){
		chomp;
		my %aln_struct = &parse_sam_line($_);
		next if $aln_struct{'mapq'} < $mqual or !$aln_struct{'unique'};
		$rep_tbl{$aln_struct{'rname'}}->{'count'}[$i]++;
	}
	$i++;
}

sub parse_sam_line{
	my($line) = @_;
	my($qname,$flag,$rname,$pos,$mapq,$cigar,$mrnm,$mpos,$isize,$seq,$qual,$opt) = split /\t/, $line;
	my @xt = split /:/, $opt;
	my $uniq = $xt[2] eq 'U'? 1 : 0;
	return ('rname' => $rname, 
			'mapq' => $mapq, 
			'unique' => $uniq);
}

# Output count table.
open CTBL, ">$tbl" or die "Open count table error:$!\n";
print CTBL join("\t", 'Name', 'Type', 'Origin', @fqe), "\n";
while(my($name,$entry) = each %rep_tbl){
	print CTBL join("\t", $name, $entry->{'type'}, $entry->{'origin'}, 
		join("\t", @{$entry->{'count'}})), "\n";
}
close CTBL;



#### Example RepBase fasta lines. ####
#>(A)n	Simple Repeat	Eukaryota
#aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa
#>(AC)n	Simple Repeat	Eukaryota
#acacacacacacacacacacacacacacacacacacacacacacacacacacacacacacacacacaca
#>(AG)n	Simple Repeat	Eukaryota
#agagagagagagagagagagagagagagagagagagagagagagagagagagagagagagagagagaga
#>(AT)n	Simple Repeat	Eukaryota
#atatatatatatatatatatatatatatatatatatatatatatatatatatatatatatatatatata
#>(C)n	Simple Repeat	Eukaryota
#ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
#######################################
