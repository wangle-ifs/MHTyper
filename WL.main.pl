#!/usr/bin/perl -w
=head1 Program Description

It works for ilmn or ion sequencing platform.

Microhaplotype Parse

=head1 Contact & Version

  Author:  wangle_02@163.com
  Version: 0.3,  Date: 2017-07-07

=head1 Command-line Option
	--bam		bam file from BWA MEM, binary format
	--sam		sam file from BWA MEM, txt format
	--rs		default= WL.rs.txt (provided by analyses tech)
	--amplicon	default= WL.amplicon.txt
	--dbSNP		default= WL.dbSNP.txt
	--outdir	directory of output files, default= WL
	--prefix	prefix of all output files, default= WL


=head1 Usage Exmples

  perl WL.main.pl --bam RM-P1.bam
		or
  perl WL.main.pl --sam RM-P1.sam

=cut

use strict;
use FindBin qw($Bin);
use Getopt::Long;
use File::Basename;
use List::Util;

my($help,$bam,$sam,$rs,$amplicon,$dbSNP,$outdir,$prefix)=("","","","$Bin/WL.rs.txt","$Bin/WL.amplicon.txt","$Bin/WL.dbSNP.txt","WL","WL");
GetOptions(
	"bam=s"		=>	\$bam,
	"sam=s"		=>	\$sam,
	"rs:s"		=>	\$rs,
	"amplicon:s"	=>	\$amplicon,
	"dbSNP:s"	=>	\$dbSNP,
	"outdir:s"	=>	\$outdir,
	"prefix:s"	=>	\$prefix,
	);

die `pod2text $0` if ($help || !$bam && !$sam);
system("mkdir -p $outdir");
my$novo="false";
my(%rs,%add);
open IN,$rs or die "$!\n";
<IN>;
while (<IN>) {
	chomp;
	my@array=split /\t/,$_;
	#mh14KK-101      ELK2B   rs28529526/rs10134526   2       14      106009477 106009572     95
	$array[4]="chr$array[4]";
	$array[5]=~s/\s+/\//g;
	my@rs=split /\//,$array[2];
	my@posi=split /\//,$array[5];
	$add{$array[0]}="$array[4]:$array[5]($array[2])\t$array[3]\t$array[1]\t$array[6]";
	for (my$i=0;$i<@rs;$i++) {
		$rs{$array[0]}{"$array[4]:$posi[$i]"}={amplicon=>$array[0],gene=>$array[1],rs=>$rs[$i],RS=>"$array[4]:$array[5]($array[2])",insert=>$array[6],};
	}
}
close IN;

my%amplicon;
open IN,$amplicon or die "$!\n";
<IN>;
while (<IN>) {
	chomp;
	my@array=split /\t/,$_;
	#mh01KK-002      1       ESRRG   216634338       216634360       216634472       216634449
	$array[1]="chr$array[1]";
	$amplicon{$array[1]}{$array[0]}={insertL=>$array[4], insertR=>$array[6], ampliconL=>$array[3], ampliconR=>$array[5],};
	foreach my$position (keys%{$rs{$array[0]}}) {
		my($chrom,$posi)=split /:/,$position;
		if (($posi-$array[4])*($posi-$array[6])>0) {
			print STDERR "$array[1]\t$posi\t$rs{$array[0]}{$position}->{'rs'}\tnot_in_amplicon\t$_\n";
		}
	}
}
close IN;

my%variants;
open IN,$dbSNP or die "$!\n";
while (<IN>) {
	chomp;
	next if /lost/;
	my@array=split;
	#22      19951207        19951207        C       G       rs4818
	$variants{"chr$array[0]:$array[1]"}={Ref=>$array[3],Alt=>$array[4],rs=>$array[5],no=>"chr$array[0]:$array[1]$array[3]=$array[3]($array[5])",};
	
}
close IN;

my%alignment;
if ($bam ne "") {
	open IN,"samtools view $bam|" or die "$!\n";
}
elsif ($sam ne "") {
	open IN,$sam or die "$!\n";
}
else {
	die "input file should be .bam or .sam!\n";
}
while (<IN>) {
	chomp;
	my@array=split;
	#M50185:14:000000000-D26EH:1:1101:10790:2694[0] 83[1] mh01KK-002[2] 98[3] 60[4] 113S114M2D24M[5] =[6] 101[7] -137[8] CAGTGTTGATTTTACACTTGTCATATACAGTGTTGTTGTTTTTTTTTTTTTTTTCAAGCAGAAGACGGCATACGAGATCGTGATGTGACTGGATCAGACGTGTGCTCTTCCGATCTTCTGGATAAGGGAGGAAGAAACTTTGCCTTTTGACCATTTCAACTATTATTATTATTACCCAGTTAGAAAGTGAATAAATGACCTAAGTGGGAAACCTGACATAGGTAGACATATTGGCTTCAGAACTAGAAGGC[9] //.;///9////////.///////99///9......9.---9---A-FDCFHGEHHHHGHHGHGGGHHHHGGHHGHHGFHBHHHHHHHHHHHHHGHGHGHHHHHHHGGGGHHHGHHGHHHHHHHHGHHHHHHHHHHHHHHHHGHHHHHHGHHHHHGHHHHHHHHHHHGHHHHHHHHHHHGHHHHHHHHHHHHHHHHGFHHHGB2FHHGFHGHHHHGGHHHHHHHHHHHHGFGGFFGGGGFFFFFFFBBBBB[10] NM:i:3[11] MD:Z:90A23^AT24[12] AS:i:125[13] XS:i:0[14]
	next unless $array[6] eq "=" && $array[5] ne "*" && $array[5]!~/H/;
	if ($array[5]=~/^(\d+)S/) {
		$array[3]+=$1;
		$array[9]=substr($array[9],$1);
		$array[5]=~s/^\d+S//;
	}
	if ($array[5]=~/(\d+)S$/) {
		$array[9]=substr($array[9],0,(length$array[9])-length($1));
		$array[5]=~s/(\d+)S$//;
	}
	@array[11..$#array]=sort@array[11..$#array];
	foreach my$amplicon (keys%{$amplicon{$array[2]}}) {
		if (($amplicon{$array[2]}{$amplicon}->{'ampliconR'}-$array[3]-0.5*(length$array[9]))*($amplicon{$array[2]}{$amplicon}->{'ampliconL'}-$array[3]-0.5*(length$array[9]))<=0) {
			$alignment{$amplicon}{$array[0]}{"$array[2] $array[3] $array[1]"}=join " ",@array[9,5,11..$#array];
			last;
		}
	}
}
close IN;

my%STAT;
foreach my$amplicon (keys%alignment) {
	foreach my$readsid (keys%{$alignment{$amplicon}}) {
		my@position=sort keys%{$alignment{$amplicon}{$readsid}};
		next unless@position==2;
		my$printout="";
		my(%varlist,%position);
		foreach my$posi (keys%{$rs{$amplicon}}) {
			$varlist{$posi}=$variants{$posi}->{'no'};
		}
		foreach my $position (@position) {
			#chrom posi flag TCTGGATAAGGGAGGAAGAAACTTTGCCTTTTGACCATTTCAACTATTATTATTATTACCCAGTTAGAAAGTGAATAAATGACCTAAATGAGAAACCTGACATAGGTAGACATATTGGCTTCAGAACTAGAAGGCAGATCGGAAGAGCACCACGTCTGAACTCCAGTCACATCACGATCTCGTATGCCGTCTTCTGCTTGAAAAAAAAAAACAATAATATCACTTTCTCTCCTAACACTTAATACCCC 
			#111M2D26M AS:i:124 MD:Z:90G20^AT26 NM:i:3 XS:i:0
			#35M1I14M AS:i:37 MD:Z:17G31 NM:i:2 XS:i:19
			my@array=split / /,"$position $alignment{$amplicon}{$readsid}{$position}";
			$alignment{$amplicon}{$readsid}{$position}=~/MD:Z:([^ ]+)/;
			my$MD=$1;
			$array[4]=~s/([MDI])(\d+)/$1 $2/g;
			my@ciga=split / /,$array[4];
			my($newQuery,$newPosi,$delPosi)=("",0,0);
			my%insertion;
			foreach my$ciga (@ciga) {
				$ciga=~/(\d+)([MID])/;
				my($num,$tag)=($1,$2);
				if ($tag eq "I") {
					$insertion{$newPosi+$delPosi+$array[1]}=substr($array[3],$newPosi,$num);
					$insertion{$newPosi+$delPosi+$array[1]}="->$insertion{$newPosi+$delPosi+$array[1]}";
					$newPosi+=$num;
				}
				elsif ($tag eq "D") {
					$delPosi+=$num;
				}
				else {
					$newQuery.=substr($array[3],$newPosi,$num);
					$newPosi+=$num;
				}
			}
			for (my$i=$array[1];$i<$array[1]+$newPosi;$i++) {
				$position{"$array[0]:$i"}++;
			}
			my(%snp,%deletion);
			($newPosi,$delPosi)=(0,0);
			while ($MD ne "") {
				if ($MD=~/^(\d+)/) {
					$newPosi+=$1;
					$MD=~s/^\d+//;
				}
				elsif ($MD=~/^([ATCG]+)/) {
					$snp{$newPosi+$delPosi+$array[1]}="$1>";
					$snp{$newPosi+$delPosi+$array[1]}.=substr($newQuery,$newPosi,length($1));
					$newPosi+=length$1;
					$MD=~s/^([ATCG]+)//;
				}
				elsif ($MD=~/^\^([ATCG]+)/) {
					$deletion{$newPosi+$delPosi+$array[1]}="$1>-";
					$delPosi+=length$1;
					$MD=~s/^\^([ATCG]+)//;
				}
			}
			my$variants="";
			foreach my$p (sort keys%snp) {
				$snp{$p}=~/([ATCG-]+)>([ATCG-]+)/;
				if (exists$variants{"$array[0]:$p"} && $variants{"$array[0]:$p"}->{'Ref'} eq $1 && $variants{"$array[0]:$p"}->{'Alt'} eq $2) {
					$varlist{"$array[0]:$p"}="$array[0]:$p$snp{$p}"."(".$variants{"$array[0]:$p"}->{'rs'}.")";
				}
				else {
					$variants.="$array[0]:$p$snp{$p},";
				}
			}
			foreach my$p (sort keys%insertion) {
				$insertion{$p}=~/([ATCG-]+)>([ATCG-]+)/;
				if (exists$variants{"$array[0]:$p"} && $variants{"$array[0]:$p"}->{'Ref'} eq $1 && $variants{"$array[0]:$p"}->{'Alt'} eq $2) {
					$varlist{"$array[0]:$p"}="$array[0]:$p$insertion{$p}"."(".$variants{"$array[0]:$p"}->{'rs'}.")";
				}
				else {
					$variants.="$array[0]:$p$insertion{$p},";
				}
			}
			foreach my$p (sort keys%deletion) {
				$deletion{$p}=~/([ATCG-]+)>([ATCG-]+)/;
				if (exists$variants{"$array[0]:$p"} && $variants{"$array[0]:$p"}->{'Ref'} eq $1 && $variants{"$array[0]:$p"}->{'Alt'} eq $2) {
					$varlist{"$array[0]:$p"}="$array[0]:$p$deletion{$p}"."(".$variants{"$array[0]:$p"}->{'rs'}.")";
				}
				else {
					$variants.="$array[0]:$p$deletion{$p},";
				}
			}
			$printout.=$variants if $variants ne "";
		}

		my@printout=split /,/,$printout;
		$printout="";

		my$microhap="";
		foreach my$posi (sort keys%varlist) {
			if (exists$position{$posi}) {
				$printout.="$varlist{$posi},";
				$microhap.=$1 if $varlist{$posi}=~/[=>]([ATCG-]+)/;
			}
			else {
				$varlist{$posi}=~s/=[ATCG-]+/>N/;
				$printout.="$varlist{$posi},";
				$microhap.="N";
			}
		}
		if ($novo eq "true") {
			my%uniqVar;
			map{$uniqVar{$_}++}@printout;
			foreach my$uniqVar (sort keys%uniqVar) {
				$printout.="$uniqVar,";
			}
		}
		$printout="NA" if $printout eq "";
		$printout=~s/,$//;
		next if $microhap=~/N/;
		$STAT{$amplicon}{$microhap}{'cnt'}++;
		$STAT{$amplicon}{$microhap}{'detail'}=$printout;
		$STAT{$amplicon}{$microhap}{'readsid'}.="$readsid," if $STAT{$amplicon}{$microhap}{'cnt'}<50;
	}
}

my%plot;
open OUT,">$outdir/$prefix.rs.txt" or die "$!\n";
print OUT "Amplicon\tReadsCnt\tMicrohap\tDetail\trsID\trsNum\tGene\tinsertSize\n";
foreach my$amplicon (sort %add) {
	if (!exists$STAT{$amplicon}) {
		print OUT "$amplicon\t0\tNA\tNA\t$add{$amplicon}\n" if exists$add{$amplicon};
		next;
	}
	my$max=0;
	foreach my$microhap (sort{$STAT{$amplicon}{$b}{'cnt'}<=>$STAT{$amplicon}{$a}{'cnt'}}keys%{$STAT{$amplicon}}) {
		$max=$STAT{$amplicon}{$microhap}{'cnt'};
		last;
	}
	foreach my$microhap (sort keys%{$STAT{$amplicon}}) {
		next unless $STAT{$amplicon}{$microhap}{'cnt'}>=$max*0.01 && $STAT{$amplicon}{$microhap}{'cnt'}>10;
		$plot{$amplicon}{$microhap}=$STAT{$amplicon}{$microhap}{'cnt'};
		print OUT "$amplicon\t$STAT{$amplicon}{$microhap}{'cnt'}\t\"$microhap\"\t$STAT{$amplicon}{$microhap}{'detail'}\t$add{$amplicon}\n";
		#print STDERR "$amplicon\t$STAT{$amplicon}{$microhap}{'cnt'}\t\"$microhap\"\t$STAT{$amplicon}{$microhap}{'detail'}\t$STAT{$amplicon}{$microhap}{'readsid'}\n";
	}
}
close OUT;

open OUT,">$outdir/$prefix.all.R" or die "$!\n";
open A,">$outdir/$prefix.amplicon.R" or die "$!\n";
my$cnt=0;
foreach my$amplicon (sort keys%plot) {
	if ($cnt%30==0) {
		my$subjpeg=int($cnt/30);
		print OUT "jpeg(filename=\"$outdir/$prefix.all-$subjpeg.jpeg\", width=3000, height=3500, quality=100)\n";
		print OUT "par(mfrow=c(6,5), oma=c(4,10,4,4), mar=c(5,10,5,3))\n";
	}
	$cnt++;
	my($value,$labels,$col,$xlim)=("","","",10);
	my@microhap=sort keys%{$plot{$amplicon}};
	my@value=values%{$plot{$amplicon}};
	@value=sort{$b<=>$a}@value;
	$xlim=2*@microhap if @microhap>5;
	foreach my$microhap (sort keys%{$plot{$amplicon}}) {
		if ($plot{$amplicon}{$microhap}<50) {
			$col.="\"red\",";
		}
		else {
			$col.="\"blue\",";
		}
		$labels.="\"$microhap\",";
		$value.="$plot{$amplicon}{$microhap},";
	}
	chop$col;
	chop$labels;
	chop$value;
	print OUT "text(barplot(c($value), xlim=c(1,$xlim), ylim=c(0,1.2*$value[0]), lwd=3, cex.axis=2, cex.names=2, cex.lab=2, names.arg=c($labels), font.lab=2, ylab=\"Reads Count\", space=1, col=c($col), main=\"$amplicon\", cex.main=4), c($value), labels=c($value), cex=c(2,2), font=10, pos=3)\n";
	print OUT "axis(side=1, labels=FALSE)\n";
	print OUT "dev.off()\n" if $cnt%29==0 && $cnt!=0;

	print A "jpeg(filename=\"$outdir/$prefix.$amplicon.jpeg\", width=500, height=500, quality=100)\n";
        print A "par(mar=c(5,7,5,3))\n";
	print A "text(barplot(c($value), xlim=c(1,$xlim), ylim=c(0,1.2*$value[0]), lwd=3, cex.axis=2, cex.names=2, cex.lab=2, names.arg=c($labels), font.lab=1, ylab=\"Reads Count\", space=1, col=c($col), main=\"$amplicon\", cex.main=2), c($value), labels=c($value), cex=c(2,2), font=10, pos=3)\n";
	print A "axis(side=1, labels=FALSE)\n";
	print A "dev.off()\n";
}
close OUT;
close A;

system("R CMD BATCH $outdir/$prefix.all.R");
system("R CMD BATCH $outdir/$prefix.amplicon.R");
__END__

