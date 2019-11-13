# MHTyper
 microhaplotype allele-calling pipeline for use with next generation sequencing data

It works for ilmn or ion sequencing platform.

Microhaplotype Parse

Contact & Version

  Author:  wangle_02@163.com
  Version: 0.3,  Date: 2017-07-07

Command-line Option
        --bam           bam file from BWA MEM, binary format
        --sam           sam file from BWA MEM, txt format
        --rs            default= WL.rs.txt (provided by analyses tech)
        --amplicon      default= WL.amplicon.txt
        --dbSNP         default= WL.dbSNP.txt
        --outdir        directory of output files, default= WL
        --prefix        prefix of all output files, default= WL


Usage Exmples

  perl WL.main.pl --bam RM-P1.bam
                or
  perl WL.main.pl --sam RM-P1.sam
