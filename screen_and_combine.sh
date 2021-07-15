#!/usr/bin/env bash

mkdir final_contigs

SAMP=$(cat)

bioawk -c fastx '{ if(length($seq)>499) {print ">" $name; print $seq }}' first_contigs/""$SAMP""_assembly/final.contigs.fa > $SAMP""_first.fasta
bioawk -c fastx '{ if(length($seq)>499) {print ">" $name; print $seq }}' second_contigs/""$SAMP""_assembly/final.contigs.fa > $SAMP""_second.fasta
cat $SAMP""_first.fasta $SAMP""_second.fasta > final_contigs/""$SAMP""_final.fasta
