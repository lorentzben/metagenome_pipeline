#!/usr/bin/env bash

mkdir final_contigs

SAMP=$(cat)

bioawk -c fastx '{ if(length($seq)>499) {print ">" $name; print $seq }}' first_contings_""$SAMP""_$
bioawk -c fastx '{ if(length($seq)>499) {print ">" $name; print $seq }}' second_contings_""$SAMP""$
cat $SAMP""_first.fasta $SAMP""_second.fasta > final_contigs/""$SAMP""_final.fasta
