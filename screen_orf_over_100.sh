#!/usr/bin/env bash

SAMP=$(cat)

bioawk -c fastx '{ if(length($seq)>99) {print ">" $name; print $seq }}' prodigal_inital/""$SAMP"".fna > orf_over_100/""$SAMP""_over_100.fasta

