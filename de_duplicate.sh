#!/usr/bin/env bash

SAMP=$(cat)

awk 'BEGIN {i = 1;} { if ($1 ~ /^>/) { tmp = h[i]; h[i] = $1; } else if (!a[$1]) { s[i] = $1; a[$1] = "1"; i++; } else { h[i] = tmp; } } END { for (j = 1; j < i; j++) { print h[j]; print s[j]; } }' <''$SAMP''_library.fasta> geneLibrary/''$SAMP''.fasta