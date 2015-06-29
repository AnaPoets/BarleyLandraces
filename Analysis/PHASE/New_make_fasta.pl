#!/usr/bin/perl -wpl

#Convert file into FASTA format
s/(\w+_2)//ig;
s/^(CIho\d+|PI\d+|WBDC\d+)/\#$1\n/ig;
s/\t//g;

