module_detection
================

Discover modules --- simple and fast

to run the program:

  python main.py motif.txt Promoter.fa


#################################################
	
it will take your Motif file, search against your Promoter file, and then all the putative modules will be listed in the output.

All result file is in the "result" folder

A sample output is like:

AAGCGCAA[1,0,0] 26...28 AAAGAAA[1,0,0]  total hits:  273
AAGCGCAA[1,0,0] 13...15 CCCTTTC[1,0,0]  total hits:  282
AAGCGCAA[1,0,0] 42...44 CCCTTTC[1,0,0]  total hits:  279
AAAGAAA[1,0,0] 1...2 AAAGAAA[1,0,0]  total hits:  490
CCCTTTC[1,0,0] 14...16 AAGCGCAA[1,0,0]  total hits:  85
CCCTTTC[1,0,0] 5...7 AAAGAAA[1,0,0]  total hits:  427
CCCTTTC[1,0,0] 26...28 CCCTTTC[1,0,0]  total hits:  355
AAGCGCAA[1,0,0] 13...15 CCCTTTC[1,0,0] 5...7 AAAGAAA[1,0,0]  total hits:  262
AAGCGCAA[1,0,0] 42...44 CCCTTTC[1,0,0] 26...28 CCCTTTC[1,0,0]  total hits:  206
AAAGAAA[1,0,0] 1...2 AAAGAAA[1,0,0] 1...2 AAAGAAA[1,0,0]  total hits:  66
CCCTTTC[1,0,0] 14...16 AAGCGCAA[1,0,0] 42...44 CCCTTTC[1,0,0]  total hits:  76
CCCTTTC[1,0,0] 14...16 AAGCGCAA[1,0,0] 42...44 CCCTTTC[1,0,0] 26...28 CCCTTTC[1,0,0]  total hits:  63

AAGCGCAA[1,0,0] is the motif, which allows for 1 mismatch
26...28 is the distance between first motif and next motif
total hits is the total number of occurence of the module
