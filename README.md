====================================================================
		    
 invMap                             2023/09/08
 ver1.0.0                           Baoji, China
		    	 	       
====================================================================

System Requirements
===================

Software Requirements
---------------------
The invMap is supported on the following operating systems 
   
   Linux (recommended)
   Windows 7
   Windows 10 

Hardware Requirements/Recommendations
-------------------------------------
   Intel Core i3 2100 or later for Windows 
   Color display capable of 1024 X 768 pixel resolution
   

Memory Requirements/Recommendations
-------------------------------------
invMap  requires approximately
800 M of available disk space on the drive or partition.

invMap requires a minimum of
16G of RAM for mapping sequence datasets to large reference genomes.


User's Guide
=================

Input Sequences Requirements
----------------------------
Input file requires fasta or fastq format

Execute Step
------------
Step 1: compile codes using make command
Step 2: ./invmap then come out help information and options explanation


A running example
-----------

./invMap genome.fa reads.fa > invMap_aligned.sam

The final mapping result is invMap_aligned.sam in SAM format.


Copyright Notice
===================
Software, documentation and related materials:
Copyright (c) 2023-2027 
Baoji University of Arts and Sciences,China
All rights reserved.
