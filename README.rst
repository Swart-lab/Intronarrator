Description
===========

Intronarrator is a set of Python scripts that predict introns in genes using RNA-seq data. The program
AUGUSTUS is used to do the actual gene prediction, using its "intronless" model. Intronarrator is driven 
by a main bash script "intronarrator.sh", 

Dependencies
============

The following programs need to be installed: AUGUSTUS, Infernal, tRNAscan-SE
version 2.0.

The following Python dependencies need to be installed: Numpy, BioPython, Pysam
if not present already.

Usage
=====

intronarrator.sh needs to be copied to the working directory. After this it
needs to be edited and properly setup before running, including setting
the paths to AUGUSTUS configuration files. You should also decide on the number
of parallel processes to run.

Helper scripts for producing stranded BAM files, and generation of "hints" for
AUGUSTUS from the BAM files are included in the directory "helper_scripts".
These should be copied and modified according to need before use.

