Description
===========
Intronarrator is a collection of Python and shell scripts that predicts introns
using RNA-seq data, then uses the program `AUGUSTUS
<https://github.com/Gaius-Augustus/Augustus>`_ to do gene prediction, using its
"intronless" model.

The motivation for this is that AUGUSTUS's intron model leads to poor prediction
of extremely short (mostly 15 bp) introns of heterotrichous ciliates, despite
having an intron length model designed to accommodate extreme variation in
intron lengths observed across eukaryotes. After trying a number of options to
alter this model, we were unable to accurately predict introns - leading, more
often than not, to predicted introns being incorrect ("faketrons"). RNA-seq
data, on the other hand, may be very extensive -- providing deep coverage of
most genes. Thus, it is possible to predict the introns directly from RNA-seq
data, and leave the rest of the gene prediction to AUGUSTUS.

The approach is illustrated below (semi-transparent green indicates
where a gene to be predicted is located): |approach| 

.. |approach| image:: images/intronarrator_approach.png

Dependencies and installation
=============================
The following programs need to be installed: AUGUSTUS, Infernal, tRNAscan-SE
version 2.0.

Code can be obtained from github:
``git clone https://github.com/Swart-lab/Intronarrator/``

An Intronarrator evironment with all the necessary depencies can be installed from the included env.yml with conda::

        cd Intronarrator
        conda env create -f env.yml
        conda activate Intronarrator

Usage
=====
intronarrator.sh needs to be copied to the working directory. After this it
needs to be edited and properly setup before running, including setting
the paths to AUGUSTUS configuration files. You should also decide on the number
of parallel processes to run.

Helper scripts for producing stranded BAM files, and generation of "hints" for
AUGUSTUS from the BAM files are included in the directory "helper_scripts".
These should be copied and modified according to need before use.

TODO and known limitations
==========================
Reads with multiple introns currently are ignored - these still need to be
properly accommodated. Therefore be warned: this software is currently
inappropriate for organisms with high densities of introns and frequent
alternative splicing (e.g. most metazoans)!
