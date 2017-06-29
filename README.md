
# This file is part of table_maker.py

# Copyright Institut Curie 2014

# This software is a computer program whose purpose is to MaxEntScan scores.

# This software is governed by the CeCILL license under French law and
# abiding by the rules of distribution of free software.  You can  use, 
# modify and/ or redistribute the software under the terms of the CeCILL
# license as circulated by CEA, CNRS and INRIA at the following URL
# "http://www.cecill.info". 

# As a counterpart to the access to the source code and  rights to copy,
# modify and redistribute granted by the license, users are provided only
# with a limited warranty  and the software's author,  the holder of the
# economic rights,  and the successive licensors  have only  limited
# liability. 

# In this respect, the user's attention is drawn to the risks associated
# with loading,  using,  modifying and/or developing or reproducing the
# software by the user in light of its specific status of free software,
# that may mean  that it is complicated to manipulate,  and  that  also
# therefore means  that it is reserved for developers  and  experienced
# professionals having in-depth computer knowledge. Users are therefore
# encouraged to load and test the software's suitability as regards their
# requirements in conditions enabling the security of their systems and/or 
# data to be ensured and,  more generally, to use and operate it in the 
# same conditions as regards security. 

# The fact that you are presently reading this means that you have had
# knowledge of the CeCILL license and that you accept its terms.

# clinTools

## Requirement

This toolkit work with python 2

To use it you must install:
  * samtools version 1.1
  * hgvs python module version 0.3.7
  * pysam python module version 0.8.3

You can use the pip module manager for example in a virtual env

```
virtualenv clintools
source clintools/bin/activate
pip install hgvs==0.3.7 pysam==0.8.3
```

Don't forget to update your PATH environnement variable to use the appropriate samtools version.

## Table maker

This script and various annotation to variant files annotated with annovar.

#### Launching

This script require 5 input :

* A configuration file which is a tabulated file with 3 field per line : sample_id, annotated_file, bam_file
* Another configuration file containing an list of selected transcript NM id (one per line)
* The gtf of refGene mRNA.
* The fasta sequence of the reference genome
* A configuration file containing the association between chromosome as reference in the fasta file and the corresponding NC id.

```
cd examples/table_maker/
../../table_maker.py -Q 20 -q 20 example.conf prefered_nm.conf hg19_mRNA.gtf hg19.fa chr_accessions_hg19_GRCh37.p13.tsv > example_table_report.tsv
```


All annotated file must contain the following mandatory fields : "Start", "End", "Chr", "Otherinfo", "Ref", "Alt", "Gene.refGene", "Func.refGene", "ExonicFunc.refGene", "AAChange.refGene"


#### Output file description

* Barcode : the sample_id from the firs configuration file
* Gene : Gene name from annovar "Gene.refGene" annotation
* Chr, start, end, ref , alt : informations extracted from vcf (annovar "Otherinfo" field)
* all not mandatory fields from annovar
* NM : the nm identifier of the impacted transcript
* cDNAchange : hgvs c. annotation
* AAchange : hgvs p. annotation
* all_hgvs : all hgvs annotation reported by annovar
* Status : "Heterozygous" if the allelic ratio is below 75 % and "Homozygous" otherwise.
* Depth, allelic ratio, number of reads supporting reference, alternative allele and all bases
* Strand bias for each read counts
* Grantham score : a score representing the physico-chemical effect of an amino acid substitution (Grantham 1974)
* can_splice_distance : the distance from the nearest canonic splice site
* MES_ref, MES_alt, MES_delta : maxentscan score (splicing) for reference sequence, variant sequence and ratio between 2 scores
* exon : the impacted or the nearest exon
* maxentscan donor and acceptor scores for reference sequence and variant sequence at variant position
* Base_around : the sequence around variant (useful to see stretch)
* minimum and maximum distances between variant position and start / end for reads


## References:

Please cite the following articles:


* Ke S, Shang S, Kalachikov SM, et al. Quantitative evaluation of all hexamers as exonic splicing elements. Genome Research. 2011;21(8):1360-1374. doi:10.1101/gr.119628.110.


* Grantham R. Amino acid difference formula to help explain protein evolution. Science. 1974 Sep 6;185(4154):862-4.

* Yeo G, Burge CB. Maximum entropy modeling of short sequence motifs with applications to RNA splicing signals. J Comput Biol. 2004;11(2-3):377-94.

