# Vcf_anotation

The script **vcf_test.py** consists of 4 functions. It is able to parse input vcf file using the **vcf_parser** function. This function takes a vcf file as an input parameter, opens using vcf module and iterates through records. In order to have a different alternative allele as a different entry in the output csv file, added an inner loop which iterates through zipped alternative allele and TR (total reads) score for each record. The output of this function is a list of dicts and consists of chrom, pos, ref, alt, depth, alt_reads, percent_alt_reads and percent_ref_reads infos for each record in vcf. 

The second function is **construct_hgvs_notation** which takes chrom, pos, ref and alt information for each variant and returns HGVS notations for different type of mutations. 

The third function **ensemble_request** takes output of vcf_parser iterates through records, constructs hgvs notation for each record and sends a GET request to VEP HGVS endpoint, parses response json and
updates each record data with gene, variant_effect, minor_allele, minor_allele_frequency, somatic, id fields.

The fourth function **response_json_parser** takes json response as an input and parses it to get specific fields. It is ran in ensemble_request function. 

Whole data which is an output of ensemble_request function is exported in a csv file.


 ## Prerequisites

In order to run the code, please make sure you have installed **PyVCF** and **requests** modules.

## Usage

You can run vcf_test.py using command: **python vcf_test.py /path/to/your/vcf/file /path/to/your/output/csv/file**

