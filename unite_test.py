import unittest
import sys
import csv
import requests
import vcf_test


# Some initialization code if needed
vcf_file = sys.argv[1]
parsed_vcf = vcf_test.vcf_parser(vcf_file)

class test_vcf_test(unittest.TestCase):
    # Check the record counts
    def test_count_records(self):
        count = len(parsed_vcf)
        self.assertEqual(count, 11801)
    # Check if snp Hgvs notation is correct
    def test_snp_hgvs_notation(self):
        chrom = parsed_vcf[0]['chrom']
        pos = parsed_vcf[0]['pos']
        ref = parsed_vcf[0]['ref']
        alt = parsed_vcf[0]['alt']
        hgvs_notation = vcf_test.construct_hgvs_notation(chrom, pos, ref, alt)
        self.assertEqual(hgvs_notation, '1:g.1158631A>G')
    # Check if delins Hgvcs notation is correct
    def test_delins_hgvs_notation(self):
        chrom = parsed_vcf[9]['chrom']
        pos = parsed_vcf[9]['pos']
        ref = parsed_vcf[9]['ref']
        alt = parsed_vcf[9]['alt']
        hgvs_notation = vcf_test.construct_hgvs_notation(chrom, pos, ref, alt)
        self.assertEqual(hgvs_notation, '1:g.1647722_1647730delinsTCTAGGATG')
    # Check if insertions Hgvs notation is correct
    def test_insertion_hgvs_notation(self):
        chrom = parsed_vcf[12]['chrom']
        pos = parsed_vcf[12]['pos']
        ref = parsed_vcf[12]['ref']
        alt = parsed_vcf[12]['alt']
        hgvs_notation = vcf_test.construct_hgvs_notation(chrom, pos, ref, alt)
        self.assertEqual(hgvs_notation, '1:g.1647893_1647894insTTTCTT')
    # Check if response gene name for the first entry is correct
    def test_response_for_gene_name(self):
        whole_data_for_first = vcf_test.ensemble_request([parsed_vcf[0]])
        self.assertEqual(whole_data_for_first[0]['gene_name'], 'SDF4')
    # Check if response minor allele frequency for the first entry is correct
    def test_response_for_minor_allele(self):
        whole_data_for_first = vcf_test.ensemble_request([parsed_vcf[1]])
        self.assertEqual(whole_data_for_first[0]['minor_allele_frequency'], 0.0481)

    # Check if response minor rs_id for the first entry is correct
    def test_response_for_minor_allele(self):
            whole_data_for_first = vcf_test.ensemble_request([parsed_vcf[3]])
            self.assertEqual(whole_data_for_first[0]['id'], 'rs307348')



if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(test_vcf_test)
    unittest.TextTestRunner().run(suite)