import vcf
import requests
import csv
import sys

# Function for parsing vcf file
def vcf_parser(file_path):
    vcf_records = vcf.Reader(open(file_path, 'r'))
    parsed_vcf_data = []
    # Iterate through records and collect info
    for record in vcf_records:
        for alt_allele, total_reads in zip(record.ALT, record.INFO['TR']):
            # Calculation of alternative and reference reads percent with total coverage (TC) and total reads scores
            total_coverage = record.INFO['TC']
            percent_alt_reads = round((total_reads / total_coverage) * 100, 2)
            percent_ref_reads = round(100 - percent_alt_reads, 2)
            entry = {
                'chrom': record.CHROM,
                'pos': record.POS,
                'ref': record.REF,
                'alt': str(alt_allele),
                'depth': total_coverage,
                'alt_reads': total_reads,
                'percent_alt_reads': percent_alt_reads,
                'percent_ref_reads': percent_ref_reads

            }

            # Append parsed_data list with parsed information from vcf file
            parsed_vcf_data.append(entry)

    return parsed_vcf_data

# Function for constructing HGVS notation
def construct_hgvs_notation(chrom, pos, ref ,alt):
    # When record is SNP
    if len(ref) == 1 and len(alt) == 1:
        return f"{chrom}:g.{pos}{ref}>{alt}"
    # When record is deletion
    elif len(ref) > 1 and len(alt) == 1:
        return f"{chrom}:g.{pos}_{pos + len(alt) - 1}del{ref[1:]}"
    # When record is insertion
    elif len(ref) == 1 and len(alt) > 1:
        return f"{chrom}:g.{pos}_{pos + 1}ins{alt[1:]}"
    # When record is duplication
    elif len(alt) > len(ref) and alt.startswith(ref):
        duplicated_seq = alt[len(ref):]
        end_pos = pos + len(duplicated_seq) - 1
        return f"{chrom}:g.{pos}_{end_pos}dup{duplicated_seq}"
    # When record is delin
    elif len(ref) > 1 and len(alt) > 1:
        end_pos = pos + len(ref) - 1
        return f"{chrom}:g.{pos}_{end_pos}delins{alt}"
    else:
        return None


# Function for sending a get request to HGVS endpoint
def ensemble_request(vcf_data):
    server = "https://grch37.rest.ensembl.org"
    # Iterate through every record in vcf_data list and sending get request for fetching info
    for record_info in vcf_data:
        hgvs_notation = construct_hgvs_notation(record_info['chrom'], record_info['pos'], record_info['ref'], record_info['alt'])
        vep_request = requests.get(f'{server}/vep/human/hgvs/{hgvs_notation}', headers={"Content-Type": "application/json"} )
        if not vep_request.ok:
            print(f"Failed to retrieve data for {hgvs_notation}")
            continue
        decoded = vep_request.json()
        record_info.update(response_json_parser(decoded))


    return vcf_data

# Function for parsing response json
def response_json_parser(response_json):
    # Iterate through response json and extract data
    for json_records in response_json:
        vep_data = {'gene_name': '',
                    'variant_effect': '',
                    'minor_allele': '',
                    'minor_allele_frequency': '',
                    'somatic': '',
                    'id': ''}
        # Finding transcipts_consequences in records and extracting gene_symbol
        if 'transcript_consequences' in json_records:
            for consequence in json_records['transcript_consequences']:
                if 'gene_symbol' in consequence:
                    vep_data['gene_name'] = consequence['gene_symbol']
                    break
        # Extracting variant effect if exists
        vep_data['variant_effect'] = json_records.get('most_severe_consequence', '')
        # Extract minor allele, frequency and ID from colocated variants list
        if 'colocated_variants' in json_records:
            for variant in json_records['colocated_variants']:
                if 'minor_allele' in variant:
                    vep_data['minor_allele'] = variant['minor_allele']
                    vep_data['minor_allele_frequency'] = variant.get('minor_allele_freq', '')

                is_somatic = variant.get('somatic', False)
                cosmic_id = variant.get('id', '') if variant.get('id', '').startswith("COS") else None
                rs_id = variant.get('id', '') if variant.get('id', '').startswith("rs") else None

                if is_somatic and cosmic_id:
                    vep_data['id'] = cosmic_id
                    vep_data['somatic'] = variant.get('somatic')
                    break
                elif not is_somatic and rs_id:
                    vep_data['id'] = rs_id
                    break
    return vep_data

# Code execution and export output in csv file
if __name__ == '__main__':
    vcf_file_path = sys.argv[1]
    csv_file_path = sys.argv[2]
    vcf_parsed_data = vcf_parser(vcf_file_path)
    full_data = ensemble_request(vcf_parsed_data)
    with open(csv_file_path, 'w', newline='', encoding='utf-8') as csvfile:
        fieldnames = ['chrom', 'pos', 'ref', 'alt', 'depth', 'alt_reads', 'percent_alt_reads', 'percent_ref_reads',
                      'gene_name', 'variant_effect', 'minor_allele', 'minor_allele_frequency', 'somatic', 'id']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

        writer.writeheader()
        for row in full_data:
            writer.writerow(row)

