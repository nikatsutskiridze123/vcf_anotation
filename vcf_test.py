import vcf
import requests
import csv
import sys
import json

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
# Function for sending a post request to HGVS endpoint in batches


# Function for sending a post request to HGVS endpoint in batches
def ensemble_request(vcf_data, batch_size=300):
    server = "http://grch37.rest.ensembl.org"
    endpoint = "/vep/human/hgvs"
    url = f'{server}{endpoint}'
    headers = {"Content-Type": "application/json", "Accept": "application/json"}

    # Split data into batches size of 300
    batches = [vcf_data[i:i + batch_size] for i in range(0, len(vcf_data), batch_size)]

    # Function to send a single batch
    def send_batch(batch):
        payload = json.dumps({"hgvs_notations": batch})
        vep_request = requests.post(url, headers=headers, data=payload)
        if not vep_request.ok:
            print("Failed to retrieve data for a batch")
            return []
        return vep_request.json()


    # Process each batch
    for batch in batches:
        hgvs_notations = [construct_hgvs_notation(record['chrom'], record['pos'], record['ref'], record['alt']) for
                          record in batch]
        batch_response = send_batch(hgvs_notations)

        for record_info, json_record in zip(batch, batch_response):
            # Print json_record for debugging


            update = response_json_parser(json_record)
            record_info.update(update)

    return vcf_data




# Function for parsing response json
def response_json_parser(response_json):
    # Initialize the dictionary with default values
    vep_data = {
        'gene_name': '',
        'variant_effect': '',
        'minor_allele': '',
        'minor_allele_frequency': '',
        'somatic': '',
        'id': ''
    }

    # Check if transcript_consequences is in response_json and process it
    if 'transcript_consequences' in response_json:
        for consequence in response_json['transcript_consequences']:
            if 'gene_symbol' in consequence:
                vep_data['gene_name'] = consequence['gene_symbol']
                break

    # Extract other fields from response_json
    vep_data['variant_effect'] = response_json.get('most_severe_consequence', '')

    # Process colocated_variants if it exists
    if 'colocated_variants' in response_json:
        for variant in response_json['colocated_variants']:
            is_somatic = variant.get('somatic', False)

            if 'minor_allele' in variant:
                vep_data['minor_allele'] = variant['minor_allele']
                vep_data['minor_allele_frequency'] = variant.get('minor_allele_freq', '')

            # Check if variant is somatic
            if is_somatic:
                vep_data['somatic'] = 1
                if variant.get('id', '').startswith("COS"):
                    vep_data['id'] = variant['id']
                    break  
            else:
                if variant.get('id', '').startswith("rs"):
                    vep_data['id'] = variant['id']


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