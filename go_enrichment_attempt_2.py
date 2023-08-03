import os
import re
import csv

def process_gff_files(directory):
    category_data = {}  # Dictionary to store data for each category

    # Loop through the files in the directory
    for filename in os.listdir(directory):
        if filename.startswith("P") and filename.endswith(".gff"):
            category = int(re.search(r'P(\d+)', filename).group(1))  # Extract the category number (P1, P2, ..., P10)
            if category not in category_data:
                category_data[category] = {}

            # Read the content of the file
            with open(os.path.join(directory, filename), "r") as file:
                content = file.read()

                # Use regular expression to find ';product=Y'
                matches = re.findall(r';product=(.*?)(?:;|\n|$)', content)

                # Count the occurrences and accumulate them in the dictionary
                for match in matches:
                    match = match.strip().lower()  # Convert product to lowercase
                    category_data[category][match] = category_data[category].get(match, 0) + 1

    return category_data

def load_manual_uniprot_mapping(file_path):
    manual_uniprot_mapping = {}
    with open(file_path, "r") as file:
        for line in file:
            product, uniprot_id = line.strip().split("\t")
            # Convert product to lowercase for case-insensitive matching
            manual_uniprot_mapping[product.lower()] = uniprot_id

    return manual_uniprot_mapping

def load_go_terms(file_path):
    go_terms = {}
    with open(file_path, "r") as file:
        reader = csv.reader(file, delimiter='\t')
        next(reader)  # Skip the header row
        for fields in reader:
            if len(fields) >= 11:
                uniprot_id = fields[1]
                bp_terms = fields[8].split(";") if fields[8] else []
                cc_terms = fields[9].split(";") if fields[9] else []
                mf_terms = fields[10].split(";") if fields[10] else []

                # Convert GO terms to lowercase for case-insensitive matching, remove whitespace, 
                bp_terms = [term.lower() for term in bp_terms]
                bp_terms = [s.strip() for s in bp_terms]
                bp_terms = [i.replace('"', '') for i in bp_terms]
                print(bp_terms)
                cc_terms = [term.lower() for term in cc_terms]
                cc_terms = [s.strip() for s in cc_terms]
                cc_terms = [i.replace('"', '') for i in cc_terms]
                mf_terms = [term.lower() for term in mf_terms]
                mf_terms = [s.strip() for s in mf_terms]
                mf_terms = [i.replace('"', '') for i in mf_terms]
                go_terms[uniprot_id] = {
                    "biological_process": bp_terms,
                    "cellular_component": cc_terms,
                    "molecular_function": mf_terms
                }

    return go_terms


def count_all_go_terms(category_data, manual_uniprot_mapping, go_terms_data):
    # Dictionary to count occurrences of all GO terms for each propagation
    go_terms_counts = {}

    for category, data in category_data.items():
        go_terms_counts[category] = {"Unknown": 0}

        for product, count in data.items():
            uniprot_id = manual_uniprot_mapping.get(product.lower())
            if uniprot_id:
                go_terms_info = go_terms_data.get(uniprot_id)
                print(product,go_terms_info)
                if go_terms_info:
                    for go_term in go_terms_info["biological_process"]:
                        go_terms_counts[category][go_term] = go_terms_counts[category].get(go_term, 0) + count
                    for go_term in go_terms_info["cellular_component"]:
                        go_terms_counts[category][go_term] = go_terms_counts[category].get(go_term, 0) + count
                    for go_term in go_terms_info["molecular_function"]:
                        go_terms_counts[category][go_term] = go_terms_counts[category].get(go_term, 0) + count
            else:
                go_terms_counts[category]["Unknown"] += count

    return go_terms_counts

# Example usage with the provided directory path
directory_path = "/Users/pimswart/P_all_separate_representatives_manualremovedbacterial_pharokka_gffs/"
result = process_gff_files(directory_path)

# Load the manual Uniprot mapping from /Users/pimswart/manual_uniprot.txt
manual_uniprot_file = "/Users/pimswart/manual_uniprot.txt"
manual_uniprot_mapping = load_manual_uniprot_mapping(manual_uniprot_file)

# Load GO terms from /Users/pimswart/Downloads/idmapping_2023_07_25_2.tsv
go_terms_file = "/Users/pimswart/Downloads/idmapping_2023_07_25_2.tsv"
go_terms_data = load_go_terms(go_terms_file)

# Count occurrences of all GO terms for each propagation
go_terms_counts = count_all_go_terms(result, manual_uniprot_mapping, go_terms_data)

# Write the data to a CSV file
output_file = "go_terms_counts.csv"
with open(output_file, "w", newline="") as csvfile:
    csv_writer = csv.writer(csvfile)
    
    # Write the header row
    header_row = ["Propagation", "GO Term", "Occurrences"]
    csv_writer.writerow(header_row)
    
    # Write the data rows
    for category, go_terms in go_terms_counts.items():
        for go_term, count in go_terms.items():
            row = [f"P{category}", go_term, count]
            csv_writer.writerow(row)

print(f"Data written to {output_file}")






def add_go_term_type_to_output(output_file, go_terms_file):
    # Load GO terms from idmapping file
    go_terms_data = load_go_terms(go_terms_file)

    # Create a dictionary to map each GO term to its corresponding type
    go_term_type_mapping = {}
    for uniprot_id, go_info in go_terms_data.items():
        for go_term in go_info["biological_process"]:
            go_term_type_mapping[go_term.lower()] = "Biological Process"
        for go_term in go_info["cellular_component"]:
            go_term_type_mapping[go_term.lower()] = "Cellular Component"
        for go_term in go_info["molecular_function"]:
            go_term_type_mapping[go_term.lower()] = "Molecular Function"

    # Update the output file with the new "GO Term Type" column
    temp_file = output_file + ".temp"
    with open(output_file, "r") as csvfile, open(temp_file, "w", newline="") as temp_csvfile:
        csv_reader = csv.reader(csvfile)
        csv_writer = csv.writer(temp_csvfile)

        # Write the header row with the additional "GO Term Type" column
        header_row = next(csv_reader)
        header_row.append("GO Term Type")
        csv_writer.writerow(header_row)

        # Write the data rows with the new "GO Term Type" information
        for row in csv_reader:
            go_term = row[1].lower()
            go_term_type = go_term_type_mapping.get(go_term, "Unknown")
            row.append(go_term_type)
            csv_writer.writerow(row)

    # Replace the original output file with the updated file
    os.replace(temp_file, output_file)

# Example usage after writing the data to the CSV file
output_file = "go_terms_counts.csv"
go_terms_file = "/Users/pimswart/Downloads/idmapping_2023_07_25_2.tsv"
add_go_term_type_to_output(output_file, go_terms_file)

print(f"GO Term Type added to {output_file}")
