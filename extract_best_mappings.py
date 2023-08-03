# #METHOD 1 - 'MAPPING QUALITY'
# import pysam

# def extract_best_mappings(sam_file):
#     best_mappings = {}  # Dictionary to store the best mappings for each query sequence

#     with pysam.AlignmentFile(sam_file, "r") as sam:
#         for alignment in sam:
#             query_sequence = alignment.query_name
#             mapping_quality = alignment.mapping_quality

#             # Check if the current alignment is better than the previous best mapping
#             if query_sequence not in best_mappings or mapping_quality > best_mappings[query_sequence]["mapping_quality"]:
#                 best_mappings[query_sequence] = {
#                     "mapping_quality": mapping_quality,
#                     "reference_name": alignment.reference_name,
#                     # Add more fields as needed
#                 }

#     return best_mappings

# # Usage example
# sam_file = "/Users/pimswart/minimap_P_all_representatives_manualremovedbacterial_fasta_catted_cat_references_namecharactersfixedforcircoletto_allprophageregionsmaskednotjustregion2_Ecolibl21DE3_LE392_and_maskedprophageregions_and_parentalphages.sam"
# output_file = "/Users/pimswart/1_extracted_best_mappings_minimap_P_all_representatives_manualremovedbacterial_fasta_catted_cat_references_namecharactersfixedforcircoletto_allprophageregionsmaskednotjustregion2_Ecolibl21DE3_LE392_and_maskedprophageregions_and_parentalphages.txt"
# best_mappings = extract_best_mappings(sam_file)

# # Write the best mappings to the output file
# with open(output_file, "w") as f:
#     for query_sequence, mapping in best_mappings.items():
#         f.write(f"{query_sequence}\t{mapping['reference_name']}\n")








#METHOD 2 - 'PRIMARY ALIGNMENT'
import pysam

def extract_primary_mappings(sam_file):
    primary_mappings = {}  # Dictionary to store the primary mappings for each query sequence

    with pysam.AlignmentFile(sam_file, "r") as sam:
        for alignment in sam:
            query_sequence = alignment.query_name
            mapping_quality = alignment.mapping_quality

            if not alignment.is_secondary and not alignment.is_supplementary:
                # Check if the current alignment is the primary alignment for the query sequence
                if query_sequence not in primary_mappings:
                    primary_mappings[query_sequence] = {
                        "mapping_quality": mapping_quality,
                        "reference_name": alignment.reference_name,
                        # Add more fields as needed
                    }
                else:
                    # Update the primary alignment if the current alignment has higher mapping quality
                    if mapping_quality > primary_mappings[query_sequence]["mapping_quality"]:
                        primary_mappings[query_sequence]["mapping_quality"] = mapping_quality
                        primary_mappings[query_sequence]["reference_name"] = alignment.reference_name

    return primary_mappings

# Usage example (THESE ARE OLD FILENAMES)
sam_file = "/Users/pimswart/minimap_P_all_representatives_manualremovedbacterial_fasta_catted_cat_references_namecharactersfixedforcircoletto_allprophageregionsmaskednotjustregion2_Ecolibl21DE3_LE392_and_maskedprophageregions_and_parentalphages.sam"
output_file = "/Users/pimswart/2_extracted_best_mappings_minimap_P_all_representatives_manualremovedbacterial_fasta_catted_cat_references_namecharactersfixedforcircoletto_allprophageregionsmaskednotjustregion2_Ecolibl21DE3_LE392_and_maskedprophageregions_and_parentalphages.txt"
primary_mappings = extract_primary_mappings(sam_file)

# Write the primary mappings to the output file
with open(output_file, "w") as f:
    for query_sequence, mapping in primary_mappings.items():
        f.write(f"{query_sequence}\t{mapping['reference_name']}\n")
