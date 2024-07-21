import json
import sys
import os


def simplify_and_sum_transcripts(data):
    gene_totals_across_conditions = {}
    simplified_data = {}

    # Sum transcript values and collect them across all conditions
    for sample_id, genes in data.items():
        simplified_data[sample_id] = {}
        for gene_id, gene_data in genes.items():
            transcripts = gene_data.get("transcripts", {})
            total_value = 0.0
            simplified_transcripts = {}
            for transcript_id, transcript_details in transcripts.items():
                transcript_value = (
                    transcript_details.get("value", 0.0)
                    if isinstance(transcript_details, dict)
                    else 0.0
                )
                simplified_transcripts[transcript_id] = transcript_value
                total_value += transcript_value

            gene_data_copy = (
                gene_data.copy()
            )  # Make a copy to avoid modifying the original
            gene_data_copy["transcripts"] = simplified_transcripts
            gene_data_copy["value"] = (
                total_value  # Replace the gene-level value with the sum of transcript values
            )
            simplified_data[sample_id][gene_id] = gene_data_copy

            if gene_id not in gene_totals_across_conditions:
                gene_totals_across_conditions[gene_id] = []
            gene_totals_across_conditions[gene_id].append(total_value)

    # Determine which genes to remove
    genes_to_remove = [
        gene_id
        for gene_id, totals in gene_totals_across_conditions.items()
        if all(total < 5 for total in totals)
    ]

    # Remove genes from the simplified data structure
    for sample_id, genes in simplified_data.items():
        for gene_id in genes_to_remove:
            if gene_id in genes:
                del genes[gene_id]

    return simplified_data


def read_json(file_path):
    with open(file_path, "r") as file:
        return json.load(file)


def write_json(data, file_path):
    with open(file_path, "w") as file:
        json.dump(data, file, indent=4)


def main():
    if len(sys.argv) != 2:
        print("Usage: python script.py <input_file_path>")
        sys.exit(1)

    input_file_path = sys.argv[1]
    base, ext = os.path.splitext(input_file_path)
    output_file_path = f"{base}_simplified{ext}"

    try:
        # Load the gene data from the specified input JSON file
        gene_dict = read_json(input_file_path)

        # Simplify the transcripts, sum their values, and remove genes under a threshold across all conditions
        modified_gene_dict = simplify_and_sum_transcripts(gene_dict)

        # Save the modified gene data to the newly named output JSON file
        write_json(modified_gene_dict, output_file_path)

        print(f"Modified gene data has been saved to {output_file_path}")

    except Exception as e:
        print(f"Error: {str(e)}")
        sys.exit(1)


if __name__ == "__main__":
    main()
