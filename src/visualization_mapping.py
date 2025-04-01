import mygene
import logging
from typing import Dict, List, Tuple, Optional

class GeneMapper:
    def __init__(self):
        self.mg = mygene.MyGeneInfo()
        self.logger = logging.getLogger('IsoQuant.visualization.mapping')

    def get_gene_info_from_mygene(self, ensembl_ids: List[str]) -> Dict[str, Dict]:
        """
        Query MyGene.info API for gene information using batch query.
        
        Args:
            ensembl_ids: List of Ensembl gene IDs
            
        Returns:
            Dict mapping query IDs to gene information
        """
        try:
            # Batch query for gene information
            results = self.mg.querymany(
                ensembl_ids,
                scopes='ensembl.gene',  # Only search for gene IDs
                fields=['symbol', 'name'],  # Only get essential fields
                species='human',
                as_dataframe=False,
                returnall=True
            )

            # Process results
            mapping = {}
            for hit in results['out']:
                query_id = hit.get('query', '')
                if 'notfound' in hit:
                    self.logger.debug(f"Gene ID not found: {query_id}")
                    continue
                
                mapping[query_id] = {
                    'symbol': hit.get('symbol', query_id),
                    'name': hit.get('name', hit.get('symbol', query_id))
                }

            # Log query statistics
            self.logger.debug(
                f"MyGene.info query stats: "
                f"Total={len(ensembl_ids)}, "
                f"Found={len(mapping)}, "
                f"Missing={len(ensembl_ids) - len(mapping)}"
            )

            return mapping

        except Exception as e:
            self.logger.error(f"Failed to fetch info from MyGene.info: {str(e)}")
            return {}

    def map_genes(self, gene_ids: List[str], updated_gene_dict: Dict) -> Dict[str, Tuple[str, str]]:
        """
        Map Ensembl gene IDs to symbols, using multiple fallback methods:
        1. IsoQuant's updated_gene_dict
        2. MyGene.info
        3. Parse symbol from Ensembl ID if possible
        
        Returns:
            Dict mapping Ensembl IDs to (symbol, gene_name) tuples
        """
        mapping = {}
        unmapped_ids = []

        # First try to map using updated_gene_dict
        for gene_id in gene_ids:
            symbol_found = False
            for gene_category, genes in updated_gene_dict.items():
                if gene_id in genes:
                    gene_info = genes[gene_id]
                    # Only use name if it's not empty and not the same as gene_id
                    if gene_info.get("name") and gene_info["name"] != gene_id:
                        mapping[gene_id] = (gene_info["name"], gene_info["name"])
                        symbol_found = True
                        break
            
            if not symbol_found:
                unmapped_ids.append(gene_id)

        # For unmapped genes, try MyGene.info batch query
        if unmapped_ids:
            self.logger.debug(f"Querying MyGene.info for {len(unmapped_ids)} unmapped genes")
            mygene_results = self.get_gene_info_from_mygene(unmapped_ids)
            
            remaining_unmapped = []
            for gene_id in unmapped_ids:
                if gene_id in mygene_results:
                    info = mygene_results[gene_id]
                    mapping[gene_id] = (info['symbol'], info['name'])
                    self.logger.debug(f"Mapped {gene_id} to {info['symbol']} using MyGene.info")
                else:
                    remaining_unmapped.append(gene_id)

            # For still unmapped genes, try to extract info from Ensembl ID
            for gene_id in remaining_unmapped:
                # Try to extract meaningful info from Ensembl ID
                if gene_id.startswith('ENSG'):
                    # For novel genes, use the last part of the ID as a temporary symbol
                    temp_symbol = f"GENE_{gene_id.split('0')[-1]}"
                    mapping[gene_id] = (temp_symbol, gene_id)
                    self.logger.warning(f"Using derived symbol {temp_symbol} for {gene_id}")
                else:
                    mapping[gene_id] = (gene_id, gene_id)
                    self.logger.warning(f"Could not map {gene_id} using any method")

        return mapping

    def map_gene_symbols(self, feature_ids: List[str], level: str, updated_gene_dict: Dict = None) -> Dict[str, Dict[str, Optional[str]]]:
        """
        Map feature IDs to gene and transcript names using updated gene dictionary.

        Args:
            feature_ids: List of feature IDs (gene symbols or transcript IDs)
            level: Analysis level ("gene" or "transcript")
            updated_gene_dict: Optional updated gene dictionary

        Returns:
            Dict[str, Dict[str, Optional[str]]]: Mapping from feature ID to a dictionary
                                                  containing 'transcript_symbol' and 'gene_name'.
                                                  'transcript_symbol' is None for gene-level analysis.
        """
        mapping: Dict[str, Dict[str, Optional[str]]] = {}
        unmapped_gene_ids_batch: List[str] = [] # Initialize list to collect unmapped gene IDs for batch query

        for feature_id in feature_ids:
            if level == "gene":
                # Gene-level mapping: Search in updated_gene_dict, fallback to batched MyGene API
                gene_name = None
                found_in_dict = False
                if updated_gene_dict:
                    for condition, condition_gene_dict in updated_gene_dict.items():
                        if feature_id in condition_gene_dict:
                            found_in_dict = True
                            gene_name = condition_gene_dict[feature_id].get("name")
                            break
                    if not found_in_dict:
                        unmapped_gene_ids_batch.append(feature_id) # Add to batch list for MyGene query
                else:
                    unmapped_gene_ids_batch.append(feature_id) # Add to batch list for MyGene query


                mapping[feature_id] = {
                    "transcript_symbol": gene_name, # For gene-level, use gene name as transcript_symbol
                    "gene_name": gene_name if gene_name else feature_id
                }

            elif level == "transcript":
                # Transcript-level mapping: Search for transcript name in updated_gene_dict across all conditions
                gene_name = None
                transcript_symbol = None
                gene_found_for_transcript = False # Flag to track if gene is found for transcript

                if updated_gene_dict:
                    for condition, condition_gene_dict in updated_gene_dict.items(): # Iterate through conditions
                        for gene_id, gene_data in condition_gene_dict.items(): # Iterate through genes in each condition
                            if "transcripts" in gene_data and feature_id in gene_data["transcripts"]:
                                gene_found_for_transcript = True
                                transcript_info = gene_data["transcripts"].get(feature_id, {})
                                transcript_symbol = transcript_info.get("name")
                                mapping[feature_id] = {
                                    "transcript_symbol": f"{transcript_symbol} ({gene_data.get('name')})" if feature_id.startswith("transcript") else transcript_symbol,
                                    "gene_name": gene_data.get("name") # Get gene_name from gene_data
                                }
                                self.logger.debug(f"Transcript-level mapping: Found transcript {feature_id}, gene_data: {gene_data}") # Debug log to inspect gene_data
                                break # Found transcript, exit inner loop (genes in condition)
                        if gene_found_for_transcript: # If transcript found in any gene in this condition, exit condition loop
                            break
                    if not gene_found_for_transcript:
                        self.logger.debug(f"Transcript-level mapping: No gene found for Transcript ID {feature_id} in updated_gene_dict across any condition")
                        mapping[feature_id] = { # Assign mapping here for not found case
                            "transcript_symbol": f"{feature_id} (No gene name)", # Indicate no gene name available
                            "gene_name": None
                        }
                else: # If updated_gene_dict is None
                    self.logger.debug("Transcript-level mapping: updated_gene_dict is None")
                    mapping[feature_id] = {
                        "transcript_symbol": f"{feature_id} (No gene name)", # Indicate no gene name available when dict is None
                        "gene_name": None # gene_name is None
                    }
                    self.logger.debug(f"Transcript-level mapping: Using feature_id as transcript_symbol, no gene name available (updated_gene_dict is None)") # Debug log

            else:
                raise ValueError(f"Invalid level: {level}. Must be 'gene' or 'transcript'.")

        # Perform batched MyGene API query for all unmapped gene IDs at once (gene-level only)
        if level == "gene" and unmapped_gene_ids_batch:
            self.logger.debug(f"Gene-level mapping: Performing batched MyGene API query for {len(unmapped_gene_ids_batch)} gene IDs")
            mygene_batch_info = self.get_gene_info_from_mygene(unmapped_gene_ids_batch) # Batched query

            if mygene_batch_info:
                for feature_id in unmapped_gene_ids_batch: # Iterate through the unmapped IDs
                    if feature_id in mygene_batch_info: # Check if MyGene returned info for this ID
                        gene_name_from_mygene = mygene_batch_info[feature_id].get('symbol')
                        if gene_name_from_mygene:
                            mapping[feature_id]["gene_name"] = gene_name_from_mygene # Update gene_name in mapping
                            mapping[feature_id]["transcript_symbol"] = gene_name_from_mygene # Update transcript_symbol
                    else:
                        self.logger.debug(f"Gene-level mapping: Batched MyGene API did not return info for Feature ID {feature_id}")
            else:
                self.logger.warning("Gene-level mapping: Batched MyGene API query failed or returned no results.")


        return mapping 