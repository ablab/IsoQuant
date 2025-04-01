import logging
import pandas as pd
from typing import Dict, List, Tuple, Optional, Union
from pathlib import Path
from rpy2 import robjects
from rpy2.robjects import r, Formula
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter
from src.visualization_plotter import ExpressionVisualizer
from src.visualization_mapping import GeneMapper
import numpy as np
from scipy.stats import gmean
from sklearn.decomposition import PCA
from rpy2.rinterface_lib import callbacks

class DifferentialAnalysis:
    def __init__(
        self,
        output_dir: Path,
        viz_output: Path,
        ref_conditions: List[str],
        target_conditions: List[str],
        updated_gene_dict: Dict[str, Dict],
        ref_only: bool = False,
        dictionary_builder: "DictionaryBuilder" = None,
    ):
        """Initialize differential expression analysis."""
        def quiet_cb(x):
            pass
        
        # Silence R stdout/stderr
        callbacks.logger.setLevel(logging.WARNING)  # Affects R's logging only
        callbacks.consolewrite_print = quiet_cb
        callbacks.consolewrite_warnerror = quiet_cb
        
        self.output_dir = Path(output_dir)
        self.deseq_dir = Path(viz_output) / "differential_expression"
        self.deseq_dir.mkdir(parents=True, exist_ok=True)
        self.ref_conditions = ref_conditions
        self.target_conditions = target_conditions
        self.ref_only = ref_only
        self.updated_gene_dict = updated_gene_dict
        self.dictionary_builder = dictionary_builder
        
        # Create a single logger for this class
        self.logger = logging.getLogger('IsoQuant.visualization.differential_exp')
        
        # Get transcript mapping if available
        self.transcript_map = {}
        if hasattr(self.dictionary_builder, 'config') and hasattr(self.dictionary_builder.config, 'transcript_map'):
            self.transcript_map = self.dictionary_builder.config.transcript_map
            if self.transcript_map:
                self.logger.info(f"Using transcript mapping from dictionary_builder with {len(self.transcript_map)} entries for DESeq2 analysis")
            else:
                # Try to load transcript mapping directly from file
                self.logger.info("Transcript mapping from dictionary_builder is empty, trying to load it directly from file")
                self._load_transcript_mapping_from_file()
        else:
            # Try to load transcript mapping directly from file
            self.logger.info("No transcript mapping available from dictionary_builder, trying to load it directly from file")
            self._load_transcript_mapping_from_file()
        
        self.transcript_to_gene = self._create_transcript_to_gene_map()
        self.visualizer = ExpressionVisualizer(self.deseq_dir)
        self.gene_mapper = GeneMapper()

    def _load_transcript_mapping_from_file(self):
        """Load transcript mapping directly from the transcript_mapping.tsv file."""
        mapping_file = self.output_dir / "transcript_mapping.tsv"
        
        if not mapping_file.exists():
            self.logger.warning(f"Transcript mapping file not found at {mapping_file}")
            return
        
        try:
            # Load the transcript mapping file
            self.logger.info(f"Loading transcript mapping from {mapping_file}")
            self.transcript_map = {}
            
            # Skip header and read the mapping
            with open(mapping_file, 'r') as f:
                header = f.readline()  # Skip header
                for line in f:
                    parts = line.strip().split('\t')
                    if len(parts) == 2:
                        transcript_id, canonical_id = parts
                        self.transcript_map[transcript_id] = canonical_id
            
            self.logger.info(f"Successfully loaded {len(self.transcript_map)} transcript mappings from file")
            
            # Log some examples for debugging
            sample_items = list(self.transcript_map.items())[:5]
            for orig, canon in sample_items:
                self.logger.debug(f"Mapping sample: {orig} → {canon}")
        except Exception as e:
            self.logger.error(f"Failed to load transcript mapping: {str(e)}")

    def _create_transcript_to_gene_map(self) -> Dict[str, str]:
        """
        Create a mapping from transcript IDs to gene names.

        Returns:
            Dict[str, str]: Mapping from transcript ID to gene name.
        """
        transcript_map = {}
        for gene_category, genes in self.updated_gene_dict.items():
            for gene_id, gene_info in genes.items():
                gene_name = gene_info.get("name", gene_id)
                transcripts = gene_info.get("transcripts", {})
                for transcript_id, transcript_info in transcripts.items():
                    transcript_name = transcript_info.get("name", gene_name)
                    transcript_map[transcript_id] = transcript_name
        return transcript_map

    def run_complete_analysis(self) -> Tuple[Path, Path, pd.DataFrame, pd.DataFrame]:
        """Run differential expression analysis for both genes and transcripts.
        
        Returns:
            Tuple containing:
                - Path to gene results file
                - Path to transcript results file
                - DataFrame of transcript counts
                - DataFrame of DESeq2 gene-level results
        """
        self.logger.info("Starting differential expression analysis")

        valid_transcripts = set()
        for condition_genes in self.updated_gene_dict.values():
            for gene_info in condition_genes.values():
                valid_transcripts.update(gene_info.get("transcripts", {}).keys())

        # --- 1. Load Count Data ---
        gene_counts = self._get_condition_data("gene_grouped_counts.tsv")
        transcript_counts = self._get_condition_data("transcript_grouped_counts.tsv")
        self.logger.debug(f"Transcript counts shape after loading: {transcript_counts.shape}")
        self.logger.debug(f"Gene counts shape after loading: {gene_counts.shape}")

        # --- 2. Novel Transcript Filtering (Transcript Level) ---
        if self.dictionary_builder:
            novel_transcript_ids = self.dictionary_builder.get_novel_feature_ids()[1]
            self.logger.debug(f"Number of novel transcripts identified: {len(novel_transcript_ids)}")

            original_transcript_count_novel_filter = transcript_counts.shape[0]
            transcript_counts = transcript_counts[~transcript_counts.index.isin(novel_transcript_ids)] # Filter out novel transcripts
            novel_filtered_count = transcript_counts.shape[0]
            removed_novel_count = original_transcript_count_novel_filter - novel_filtered_count
            self.logger.info(f"Novel transcript filtering: Removed {removed_novel_count} transcripts from novel genes ({removed_novel_count / original_transcript_count_novel_filter * 100:.1f}%)")
            self.logger.debug(f"Transcript counts shape after novel gene filtering: {transcript_counts.shape}")
        else:
            self.logger.info("Novel transcript filtering: Skipped (no dictionary builder)")

        # --- 3. Valid Transcript Filtering (Transcript Level) ---
        original_transcript_count_valid_filter = transcript_counts.shape[0]
        transcript_counts = transcript_counts[transcript_counts.index.isin(valid_transcripts)] # Filter to valid transcripts
        valid_transcript_filtered_count = transcript_counts.shape[0]
        removed_valid_transcript_count = original_transcript_count_valid_filter - valid_transcript_filtered_count
        self.logger.info(f"Valid transcript filtering: Removed {removed_valid_transcript_count} transcripts not in updated_gene_dict ({removed_valid_transcript_count / original_transcript_count_valid_filter * 100:.1f}%)")
        self.logger.debug(f"Transcript counts shape after valid transcript filtering: {transcript_counts.shape}")

        if transcript_counts.empty:
            self.logger.error("No valid transcripts found after filtering.")
            raise ValueError("No valid transcripts found after filtering.")

        # --- 4. Count-based Filtering (Gene and Transcript Levels) ---
        gene_counts_filtered = self._filter_counts(gene_counts, level="gene")
        transcript_counts_filtered = self._filter_counts(transcript_counts, level="transcript") # Filter transcript counts AFTER novel and valid transcript filtering

        self.transcript_count_data = transcript_counts_filtered # Store filtered transcript counts

        # --- 5. Run DESeq2 Analysis ---
    
        deseq2_results_gene_file, gene_normalized_counts = self._run_level_analysis(
            level="gene",
            pattern="gene_grouped_counts.tsv",
            count_data=gene_counts_filtered,
            coldata=self._build_design_matrix(gene_counts_filtered)
        )
        deseq2_results_transcript_file, transcript_normalized_counts = self._run_level_analysis(
            level="transcript",
            pattern="transcript_model_grouped_counts.tsv" if not self.ref_only else "transcript_grouped_counts.tsv",
            count_data=transcript_counts_filtered
        )

        # Load the gene-level results for GSEA
        deseq2_results_df = pd.read_csv(deseq2_results_gene_file)

        # Update how we create the labels
        target_label = "+".join(self.target_conditions)
        reference_label = "+".join(self.ref_conditions)

        # --- Visualize Gene-Level Results ---
        self.visualizer.visualize_results(
            results=deseq2_results_df,
            target_label=target_label,
            reference_label=reference_label,
            min_count=10,
            feature_type="genes",
        )
        self.logger.info(f"Gene-level visualizations saved to {self.deseq_dir}")

        # Run PCA with correct labels
        normalized_gene_counts = self._median_ratio_normalization(gene_counts_filtered)
        self._run_pca(
            normalized_gene_counts, 
            "gene", 
            coldata=self._build_design_matrix(gene_counts_filtered),
            target_label=target_label,
            reference_label=reference_label
        )

        # --- Visualize Transcript-Level Results ---
        self.visualizer.visualize_results(
            results=pd.read_csv(deseq2_results_transcript_file),
            target_label=target_label,
            reference_label=reference_label,
            min_count=10,
            feature_type="transcripts",
        )
        self.logger.info(f"Transcript-level visualizations saved to {self.deseq_dir}")

        # Run PCA with correct labels for transcript level
        normalized_transcript_counts = self._median_ratio_normalization(transcript_counts_filtered)
        self._run_pca(
            normalized_transcript_counts, 
            "transcript", 
            coldata=self._build_design_matrix(transcript_counts_filtered),
            target_label=target_label,
            reference_label=reference_label
        )

        return deseq2_results_gene_file, deseq2_results_transcript_file, transcript_counts_filtered, deseq2_results_df

    def _run_level_analysis(
        self, level: str, count_data: pd.DataFrame, pattern: Optional[str] = None, coldata=None
    ) -> Tuple[Path, pd.DataFrame]:
        """
        Run DESeq2 analysis for a specific level and return results.

        Args:
            level: Analysis level ("gene" or "transcript")
            pattern: Optional pattern for output file naming (not used for data loading anymore)
            count_data: PRE-FILTERED count data DataFrame

        Returns:
            Tuple containing: (results_path, results_df)
        """
        # --- SIMPLIFIED: _run_level_analysis now assumes count_data is already loaded and filtered ---

        if count_data.empty:
            self.logger.error(f"Input count data is empty for level: {level}")
            raise ValueError(f"Input count data is empty for level: {level}")

        filtered_data = count_data.copy() # Work with a copy to avoid modifying original

        # Create design matrix and run DESeq2
        coldata = self._build_design_matrix(filtered_data)
        results, normalized_counts_r = self._run_deseq2(filtered_data, coldata, level)

        # Process results
        results.index.name = "feature_id"
        results.reset_index(inplace=True)
        mapping = self._map_gene_symbols(results["feature_id"].unique(), level)
        
        # Add transcript_symbol and gene_name columns
        results["transcript_symbol"] = results["feature_id"].map(lambda x: mapping[x]["transcript_symbol"])
        results["gene_name"] = results["feature_id"].map(lambda x: mapping[x]["gene_name"])

        # Drop transcript_symbol column for gene-level analysis as it's redundant
        if level == "gene":
            results = results.drop(columns=["transcript_symbol"])

        # Save results
        target_label = "+".join(self.target_conditions)
        reference_label = "+".join(self.ref_conditions)
        outfile = self.deseq_dir / f"DE_{level}_{target_label}_vs_{reference_label}.csv"
        results.to_csv(outfile, index=False)
        self.logger.info(f"Saved DESeq2 results to {outfile}")

        # Write top genes
        self._write_top_genes(results, level)

        # No normalized counts returned from _run_deseq2 anymore
        return outfile, pd.DataFrame() # Return empty DataFrame for normalized counts

    def _get_merged_transcript_counts(self, pattern: str) -> pd.DataFrame:
        """
        Get transcript count data and apply transcript mapping to create a merged grouped dataframe.
        This preserves the individual sample columns needed for DESeq2, but merges identical transcripts.
        """
        self.logger.info(f"Creating merged transcript count matrix with pattern: {pattern}")
        
        # Adjust pattern if needed
        adjusted_pattern = pattern
        if not self.ref_only and pattern == "transcript_grouped_counts.tsv":
            adjusted_pattern = "transcript_model_grouped_counts.tsv"
        
        self.logger.info(f"Using file pattern: {adjusted_pattern}")
        
        # Store sample dataframes
        all_sample_dfs = []
        
        # Process each condition directory
        for condition in self.ref_conditions + self.target_conditions:
            condition_dir = Path(self.output_dir) / condition
            count_files = list(condition_dir.glob(f"*{adjusted_pattern}"))
            
            if not count_files:
                self.logger.error(f"No count files found for condition: {condition}")
                raise FileNotFoundError(f"No count files matching {adjusted_pattern} found in {condition_dir}")
            
            # Load each count file
            for file_path in count_files:
                self.logger.info(f"Reading count data from: {file_path}")
                
                # Load the file
                df = pd.read_csv(file_path, sep="\t")
                if "#feature_id" not in df.columns and df.columns[0].startswith("#"):
                    # Rename first column if it's the feature ID column but named differently
                    df.rename(columns={df.columns[0]: "#feature_id"}, inplace=True)
                
                # Set feature_id as index
                df.set_index("#feature_id", inplace=True)
                
                # For multi-condition data (typical in sample files)
                # We need to prefix each column with the condition name
                for col in df.columns:
                    df.rename(columns={col: f"{condition}_{col}"}, inplace=True)
                
                all_sample_dfs.append(df)
        
        # Concatenate all dataframes to get the full matrix
        if not all_sample_dfs:
            self.logger.error("No sample data frames found")
            raise ValueError("No sample data found")
        
        # Combine all sample dataframes
        combined_df = pd.concat(all_sample_dfs, axis=1)
        self.logger.info(f"Combined count data shape before mapping: {combined_df.shape}")
        
        # Apply transcript mapping if available
        if not hasattr(self, 'transcript_map') or not self.transcript_map:
            self.logger.info("No transcript mapping available, using raw counts")
            return combined_df
        
        # Log transcript mapping info
        self.logger.info(f"Applying transcript mapping with {len(self.transcript_map)} mappings")
        
        # Get unique transcript IDs and create mapping dictionary
        unique_transcripts = combined_df.index.unique()
        transcript_groups = {}
        
        # Group transcripts by their canonical ID
        for transcript_id in unique_transcripts:
            canonical_id = self.transcript_map.get(transcript_id, transcript_id)
            if canonical_id not in transcript_groups:
                transcript_groups[canonical_id] = []
            transcript_groups[canonical_id].append(transcript_id)
        
        # Create the merged dataframe
        merged_df = pd.DataFrame(index=list(transcript_groups.keys()), columns=combined_df.columns)
        
        # Track merge statistics
        total_transcripts = len(unique_transcripts)
        merged_groups = 0
        merged_transcripts = 0
        
        # For each canonical transcript ID, sum the counts from all transcripts that map to it
        for canonical_id, transcript_ids in transcript_groups.items():
            if len(transcript_ids) == 1:
                # Just one transcript, copy the row directly
                merged_df.loc[canonical_id] = combined_df.loc[transcript_ids[0]]
            else:
                # Multiple transcripts map to this canonical ID, sum their counts
                merged_df.loc[canonical_id] = combined_df.loc[transcript_ids].sum()
                merged_groups += 1
                merged_transcripts += len(transcript_ids) - 1  # Count transcripts beyond the first one
                
                # Log details of significant merges (more than 2 transcripts or interesting transcripts)
                if len(transcript_ids) > 2 or any("ENST" in t for t in transcript_ids):
                    self.logger.info(f"Merged transcript group for {canonical_id}: {transcript_ids}")
        
        # Log merge statistics
        self.logger.info(f"Transcript merging complete: {merged_groups} canonical IDs had multiple transcripts")
        self.logger.info(f"Merged {merged_transcripts} transcripts into canonical IDs ({merged_transcripts/total_transcripts:.1%} of total)")
        self.logger.info(f"Final merged count matrix shape: {merged_df.shape}")
        
        return merged_df

    def _get_condition_data(self, pattern: str) -> pd.DataFrame:
        """Get count data for differential expression analysis."""
        if pattern == "transcript_grouped_counts.tsv":
            # For transcript data, use our merged function
            return self._get_merged_transcript_counts(pattern)
        elif pattern == "gene_grouped_counts.tsv":
            # For gene data, use a simpler approach (no merging needed)
            self.logger.info(f"Loading gene count data with pattern: {pattern}")
            
            # Store sample dataframes
            all_sample_dfs = []
            
            # Process each condition directory
            for condition in self.ref_conditions + self.target_conditions:
                condition_dir = Path(self.output_dir) / condition
                count_files = list(condition_dir.glob(f"*{pattern}"))
                
                if not count_files:
                    self.logger.error(f"No gene count files found for condition: {condition}")
                    raise FileNotFoundError(f"No count files matching {pattern} found in {condition_dir}")
                
                # Load each count file
                for file_path in count_files:
                    self.logger.info(f"Reading gene count data from: {file_path}")
                    
                    # Load the file
                    df = pd.read_csv(file_path, sep="\t")
                    if "#feature_id" not in df.columns and df.columns[0].startswith("#"):
                        # Rename first column if it's the feature ID column but named differently
                        df.rename(columns={df.columns[0]: "#feature_id"}, inplace=True)
                    
                    # Set feature_id as index
                    df.set_index("#feature_id", inplace=True)
                    
                    # For multi-condition data (typical in sample files)
                    # We need to prefix each column with the condition name
                    for col in df.columns:
                        df.rename(columns={col: f"{condition}_{col}"}, inplace=True)
                    
                    all_sample_dfs.append(df)
            
            # Concatenate all dataframes to get the full matrix
            if not all_sample_dfs:
                self.logger.error("No gene sample data frames found")
                raise ValueError("No gene sample data found")
            
            # Combine all sample dataframes
            combined_df = pd.concat(all_sample_dfs, axis=1)
            self.logger.info(f"Combined gene count data shape: {combined_df.shape}")
            return combined_df
        else:
            self.logger.error(f"Unsupported count pattern: {pattern}")
            raise ValueError(f"Unsupported count pattern: {pattern}")

    def _filter_counts(self, count_data: pd.DataFrame, min_count: int = 10, level: str = "gene") -> pd.DataFrame:
        """
        Filter features based on counts.
        
        For genes: Keep if mean count >= min_count in either condition group
        For transcripts: Keep if count >= min_count in at least half of all samples
        """
        if level == "transcript":
            total_samples = len(count_data.columns)
            min_samples_required = total_samples // 2
            samples_passing = (count_data >= min_count).sum(axis=1)
            keep_features = samples_passing >= min_samples_required
            
            self.logger.info(
                f"Transcript filtering: Keeping transcripts with counts >= {min_count} "
                f"in at least {min_samples_required}/{total_samples} samples"
            )
        else:  # gene level
            ref_cols = [
                col for col in count_data.columns
                if any(col.startswith(f"{cond}_") for cond in self.ref_conditions)
            ]
            tgt_cols = [
                col for col in count_data.columns
                if any(col.startswith(f"{cond}_") for cond in self.target_conditions)
            ]

            ref_means = count_data[ref_cols].mean(axis=1)
            tgt_means = count_data[tgt_cols].mean(axis=1)
            keep_features = (ref_means >= min_count) | (tgt_means >= min_count)
            
            self.logger.info(
                f"Gene filtering: Keeping genes with mean count >= {min_count} "
                f"in either condition group"
            )

        filtered_data = count_data[keep_features]
        self.logger.info(
            f"After filtering: Retained {filtered_data.shape[0]} features"
        )
        
        return filtered_data

    def _build_design_matrix(self, count_data: pd.DataFrame) -> pd.DataFrame:
        """Create experimental design matrix for DESeq2.
        
        Each column in the count data (sample) needs to be assigned to either 
        the reference or target group for differential expression analysis.
        """
        groups = []
        condition_assignments = []
        sample_ids = []
        
        self.logger.info("Building experimental design matrix")
        
        for sample in count_data.columns:
            # Extract the condition from the sample name
            # Matches pattern: conditionname_sampleid
            # The column name should start with the condition name followed by an underscore
            condition = None
            for cond in self.ref_conditions + self.target_conditions:
                if sample.startswith(f"{cond}_"):
                    condition = cond
                    # Extract the sample ID (everything after the condition name and underscore)
                    sample_id = sample[len(condition)+1:]
                    break
            
            if condition is None:
                self.logger.error(f"Could not determine condition for sample: {sample}")
                raise ValueError(f"Sample column '{sample}' does not match any specified condition")
            
            # Assign to reference or target group
            if condition in self.ref_conditions:
                groups.append("Reference")
            else:
                groups.append("Target")
                
            # Store the condition and sample ID for additional information
            condition_assignments.append(condition)
            sample_ids.append(sample)
        
        # Create the design matrix DataFrame
        design_matrix = pd.DataFrame({
            "group": groups,
            "condition": condition_assignments,
            "sample_id": sample_ids
        }, index=count_data.columns)
        
        # Log the design matrix for debugging
        self.logger.info(f"Design matrix:\n{design_matrix}")
        
        return design_matrix

    def _run_deseq2(
        self, count_data: pd.DataFrame, coldata: pd.DataFrame, level: str
    ) -> Tuple[pd.DataFrame, pd.DataFrame]:
        """Run DESeq2 analysis."""
        deseq2 = importr("DESeq2")
        count_data = count_data.fillna(0).round().astype(int)

        with localconverter(robjects.default_converter + pandas2ri.converter):
            # Convert count_data and coldata to R DataFrames explicitly before creating DESeqDataSet
            count_data_r = pandas2ri.py2rpy(count_data)
            coldata_r = pandas2ri.py2rpy(coldata)

            dds = deseq2.DESeqDataSetFromMatrix(
                countData=count_data_r, colData=coldata_r, design=Formula("~ group")
            )
            dds = deseq2.DESeq(dds)
            res = deseq2.results(
                dds, contrast=robjects.StrVector(["group", "Target", "Reference"])
            )
            # No normalized counts from DESeq2 anymore

            return pd.DataFrame(
                robjects.conversion.rpy2py(r("data.frame")(res)), index=count_data.index
            ), pd.DataFrame() # Return empty DataFrame for normalized counts

    def _map_gene_symbols(self, feature_ids: List[str], level: str) -> Dict[str, Dict[str, Optional[str]]]:
        """
        Map feature IDs to gene and transcript names using GeneMapper class.
        
        For transcripts that have been mapped to canonical IDs, ensure we properly handle the mapping.

        Args:
            feature_ids: List of feature IDs (gene symbols or transcript IDs)
            level: Analysis level ("gene" or "transcript")

        Returns:
            Dict[str, Dict[str, Optional[str]]]: Mapping from feature ID to a dictionary
                                                  containing 'transcript_symbol' and 'gene_name'.
                                                  'transcript_symbol' is None for gene-level analysis.
        """
        # Check if we need to handle canonical transcript IDs
        if level == "transcript" and self.transcript_map:
            # Create a mapping from canonical IDs to original IDs for reverse lookup
            canonical_to_original = {}
            for original, canonical in self.transcript_map.items():
                if canonical not in canonical_to_original:
                    canonical_to_original[canonical] = []
                canonical_to_original[canonical].append(original)
            
            # Process feature_ids that may include canonical IDs
            result = {}
            for feature_id in feature_ids:
                # First try to map directly
                direct_map = self.gene_mapper.map_gene_symbols([feature_id], level, self.updated_gene_dict)
                
                # If direct mapping worked, use it
                if feature_id in direct_map and direct_map[feature_id]["gene_name"]:
                    result[feature_id] = direct_map[feature_id]
                    continue
                
                # If this is a canonical ID, try to map using one of its original IDs
                if feature_id in canonical_to_original:
                    for original_id in canonical_to_original[feature_id]:
                        original_map = self.gene_mapper.map_gene_symbols([original_id], level, self.updated_gene_dict)
                        if original_id in original_map and original_map[original_id]["gene_name"]:
                            # Use the original ID's mapping but keep the canonical ID as the transcript symbol
                            result[feature_id] = {
                                "transcript_symbol": feature_id,
                                "gene_name": original_map[original_id]["gene_name"]
                            }
                            self.logger.debug(f"Mapped canonical ID {feature_id} using original ID {original_id}")
                            break
                
                # If still not mapped, use a default mapping
                if feature_id not in result:
                    result[feature_id] = {
                        "transcript_symbol": feature_id,
                        "gene_name": feature_id.split('.')[0] if '.' in feature_id else feature_id
                    }
            
            return result
        
        # For gene level or when no transcript mapping is available, use the original method
        return self.gene_mapper.map_gene_symbols(feature_ids, level, self.updated_gene_dict)

    def _write_top_genes(self, results: pd.DataFrame, level: str) -> None:
        """Write genes associated with top 100 transcripts by absolute fold change to file."""
        results["abs_stat"] = abs(results["stat"])

        if level == "transcript":
            # where baseMean is greater than 500
            top_transcripts = results[results["baseMean"] > 500].nlargest(len(results), "abs_stat") 

            unique_genes = set()
            top_unique_gene_transcripts = []
            transcript_count = 0
            unique_gene_count = 0

            for _, transcript_row in top_transcripts.iterrows():
                gene_name = transcript_row["gene_name"]
                if gene_name not in unique_genes:
                    unique_genes.add(gene_name)
                    top_unique_gene_transcripts.append(transcript_row)
                    unique_gene_count += 1
                    if unique_gene_count >= 100: # Stop when we reach 500 unique genes
                        break
                transcript_count += 1 # Keep track of total transcripts considered

            top_genes = [row["gene_name"] for row in top_unique_gene_transcripts] # Extract gene names from selected transcripts

            # Write to file
            top_genes_file = self.deseq_dir / "genes_of_top_100_DE_transcripts.txt"
            pd.Series(top_genes).to_csv(top_genes_file, index=False, header=False)
            self.logger.debug(f"Wrote genes of top 100 DE transcripts to {top_genes_file}")
        else:
            # For gene-level analysis, keep original behavior
            # top_genes = results.nlargest(100, "abs_stat")["symbol"] # OLD: was writing symbols (gene IDs)
            top_genes = results.nlargest(100, "abs_stat")["gene_name"]
            top_genes_file = self.deseq_dir / "top_100_DE_genes.txt"
            top_genes.to_csv(top_genes_file, index=False, header=False)
            self.logger.debug(f"Wrote top 100 DE genes to {top_genes_file}")

    def _run_pca(self, normalized_counts, level, coldata, target_label, reference_label):
        """Run PCA analysis and create visualization."""
        self.logger.info(f"Running PCA for {level} level...")

        # Run PCA - Calculate 10 components
        pca = PCA(n_components=10) # Keep n_components=10 to generate scree plot with 10 components
        log_normalized_counts = np.log2(normalized_counts + 1)
        pca_result = pca.fit_transform(log_normalized_counts.transpose())
        #map the feature names to gene names using the gene_mapper
        feature_names = normalized_counts.index.tolist()
        gene_names = self.gene_mapper.map_gene_symbols(feature_names, level, self.updated_gene_dict)

        # Get explained variance ratio and loadings
        explained_variance = pca.explained_variance_ratio_
        loadings = pca.components_  # Loadings are in pca.components_

        # Create DataFrame with columns for all 10 components
        pc_columns = [f'PC{i+1}' for i in range(10)] # Generate column names: PC1, PC2, ..., PC10
        pca_df = pd.DataFrame(data=pca_result, columns=pc_columns, index=log_normalized_counts.columns) # Use all 10 column names
        pca_df['group'] = coldata['group'].values

        title = f"{level.capitalize()} Level PCA: {target_label} vs {reference_label}\nPC1 ({100*explained_variance[0]:.2f}%) / PC2 ({100*explained_variance[1]:.2f}%)"

        # Use the plotter's PCA method, passing explained variance and loadings
        self.visualizer.plot_pca(
            pca_df=pca_df, # pca_df now contains 10 components
            title=title,
            output_prefix=f"pca_{level}",
            explained_variance=explained_variance, # Pass explained variance (for scree plot)
            loadings=loadings, # Pass loadings (for loadings output)
            feature_names=gene_names # Pass feature names (gene names)
        )

    def _median_ratio_normalization(self, count_data: pd.DataFrame) -> pd.DataFrame:
        """
        Perform median-by-ratio normalization on count data.
        This is similar to the normalization used in DESeq2.
        Handles zeros and potential data type issues safely.
        """
        try:
            # Convert to numeric and handle any non-numeric values
            count_data_numeric = count_data.apply(pd.to_numeric, errors='coerce').fillna(0)
            
            # Ensure all values are positive or zero
            count_data_nonneg = count_data_numeric.clip(lower=0)
            
            # Add pseudocount to avoid zeros (1 is a common choice)
            count_data_safe = count_data_nonneg + 1
            
            # Check data types and values
            self.logger.debug(f"Count data shape: {count_data_safe.shape}")
            self.logger.debug(f"Count data dtype: {count_data_safe.values.dtype}")
            self.logger.debug(f"Min value: {count_data_safe.values.min()}, Max value: {count_data_safe.values.max()}")
            
            # Convert to numpy array
            counts_numpy = count_data_safe.values.astype(float)
            
            # Alternative geometric mean calculation
            # Use log1p which is log(1+x) to handle zeros more safely
            log_counts = np.log(counts_numpy)
            row_means = np.mean(log_counts, axis=1)
            geometric_means = np.exp(row_means)
            
            # Check for any invalid values in geometric means
            if np.any(~np.isfinite(geometric_means)):
                self.logger.warning("Found non-finite values in geometric means, replacing with 1.0")
                geometric_means[~np.isfinite(geometric_means)] = 1.0
            
            # Calculate ratio of each count to the geometric mean
            # Reshape geometric_means for broadcasting
            geo_means_col = geometric_means.reshape(-1, 1)
            ratios = counts_numpy / geo_means_col
            
            # Calculate size factor for each sample (median of ratios)
            size_factors = np.median(ratios, axis=0)
            
            # Check for any invalid values in size factors
            if np.any(size_factors <= 0) or np.any(~np.isfinite(size_factors)):
                self.logger.warning("Found invalid size factors, replacing with 1.0")
                size_factors[~np.isfinite(size_factors) | (size_factors <= 0)] = 1.0
            
            # Log size factors
            self.logger.info(f"Size factors: {size_factors}")
            
            # Normalize counts by dividing by size factors
            # Use original count data (without pseudocount) for the final normalization
            normalized_counts = pd.DataFrame(
                count_data_nonneg.values / size_factors,
                index=count_data.index,
                columns=count_data.columns
            )
            
            # Fill any NaN values with 0
            normalized_counts = normalized_counts.fillna(0)
            
            return normalized_counts
            
        except Exception as e:
            self.logger.error(f"Error in median ratio normalization: {str(e)}")
            self.logger.error("Falling back to simple TPM-like normalization")
            
            # Fallback normalization (similar to TPM)
            # Sum each column and divide counts by the sum
            col_sums = count_data.sum(axis=0)
            col_sums = col_sums.replace(0, 1)  # Avoid division by zero
            
            # Normalize each column by its sum and multiply by 1e6 (similar to TPM scaling)
            normalized_counts = count_data.div(col_sums, axis=1) * 1e6
            
            return normalized_counts