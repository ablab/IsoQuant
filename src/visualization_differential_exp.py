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
        filter_min_count: int = 10,
        pca_n_components: int = 10,
        top_transcripts_base_mean: int = 500,
        top_n_genes: int = 100,
        log_level: int = logging.INFO, # Allow configuring log level
        tech_rep_dict: Dict[str, str] = None,
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
        
        # Configurable parameters
        self.filter_min_count = filter_min_count
        self.pca_n_components = pca_n_components
        self.top_transcripts_base_mean = top_transcripts_base_mean
        self.top_n_genes = top_n_genes # Used for both gene and transcript top list size
        
        # Create a single logger for this class
        self.logger = logging.getLogger('IsoQuant.visualization.differential_exp')
        self.logger.setLevel(log_level) # Set logger level
        
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
        self.tech_rep_dict = tech_rep_dict

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
        """
        Run differential expression analysis for both genes and transcripts.
        Orchestrates loading, filtering, DESeq2 execution, and visualization.

        Returns:
            Tuple containing:
                - Path to gene results file
                - Path to transcript results file
                - DataFrame of transcript counts (filtered but not normalized)
                - DataFrame of DESeq2 gene-level results (unfiltered by significance)
        """
        self.logger.info("Starting differential expression analysis workflow.")

        # --- 1. Load and Filter Data ---
        gene_counts_filtered, transcript_counts_filtered = self._load_and_filter_data()

        # Store filtered transcript counts (as required by original return signature)
        self.transcript_count_data = transcript_counts_filtered

        # --- 2. Run DESeq2 Analysis (Gene Level) ---
        (deseq2_results_gene_file,
         deseq2_results_df_gene,
         gene_normalized_counts) = self._perform_level_analysis("gene", gene_counts_filtered)

        # --- 3. Run DESeq2 Analysis (Transcript Level) ---
        (deseq2_results_transcript_file,
         deseq2_results_df_transcript,
         transcript_normalized_counts) = self._perform_level_analysis("transcript", transcript_counts_filtered)

        # --- 4. Generate Visualizations ---
        self._generate_visualizations(
            gene_counts_filtered=gene_counts_filtered, # Pass filtered counts for coldata generation
            transcript_counts_filtered=transcript_counts_filtered, # Pass filtered counts for coldata generation
            gene_normalized_counts=gene_normalized_counts,
            transcript_normalized_counts=transcript_normalized_counts,
            deseq2_results_df_gene=deseq2_results_df_gene,
            deseq2_results_df_transcript=deseq2_results_df_transcript
        )

        self.logger.info("Differential expression analysis workflow complete.")
        # Return signature matches original: results files, filtered transcript counts, gene results df
        return deseq2_results_gene_file, deseq2_results_transcript_file, transcript_counts_filtered, deseq2_results_df_gene

    def _load_and_filter_data(self) -> Tuple[pd.DataFrame, pd.DataFrame]:
        """Loads, filters (novelty, validity, counts), and returns gene and transcript count data."""
        self.logger.info("Loading and filtering count data...")

        # --- Load Count Data ---
        gene_counts = self._get_condition_data("gene_grouped_counts.tsv")
        transcript_counts = self._get_condition_data("transcript_grouped_counts.tsv")
        self.logger.debug(f"Raw transcript counts shape: {transcript_counts.shape}")
        self.logger.debug(f"Raw gene counts shape: {gene_counts.shape}")

        # --- Apply Transcript-Specific Filters ---
        transcript_counts_filtered = self._apply_transcript_filters(transcript_counts)

        # --- Apply Count-based Filtering (Gene Level) ---
        gene_counts_filtered = self._filter_counts(gene_counts, level="gene")

        if gene_counts_filtered.empty:
             self.logger.error("No genes remaining after count filtering.")
             raise ValueError("No genes remaining after count filtering.")
        if transcript_counts_filtered.empty:
             self.logger.error("No transcripts remaining after count filtering.")
             raise ValueError("No transcripts remaining after count filtering.")

        self.logger.info("Data loading and filtering complete.")
        self.logger.info(f"Final gene counts shape: {gene_counts_filtered.shape}")
        self.logger.info(f"Final transcript counts shape: {transcript_counts_filtered.shape}")

        return gene_counts_filtered, transcript_counts_filtered

    def _apply_transcript_filters(self, transcript_counts: pd.DataFrame) -> pd.DataFrame:
        """Applies novel, valid, and count-based filters specifically to transcript data."""
        self.logger.debug(f"Applying filters to transcript data (initial shape: {transcript_counts.shape})")

        # --- Valid Transcript Set ---
        # Determine the set of transcripts considered valid based on the updated_gene_dict
        valid_transcripts = set()
        for condition_genes in self.updated_gene_dict.values():
            for gene_info in condition_genes.values():
                valid_transcripts.update(gene_info.get("transcripts", {}).keys())
        if not valid_transcripts:
             self.logger.warning("No valid transcripts found in updated_gene_dict. Skipping validity filter.")
        self.logger.debug(f"Found {len(valid_transcripts)} valid transcript IDs in updated_gene_dict.")

        # --- Novel Transcript Filtering ---
        if self.dictionary_builder:
            novel_transcript_ids = self.dictionary_builder.get_novel_feature_ids()[1] # Assuming index 1 is transcripts
            self.logger.debug(f"Number of novel transcripts identified: {len(novel_transcript_ids)}")
            original_count = transcript_counts.shape[0]
            transcript_counts = transcript_counts[~transcript_counts.index.isin(novel_transcript_ids)]
            removed_count = original_count - transcript_counts.shape[0]
            perc_removed = (removed_count / original_count * 100) if original_count > 0 else 0
            self.logger.info(f"Novel Gene filtering: Removed {removed_count} transcripts ({perc_removed:.1f}%)")
            self.logger.debug(f"Shape after novel filtering: {transcript_counts.shape}")
        else:
            self.logger.info("Novel transcript filtering: Skipped (no dictionary builder).")



        if transcript_counts.empty:
            self.logger.warning("No transcripts remaining after novel gene filtering. Count filtering will be skipped.")
            return transcript_counts # Return empty dataframe

        # --- Count-based Filtering (Transcript Level) ---
        transcript_counts_filtered = self._filter_counts(transcript_counts, level="transcript")

        self.logger.debug(f"Final transcript counts shape after all filters: {transcript_counts_filtered.shape}")
        return transcript_counts_filtered

    def _perform_level_analysis(
        self, level: str, count_data: pd.DataFrame
    ) -> Tuple[Path, pd.DataFrame, pd.DataFrame]:
        """
        Runs DESeq2 analysis for a specific level (gene or transcript).

        Args:
            level: Analysis level ("gene" or "transcript").
            count_data: PRE-FILTERED count data DataFrame for the level.

        Returns:
            Tuple containing:
                - Path to the saved DESeq2 results CSV file.
                - DataFrame of the DESeq2 results.
                - DataFrame of the DESeq2 normalized counts.
        """
        self.logger.info(f"Performing DESeq2 analysis for level: {level}")

        if count_data.empty:
            self.logger.error(f"Input count data is empty for level: {level}")
            raise ValueError(f"Input count data is empty for level: {level}")

        # Create design matrix
        coldata = self._build_design_matrix(count_data)

        # Run DESeq2 - Now returns results and normalized counts
        results_df, normalized_counts_df = self._run_deseq2(count_data, coldata, level)

        # --- Process DESeq2 Results ---
        results_df.index.name = "feature_id"
        results_df.reset_index(inplace=True) # Keep feature_id as a column

        # Map gene symbols/names
        mapping = self._map_gene_symbols(results_df["feature_id"].unique(), level)

        # Add transcript_symbol and gene_name columns safely using .get
        results_df["transcript_symbol"] = results_df["feature_id"].map(
            lambda x: mapping.get(x, {}).get("transcript_symbol", x) # Default to feature_id if not found
        )
        results_df["gene_name"] = results_df["feature_id"].map(
            lambda x: mapping.get(x, {}).get("gene_name", x.split('.')[0] if '.' in x else x) # Default to feature_id logic if not found
        )


        # Drop transcript_symbol column for gene-level analysis as it's redundant
        if level == "gene":
            results_df = results_df.drop(columns=["transcript_symbol"], errors='ignore') # Use errors='ignore'

        # --- Save Results ---
        target_label = "+".join(self.target_conditions)
        reference_label = "+".join(self.ref_conditions)
        # Use the pattern argument passed to _get_condition_data if needed, or derive filename like this
        outfile = self.deseq_dir / f"DE_{level}_{target_label}_vs_{reference_label}.csv"
        results_df.to_csv(outfile, index=False)
        self.logger.info(f"Saved DESeq2 results ({results_df.shape[0]} features) to {outfile}")

        # --- Write Top Genes/Transcripts ---
        self._write_top_genes(results_df, level)

        self.logger.info(f"DESeq2 analysis complete for level: {level}")
        return outfile, results_df, normalized_counts_df

    def _generate_visualizations(
        self,
        gene_counts_filtered: pd.DataFrame,
        transcript_counts_filtered: pd.DataFrame,
        gene_normalized_counts: pd.DataFrame,
        transcript_normalized_counts: pd.DataFrame,
        deseq2_results_df_gene: pd.DataFrame,
        deseq2_results_df_transcript: pd.DataFrame,
    ):
        """Generates PCA plots and other visualizations based on DESeq2 results and normalized counts."""
        self.logger.info("Generating visualizations...")
        target_label = "+".join(self.target_conditions)
        reference_label = "+".join(self.ref_conditions)

        # --- Visualize Gene-Level DE Results ---
        self.visualizer.visualize_results(
            results=deseq2_results_df_gene, # Use DataFrame directly
            target_label=target_label,
            reference_label=reference_label,
            min_count=self.filter_min_count, # Use configured value
            feature_type="genes",
        )
        self.logger.info(f"Gene-level DE summary visualizations saved to {self.deseq_dir}")

        # --- Run PCA (Gene Level) ---
        if not gene_normalized_counts.empty:
            gene_coldata = self._build_design_matrix(gene_counts_filtered) # Need coldata matching the counts used
            self._run_pca(
                normalized_counts=gene_normalized_counts,
                level="gene",
                coldata=gene_coldata,
                target_label=target_label,
                reference_label=reference_label
            )
        else:
            self.logger.warning("Skipping gene-level PCA: Normalized counts are empty.")

        # --- Visualize Transcript-Level DE Results ---
        self.visualizer.visualize_results(
            results=deseq2_results_df_transcript, # Use DataFrame directly
            target_label=target_label,
            reference_label=reference_label,
            min_count=self.filter_min_count, # Use configured value
            feature_type="transcripts",
        )
        self.logger.info(f"Transcript-level DE summary visualizations saved to {self.deseq_dir}")

        # --- Run PCA (Transcript Level) ---
        if not transcript_normalized_counts.empty:
            transcript_coldata = self._build_design_matrix(transcript_counts_filtered) # Need coldata matching the counts used
            self._run_pca(
                normalized_counts=transcript_normalized_counts,
                level="transcript",
                coldata=transcript_coldata,
                target_label=target_label,
                reference_label=reference_label
            )
        else:
             self.logger.warning("Skipping transcript-level PCA: Normalized counts are empty.")

        self.logger.info("Visualizations generated.")

    def _get_merged_transcript_counts(self, pattern: str) -> pd.DataFrame:
        """
        Get transcript count data and apply transcript mapping to create a merged grouped dataframe.
        This preserves the individual sample columns needed for DESeq2, but merges identical transcripts.
        """
        self.logger.debug(f"Creating merged transcript count matrix with pattern: {pattern}")
        
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
        
        # Apply technical replicate merging before transcript mapping
        combined_df = self._merge_technical_replicates(combined_df)
        
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
                    self.logger.debug(f"Merged transcript group for {canonical_id}: {transcript_ids}")
        
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
            
            # Apply technical replicate merging
            combined_df = self._merge_technical_replicates(combined_df)
            
            return combined_df
        else:
            self.logger.error(f"Unsupported count pattern: {pattern}")
            raise ValueError(f"Unsupported count pattern: {pattern}")

    def _filter_counts(self, count_data: pd.DataFrame, level: str = "gene") -> pd.DataFrame:
        """
        Filter features based on counts using the configured threshold.

        For genes: Keep if mean count >= configured min_count in either condition group.
        For transcripts: Keep if count >= configured min_count in at least half of all samples.
        """
        if count_data.empty:
            self.logger.warning(f"Input count data for filtering ({level}) is empty. Returning empty DataFrame.")
            return count_data

        # Use the configured minimum count threshold
        min_count_threshold = self.filter_min_count

        if level == "transcript":
            total_samples = len(count_data.columns)
            min_samples_required = max(1, total_samples // 2) # Ensure at least 1 sample is required
            samples_passing = (count_data >= min_count_threshold).sum(axis=1)
            keep_features = samples_passing >= min_samples_required

            self.logger.info(
                f"Transcript filtering: Keeping transcripts with counts >= {min_count_threshold} "
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

            # Handle cases where one condition might have no samples after potential upstream filtering
            if not ref_cols:
                self.logger.warning("No reference columns found for gene count filtering.")
                ref_means = pd.Series(0, index=count_data.index) # Assign 0 mean if no ref samples
            else:
                ref_means = count_data[ref_cols].mean(axis=1)

            if not tgt_cols:
                 self.logger.warning("No target columns found for gene count filtering.")
                 tgt_means = pd.Series(0, index=count_data.index) # Assign 0 mean if no target samples
            else:
                tgt_means = count_data[tgt_cols].mean(axis=1)

            keep_features = (ref_means >= min_count_threshold) | (tgt_means >= min_count_threshold)

            self.logger.info(
                f"Gene filtering: Keeping genes with mean count >= {min_count_threshold} "
                f"in either reference or target condition group"
            )

        filtered_data = count_data[keep_features]
        removed_count = count_data.shape[0] - filtered_data.shape[0]
        self.logger.info(
            f"After count filtering ({level}): Retained {filtered_data.shape[0]} / {count_data.shape[0]} features "
            f"(Removed {removed_count})"
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
        self.logger.debug(f"Design matrix:\n{design_matrix}")
        
        return design_matrix

    def _run_deseq2(
        self, count_data: pd.DataFrame, coldata: pd.DataFrame, level: str
    ) -> Tuple[pd.DataFrame, pd.DataFrame]:
        """
        Run DESeq2 analysis and return results and normalized counts.

        Args:
            count_data: Raw count data (filtered).
            coldata: Design matrix.
            level: Analysis level (gene/transcript).

        Returns:
            Tuple[pd.DataFrame, pd.DataFrame]: DESeq2 results, DESeq2 normalized counts.
        """
        self.logger.info(f"Running DESeq2 for {level} level...")
        deseq2 = importr("DESeq2")
        # Ensure counts are integers for DESeq2
        count_data = count_data.fillna(0).round().astype(int)

        # Ensure count data has no negative values before passing to R
        if (count_data < 0).any().any():
             self.logger.warning(f"Negative values found in count data for {level}. Clamping to 0.")
             count_data = count_data.clip(lower=0)

        if count_data.empty:
            self.logger.error(f"Count data is empty before running DESeq2 for {level}.")
            # Return empty dataframes if counts are empty
            return pd.DataFrame(), pd.DataFrame(index=count_data.index, columns=count_data.columns)

        try:
            with localconverter(robjects.default_converter + pandas2ri.converter):
                # Convert count_data and coldata to R DataFrames
                count_data_r = robjects.conversion.py2rpy(count_data)
                coldata_r = robjects.conversion.py2rpy(coldata)

                # Create DESeqDataSet
                self.logger.debug("Creating DESeqDataSet...")
                dds = deseq2.DESeqDataSetFromMatrix(
                    countData=count_data_r, colData=coldata_r, design=Formula("~ group")
                )

                # Run DESeq analysis
                self.logger.debug("Running DESeq()...")
                dds = deseq2.DESeq(dds)

                # Get results
                self.logger.debug("Extracting results()...")
                res = deseq2.results(
                    dds, contrast=robjects.StrVector(["group", "Target", "Reference"])
                )
                res_df = robjects.conversion.rpy2py(r("as.data.frame")(res)) # Convert to R dataframe first for stability
                res_df.index = count_data.index # Assign original feature IDs as index

                # Extract dispersion estimates
                self.logger.debug("Extracting dispersion estimates...")
                dispersions_r = r['dispersions'](dds)
                dispersions_py = robjects.conversion.rpy2py(dispersions_r)
                
                # Add dispersion estimates to results DataFrame
                res_df['dispersion'] = dispersions_py

                # Extract size factors
                self.logger.debug("Extracting size factors...")
                size_factors_r = r['sizeFactors'](dds)
                size_factors_py = robjects.conversion.rpy2py(size_factors_r)
                
                # Create size factors DataFrame with sample names
                size_factors_df = pd.DataFrame({
                    'sample': count_data.columns,
                    'size_factor': size_factors_py
                })
                
                # Save size factors to file
                target_label = "+".join(self.target_conditions)
                reference_label = "+".join(self.ref_conditions)
                size_factors_file = self.deseq_dir / f"size_factors_{level}_{target_label}_vs_{reference_label}.csv"
                size_factors_df.to_csv(size_factors_file, index=False)
                self.logger.info(f"Size factors saved to {size_factors_file}")

                # Correct way to call the R 'counts' function on the dds object
                # Ensure 'r' is imported: from rpy2.robjects import r
                normalized_counts_r = r['counts'](dds, normalized=True)

                # Convert R matrix to pandas DataFrame
                normalized_counts_py = robjects.conversion.rpy2py(normalized_counts_r)
                # Ensure DataFrame structure matches original count_data (features x samples)
                normalized_counts_df = pd.DataFrame(normalized_counts_py, index=count_data.index, columns=count_data.columns)

                # Generate dispersion and count summaries
                self._generate_dispersion_summary(res_df, level)

                self.logger.info(f"DESeq2 run completed for {level}. Results shape: {res_df.shape}, Normalized counts shape: {normalized_counts_df.shape}")
                return res_df, normalized_counts_df

        except Exception as e:
            self.logger.error(f"Error running DESeq2 for {level}: {str(e)}")
            # Return empty DataFrames on error to avoid downstream issues
            return pd.DataFrame(), pd.DataFrame(index=count_data.index, columns=count_data.columns)

    def _generate_dispersion_summary(self, results_df: pd.DataFrame, level: str) -> None:
        """
        Generate summary statistics for average read counts and dispersion estimates.
        Saves summary to a file and logs key statistics.

        Args:
            results_df: DESeq2 results DataFrame with baseMean and dispersion columns
            level: Analysis level (gene/transcript)
        """
        if results_df.empty:
            self.logger.warning(f"Cannot generate dispersion summary for {level}: Results DataFrame is empty.")
            return

        self.logger.info(f"Generating dispersion and count summary for {level} level...")

        # Check if required columns exist
        required_cols = ['baseMean', 'dispersion']
        missing_cols = [col for col in required_cols if col not in results_df.columns]
        if missing_cols:
            self.logger.warning(f"Cannot generate complete summary for {level}: Missing columns {missing_cols}")
            return

        # Remove NaN values for summary statistics
        clean_data = results_df[['baseMean', 'dispersion']].dropna()
        
        if clean_data.empty:
            self.logger.warning(f"No valid data for dispersion summary for {level} after removing NaN values.")
            return

        # Calculate summary statistics
        summary_stats = {
            'level': level,
            'total_features': len(results_df),
            'features_with_valid_data': len(clean_data),
            
            # Average read count (baseMean) statistics
            'baseMean_mean': clean_data['baseMean'].mean(),
            'baseMean_median': clean_data['baseMean'].median(),
            'baseMean_std': clean_data['baseMean'].std(),
            'baseMean_min': clean_data['baseMean'].min(),
            'baseMean_max': clean_data['baseMean'].max(),
            'baseMean_q25': clean_data['baseMean'].quantile(0.25),
            'baseMean_q75': clean_data['baseMean'].quantile(0.75),
            
            # Dispersion statistics
            'dispersion_mean': clean_data['dispersion'].mean(),
            'dispersion_median': clean_data['dispersion'].median(),
            'dispersion_std': clean_data['dispersion'].std(),
            'dispersion_min': clean_data['dispersion'].min(),
            'dispersion_max': clean_data['dispersion'].max(),
            'dispersion_q25': clean_data['dispersion'].quantile(0.25),
            'dispersion_q75': clean_data['dispersion'].quantile(0.75),
        }

        # Add size factor statistics if available
        target_label = "+".join(self.target_conditions)
        reference_label = "+".join(self.ref_conditions)
        size_factors_file = self.deseq_dir / f"size_factors_{level}_{target_label}_vs_{reference_label}.csv"
        
        if size_factors_file.exists():
            try:
                size_factors_df = pd.read_csv(size_factors_file)
                if 'size_factor' in size_factors_df.columns:
                    sf_data = size_factors_df['size_factor'].dropna()
                    summary_stats.update({
                        'size_factor_mean': sf_data.mean(),
                        'size_factor_median': sf_data.median(),
                        'size_factor_std': sf_data.std(),
                        'size_factor_min': sf_data.min(),
                        'size_factor_max': sf_data.max(),
                        'size_factor_q25': sf_data.quantile(0.25),
                        'size_factor_q75': sf_data.quantile(0.75),
                    })
                    self.logger.info(f"  Size factors: mean={summary_stats['size_factor_mean']:.4f}, median={summary_stats['size_factor_median']:.4f}, range={summary_stats['size_factor_min']:.4f}-{summary_stats['size_factor_max']:.4f}")
            except Exception as e:
                self.logger.warning(f"Could not read size factors file: {e}")

        # Log key statistics
        self.logger.info(f"{level.capitalize()} level summary:")
        self.logger.info(f"  Total features: {summary_stats['total_features']}")
        self.logger.info(f"  Features with valid data: {summary_stats['features_with_valid_data']}")
        self.logger.info(f"  Average read count (baseMean): mean={summary_stats['baseMean_mean']:.2f}, median={summary_stats['baseMean_median']:.2f}")
        self.logger.info(f"  Dispersion: mean={summary_stats['dispersion_mean']:.4f}, median={summary_stats['dispersion_median']:.4f}")

        # Create a more detailed summary for significant DE genes/transcripts
        if 'padj' in results_df.columns:
            significant_features = results_df[results_df['padj'] < 0.05].dropna(subset=['baseMean', 'dispersion'])
            if not significant_features.empty:
                summary_stats.update({
                    'significant_features_count': len(significant_features),
                    'significant_baseMean_mean': significant_features['baseMean'].mean(),
                    'significant_baseMean_median': significant_features['baseMean'].median(),
                    'significant_dispersion_mean': significant_features['dispersion'].mean(),
                    'significant_dispersion_median': significant_features['dispersion'].median(),
                })
                
                self.logger.info(f"  Significant DE features (padj < 0.05): {summary_stats['significant_features_count']}")
                self.logger.info(f"  Significant features - Average read count: mean={summary_stats['significant_baseMean_mean']:.2f}, median={summary_stats['significant_baseMean_median']:.2f}")
                self.logger.info(f"  Significant features - Dispersion: mean={summary_stats['significant_dispersion_mean']:.4f}, median={summary_stats['significant_dispersion_median']:.4f}")

        # Save summary to file
        summary_file = self.deseq_dir / f"dispersion_count_summary_{level}_{target_label}_vs_{reference_label}.txt"
        
        with open(summary_file, 'w') as f:
            f.write(f"Dispersion and Count Summary for {level.capitalize()} Level Analysis\n")
            f.write(f"Comparison: {target_label} vs {reference_label}\n")
            f.write("=" * 60 + "\n\n")
            
            for key, value in summary_stats.items():
                if isinstance(value, float):
                    f.write(f"{key}: {value:.6f}\n")
                else:
                    f.write(f"{key}: {value}\n")

        self.logger.info(f"Dispersion and count summary saved to {summary_file}")

        # Also save detailed data for further analysis
        detailed_file = self.deseq_dir / f"detailed_dispersion_data_{level}_{target_label}_vs_{reference_label}.csv"
        
        # Include feature mapping information if available
        detailed_data = results_df[['baseMean', 'dispersion', 'log2FoldChange', 'pvalue', 'padj']].copy()
        if 'gene_name' in results_df.columns:
            detailed_data['gene_name'] = results_df['gene_name']
        if 'transcript_symbol' in results_df.columns:
            detailed_data['transcript_symbol'] = results_df['transcript_symbol']
        
        detailed_data.to_csv(detailed_file)
        self.logger.info(f"Detailed dispersion data saved to {detailed_file}")

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
        """Write top genes/transcripts based on absolute statistic value to file."""
        if results.empty or 'stat' not in results.columns:
             self.logger.warning(f"Cannot write top genes for {level}: Results DataFrame is empty or missing 'stat' column.")
             return

        # Ensure 'stat' column is numeric, fill NaNs that might cause issues
        results['stat'] = pd.to_numeric(results['stat'], errors='coerce').fillna(0)
        results["abs_stat"] = abs(results["stat"])

        # Use configured number of top genes/transcripts
        top_n = self.top_n_genes

        if level == "transcript":
            # Use configured base mean threshold
            base_mean_threshold = self.top_transcripts_base_mean
             # Ensure 'baseMean' column is numeric, fill NaNs
            if 'baseMean' not in results.columns:
                self.logger.warning(f"Cannot apply baseMean filter for {level}: 'baseMean' column missing. Considering all transcripts.")
                filtered_results = results
            else:
                results['baseMean'] = pd.to_numeric(results['baseMean'], errors='coerce').fillna(0)
                filtered_results = results[results["baseMean"] > base_mean_threshold]

            if filtered_results.empty:
                 self.logger.warning(f"No transcripts found with baseMean > {base_mean_threshold}. Top genes file will be empty.")
                 top_unique_gene_transcripts_df = pd.DataFrame() # Empty dataframe
            else:
                # Sort by absolute statistic value
                top_transcripts = filtered_results.sort_values("abs_stat", ascending=False)

                # Ensure 'gene_name' column exists
                if 'gene_name' not in top_transcripts.columns:
                     self.logger.error(f"Cannot extract top unique genes for {level}: 'gene_name' column missing.")
                     return

                # Get top N unique genes based on the highest ranked transcript for each gene
                top_unique_gene_transcripts_df = top_transcripts.drop_duplicates(subset=['gene_name'], keep='first').head(top_n)
                self.logger.info(f"Highest adjusted p-value in top {top_n} unique genes: {top_unique_gene_transcripts_df['padj'].max()}")

            top_genes_list = top_unique_gene_transcripts_df["gene_name"].tolist() if not top_unique_gene_transcripts_df.empty else []

            # Write to file
            target_label = "+".join(self.target_conditions)
            reference_label = "+".join(self.ref_conditions)
            top_genes_file = self.deseq_dir / f"genes_of_top_{top_n}_DE_transcripts_{target_label}_vs_{reference_label}.txt"

            pd.Series(top_genes_list).to_csv(top_genes_file, index=False, header=False)
            self.logger.info(f"Wrote {len(top_genes_list)} unique genes (from top {top_n} DE transcripts with baseMean > {base_mean_threshold}) to {top_genes_file}")

        else: # Gene level
             # Ensure 'gene_name' column exists for gene level as well
            if 'gene_name' not in results.columns:
                 self.logger.error(f"Cannot extract top genes for {level}: 'gene_name' column missing.")
                 return

            # Get top N genes directly by absolute statistic
            top_genes_df = results.nlargest(top_n, "abs_stat")
            top_genes_list = top_genes_df["gene_name"].tolist()

            # Write to file
            target_label = "+".join(self.target_conditions)
            reference_label = "+".join(self.ref_conditions)
            top_genes_file = self.deseq_dir / f"top_{top_n}_DE_genes_{target_label}_vs_{reference_label}.txt"

            pd.Series(top_genes_list).to_csv(top_genes_file, index=False, header=False)
            self.logger.info(f"Wrote top {len(top_genes_list)} DE genes to {top_genes_file}")

    def _run_pca(self, normalized_counts, level, coldata, target_label, reference_label):
        """Run PCA analysis and create visualization using DESeq2 normalized counts."""
        self.logger.info(f"Running PCA for {level} level using DESeq2 normalized counts...")

        if normalized_counts.empty:
            self.logger.warning(f"Skipping PCA for {level}: Normalized counts data is empty.")
            return

        # Basic check for variance - PCA fails if variance is zero
        if normalized_counts.var().sum() == 0:
            self.logger.warning(f"Skipping PCA for {level}: Data has zero variance.")
            return

        # Use configured number of components
        n_components = min(self.pca_n_components, normalized_counts.shape[0], normalized_counts.shape[1]) # Cannot exceed number of features or samples
        if n_components < 2:
             self.logger.warning(f"Skipping PCA for {level}: Not enough features/samples ({normalized_counts.shape}) for {n_components} components.")
             return
        if n_components != self.pca_n_components:
             self.logger.warning(f"Reducing number of PCA components to {n_components} due to data dimensions.")


        # Log transform the DESeq2 normalized counts (add 1 to handle zeros)
        # Ensure data is numeric before transformation
        log_normalized_counts = np.log2(normalized_counts.apply(pd.to_numeric, errors='coerce').fillna(0) + 1)


        # Check for NaNs/Infs after log transform which can happen if counts were negative (though clamped earlier) or exactly -1
        if np.isinf(log_normalized_counts).any().any() or np.isnan(log_normalized_counts).any().any():
            self.logger.warning(f"NaNs or Infs found in log-transformed counts for {level}. Replacing with 0. This might indicate issues with count data.")
            log_normalized_counts = log_normalized_counts.replace([np.inf, -np.inf], 0).fillna(0)


        try:
            pca = PCA(n_components=n_components)
            # Transpose because PCA expects samples as rows, features as columns
            pca_result = pca.fit_transform(log_normalized_counts.transpose())

            # Map feature IDs (index of normalized_counts) to gene names
            feature_ids = normalized_counts.index.tolist()
            # Use the mapping function - ensure it handles potential errors/missing keys
            gene_mapping_dict = self._map_gene_symbols(feature_ids, level)
            # Create a list of gene names in the same order as features
            feature_names_mapped = [gene_mapping_dict.get(fid, {}).get('gene_name', fid) for fid in feature_ids]


            # Get explained variance ratio and loadings
            explained_variance = pca.explained_variance_ratio_
            loadings = pca.components_  # Loadings are in pca.components_

            # Create DataFrame with columns for all calculated components
            pc_columns = [f'PC{i+1}' for i in range(n_components)]
            pca_df = pd.DataFrame(data=pca_result[:, :n_components], columns=pc_columns, index=log_normalized_counts.columns) # Use sample names as index

            # Add group information from coldata, ensuring index alignment
            # It's safer to reset index on coldata if it uses sample names as index too
            if coldata.index.equals(pca_df.index):
                 pca_df['group'] = coldata['group'].values
            else:
                 self.logger.warning(f"Index mismatch between PCA results and coldata for {level}. Group information might be incorrect.")
                 # Attempt to merge or handle, here just assigning potentially misaligned
                 pca_df['group'] = coldata['group'].values[:len(pca_df)]


            # Title focuses on PC1/PC2 for the scatter plot, even if more components were calculated
            pc1_var = explained_variance[0] * 100 if len(explained_variance) > 0 else 0
            pc2_var = explained_variance[1] * 100 if len(explained_variance) > 1 else 0
            title = f"{level.capitalize()} Level PCA: {target_label} vs {reference_label}\nPC1 ({pc1_var:.2f}%) / PC2 ({pc2_var:.2f}%)"


            # Use the plotter's PCA method, passing explained variance and loadings
            self.visualizer.plot_pca(
                pca_df=pca_df, # pca_df contains n_components columns
                title=title,
                output_prefix=f"pca_{level}",
                explained_variance=explained_variance, # Pass full explained variance for scree plot
                loadings=loadings, # Pass loadings
                # Pass the mapped gene names corresponding to the features (rows of normalized_counts)
                feature_names=feature_names_mapped
            )
            self.logger.info(f"PCA plots saved for {level} level.")

        except Exception as e:
            self.logger.error(f"Error during PCA calculation or plotting for {level}: {str(e)}")

    def _merge_technical_replicates(self, count_data: pd.DataFrame) -> pd.DataFrame:
        """
        Merge technical replicates by summing counts for samples in the same replicate group.
        
        Args:
            count_data: DataFrame with samples as columns and features as rows
            
        Returns:
            DataFrame with technical replicates merged
        """
        if not self.tech_rep_dict:
            self.logger.info("No technical replicates specified, returning original data")
            return count_data
        
        self.logger.info(f"Merging technical replicates using {len(self.tech_rep_dict)} mappings")
        
        # Create a mapping from sample columns to replicate groups
        sample_to_group = {}
        for col in count_data.columns:
            # Extract the base sample name (remove condition prefix if present)
            base_sample = col
            for condition in self.ref_conditions + self.target_conditions:
                if col.startswith(f"{condition}_"):
                    base_sample = col[len(condition)+1:]
                    break
            
            # Check if this sample is in the technical replicates mapping
            if base_sample in self.tech_rep_dict:
                group_name = self.tech_rep_dict[base_sample]
                # Reconstruct the group name with condition prefix
                condition_prefix = col.replace(base_sample, "").rstrip("_")
                if condition_prefix:
                    full_group_name = f"{condition_prefix}_{group_name}"
                else:
                    full_group_name = group_name
                sample_to_group[col] = full_group_name
            else:
                # Keep original sample name if not in technical replicates
                sample_to_group[col] = col
        
        # Group samples by their replicate groups
        group_to_samples = {}
        for sample, group in sample_to_group.items():
            if group not in group_to_samples:
                group_to_samples[group] = []
            group_to_samples[group].append(sample)
        
        # Create merged DataFrame
        merged_data = pd.DataFrame(index=count_data.index)
        
        merge_stats = {"merged_groups": 0, "original_samples": len(count_data.columns)}
        
        for group_name, samples in group_to_samples.items():
            if len(samples) == 1:
                # No merging needed, just rename
                merged_data[group_name] = count_data[samples[0]]
            else:
                # Sum technical replicates
                merged_data[group_name] = count_data[samples].sum(axis=1)
                merge_stats["merged_groups"] += 1
                self.logger.debug(f"Merged technical replicates for group {group_name}: {samples}")
        
        merge_stats["final_samples"] = len(merged_data.columns)
        self.logger.info(
            f"Technical replicate merging complete: "
            f"{merge_stats['original_samples']} samples -> {merge_stats['final_samples']} samples "
            f"({merge_stats['merged_groups']} groups had multiple replicates)"
        )
        
        return merged_data