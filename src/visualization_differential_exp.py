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
        
        self.transcript_to_gene = self._create_transcript_to_gene_map()
        self.visualizer = ExpressionVisualizer(self.deseq_dir)
        self.gene_mapper = GeneMapper()

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

    def _get_condition_data(self, pattern: str) -> pd.DataFrame:
        """Combine count data from all conditions."""
        all_counts = []

        # Modify pattern if ref_only is False and pattern is transcript_grouped_counts - CORRECTED LOGIC
        adjusted_pattern = pattern # Initialize adjusted_pattern to the original pattern
        if not self.ref_only and pattern == "transcript_grouped_counts.tsv": # Only adjust if ref_only is False AND base pattern is transcript_grouped_counts
            adjusted_pattern = "transcript_model_grouped_counts.tsv"

        for condition in self.ref_conditions + self.target_conditions:
            condition_dir = self.output_dir / condition
            count_files = list(condition_dir.glob(f"*{adjusted_pattern}")) # Use adjusted pattern

            for file_path in count_files:
                self.logger.debug(f"Reading count data from: {file_path}")
                df = pd.read_csv(file_path, sep="\t", dtype={"#feature_id": str})
                df.set_index("#feature_id", inplace=True)

                # Rename columns to include condition
                for col in df.columns:
                    if col.lower() != "count":
                        df = df.rename(columns={col: f"{condition}_{col}"})

                all_counts.append(df)

        return pd.concat(all_counts, axis=1)

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
        """Create experimental design matrix."""
        groups = []
        for sample in count_data.columns:
            if any(sample.startswith(f"{cond}_") for cond in self.ref_conditions):
                groups.append("Reference")
            else:
                groups.append("Target")

        return pd.DataFrame({"group": groups}, index=count_data.columns)

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

        Args:
            feature_ids: List of feature IDs (gene symbols or transcript IDs)
            level: Analysis level ("gene" or "transcript")

        Returns:
            Dict[str, Dict[str, Optional[str]]]: Mapping from feature ID to a dictionary
                                                  containing 'transcript_symbol' and 'gene_name'.
                                                  'transcript_symbol' is None for gene-level analysis.
        """
        return self.gene_mapper.map_gene_symbols(feature_ids, level, self.updated_gene_dict)

    def _write_top_genes(self, results: pd.DataFrame, level: str) -> None:
        """Write genes associated with top 100 transcripts by absolute fold change to file."""
        results["abs_stat"] = abs(results["stat"])

        if level == "transcript":
            top_transcripts = results.nlargest(len(results), "abs_stat") # Get ALL transcripts ranked by abs_stat

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
                    if unique_gene_count >= 100: # Stop when we reach 100 unique genes
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
        """
        # 1. Calculate geometric mean for each feature (row)
        geometric_means = count_data.apply(gmean, axis=1)

        # 2. Handle rows with zero geometric mean (replace with NaN to avoid division by zero)
        geometric_means[geometric_means == 0] = np.nan

        # 3. Calculate ratio of each count to the geometric mean
        count_ratios = count_data.divide(geometric_means, axis=0)

        # 4. Calculate size factor for each sample (column) as the median of ratios
        size_factors = count_ratios.median(axis=0)

        # 5. Normalize counts by dividing by size factors
        normalized_counts = count_data.divide(size_factors, axis=1)

        # 6. Fill NaN values with 0 after normalization
        normalized_counts = normalized_counts.fillna(0)

        return normalized_counts