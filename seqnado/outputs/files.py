from typing import Protocol, List, Any
from pathlib import Path

from pydantic import BaseModel, computed_field, Field, field_validator
from typing import Literal
from snakemake.io import expand
from loguru import logger

from seqnado import (
    Assay,
    PileupMethod,
    DataScalingTechnique,
    PeakCallingMethod,
    MethylationMethod,
    QuantificationMethod,
    SpikeInMethod,
)
from seqnado.core import AssaysWithHeatmaps, AssaysWithSpikein, AssaysWithPeakCalling
from seqnado.inputs import (
    FastqCollection,
    FastqCollectionForIP,
    SampleGroups,
    BamCollection,
    BigWigCollection,
)
from seqnado.config.configs import QCConfig


class FileCollection(Protocol):
    @property
    def files(self) -> List[str]:
        """Return a list of file paths."""
        pass


class BasicFileCollection(BaseModel):
    files: List[str] = Field(default_factory=list)


class QCFiles(BaseModel):
    assay: Assay
    samples: FastqCollection | FastqCollectionForIP | BamCollection | BigWigCollection
    config: QCConfig = QCConfig()
    output_dir: str = "seqnado_output"

    @property
    def fastqc_files(self) -> list[str]:
        if isinstance(self.samples, (FastqCollection, FastqCollectionForIP)):
            files = []
            for sample in self.samples.sample_names:
                if self.samples.is_paired_end(sample):
                    # Paired-end: generate files for both reads
                    files.extend(
                        [
                            f"{self.output_dir}/qc/fastqc_raw/{sample}_1_fastqc.html",
                            f"{self.output_dir}/qc/fastqc_raw/{sample}_2_fastqc.html",
                        ]
                    )
                else:
                    # Single-end: generate file without read number
                    files.append(
                        f"{self.output_dir}/qc/fastqc_raw/{sample}_fastqc.html"
                    )
            return files
        return []

    @property
    def fastqscreen_files(self) -> list[str]:
        if (
            isinstance(self.samples, (FastqCollection, FastqCollectionForIP))
            and self.config.run_fastq_screen
        ):
            files = []
            for sample in self.samples.sample_names:
                if self.samples.is_paired_end(sample):
                    # Paired-end: generate files for both reads
                    files.extend(
                        [
                            f"{self.output_dir}/qc/fastq_screen/{sample}_1_screen.html",
                            f"{self.output_dir}/qc/fastq_screen/{sample}_2_screen.html",
                        ]
                    )
                else:
                    # Single-end: generate file without read number
                    files.append(
                        f"{self.output_dir}/qc/fastq_screen/{sample}_screen.html"
                    )
            return files
        return []

    @property
    def qualimap_files(self) -> list[str]:
        if not isinstance(self.samples, (BigWigCollection)):
            match self.assay:
                case Assay.RNA:
                    return expand(
                        f"{self.output_dir}/qc/qualimap_rnaseq/{{sample}}/qualimapReport.html",
                        sample=self.samples.sample_names,
                    )
                case _:
                    return expand(
                        f"{self.output_dir}/qc/qualimap_bamqc/{{sample}}/qualimapReport.html",
                        sample=self.samples.sample_names,
                    )
        return []

    @computed_field
    @property
    def files(self) -> List[str]:
        return [*self.fastqc_files, *self.fastqscreen_files, *self.qualimap_files]


class SeqNadoReportFile(BaseModel):
    output_dir: str = "seqnado_output"

    @property
    def report_file(self) -> list[str]:
        return [
            f"{self.output_dir}/seqnado_report.html",
        ]

    @property
    def protocol_file(self) -> list[str]:
        return [
            f"{self.output_dir}/protocol.txt",
        ]

    @computed_field
    @property
    def files(self) -> list[str]:
        return [*self.report_file, *self.protocol_file]


class BigWigFiles(BaseModel):
    assay: Assay
    names: list[str] = Field(default_factory=list)
    pileup_methods: list[PileupMethod]
    scale_methods: list[DataScalingTechnique] = [DataScalingTechnique.UNSCALED]
    spikein_methods: list[SpikeInMethod] = Field(default_factory=list)
    output_dir: str = "seqnado_output"

    @property
    def prefix(self) -> str:
        return f"{self.output_dir}/bigwigs/"

    @property
    def is_rna(self) -> bool:
        return self.assay == Assay.RNA

    @property
    def incompatible_methods(self) -> dict[PileupMethod, list[DataScalingTechnique]]:
        return {
            PileupMethod.HOMER: [
                DataScalingTechnique.CSAW,
                DataScalingTechnique.SPIKEIN,
            ],
            PileupMethod.BAMNADO: [
                DataScalingTechnique.CSAW,
                DataScalingTechnique.SPIKEIN,
            ],
            PileupMethod.DEEPTOOLS: [Assay.MCC],
        }

    def _is_compatible(
        self, method: PileupMethod, scale: DataScalingTechnique, assay: Assay
    ) -> bool:
        if scale in self.incompatible_methods.get(method, []):
            return False
        if assay in self.incompatible_methods.get(method, []):
            return False
        return True

    def generate_bigwig_paths(self) -> list[str]:
        paths = []

        for method in self.pileup_methods:
            for scale in self.scale_methods:
                if not self._is_compatible(method, scale, self.assay):
                    continue

                # For spike-in scaling, iterate over all spikein methods
                if scale == DataScalingTechnique.SPIKEIN and self.spikein_methods:
                    for spikein_method in self.spikein_methods:
                        scale_path = f"{scale.value}/{spikein_method.value}"
                        if self.is_rna:
                            for name in self.names:
                                for strand in ["plus", "minus"]:
                                    path = f"{self.prefix}{method.value}/{scale_path}/{name}_{strand}.bigWig"
                                    paths.append(path)
                        else:
                            for name in self.names:
                                path = f"{self.prefix}{method.value}/{scale_path}/{name}.bigWig"
                                paths.append(path)
                else:
                    scale_path = scale.value
                    if self.is_rna:
                        for name in self.names:
                            for strand in ["plus", "minus"]:
                                path = f"{self.prefix}{method.value}/{scale_path}/{name}_{strand}.bigWig"
                                paths.append(path)
                    else:
                        for name in self.names:
                            path = (
                                f"{self.prefix}{method.value}/{scale_path}/{name}.bigWig"
                            )
                            paths.append(path)

        return paths

    @computed_field
    @property
    def files(self) -> list[str]:
        return self.generate_bigwig_paths()


class PeakCallingFiles(BaseModel):
    assay: Assay
    names: list[str]
    peak_calling_method: list[PeakCallingMethod]
    output_dir: str = "seqnado_output"
    is_merged: bool = False

    @property
    def prefix(self) -> str:
        return f"{self.output_dir}/peaks/"

    @field_validator("assay")
    def validate_assay(cls, value):
        if value not in AssaysWithPeakCalling:
            raise ValueError(f"Invalid assay for peak calling: {value}")
        return value

    @property
    def peak_files(self) -> list[str]:
        return expand(
            self.prefix + "{method}/{sample}.bed"
            if not self.is_merged
            else self.prefix + "{method}/merged/{sample}.bed",
            sample=self.names,
            method=[m.value for m in self.peak_calling_method],
        )

    @computed_field
    @property
    def files(self) -> list[str]:
        """Return a list of peak calling files."""
        return self.peak_files


class MotifFiles(BaseModel):
    """Motif analysis output files for ChIP-seq, ATAC-seq, and CUT&Tag assays."""

    assay: Assay
    names: list[str]
    peak_calling_method: list[PeakCallingMethod]
    output_dir: str = "seqnado_output"
    run_homer: bool = True
    run_meme: bool = False

    @field_validator("assay")
    def validate_assay(cls, value):
        if value not in AssaysWithPeakCalling:
            raise ValueError(f"Invalid assay for motif analysis: {value}")
        return value

    @property
    def homer_files(self) -> list[str]:
        """Return HOMER motif analysis output files."""
        if not self.run_homer:
            return []
        return expand(
            f"{self.output_dir}/motifs/homer/{{method}}/{{sample}}/.completed",
            sample=self.names,
            method=[m.value for m in self.peak_calling_method],
        )

    @property
    def meme_files(self) -> list[str]:
        """Return MEME-ChIP motif analysis output files."""
        if not self.run_meme:
            return []
        return expand(
            f"{self.output_dir}/motifs/meme/{{method}}/{{sample}}/.completed",
            sample=self.names,
            method=[m.value for m in self.peak_calling_method],
        )

    @computed_field
    @property
    def files(self) -> list[str]:
        """Return all motif analysis files."""
        return [*self.homer_files, *self.meme_files]


class HeatmapFiles(BaseModel):
    assay: Assay

    @field_validator("assay")
    def validate_assay(cls, value):
        if value not in AssaysWithHeatmaps:
            raise ValueError(f"Invalid assay for heatmap: {value}")
        return value

    output_dir: str = "seqnado_output"

    @property
    def heatmap_files(self) -> list[str]:
        return [
            f"{self.output_dir}/heatmap/heatmap.pdf",
            f"{self.output_dir}/heatmap/metaplot.pdf",
        ]

    @computed_field
    @property
    def files(self) -> list[str]:
        return self.heatmap_files


class HubFiles(BaseModel):
    hub_dir: Path
    hub_name: str

    @property
    def hub_txt(self) -> Path:
        return self.hub_dir / f"{self.hub_name}.hub.txt"

    @computed_field
    @property
    def files(self) -> List[str]:
        return [str(self.hub_txt)]


class SpikeInFiles(BaseModel):
    assay: Assay
    names: list[str]
    method: list[SpikeInMethod]
    output_dir: str = "seqnado_output"

    @field_validator("assay")
    def validate_assay(cls, value):
        if value not in AssaysWithSpikein:
            raise ValueError(f"Invalid assay for spike-in: {value}")
        return value

    @property
    def norm_factors(self) -> list[str]:
        return expand(
            f"{self.output_dir}/resources/{{method}}/normalisation_factors.tsv",
            method=[m.value for m in self.method],
        )

    @computed_field
    @property
    def files(self) -> list[str]:
        return self.norm_factors


class PlotFiles(BaseModel):
    coordinates: Path
    file_format: Literal["svg", "png", "pdf"] = "svg"
    output_dir: str = "seqnado_output"

    @property
    def plot_names(self):
        import pandas as pd

        plots = []

        try:
            # Read BED file using pandas (pyranges replacement)
            bed_columns = ["Chromosome", "Start", "End", "Name", "Score", "Strand"]
            coords_df = pd.read_csv(
                str(self.coordinates), sep="\t", header=None, comment="#"
            )
            coords_df.columns = bed_columns[: len(coords_df.columns)]

            outdir = Path(f"{self.output_dir}/genome_browser_plots/")
            for region in coords_df.itertuples():
                fig_name = (
                    f"{region.Chromosome}-{region.Start}-{region.End}"
                    if not hasattr(region, "Name") or not region.Name
                    else region.Name
                )
                plots.append(outdir / f"{fig_name}.{self.file_format}")

        except FileNotFoundError:
            logger.warning(
                f"Could not find plotting coordinates file: {self.coordinates}"
            )

        return plots

    @computed_field
    @property
    def files(self) -> List[str]:
        """Return a list of plot files."""
        return [str(p) for p in self.plot_names]


class SNPFilesRaw(BaseModel):
    assay: Assay
    names: list[str]
    output_dir: str = "seqnado_output"

    @property
    def snp_files(self) -> list[str]:
        return expand(
            f"{self.output_dir}/variant/{{sample}}.vcf.gz",
            sample=self.names,
        )

    @computed_field
    @property
    def files(self) -> list[str]:
        return self.snp_files


class SNPFilesAnnotated(BaseModel):
    assay: Assay
    names: list[str]
    output_dir: str = "seqnado_output"

    @property
    def anno_snp_files(self) -> list[str]:
        return expand(
            f"{self.output_dir}/variant/{{sample}}.anno.vcf.gz",
            sample=self.names,
        )

    @computed_field
    @property
    def files(self) -> list[str]:
        return self.anno_snp_files


class MethylationFiles(BaseModel):
    assay: Assay
    names: list[str]
    genomes: List[str]
    method: MethylationMethod
    output_dir: str = "seqnado_output"

    @property
    def prefix(self) -> str:
        return f"{self.output_dir}/methylation/"

    @property
    def split_bams_files(self) -> List[str]:
        return expand(
            f"{self.output_dir}/aligned/spikein/{{sample}}_{{genome}}.bam",
            sample=self.names,
            genome=self.genomes,
        )

    @property
    def methyldackel_files(self) -> List[str]:
        """
        For TAPS, return the inverted file ({sample}_{genome}_CpG_inverted.bedGraph).
        For other methods, return the direct file ({sample}_{genome}_CpG.bedGraph).
        """
        if self.method.value == "taps":
            file_pattern = f"{self.output_dir}/methylation/methyldackel/{{sample}}_{{genome}}_CpG_inverted.bedGraph"
        else:
            file_pattern = f"{self.output_dir}/methylation/methyldackel/{{sample}}_{{genome}}_CpG.bedGraph"
        return expand(
            file_pattern,
            sample=self.names,
            genome=self.genomes,
        )

    @property
    def methylation_bias(self) -> List[str]:
        """Return the methylation bias files."""
        files = []
        files.append(f"{self.output_dir}/methylation/methylation_conversion.tsv")
        files.extend(
            expand(
                f"{self.output_dir}/methylation/methyldackel/bias/{{sample}}_{{genome}}.txt",
                sample=self.names,
                genome=self.genomes,
            )
        )
        return files

    @computed_field
    @property
    def files(self) -> list[str]:
        """Return a list of methylation files."""
        return [
            *self.split_bams_files,
            *self.methyldackel_files,
            *self.methylation_bias,
        ]


class CRISPRFiles(BaseModel):
    """Returns a list of CRISPR-related files including MAGeCK outputs."""

    use_mageck: bool = False
    output_dir: str = "seqnado_output"

    @property
    def feature_counts_files(self) -> list[str]:
        """Return featureCounts files (always generated)."""
        return [
            f"{self.output_dir}/readcounts/feature_counts/read_counts.tsv",
        ]

    @property
    def mageck_files(self) -> list[str]:
        """Return MAGeCK files if enabled for CRISPR assays."""
        if not self.use_mageck:
            return []

        return [
            f"{self.output_dir}/readcounts/mageck/mageck_count.count.txt",
            f"{self.output_dir}/readcounts/mageck/mageck_mle.gene_summary.txt",
            f"{self.output_dir}/readcounts/mageck/mageck_mle.sgrna_summary.txt",
        ]

    @computed_field
    @property
    def files(self) -> list[str]:
        """Return a list of CRISPR-specific files."""
        return [*self.feature_counts_files, *self.mageck_files]


class BigBedFiles(BaseModel):
    bed_files: list[Path] = Field(default_factory=list)

    @field_validator("bed_files", mode="before")
    def validate_bed_files(cls, v: list[Path | str]) -> list[Path]:
        return [Path(f) for f in v]

    @computed_field
    @property
    def files(self) -> list[str]:
        """Return a list of bigBed files."""
        return [str(f.with_suffix(".bb")) for f in self.bed_files if f.suffix == ".bed"]


class ContactFiles(BaseModel):
    assay: Assay
    names: list[str]
    output_dir: str = "seqnado_output"

    @computed_field
    @property
    def files(self) -> List[str]:
        """Return a list of contact files."""
        # Simplified to just return basic files based on sample names
        # since design_dataframe and viewpoints_grouped are not provided
        return [
            f"{self.output_dir}/contacts/{name}/{name}.mcool" for name in self.names
        ]


class QuantificationFiles(BaseModel):
    """Base class for quantification files."""

    assay: Assay
    methods: list[QuantificationMethod] = Field(default_factory=list)
    names: list[str]
    groups: SampleGroups
    output_dir: str = "seqnado_output"

    class Config:
        arbitrary_types_allowed = True

    @property
    def prefix(self) -> str:
        return f"{self.output_dir}/readcounts"

    @field_validator("methods", mode="before")
    def validate_methods_and_assays(
        cls, v: list[QuantificationMethod], info
    ) -> list[QuantificationMethod]:
        # In Pydantic v2, use info.data to access other fields
        assay = info.data.get("assay") if hasattr(info, "data") else None
        if assay == Assay.RNA:
            return [
                m
                for m in v
                if m
                in [QuantificationMethod.FEATURE_COUNTS, QuantificationMethod.SALMON]
            ]
        return [m for m in v if m == QuantificationMethod.FEATURE_COUNTS]

    @property
    def combined_counts_file(self) -> list[str]:
        """Return the combined read counts file."""
        files = []
        for m in self.methods:
            if m == QuantificationMethod.FEATURE_COUNTS:
                files.append(f"{self.prefix}/feature_counts/read_counts.tsv")
            elif m == QuantificationMethod.SALMON:
                files.append(f"{self.prefix}/salmon/salmon_counts.csv")
        return files

    @property
    def grouped_counts_files(self) -> list[str]:
        """Return the grouped read counts files."""
        files = []
        for group in self.groups.groups:
            files.extend(
                expand(
                    f"{self.prefix}{group}/read_counts.tsv",
                    group=group,
                    methods=[m.value for m in self.methods],
                )
            )
        return files

    @computed_field
    @property
    def files(self) -> list[str]:
        """Return a list of quantification files."""
        return [*self.combined_counts_file, *self.grouped_counts_files]


class GeoSubmissionFiles(BaseModel):
    """Class to handle files for GEO submission."""

    assay: Assay
    names: list[str]
    seqnado_files: list[str] = Field(
        default_factory=list,
        description="Unfiltered list of files. Will be filtered based on allowed extensions and added.",
    )
    allowed_extensions: list[str] = Field(
        default_factory=lambda: [".bigWig", ".bed", ".tsv", ".vcf.gz"]
    )
    output_dir: str = "seqnado_output"
    samples: Any = (
        None  # Optional: FastqCollection or FastqCollectionForIP to check if paired-end
    )

    class Config:
        arbitrary_types_allowed = True

    @property
    def default_files(self):
        return [
            f"{self.output_dir}/geo_submission/md5sums.txt",
            f"{self.output_dir}/geo_submission/raw_data_checksums.txt",
            f"{self.output_dir}/geo_submission/processed_data_checksums.txt",
            f"{self.output_dir}/geo_submission/samples_table.txt",
        ]

    @property
    def upload_directory(self) -> Path:
        return Path(f"{self.output_dir}/geo_submission") / self.assay.clean_name

    @property
    def upload_instructions(self) -> Path:
        return Path(f"{self.output_dir}/geo_submission") / "upload_instructions.txt"

    @property
    def raw_files(self) -> list[str]:
        """Return FASTQ files for raw data."""
        base_dir = Path(f"{self.output_dir}/geo_submission")
        files = []

        for sample in self.names:
            # Check if this sample is paired-end
            is_paired = True  # Default to paired-end for backwards compatibility
            if self.samples is not None and hasattr(self.samples, "is_paired_end"):
                try:
                    is_paired = self.samples.is_paired_end(sample)
                except Exception:
                    # If we can't determine, assume paired-end
                    is_paired = True

            if is_paired:
                # Paired-end: add both R1 and R2
                files.extend(
                    [
                        str(base_dir / f"{sample}_1.fastq.gz"),
                        str(base_dir / f"{sample}_2.fastq.gz"),
                    ]
                )
            else:
                # Single-end: add only the sample file
                files.append(str(base_dir / f"{sample}.fastq.gz"))

        return files

    @property
    def processed_data_files(self) -> list[str]:
        """Return processed files for GEO submission."""
        base_dir = Path(f"{self.output_dir}/geo_submission")
        files = []
        for file in self.seqnado_files:
            file_path = Path(file)
            # Check both single suffix and compound suffix (e.g., .vcf.gz)
            ext = "".join(file_path.suffixes)
            if (
                ext in self.allowed_extensions
                or file_path.suffix in self.allowed_extensions
            ):
                # Need to flatten the file un-nest the directory structure
                # e.g. seqnado_output/bigwigs/METHOD/SCALE/NAME.bigWig
                basename = file_path.stem
                scale_method = file_path.parent.name
                method = file_path.parent.parent.name

                files.append(str(base_dir / f"{basename}_{method}_{scale_method}{ext}"))

        return files

    @computed_field
    @property
    def files(self) -> list[str]:
        """Return a list of all files for GEO submission."""
        return [*self.default_files, *self.raw_files, *self.processed_data_files]

    @property
    def metadata_files(self) -> list[str]:
        """Return only the metadata files that are declared as explicit rule outputs.
        
        Note: The individual symlinked FASTQ and processed files are created by geo_symlink rule
        but are not included here as they cannot be declared as dynamic rule outputs in Snakemake.
        """
        return self.default_files


class GEOFiles(BaseModel):
    """
    Class to manage GEO submission files including raw FASTQ files and processed data files.

    This class handles the renaming and organization of files for GEO submission by:
    - Flattening nested directory structures
    - Extracting method and scale information from file paths
    - Creating appropriate symlink names
    """

    make_geo_submission_files: bool
    assay: Assay
    design: Any  # pd.DataFrame, but use Any to avoid circular import
    sample_names: List[str]
    config: Any  # SeqnadoConfig, but use Any to avoid issues
    processed_files: List[str] = Field(default_factory=list)

    class Config:
        arbitrary_types_allowed = True

    @property
    def raw_files(self) -> dict[str, List[str]]:
        """
        Get raw FASTQ files for each sample.

        Returns:
            dict[str, List[str]]: Dictionary mapping sample names to lists of FASTQ filenames
        """
        raw_files_dict = {}

        for sample in self.sample_names:
            # Check if sample is paired-end by looking at the design DataFrame
            is_paired = False
            if self.design is not None and not (
                hasattr(self.design, "empty") and self.design.empty
            ):
                try:
                    import pandas as pd

                    # For IP-based assays (ChIP, CUT&TAG), sample names include IP/control suffix
                    # e.g., "chip-rx_MLL" where sample_id is "chip-rx" and ip is "MLL"
                    # Try to find the matching row by checking if sample contains sample_id
                    sample_row = None

                    # First try exact match with sample_id
                    if "sample_id" in self.design.columns:
                        sample_row = self.design[self.design["sample_id"] == sample]

                        # If no exact match, try to match the base part (for IP assays)
                        if sample_row.empty:
                            # Try to find sample_id that is a prefix of the sample name
                            for _, row in self.design.iterrows():
                                sid = str(row["sample_id"])
                                if sample.startswith(sid + "_") or sample == sid:
                                    sample_row = pd.DataFrame([row])
                                    break

                    if sample_row is not None and not sample_row.empty:
                        # Check r2 column to determine if paired-end
                        if "r2" in self.design.columns:
                            r2_value = sample_row["r2"].iloc[0]
                            is_paired = (
                                pd.notna(r2_value) and str(r2_value).strip() != ""
                            )
                        else:
                            # No r2 column, assume paired-end for safety
                            is_paired = True
                    else:
                        # Can't find sample in design, assume paired-end for backwards compatibility
                        is_paired = True
                except Exception:
                    # If we can't determine, assume paired-end for backwards compatibility
                    is_paired = True
            else:
                # If no design available, assume paired-end for backwards compatibility
                is_paired = True

            if is_paired:
                # Paired-end: return both R1 and R2
                raw_files_dict[sample] = [
                    f"{sample}_1.fastq.gz",
                    f"{sample}_2.fastq.gz",
                ]
            else:
                # Single-end: return only the sample file
                raw_files_dict[sample] = [
                    f"{sample}.fastq.gz",
                ]

        return raw_files_dict

    @property
    def processed_data_files(self):
        """
        Get processed data files with renamed paths for GEO submission.

        This method:
        1. Filters files by allowed extensions (.bigWig, .bed, .tsv, .vcf.gz)
        2. Extracts method and scale information from the directory structure
        3. Creates flat filenames incorporating the metadata

        Returns:
            pd.DataFrame: DataFrame with columns 'path' (original) and 'output_file_name' (for GEO)
        """
        import pandas as pd

        if not self.make_geo_submission_files:
            return pd.DataFrame(columns=["path", "output_file_name"])

        allowed_extensions = [".bigWig", ".bed", ".tsv", ".vcf.gz"]
        processed_data = []

        for file_path_str in self.processed_files:
            file_path = Path(file_path_str)

            # Check if file has an allowed extension
            # Handle both single suffix (.bed) and compound suffix (.vcf.gz)
            ext = "".join(file_path.suffixes)
            if (
                ext not in allowed_extensions
                and file_path.suffix not in allowed_extensions
            ):
                continue

            # Extract metadata from path structure
            parts = file_path.parts
            basename = file_path.stem

            # Handle compound extensions like .vcf.gz
            if len(file_path.suffixes) > 1:
                ext = "".join(file_path.suffixes)
            else:
                ext = file_path.suffix

            # Try to extract method and scale from directory structure
            try:
                # Find the index of known directories
                if "bigwigs" in parts:
                    bigwigs_idx = parts.index("bigwigs")
                    method = (
                        parts[bigwigs_idx + 1] if len(parts) > bigwigs_idx + 1 else None
                    )
                    scale = (
                        parts[bigwigs_idx + 2] if len(parts) > bigwigs_idx + 2 else None
                    )

                    if method and scale:
                        output_name = f"{basename}_{method}_{scale}{ext}"
                    elif method:
                        output_name = f"{basename}_{method}{ext}"
                    else:
                        output_name = f"{basename}{ext}"

                elif "peaks" in parts:
                    peaks_idx = parts.index("peaks")
                    method = (
                        parts[peaks_idx + 1] if len(parts) > peaks_idx + 1 else None
                    )

                    # Check if this is a merged peak file
                    if "merged" in parts:
                        if method:
                            output_name = f"{basename}_{method}_merged{ext}"
                        else:
                            output_name = f"{basename}_merged{ext}"
                    else:
                        if method and len(parts) > peaks_idx + 2:
                            output_name = f"{basename}_peaks_{method}{ext}"
                        elif method:
                            output_name = f"{basename}_{method}{ext}"
                        else:
                            output_name = f"{basename}{ext}"

                elif "quantification" in parts:
                    quant_idx = parts.index("quantification")
                    method = (
                        parts[quant_idx + 1] if len(parts) > quant_idx + 1 else None
                    )

                    if method:
                        output_name = f"{basename}_{method}{ext}"
                    else:
                        output_name = f"{basename}{ext}"

                elif "variant" in parts:
                    # Variant files are typically just {NAME}.vcf.gz
                    output_name = f"{basename}{ext}"

                else:
                    # Fallback: just use the basename
                    output_name = f"{basename}{ext}"

            except (ValueError, IndexError):
                # If we can't parse the path, just use the basename
                output_name = f"{basename}{ext}"

            processed_data.append(
                {"path": str(file_path), "output_file_name": output_name}
            )

        return pd.DataFrame(processed_data)

    @property
    def metadata(self):
        """
        Get GEO submission metadata table.

        Returns:
            pd.DataFrame: Metadata table formatted for GEO submission
        """
        import pandas as pd

        if self.design is None or (hasattr(self.design, "empty") and self.design.empty):
            # Create minimal metadata if design is not available
            metadata = pd.DataFrame(
                {
                    "Sample name": self.sample_names,
                    "title": self.sample_names,
                    "source name": [f"Sample {name}" for name in self.sample_names],
                    "organism": [""] * len(self.sample_names),
                }
            )
        else:
            # Use the existing design dataframe as base
            metadata = self.design.copy()

            # Ensure required columns exist
            if (
                "Sample name" not in metadata.columns
                and "samplename" in metadata.columns
            ):
                metadata["Sample name"] = metadata["samplename"]
            elif "Sample name" not in metadata.columns:
                # If the design has fewer rows than sample_names (e.g., due to merging/consensus),
                # we need to create a row for each sample
                if len(metadata) != len(self.sample_names):
                    # Create a new metadata DataFrame with one row per sample
                    metadata = pd.DataFrame(
                        {
                            "Sample name": self.sample_names,
                        }
                    )
                else:
                    metadata["Sample name"] = self.sample_names

        return metadata
