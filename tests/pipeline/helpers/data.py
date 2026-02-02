"""Test data management for SeqNado pipeline tests."""

from __future__ import annotations

import json
from pathlib import Path
from typing import ClassVar

from pydantic import BaseModel, model_validator

from .utils import download_with_retry, extract_tar, get_fastq_pattern


class GenomeResources(BaseModel):
    """Validated genome resources for an assay."""

    assay: str | None = None  # Used for validation
    star_index: Path | None = None
    bt2_index: Path | None = None
    chromosome_sizes: Path | None = None
    gtf: Path | None = None
    blacklist: Path | None = None
    genes: Path | None = None
    fasta: Path | None = None
    viewpoints: Path | None = None
    plot_coords: Path | None = None

    # Configuration constants
    REFERENCE_URL: ClassVar[str] = (
        "https://userweb.molbiol.ox.ac.uk/public/project/milne_group/seqnado/test_data/reference"
    )

    # Base resources shared by multiple assays

    _DEFAULT_RESOURCES: ClassVar[dict] = {
        "bt2_index": ("bt2_chr21", "bt2_chr21", "hg38"),
        "chromosome_sizes": "chr21.chrom.sizes",
        "gtf": "chr21.gtf",
        "blacklist": "hg38_chr21-blacklist.bed",
        "genes": "hg38_genes.bed",
        "fasta": "chr21.fa",
        "plot_coords": "plotting_coordinates.bed",
    }

    # Resources keyed directly by assay name
    RESOURCES: ClassVar[dict] = {
        "atac": _DEFAULT_RESOURCES,
        "cat": _DEFAULT_RESOURCES,
        "chip": _DEFAULT_RESOURCES,
        "chip-rx": {
            **_DEFAULT_RESOURCES,
            "bt2_index": ("bt2_chr21_dm6_chr2L", "bt2_chr21_dm6_chr2L", "hg38"),
            "chromosome_sizes": "chr21_dm6.fa.fai",
        },
        "crispr": {
            **_DEFAULT_RESOURCES, 
            "bt2_index": (
                "bt2_TKOv3_guides_chr21",
                "TKOv3_guides_chr21",
                "TKOv3_guides",
            ),
            "chromosome_sizes": "TKOv3_guides_chr21.fasta.fai",
            "gtf": "TKOv3_guides_chr21.saf",
        },
        "mcc": {
            **_DEFAULT_RESOURCES, 
            "viewpoints": "mcc_viewpoints.bed",
        },
        "meth": {
            **_DEFAULT_RESOURCES,
            "bt2_index": ("bt2_chr21_meth", "chr21_meth", "hg38_meth"),
            "chromosome_sizes": "chr21_meth.fa.fai",
            "fasta": "chr21_meth.fa",
        },
        "rna": {
            **_DEFAULT_RESOURCES,
            "star_index": "STAR_chr21_rna_spikein",
            "gtf": "chr21_rna_spikein.gtf",
        },
        "snp": _DEFAULT_RESOURCES,
    }

    @model_validator(mode="after")
    def validate_required_fields(self):
        """Validate that required fields for this assay are present based on its template."""
        if not self.assay:
            return self

        # Prefer exact match, then longest substring match (e.g., 'chip-rx' before 'chip')
        assay_lower = self.assay.lower()
        assay_key = next(
            (k for k in sorted(self.RESOURCES, key=len, reverse=True) if k in assay_lower),
            None
        )
        if not assay_key:
            return self

        for field in self.RESOURCES[assay_key]:
            if getattr(self, field, None) is None:
                raise ValueError(f"{field} required for {self.assay} but not found")
        return self

    def write_config(self, config_file: Path) -> None:
        """Write genome_config.json from this instance."""
        config_file.parent.mkdir(parents=True, exist_ok=True)

        existing = {}
        if config_file.exists():
            with open(config_file) as f:
                existing = json.load(f)

        existing[self.assay] = {
            k: str(v) if v else None
            for k, v in self.model_dump(exclude={"assay"}).items()
        }

        with open(config_file, "w") as f:
            json.dump(existing, f, indent=2)

    @classmethod
    def download_resources(cls, genome_path: Path, assay: str) -> GenomeResources:
        """Download genome resources for an assay.

        Args:
            genome_path: Path to store genome files
            assay: Assay type (e.g., 'rna', 'chip', 'atac')

        Returns:
            GenomeResources: Validated genome resources

        Raises:
            FileNotFoundError: If required resources cannot be obtained
        """
        genome_path.mkdir(parents=True, exist_ok=True)

        # Look up resources for this assay
        # Prefer exact match, then longest substring match (e.g., 'chip-rx' before 'chip')
        assay_lower = assay.lower()
        assay_key = next(
            (k for k in sorted(cls.RESOURCES, key=len, reverse=True) if k in assay_lower),
            "chip"
        )
        template = cls.RESOURCES[assay_key]
        ref_url = cls.REFERENCE_URL
        config: dict[str, Path | None] = {}

        for key, value in template.items():
            if key == "bt2_index":
                bt2_dir, bt2_prefix, db_name = value
                config[key] = cls._ensure_index(
                    genome_path, ref_url, bt2_dir, prefix=bt2_prefix
                )
                if config[key]:
                    cls._update_fastq_screen_config(genome_path, config[key], db_name)
            elif key == "star_index":
                config[key] = cls._ensure_index(genome_path, ref_url, value)
            else:
                dest = genome_path / value
                if not dest.exists():
                    download_with_retry(f"{ref_url}/{value}", dest)
                config[key] = dest

        # Create and validate instance (validation happens automatically via model_validator)
        return cls(assay=assay, **config)

    @staticmethod
    def _ensure_index(
        genome_path: Path,
        ref_url: str,
        index_name: str,
        prefix: str | None = None,
    ) -> Path | None:
        """Download and extract an index (STAR or Bowtie2)."""
        is_bt2 = prefix is not None

        if is_bt2:
            index_dir = genome_path / index_name
            dest = index_dir / prefix
            if list(index_dir.glob(f"{prefix}.*.bt2*")):
                return dest
        else:
            dest = genome_path / index_name
            index_dir = dest
            if dest.exists():
                return dest

        tar_path = genome_path / f"{index_name}.tar.gz"
        try:
            download_with_retry(f"{ref_url}/{tar_path.name}", tar_path)
            if is_bt2:
                index_dir.mkdir(parents=True, exist_ok=True)
                extract_tar(tar_path, index_dir)
            else:
                extract_tar(tar_path, genome_path)
            tar_path.unlink(missing_ok=True)
        except Exception as e:
            print(f"[WARNING] Could not download index {index_name}: {e}")
            tar_path.unlink(missing_ok=True)
            return None

        # Verify and return
        if is_bt2:
            return dest if list(index_dir.glob(f"{prefix}.*.bt2*")) else None
        return dest if dest.exists() else None

    @staticmethod
    def _update_fastq_screen_config(
        genome_path: Path, bt2_index_path: Path, db_name: str
    ) -> None:
        """Append genome entry to fastq_screen.conf."""
        config_path = genome_path / "fastq_screen.conf"

        # Check if this database entry already exists
        existing_entries = set()
        if config_path.exists():
            with open(config_path) as f:
                for line in f:
                    if line.startswith("DATABASE"):
                        parts = line.strip().split("\t")
                        if len(parts) >= 2:
                            existing_entries.add(parts[1])

        if db_name not in existing_entries:
            with open(config_path, "a") as f:
                f.write(f"DATABASE\t{db_name}\t{bt2_index_path}\n")


class FastqFiles(BaseModel):
    """Manages FASTQ file downloads and organization by assay."""

    fastq_dir: Path
    files: dict[str, list[Path]] = {}

    # Configuration constants
    FQ_URL: ClassVar[str] = (
        "https://userweb.molbiol.ox.ac.uk/public/project/milne_group/seqnado/test_data/fastq.tar.gz"
    )

    @classmethod
    def download(cls, fastq_dir: Path, selected_assays: list[str]) -> FastqFiles:
        """Download FASTQ files and organize them by assay.

        Args:
            fastq_dir: Directory to store FASTQ files
            selected_assays: List of assay types to get FASTQs for

        Returns:
            FastqFiles: Instance with files organized by assay

        Raises:
            RuntimeError: If FASTQ extraction fails
        """
        fastq_dir.mkdir(parents=True, exist_ok=True)

        if not list(fastq_dir.glob("*.fastq.gz")):
            tar_path = fastq_dir.parent / "fastq.tar.gz"
            download_with_retry(
                cls.FQ_URL, tar_path
            )
            extract_tar(tar_path, fastq_dir.parent, flatten=False)

            nested = fastq_dir / "fastq"
            if nested.is_dir():
                for f in nested.glob("*.fastq.gz"):
                    f.rename(fastq_dir / f.name)
                nested.rmdir()

            tar_path.unlink(missing_ok=True)

            if not list(fastq_dir.glob("*.fastq.gz")):
                raise RuntimeError(f"FASTQ extraction failed - no files in {fastq_dir}")

        files = {
            a: sorted(fastq_dir.glob(get_fastq_pattern(a))) for a in selected_assays
        }
        return cls(fastq_dir=fastq_dir, files=files)
