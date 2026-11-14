"""Pinned container images used by the SeqNado workflows.

Every container referenced by a Snakemake rule, ``tools.json`` or a workflow
script must be defined here, pinned to an immutable reference: a version tag
or an ``@sha256:`` digest. Floating tags such as ``:latest`` are not allowed.

Why: Snakemake caches each image as a ``.sif`` under ``.snakemake/singularity/``
keyed on the image URI. With a floating tag the URI never changes, so after a
new image is published Snakemake keeps running the stale cached copy against
the updated workflow code. Pinning means each SeqNado release names the exact
images it was tested with, and a new pin gets a new cache entry automatically.

To update an image, change its reference below. For images that only publish
``latest`` (seqnado_pipeline, mccnado, seqnado_ml_cpu), look up the digest of
the newly published image, e.g.::

    crane digest ghcr.io/alsmith151/seqnado_pipeline:latest

and paste it in. ``tests/unit/test_containers.py`` checks that no rule, script
or ``tools.json`` entry uses an image that is not listed here.
"""

from __future__ import annotations

from typing import Final

# General-purpose pipeline image (most rules run in this one).
SEQNADO_PIPELINE: Final[str] = (
    "oras://ghcr.io/alsmith151/seqnado_pipeline"
    "@sha256:af68bf364e6a507d81c1e98c2b159704404e3aa955f65bd1fecba87888bd74ae"
)

# Machine-learning image (LanceOtron peak calling, CPU build).
SEQNADO_ML_CPU: Final[str] = (
    "oras://ghcr.io/alsmith151/seqnado_ml_cpu"
    "@sha256:065f77bba773f15130aef1c30f020b5fdadbfd58decd7dbbf92c81ef2fe4857c"
)

# Independently versioned tool images.
BAMNADO: Final[str] = "docker://ghcr.io/alsmith151/bamnado:0.9.0"
MCCNADO: Final[str] = (
    "docker://ghcr.io/alsmith151/mccnado"
    "@sha256:bcc9deb2c21a912e8e5afd1ebd0b7719a140f50ed850cd072bc7b34c4a5aa8ba"
)
# Pinned by digest (= v0.4.9): the image reports `quantnado --version` as
# 0.0.0, so tag-based version detection is unreliable for it anyway.
QUANTNADO: Final[str] = (
    "docker://ghcr.io/milne-group/quantnado"
    "@sha256:b07a8d5d26d49fbe4baa4072735aae9bf30884b854f7735d2b15a8531e39c7f9"
)

PLOTNADO: Final[str] = "docker://ghcr.io/alsmith151/plotnado:0.3.3"

# TODO: pin lanceotron-mcc (GPU MCC peak calling, run via `apptainer exec --nv`
# in the call_mcc_peaks rule). Also Sylabs-hosted with only `latest`.
LANCEOTRON_MCC: Final[str] = "library://asmith151/lanceotron/lanceotron-mcc:latest"

# Third-party biocontainers (already immutable build tags).
MAGECK: Final[str] = "docker://quay.io/biocontainers/mageck:0.5.9.5--py310h184ae93_8"
MEME: Final[str] = "docker://quay.io/biocontainers/meme:5.5.9--pl5321h1ca524f_0"
SRA_TOOLS: Final[str] = "docker://quay.io/biocontainers/sra-tools:3.0.10--h9f5acd7_0"
MACS2: Final[str] = (
    "docker://quay.io/biocontainers/macs2:2.1.1.20160309--py27r3.3.1_1"
)
MACS3: Final[str] = "docker://quay.io/biocontainers/macs3:3.0.3--py39h0699b22_0"

#: All images SeqNado may run, keyed by constant name.
ALL_CONTAINERS: Final[dict[str, str]] = {
    name: value
    for name, value in dict(globals()).items()
    if name.isupper() and isinstance(value, str) and "://" in value
}

#: Images still on a floating tag. Must stay empty apart from known exceptions.
UNPINNED_ALLOWLIST: Final[frozenset[str]] = frozenset({LANCEOTRON_MCC})
