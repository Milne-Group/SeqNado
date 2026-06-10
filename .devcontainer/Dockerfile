FROM mambaorg/micromamba:2-noble

# Build context is repo root (set in devcontainer.json) so we can share
# the pipeline environment directly instead of duplicating it.
COPY containers/pipeline/environment_pipeline.yml /opt/environment_pipeline.yml
COPY containers/pipeline/multiqc_config.yaml /opt/seqnado/multiqc_config.yaml
COPY containers/pipeline/seqnado.png /opt/seqnado/seqnado.png

ENV PATH="/opt/conda/bin:$PATH" \
    CONDA_DEFAULT_ENV=base \
    CONDA_PREFIX="/opt/conda" \
    LD_LIBRARY_PATH="/usr/lib:/lib" \
    MAMBA_ROOT_PREFIX="/opt/conda"

USER root

RUN set -eux && \
    apt-get update && DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
      ca-certificates curl gcc libc6-dev pkg-config libssl-dev && \
    rm -rf /var/lib/apt/lists/* && \
    micromamba install -y -n base -f /opt/environment_pipeline.yml && \
    apt-get purge -y gcc libc6-dev pkg-config libssl-dev && \
    apt-get autoremove -y && \
    apt-get clean && \
    find /opt/conda/ -type f -name '*.a' -delete && \
    find /opt/conda/ -type f -name '*.pyc' -delete && \
    find /opt/conda/ -name '__pycache__' -type d -exec rm -rf {} + && \
    rm -rf /root/.cache && \
    micromamba clean -ay

# Install uv
COPY --from=ghcr.io/astral-sh/uv:latest /uv /uvx /usr/local/bin/

WORKDIR /workspace
CMD ["/bin/bash"]
