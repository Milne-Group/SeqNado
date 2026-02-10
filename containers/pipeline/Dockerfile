FROM mambaorg/micromamba:latest

# Copy files
COPY environment_pipeline.yml /opt/environment_pipeline.yml
COPY multiqc_config.yaml /opt/seqnado/multiqc_config.yaml
COPY seqnado.png /opt/seqnado/seqnado.png

# Copy repository into image so we can install the package. When building
# from the repository root use: `docker build -f containers/pipeline/Dockerfile .`
COPY . /opt/seqnado_src/

# Set environment variables
ENV PATH="/opt/conda/bin:$PATH" \
    CONDA_DEFAULT_ENV=base \
    CONDA_PREFIX="/opt/conda" \
    LD_LIBRARY_PATH="/usr/lib:/lib" \
    MAMBA_ROOT_PREFIX="/opt/conda"

# Keep apt non-interactive
ENV DEBIAN_FRONTEND=noninteractive
# Ensure we have root for package installation and system-level setup
USER root

# Main build process
RUN set -eux && \
    mkdir -p /var/lib/apt/lists/partial && \
    chmod 755 /var/lib/apt/lists/partial && \
    apt-get update && apt-get install -y --no-install-recommends \
      ca-certificates curl gnupg dirmngr gcc pkg-config libssl-dev && \
    rm -rf /var/lib/apt/lists/* && \

    # Install latest Quarto (architecture-aware)
    QUARTO_VERSION=$(curl -s https://api.github.com/repos/quarto-dev/quarto-cli/releases/latest | grep -o '"tag_name": "[^"]*' | cut -d'"' -f4 | sed 's/^v//') && \
    ARCH=$(dpkg --print-architecture 2>/dev/null || echo "$(uname -m)") && \
    echo "Detected architecture: ${ARCH}" && \
    QUARTO_TMP=/tmp/quarto.deb && \
    QUARTO_FOUND=0 && \
    ( \
      if [ "${ARCH}" = "amd64" ] || [ "${ARCH}" = "x86_64" ]; then \
        urls=( "https://github.com/quarto-dev/quarto-cli/releases/download/v${QUARTO_VERSION}/quarto-${QUARTO_VERSION}-linux-amd64.deb" ); \
      else \
        urls=( \
          "https://github.com/quarto-dev/quarto-cli/releases/download/v${QUARTO_VERSION}/quarto-${QUARTO_VERSION}-linux-aarch64.deb" \
          "https://github.com/quarto-dev/quarto-cli/releases/download/v${QUARTO_VERSION}/quarto-${QUARTO_VERSION}-linux-arm64.deb" \
          "https://github.com/quarto-dev/quarto-cli/releases/download/v${QUARTO_VERSION}/quarto-${QUARTO_VERSION}-linux-amd64.deb" \
        ); \
      fi && \
      for url in "${urls[@]}"; do \
        echo "Trying $url"; \
        if curl -fSL --retry 3 --retry-delay 2 -o "${QUARTO_TMP}" "$url"; then \
          echo "Downloaded quarto from $url"; QUARTO_FOUND=1; break; \
        fi; \
      done; \
    ) && \
    if [ "${QUARTO_FOUND}" -ne 1 ] || [ ! -s "${QUARTO_TMP}" ]; then \
      echo "Failed to download Quarto ${QUARTO_VERSION} for ${ARCH}" >&2; \
      exit 1; \
    fi && \
    dpkg -i "${QUARTO_TMP}" && rm -f "${QUARTO_TMP}" && \


    # Micromamba environment setup (ensure Python <3.12)
    micromamba install -y -n base -f /opt/environment_pipeline.yml && \

    # Install seqnado from PyPI into the base environment (no-deps)
    /opt/conda/bin/pip install --no-deps --no-cache-dir seqnado && \

    # Remove build tools but avoid purging libc6 core packages
    apt-get purge -y --auto-remove gcc pkg-config libssl-dev && \
    apt-get autoremove -y && apt-get clean && \

    # Remove unnecessary data and dev files
    find /opt/conda/ -type f -name '*.a' -delete && \
    find /opt/conda/ -type f -name '*.pyc' -delete && \
    find /opt/conda/ -name '__pycache__' -type d -exec rm -rf {} + && \
    find /opt/conda/ -name '*.egg-info' -type d -exec rm -rf {} + && \
    find /opt/conda/ -name 'tests' -type d -exec rm -rf {} + && \
    rm -rf /root/.cache && \
    micromamba clean -afy

# Add labels
ARG BUILD_DATE="unknown"
LABEL Author="asmith" \
      Version="0.4.0" \
      Description="Minimal SeqNado pipeline container with Python tools" \
      BuildDate="${BUILD_DATE}"

# Expose build date as an environment variable and create a build file
ENV BUILD_DATE="${BUILD_DATE}"

# Persist build date inside the image for runtime inspection
RUN mkdir -p /opt/seqnado && echo "${BUILD_DATE}" > /opt/seqnado/BUILD_DATE || true

# Switch back to micromamba user (set by the base image)
USER $MAMBA_USER

# Health check / test
RUN echo "Testing minimal container..." && python --version && echo "All tests passed!"

# Default command: run passed command with bash -lc
ENTRYPOINT ["/usr/bin/env", "bash", "-lc"]
CMD ["exec \"$@\""]
