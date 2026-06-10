# Dockerfile for SeqNado ML CPU
FROM python:3.12-slim

LABEL author="asmith"
LABEL software="SeqNado"
LABEL version="1.0"

# Copy repositories
COPY lanceotron-mcc/ /opt/repos/lanceotron-mcc/

# Upgrade pip, setuptools, and wheel
RUN apt-get update && apt-get install -y build-essential zlib1g-dev curl && rm -rf /var/lib/apt/lists/* \
    && python -m pip install --no-cache-dir --upgrade pip setuptools wheel

# Install PyTorch and related packages from PyTorch index
RUN python -m pip install --no-cache-dir --upgrade --no-input \
    torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cpu

# Install other Python dependencies
RUN python -m pip install --no-cache-dir \
    loguru \
    pandas \
    polars \
    pydantic \
    pyBigWig \
    tqdm \
    typer \
    cyclopts \
    joblib \
    pyranges \
    crested \
    modisco-lite \
    natsort \
    tensorflow

# Install local packages
RUN python -m pip install --no-cache-dir --no-deps \
    git+https://github.com/alsmith151/lanceotron-mcc \
    lanceotron

# Clean up to reduce image size
RUN rm -rf /opt/repos && \
    python -m pip cache purge && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Default runscript
CMD ["/bin/bash"]

# Healthcheck (optional, basic import test)
HEALTHCHECK --interval=1m --timeout=10s --start-period=30s --retries=3 \
  CMD python -c "import lanceotron, quantnado, torch, tensorflow, lanceotron_mcc" || exit 1
