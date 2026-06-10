FROM python:3.12-slim

# Copy uv binary from official image
COPY --from=ghcr.io/astral-sh/uv:latest /uv /usr/local/bin/uv

# Install seqnado
RUN uv pip install --system seqnado

# Labels
LABEL Author="asmith" \
      Description="Minimal SeqNado CLI container"

ENTRYPOINT ["seqnado"]
CMD ["--help"]
