FROM pytorch/pytorch:2.6.0-cuda12.4-cudnn9-runtime

LABEL Author="asmith" \
      Software="SeqNado" \
      Version="1.0"

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && \
    apt-get install -y --no-install-recommends python3-pip && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

RUN python -m pip install --no-cache-dir \
    git+https://github.com/alsmith151/lanceotron-mcc \
    lanceotron \
    quantnado

HEALTHCHECK --interval=1m --timeout=10s --start-period=30s --retries=3 \
  CMD python -c "import torch; assert torch.cuda.is_available(); import lanceotron" || exit 1

CMD ["/bin/bash"]
