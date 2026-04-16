# One-shot image: install PHENOMS + build Rust extension, run kernel benchmark.
# Usage (250 frames, default 4 CPU, etc.): see docker/README.md

FROM python:3.12-bookworm

RUN apt-get update && apt-get install -y --no-install-recommends \
    curl \
    build-essential \
    pkg-config \
    && rm -rf /var/lib/apt/lists/*

ENV RUSTUP_HOME=/usr/local/rustup \
    CARGO_HOME=/usr/local/cargo
ENV PATH="/usr/local/cargo/bin:${PATH}"
RUN curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y --default-toolchain stable

WORKDIR /app
COPY . .

# maturin develop requires an activated virtualenv (or .venv); system pip alone is not enough.
RUN python -m venv /opt/venv
ENV VIRTUAL_ENV=/opt/venv
ENV PATH="/opt/venv/bin:${PATH}"

RUN pip install --no-cache-dir --upgrade pip maturin \
    && pip install --no-cache-dir -e ".[benchmark]" \
    && maturin develop --release

# Belt-and-suspenders (script also setdefaults these before imports)
ENV OMP_NUM_THREADS=4 \
    OPENBLAS_NUM_THREADS=4 \
    MKL_NUM_THREADS=4 \
    NUMEXPR_NUM_THREADS=4 \
    VECLIB_MAXIMUM_THREADS=4

CMD ["python", "scripts/benchmark_kernel.py", "--sub-frames", "15", "--out-dir", "/out"]
