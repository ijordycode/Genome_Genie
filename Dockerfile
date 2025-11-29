FROM continuumio/miniconda3

COPY src/base.yaml /BIP/base.yaml
COPY src/rbio.yaml /BIP/rbio.yaml

# Create and clean the conda environment
RUN conda env create --name seq -f /BIP/base.yaml && \
    conda clean --all --yes

RUN conda env create --name rbio -f /BIP/rbio.yaml && \
    conda clean --all --yes

COPY Snakefile /BIP/Snakefile
COPY config.yaml /BIP/config.yaml

# Activate environment when container starts
SHELL ["conda", "run", "-n", "sqe", "/bin/bash", "-c"]

# Set working directory
WORKDIR /BIP

# Ensure conda environment is activated in interactive shells
SHELL ["/bin/bash", "-c"]
RUN echo "conda activate seq" >> ~/.bashrc

# Working CMD to keep container open
CMD ["tail", "-f", "/dev/null"]

# Run snakemake and keep container alive for debugging
# CMD ["bash", "-c", "conda activate bio && snakemake --cores 1; tail -f /dev/null"]
