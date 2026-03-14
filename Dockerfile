FROM continuumio/miniconda3

COPY src/ /BIP/

# Create and clean the conda environment
RUN conda env create --name seq -f /BIP/base.yaml && \
    conda clean --all --yes

RUN conda env create --name rbio -f /BIP/rbio.yaml && \
    conda clean --all --yes

# Activate environment when container starts
SHELL ["conda", "run", "-n", "seq", "/bin/bash", "-c"]

# Set working directory
WORKDIR /BIP

# Ensure conda environment is activated in interactive shells
SHELL ["/bin/bash", "-c"]
RUN echo "conda activate seq" >> ~/.bashrc

# Working CMD to keep container open
CMD ["tail", "-f", "/dev/null"]

# Run snakemake and keep container alive for debugging
# CMD ["bash", "-c", "conda activate bio && snakemake --cores 1; tail -f /dev/null"]
