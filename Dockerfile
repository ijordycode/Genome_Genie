FROM continuumio/miniconda3

COPY env.yaml /BIP/env.yaml
COPY Snakefile /BIP/Snakefile
COPY config.yaml /BIP/config.yaml

# Create and clean the conda environment
RUN conda env create --name bio -f /BIP/env.yaml && \
    conda clean --all --yes

# Activate environment when container starts
SHELL ["conda", "run", "-n", "bio", "/bin/bash", "-c"]

# Set working directory
WORKDIR /BIP

# Ensure conda environment is activated in interactive shells
SHELL ["/bin/bash", "-c"]
RUN echo "conda activate bio" >> ~/.bashrc

# Working CMD to keep container open
CMD ["tail", "-f", "/dev/null"]

# Run snakemake and keep container alive for debugging
# CMD ["bash", "-c", "conda activate bio && snakemake --cores 1; tail -f /dev/null"]
