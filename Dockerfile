# Pipeline image: all tools from environment.yml plus the pipeline code.
# Reference data and the genotype DB are NOT in the image — mount them at /shared.
#
#   docker build --build-arg PIPELINE_VERSION=$(git describe --always --dirty) -t scrna-pipeline .
#
#   docker run --rm --user $(id -u):$(id -g) \
#       -v /groups/nils/data/scrna_pipeline:/shared \
#       -v "$PWD":/work \
#       scrna-pipeline --mode bulk --sample SRR5071686 --run-id SRR5071686_hg38
FROM mambaorg/micromamba:1.5.10

COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml /tmp/environment.yml
RUN micromamba install -y -n base -f /tmp/environment.yml && \
    micromamba clean --all --yes

COPY --chown=$MAMBA_USER:$MAMBA_USER run_pipeline.py /opt/scrna_pipeline/
COPY --chown=$MAMBA_USER:$MAMBA_USER workflow /opt/scrna_pipeline/workflow
COPY --chown=$MAMBA_USER:$MAMBA_USER config   /opt/scrna_pipeline/config
COPY --chown=$MAMBA_USER:$MAMBA_USER scripts  /opt/scrna_pipeline/scripts

ARG PIPELINE_VERSION=unknown
# HOME=/tmp: the container runs as the host user, who has no home directory inside it.
ENV SCRNA_PIPELINE_VERSION=$PIPELINE_VERSION \
    SCRNA_SHARED_ROOT=/shared \
    HOME=/tmp

WORKDIR /work
# --foreground: a detached run would be killed when the container's main process exits.
ENTRYPOINT ["/usr/local/bin/_entrypoint.sh", "python", "/opt/scrna_pipeline/run_pipeline.py", "--foreground"]
