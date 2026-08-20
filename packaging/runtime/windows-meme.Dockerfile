FROM mambaorg/micromamba:2.0.5

LABEL org.opencontainers.image.source="https://github.com/oncologylab/fp-tools"

WORKDIR /opt/fp-tools
ENV FP_TOOLS_RUNTIME=system
USER root
RUN chown "$MAMBA_USER:$MAMBA_USER" /opt/fp-tools
USER $MAMBA_USER

COPY --chown=$MAMBA_USER:$MAMBA_USER packaging/runtime/meme.yml /tmp/fp-tools-runtime-meme.yml
RUN micromamba install -y -n base -f /tmp/fp-tools-runtime-meme.yml \
    && micromamba clean --all --yes

COPY --chown=$MAMBA_USER:$MAMBA_USER pyproject.toml setup.py README.md LICENSE ./
COPY --chown=$MAMBA_USER:$MAMBA_USER src ./src
USER root
RUN apt-get update \
    && apt-get install -y --no-install-recommends gcc libc6-dev \
    && micromamba run -n base python -m pip install --no-cache-dir . \
    && micromamba run -n base python -m pip check \
    && apt-get purge -y --auto-remove gcc libc6-dev \
    && rm -rf /var/lib/apt/lists/* /root/.cache/pip
USER $MAMBA_USER

WORKDIR /work
