FROM mambaorg/micromamba:2.0.5

WORKDIR /opt/fp-tools
COPY --chown=$MAMBA_USER:$MAMBA_USER . .
RUN micromamba install -y -n base -f environment.yml && micromamba clean --all --yes

CMD ["bash"]
