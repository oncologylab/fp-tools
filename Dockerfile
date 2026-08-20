FROM mambaorg/micromamba:2.0.5

WORKDIR /opt/fp-tools
USER root
RUN chown "$MAMBA_USER:$MAMBA_USER" /opt/fp-tools
USER $MAMBA_USER
COPY --chown=$MAMBA_USER:$MAMBA_USER packaging/container/environment.yml /tmp/fp-tools-container.yml
RUN micromamba install -y -n base -f /tmp/fp-tools-container.yml && micromamba clean --all --yes

COPY --chown=$MAMBA_USER:$MAMBA_USER pyproject.toml setup.py README.md LICENSE ./
COPY --chown=$MAMBA_USER:$MAMBA_USER src ./src
RUN micromamba run -n base python -m pip install --no-cache-dir . \
    && micromamba run -n base python -m pip check

USER root
RUN mkdir -p /work && chown "$MAMBA_USER:$MAMBA_USER" /work
USER $MAMBA_USER
WORKDIR /work

EXPOSE 8891
CMD ["fp-tools-gui", "--host", "0.0.0.0", "--port", "8891", "--no-browser", "--run-dir", "/work/fp-tools-gui-runs"]
