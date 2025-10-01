FROM dolfinx/lab:stable

# Binder requirement: ensure Jupyter is present and recent
RUN python3 -m pip install --no-cache-dir notebook jupyterlab

# Optional: geometry/IO tools for the scripts
RUN python3 -m pip install --no-cache-dir gmsh meshio

# Binder requirement: run as a non-root user and make files writable
ARG NB_USER=jovyan
ARG NB_UID=1000
ENV USER=${NB_USER} NB_UID=${NB_UID} HOME=/home/${NB_USER}
USER root
RUN adduser --disabled-password --gecos "Default user" --uid ${NB_UID} ${NB_USER}

# Copy repo into home and set ownership
COPY . ${HOME}
RUN chown -R ${NB_UID} ${HOME}
USER ${NB_USER}
WORKDIR ${HOME}

# Start JupyterLab when Binder launches
CMD ["jupyter", "lab", "--ip=0.0.0.0", "--port=8888", "--no-browser", "--NotebookApp.default_url=/lab"]
