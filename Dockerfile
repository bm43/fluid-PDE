FROM dolfinx/lab:stable

# elevate, add any extra packages
USER root
RUN python3 -m pip install --no-cache-dir gmsh meshio

# copy project and ensure jovyan owns it
ARG NB_USER=jovyan
ARG NB_UID=1000
ENV HOME=/home/${NB_USER}
WORKDIR ${HOME}
COPY --chown=${NB_UID}:${NB_UID} . ${HOME}
USER ${NB_USER}

# Jupyter already provided by base image; Binder will start it
