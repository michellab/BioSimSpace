# Start from the base jupyterhub image
ARG BASE_CONTAINER=jupyter/base-notebook
FROM $BASE_CONTAINER

LABEL maintainer="Christopher Woods <Christopher.Woods@bristol.ac.uk>"

# Variables holding links to download dependencies
ENV AMBER_DOWNLOAD=https://objectstorage.eu-frankfurt-1.oraclecloud.com/p/7ODpDhtg0WbT9Dl8JHW_4pJDB1Eka_EfCL3L600gIXo/n/chryswoods/b/downloads/o/amber16.tar.bz2
ENV GROMACS_DOWNLOAD=https://objectstorage.eu-frankfurt-1.oraclecloud.com/p/2jeHceFVZcIMeyRm7T3TpaUY8JajKmE3QkZ1iQ7vE3I/n/chryswoods/b/downloads/o/gromacs.tar.bz2

USER root
WORKDIR /home

# Configure environment
ENV SIRE_SILENT_PHONEHOME=1 \
    SIRE=/opt/conda \
    AMBERHOME=/home/amber16

# Download and install amber16 into /home using wget
RUN wget $AMBER_DOWNLOAD -O amber16.tar.bz2 && \
    tar -jxvf amber16.tar.bz2 && \
    rm amber16.tar.bz2

# Download and install gromacs into /home using wget
RUN wget $GROMACS_DOWNLOAD -O gromacs.tar.bz2 && \
    tar -jxvf gromacs.tar.bz2 && \
    rm gromacs.tar.bz2

# pin the version of conda so that we don't get unexpected
# changes in python
RUN echo "python 3.7.*" >> /opt/conda/conda-meta/pinned

# Do all conda work as $NB_USER
USER $NB_USER
WORKDIR $HOME

# Install Sire via conda
RUN conda install -c conda-forge -c omnia -c michellab/label/dev sire

# Install BioSimSpace via conda
RUN conda install -c conda-forge -c omnia -c michellab/label/dev biosimspace

RUN conda install watchdog
RUN pip install --upgrade pip && \
    pip install pygtail && \
    pip install matplotlib && \
    pip install jupyterhub-tmpauthenticator

# Install and enable nglview (including patching some issues) and fileupload
COPY nglview_patch.sh $HOME
RUN pip --no-cache-dir install nglview && \
    jupyter-nbextension install nglview --py --sys-prefix && \
    jupyter-nbextension enable nglview --py --sys-prefix && \
    $HOME/nglview_patch.sh /opt/conda && \
    rm $HOME/nglview_patch.sh && \
    pip install fileupload && \
    jupyter-nbextension install fileupload --py --sys-prefix && \
    jupyter-nbextension enable fileupload --py --sys-prefix

# clean up after conda, including clearing the cache
RUN conda clean -tipsy && \
    npm cache clean --force && \
    rm -rf $CONDA_DIR/share/jupyter/lab/staging && \
    rm -rf $HOME/.cache $HOME/.jupyter $HOME/.local/share/jupyter

# add our config file
# COPY jupyter_notebook_config.py /etc/jupyter/

# End as the User to make sure that we don't
# accidentally run the container as root
USER $NB_USER
