# Start from the base jupyterhub image.
ARG BASE_CONTAINER=jupyter/base-notebook:python-3.9.12
FROM $BASE_CONTAINER

LABEL maintainer="Christopher Woods <Christopher.Woods@bristol.ac.uk>, Lester Hedges <lester.hedges@gmail.com>"

# Variables holding links to download dependencies
ENV 
GROMACS_DOWNLOAD=https://objectstorage.eu-frankfurt-1.oraclecloud.com/p/kjTvvBQ__y2d_svdfejffT6NnViSipXPA-mwY0Ve57aHk3OHDpXwQQz0uuBZC_fa/n/hugs/b/notebook/o/gromacs.tar.bz2
ENV 
PLUMED_DOWNLOAD=https://objectstorage.eu-frankfurt-1.oraclecloud.com/p/hjUMKuIo2FNNiW3gcGtvcLyN6IfZEmMMMO7RP4VvR2qfClOgeGHxV8Di1Cv6s7zl/n/hugs/b/notebook/o/plumed.tar.bz2

USER root
WORKDIR /home

# Fix symbolic link for libtinfo to avoid annoying Bash warning.
RUN rm /opt/conda/lib/libtinfo.so && \
    ln -sf /opt/conda/lib/libtinfo.so /opt/conda/lib/libtinfo.so.6

# Pull in the required packages.
RUN apt-get update && \
    apt-get -yq --no-install-recommends install git nano vim && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Configure environment.
ENV SIRE_SILENT_PHONEHOME=1 \
    SIRE=/opt/conda \
    AMBERHOME=/opt/conda \
    PATH=/home/amber/bin:/home/gromacs/bin:/home/plumed/bin:$PATH \
    LD_LIBRARY_PATH=/home/plumed/lib:/opt/conda/lib:$LD_LIBRARY_PATH \
    PLUMED_KERNEL=/home/plumed/lib/libplumedKernel.so

# Download and install GROMACS into /home using wget.
RUN wget $GROMACS_DOWNLOAD -O gromacs.tar.bz2 && \
    tar -jxvf gromacs.tar.bz2 && \
    rm gromacs.tar.bz2

# Download and install PLUMED into /home using wget.
RUN wget $PLUMED_DOWNLOAD -O plumed.tar.bz2 && \
    tar -jxvf plumed.tar.bz2 && \
    rm plumed.tar.bz2

# Do all conda work as $NB_USER.
USER $NB_USER
WORKDIR $HOME

# Install BioSimSpace via Conda and dependencies for PLUMED.
RUN conda install -y mamba

RUN mamba install -y -c conda-forge -c michellab/label/dev "biosimspace<2023.0" && \
    mamba install -y -c conda-forge ambertools libgfortran=3 nodejs openmpi-mpicxx openssh patch

RUN pip install --upgrade pip && \
    pip install jupyterhub-tmpauthenticator

# Install and enable nglview.
RUN jupyter-nbextension install nglview --py --sys-prefix && \
    jupyter-nbextension enable nglview --py --sys-prefix && \
    jupyter serverextension enable jupyterlab

# Clean up after conda, including clearing the cache.
RUN conda remove cudatoolkit --force && \
    conda clean -tipsy && \
    npm cache clean --force && \
    rm -rf $CONDA_DIR/share/jupyter/lab/staging && \
    rm -rf $HOME/.cache $HOME/.jupyter $HOME/.local/share/jupyter

# Clone source and workshop repositories.
RUN git clone https://github.com/michellab/biosimspace

# Link the right directories into the right places.
RUN rmdir work && \
    ln -s $HOME/biosimspace/demo $HOME/demo

# Add any BioSimSpace patches.
ADD patches .patches

# Patch the BioSimSpace install so that we restrict GROMACS to use a single
# MPI thread and four OpenMP thread per MPI thread.
RUN patch /opt/conda/lib/python3.9/site-packages/BioSimSpace/Process/_gromacs.py .patches/_gromacs.py.patch
RUN patch /opt/conda/lib/python3.9/site-packages/BioSimSpace/Metadynamics/_metadynamics.py .patches/_metadynamics.py.patch

# Copy the example nodes into the Python library.
RUN mkdir /opt/conda/lib/python3.9/site-packages/BioSimSpace/Node/_nodes && \
    cp $HOME/demo/*.py /opt/conda/lib/python3.9/site-packages/BioSimSpace/Node/_nodes/

# Change the executable name so that BioSimSpace workshop usage statistics
# can be tracked separately.
RUN sed -i 's/= "BioSimSpace"/= "BioSimSpace (workshop)"/g' /opt/conda/lib/python3.9/site-packages/Sire/__init__.py

# Add in our custom workshop config.
USER root

# Copy acrross JupyterHub configuration.
COPY jupyterhub_config.py /etc/jupyter/

# Add in the 'update_biosimspace' command to make things easy.
COPY update_biosimspace /usr/bin

# Run update_biosimspace in case there have been any updates that haven't been
# pushed to Conda.
RUN /usr/bin/update_biosimspace

# Allow anyone to run `sudo update_biosimspace`
RUN echo "$NB_USER ALL=(ALL) NOPASSWD:ALL" >> /etc/sudoers

COPY install_workshop /usr/bin
COPY update_workshop /usr/bin

# Install the workshop as it is very large, and too large
# for every attendee to download and install their own copy...
USER $NB_USER
RUN /usr/bin/install_workshop

USER root

# Install extra packages needed by the workshop
# (doing it here so that we don't have to relayer the image
#  - especially the above workshop data!)
RUN mamba install -y alchemlyb

# We need to downgrade to ipywidgets < 8
RUN mamba install -y "ipywidgets<8"

# We need to install the Exscientia version of MDRestraintsGenerator
RUN python -m pip install git+https://github.com/Exscientia/MDRestraintsGenerator.git

# We need to install the feature_workshop branch of Sire
COPY sire-2022.3.0-py39h3fd9d12_73.tar.bz2 /home/jovyan/sire-2022.3.0-py39h3fd9d12_73.tar.bz2

RUN mamba install -y /home/jovyan/sire-2022.3.0-py39h3fd9d12_73.tar.bz2 && \
    rm /home/jovyan/sire-2022.3.0-py39h3fd9d12_73.tar.bz2

# This is the only change since feature_workshop was last compiled and installed
COPY OpenMMMD.py /opt/conda/lib/python3.9/site-packages/Sire/Tools/OpenMMMD.py 

RUN /usr/bin/update_biosimspace

USER $NB_USER

RUN /usr/bin/update_workshop

USER root
RUN mamba install -y -c conda-forge tqdm
USER $NB_USER

