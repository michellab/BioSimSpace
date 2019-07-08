# Start from the base jupyterhub image.
ARG BASE_CONTAINER=jupyter/base-notebook
FROM $BASE_CONTAINER

LABEL maintainer="Christopher Woods <Christopher.Woods@bristol.ac.uk>, Lester Hedges <lester.hedges@gmail.com>"

# Variables holding links to download dependencies.
ENV AMBER_DOWNLOAD=https://objectstorage.eu-frankfurt-1.oraclecloud.com/p/RHIfl-cHK9QqvzCsqFRMGpijRoHRPC-ZiATsjWzUgD0/n/chryswoods/b/downloads/o/amber.tar.bz2
ENV GROMACS_DOWNLOAD=https://objectstorage.eu-frankfurt-1.oraclecloud.com/p/L_wxiMnhIDBO-CeO0M6toS8egsIvZXWQ7AGJiVEpYk4/n/chryswoods/b/downloads/o/gromacs.tar.bz2
ENV PLUMED_DOWNLOAD=https://objectstorage.eu-frankfurt-1.oraclecloud.com/p/em4QTi3Q3cH3lDhrksODgX7oZ9CCn-k1J44GzI98G4w/n/chryswoods/b/downloads/o/plumed.tar.bz2

USER root
WORKDIR /home

# Pull in the required packages.
RUN apt-get update && \
    apt-get -yq --no-install-recommends install git libgfortran3 nano vim && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Configure environment.
ENV SIRE_SILENT_PHONEHOME=1 \
    SIRE=/opt/conda \
    AMBERHOME=/home/amber \
    PATH=/home/gromacs/bin:/home/plumed/bin:$PATH \
    LD_LIBRARY_PATH=/home/plumed/lib:/opt/conda/lib:$LD_LIBRARY_PATH \
    PLUMED_KERNEL=/home/plumed/lib/libplumedKernel.so

# Download and install AmberTools into /home using wget.
RUN wget $AMBER_DOWNLOAD -O amber.tar.bz2 && \
    tar -jxvf amber.tar.bz2 && \
    rm amber.tar.bz2

# Download and install GROMACS into /home using wget.
RUN wget $GROMACS_DOWNLOAD -O gromacs.tar.bz2 && \
    tar -jxvf gromacs.tar.bz2 && \
    rm gromacs.tar.bz2

# Download and install PLUMED into /home using wget.
RUN wget $PLUMED_DOWNLOAD -O plumed.tar.bz2 && \
    tar -jxvf plumed.tar.bz2 && \
    rm plumed.tar.bz2

# Pin the version of Conda so that we don't get unexpected
# changes in Python.
RUN echo "python 3.7.*" >> /opt/conda/conda-meta/pinned

# Do all conda work as $NB_USER.
USER $NB_USER

# Install BioSimSpace via Conda and dependencies for PLUMED.
RUN conda install -y -c rdkit -c conda-forge -c omnia -c michellab/label/dev biosimspace && \
    conda install -y -c conda-forge openmpi-mpicc openssh

RUN pip install --upgrade pip && \
    pip install jupyterhub-tmpauthenticator

# Install and enable nglview and fileupload.
RUN jupyter-nbextension install nglview --py --sys-prefix && \
    jupyter-nbextension enable nglview --py --sys-prefix && \
    pip install fileupload && \
    jupyter-nbextension install fileupload --py --sys-prefix && \
    jupyter-nbextension enable fileupload --py --sys-prefix

# Clean up after conda, including clearing the cache.
RUN conda clean -tipsy && \
    npm cache clean --force && \
    rm -rf $CONDA_DIR/share/jupyter/lab/staging && \
    rm -rf $HOME/.cache $HOME/.jupyter $HOME/.local/share/jupyter

# Clone source and workshop repositories.
RUN git clone https://github.com/michellab/biosimspace && \
    git clone https://github.com/CCPBioSim/python_and_data_workshop && \
    git clone https://github.com/CCPBioSim/biosimspace_workshop

# Link the right directories into the right places.
RUN rmdir work && \
    ln -s $HOME/biosimspace/demo $HOME/demo

# Copy the example nodes into the Python library.
RUN mkdir /opt/conda/lib/python3.7/site-packages/BioSimSpace/Node/_nodes && \
    cp $HOME/demo/*.py /opt/conda/lib/python3.7/site-packages/BioSimSpace/Node/_nodes/

# Add in our custom notebook config.
USER root

# Comment out the below line if you need to test locally.
# YOU MUST UNCOMMENT THIS LINE FOR REMOTE DEPLOYMENT!
COPY jupyter_notebook_config.py /etc/jupyter/

# Add in the 'update_biosimspace' command to make things easy.
COPY update_biosimspace /usr/bin

# End as the User to make sure that we don't
# accidentally run the container as root.
USER $NB_USER
