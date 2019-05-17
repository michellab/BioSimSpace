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

# Install BioSimSpace via conda
RUN conda install -c rdkit -c conda-forge -c omnia -c michellab/label/dev biosimspace

RUN pip install --upgrade pip && \
    pip install pygtail && \
    pip install pypdb && \
    pip install jupyterhub-tmpauthenticator

# Install and enable nglview and fileupload
RUN jupyter-nbextension install nglview --py --sys-prefix && \
    jupyter-nbextension enable nglview --py --sys-prefix && \
    pip install fileupload && \
    jupyter-nbextension install fileupload --py --sys-prefix && \
    jupyter-nbextension enable fileupload --py --sys-prefix

# clean up after conda, including clearing the cache
RUN conda clean -tipsy && \
    npm cache clean --force && \
    rm -rf $CONDA_DIR/share/jupyter/lab/staging && \
    rm -rf $HOME/.cache $HOME/.jupyter $HOME/.local/share/jupyter

# Now pull in git and download the workshop
USER root
RUN apt-get update && \
    apt-get -yq --no-install-recommends install git libgfortran3 nano && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

USER $NB_USER
RUN git clone https://github.com/michellab/biosimspace && \
    git clone https://github.com/CCPBioSim/python_and_data_workshop

RUN git clone https://github.com/CCPBioSim/biosimspace_workshop

# Link the right directories into the right places...
RUN rmdir work && \
    ln -s $HOME/biosimspace/demo $HOME/demo

# Copy the example nodes into the Python library.
RUN mkdir /opt/conda/lib/python3.7/site-packages/BioSimSpace/Node/_nodes && \
    cp $HOME/demo/*.py /opt/conda/lib/python3.7/site-packages/BioSimSpace/Node/_nodes/

# Ensure that Gromacs is in our path, and /opt/conda is
# in the library path for gmx dependency of libgomp
ENV PATH=$PATH:/home/gromacs/bin \
    LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/opt/conda/lib

# Add in our custom notebook config
USER root

# Comment out the below line if you need to test locally. YOU MUST
# UNCOMMENT THIS LINE FOR REMOTE DEPLOYMENT
COPY jupyter_notebook_config.py /etc/jupyter/

# Add in the 'update_biosimspace' command to make things easy
COPY update_biosimspace /usr/bin

# End as the User to make sure that we don't
# accidentally run the container as root
USER $NB_USER
