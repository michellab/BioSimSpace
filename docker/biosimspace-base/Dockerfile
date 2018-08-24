# Based on the JupyterHub docker image that is
# Copyright (c) Jupyter Development Team.
# Distributed under the terms of the Modified BSD License.
# Ubuntu 16.04 (xenial) from 2017-07-23
# https://github.com/docker-library/official-images/commit/0ea9b38b835ffb656c497783321632ec7f87b60c

# This has been extensively modified by C.Woods, so blame him
# for any errors or problems

# Start from ubuntu :-)
FROM ubuntu:bionic-20180526@sha256:c8c275751219dadad8fa56b3ac41ca6cb22219ff117ca98fe82b42f24e1ba64e

LABEL maintainer="Christopher Woods <Christopher.Woods@bristol.ac.uk>"

USER root

# Install all OS dependencies for notebook server that starts but lacks all
# features (e.g., download as all possible file formats)
ENV DEBIAN_FRONTEND noninteractive
RUN apt-get update && apt-get -yq dist-upgrade \
 && apt-get install -yq --no-install-recommends \
    wget \
    bzip2 \
    ca-certificates \
    sudo \
    locales \
    git \
    bc \
    less \
    nano \
    libgfortran3 \
    fonts-liberation \
 && apt-get clean \
 && rm -rf /var/lib/apt/lists/*

RUN echo "en_US.UTF-8 UTF-8" > /etc/locale.gen && \
    locale-gen

# Configure environment
ENV CONDA_DIR=/opt/conda \
    SHELL=/bin/bash \
    SIRE_SILENT_PHONEHOME=1 \
    SIRE=/opt/conda \
    AMBERHOME=/home/amber16 \
    NB_USER=jovyan \
    NB_UID=1000 \
    NB_GID=100 \
    LC_ALL=en_US.UTF-8 \
    LANG=en_US.UTF-8 \
    LANGUAGE=en_US.UTF-8

ENV PATH=$CONDA_DIR/bin:$PATH \
    LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$AMBERHOME/lib64 \
    HOME=/home/$NB_USER

ADD fix-permissions /usr/local/bin/fix-permissions

# Create jovyan user with UID=1000 and in the 'users' group
# and make sure these dirs are writable by the `users` group.
RUN groupadd wheel -g 11 && \
    echo "auth required pam_wheel.so use_uid" >> /etc/pam.d/su && \
    useradd -m -s /bin/bash -N -u $NB_UID $NB_USER && \
    chmod g+w /etc/passwd && \
    fix-permissions $HOME

# Download and install amber16 into /home using wget
RUN wget https://biosimspace.blob.core.windows.net/dependencies/amber16.tar.bz2 -O amber16.tar.bz2 && \
    tar -jxvf amber16.tar.bz2 && \
    rm amber16.tar.bz2

# Download and install gromacs into /home using wget
RUN wget https://biosimspace.blob.core.windows.net/dependencies/gromacs.tar.bz2 -O gromacs.tar.bz2 && \
    tar -jxvf gromacs.tar.bz2 && \
    rm gromacs.tar.bz2

# Install the latest sire-devel package...
RUN wget https://sire.blob.core.windows.net/development/sire_devel_latest_linux.run -O sire.run && \
    chmod a+x sire.run && \
    echo "/opt/conda" | ./sire.run && \
    chown -R $NB_USER:$NB_GID /opt/conda && \
    fix-permissions /opt/conda && \
    rm sire.run

# Install biosimspace into /home using git
WORKDIR /home
RUN git clone https://github.com/michellab/biosimspace && \
    chown -R $NB_USER:$NB_GID /home/biosimspace && \
    fix-permissions /home/biosimspace

# Do all conda work as $NB_USER
USER $NB_USER
WORKDIR $HOME

RUN conda config --system --prepend channels conda-forge && \
    conda config --system --set auto_update_conda false && \
    conda config --system --set show_channel_urls true

# Install Tini
RUN conda install --quiet --yes 'tini=0.18.0' && \
    conda list tini | grep tini | tr -s ' ' | cut -d ' ' -f 1,2 >> $CONDA_DIR/conda-meta/pinned && \
    conda clean -tipsy && \
    fix-permissions $CONDA_DIR && \
    fix-permissions /home/$NB_USER

# Install Jupyter Notebook, Lab, and Hub
# Generate a notebook server config
# Cleanup temporary files
# Correct permissions
# Do all this in a single RUN command to avoid duplicating all of the
# files across image layers when the permissions change
RUN conda install --quiet --yes \
    'notebook=5.6.*' \
    'jupyterhub=0.9.*' \
    'jupyterlab=0.34.*' && \
    conda clean -tipsy && \
    jupyter labextension install @jupyterlab/hub-extension@^0.11.0 && \
    npm cache clean --force && \
    jupyter notebook --generate-config && \
    rm -rf $CONDA_DIR/share/jupyter/lab/staging && \
    rm -rf /home/$NB_USER/.cache/yarn && \
    fix-permissions $CONDA_DIR && \
    fix-permissions /home/$NB_USER

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

RUN conda config --add channels conda-forge && \
    conda install mdanalysis
RUN conda install -c omnia mdtraj
RUN conda install -c omnia pymbar

# Install BioSimSpace - we install at the end to prevent
# everything from being reinstalled every time we update
# BioSimSpace
RUN cd /home/biosimspace/python && \
    git checkout devel && \
    git pull && \
    python setup.py install

# Clean up after conda, including clearing the cache
RUN conda clean -tipsy && \
    npm cache clean --force && \
    rm -rf $CONDA_DIR/share/jupyter/lab/staging && \
    fix-permissions $CONDA_DIR && \
    rm -rf $HOME/.cache $HOME/.jupyter $HOME/.local/share/jupyter

# Link the BioSimSpace directory into $HOME so that we can access
# this from the notebook
RUN ln -s /home/biosimspace $HOME/

USER root

EXPOSE 8888
WORKDIR $HOME

# Configure container startup
ENTRYPOINT ["tini", "-g", "--"]
CMD ["start-notebook.sh"]

# Add local files as late as possible to avoid cache busting
COPY start.sh /usr/local/bin/
COPY start-notebook.sh /usr/local/bin/
COPY start-singleuser.sh /usr/local/bin/
COPY jupyter_notebook_config.py /etc/jupyter/
RUN fix-permissions /etc/jupyter/

# End as the User to make sure that we don't
# accidentally run the container as root
USER $NB_USER
