# This image is used to build the latest version of the devel
# branch of BioSimSpace.

FROM siremol/sire-devel:latest

WORKDIR $HOME

USER $FN_USER

# Disable Sire analytics.
ENV SIRE_DONT_PHONEHOME=1 \
    SIRE_SILENT_PHONEHOME=1

# Update to the latest version
RUN git clone https://github.com/michellab/BioSimSpace

RUN cd BioSimSpace && git checkout devel && git pull

RUN cd BioSimSpace/python && \
    $HOME/sire.app/bin/python setup.py install

# Clean up so that the Docker image size is minimised
RUN $HOME/sire.app/bin/conda clean -tipy

ENTRYPOINT ["bash"]
