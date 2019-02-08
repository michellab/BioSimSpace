
# This image is used to build, package and test the latest version
# of the devel branch of BioSimSpace

FROM siremol/sire-devel:latest

WORKDIR $HOME

USER $FN_USER

# Update to the latest version
RUN git clone https://github.com/michellab/BioSimSpace

WORKDIR $HOME/BioSimSpace

RUN git checkout devel && git pull

RUN cd BioSimSpace/python && \
    $HOME/sire.app/bin/python setup.py install
 
# Clean up so that the Docker image size is minimised
RUN $HOME/sire.app/bin/conda clean -tipy

ENTRYPOINT ["bash"]
