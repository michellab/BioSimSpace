
# This image is used to build, package and test the latest version
# of the devel branch of BioSimSpace

FROM siremol/sire-devel:latest

WORKDIR $HOME

USER $FN_USER

# Clone BioSimSpace and update to the latest devel.
RUN git clone https://github.com/michellab/BioSimSpace
RUN cd BioSimSpace && git checkout devel && git pull
 
# Build Sire (this should take less time as it is a small update)
RUN scl enable devtoolset-6 ./build_biosimspace.sh

ENTRYPOINT ["scl", "enable", "devtoolset-6", "bash"]
