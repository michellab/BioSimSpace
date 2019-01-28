# This image is used to build the latest version of the devel
# branch of BioSimSpace.

# Build on top of latest Sire development build.
FROM siremol/sire-devel:latest

WORKDIR $HOME

USER $FN_USER

# Clone BioSimSpace and update to the latest devel.
RUN git clone https://github.com/michellab/BioSimSpace
RUN cd BioSimSpace && git checkout devel && git pull

# Build BioSimSpace into existing Sire application.
ADD build_biosimspace.sh $HOME/build_biosimspace.sh
RUN scl enable devtoolset-6 ./build_biosimspace.sh

ENTRYPOINT ["scl", "enable", "devtoolset-6", "bash"]
