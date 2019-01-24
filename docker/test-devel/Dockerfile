# This image is used to test the latest version of the devel
# branch of BioSimSpace.

FROM biosimspace/package-devel:latest

WORKDIR $HOME

# Set installation directory.
ENV BIOSIMSPACE_INSTALL_DIR biosimspace.app

# Extract binary.
RUN ./biosimspace_devel_latest_linux.run

# Move to BioSimSpace directory, pull the latest updates and run tests.
RUN cd BioSimSpace && git checkout devel && git pull && \
    $HOME/biosimspace.app/bin/pytest -v test

ENTRYPOINT ["scl", "enable", "devtoolset-6", "bash"]
