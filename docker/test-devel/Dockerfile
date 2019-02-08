# This image is used to test the latest version of the devel
# branch of BioSimSpace.

FROM biosimspace/biosimspace-devel:latest

WORKDIR $HOME

# Move to BioSimSpace directory, pull the latest updates and run tests.
RUN cd BioSimSpace && \
    $HOME/sire.app/bin/pytest -v test

ENTRYPOINT ["bash"]
