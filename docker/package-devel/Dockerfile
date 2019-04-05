# This image is used to package the latest version of the devel
# branch of BioSimSpace.

FROM biosimspace/biosimspace-devel:latest

WORKDIR $HOME

# Add the modified Sire installation script.
ADD install_sire.sh .

# Move modified installation script into Sire application.
RUN mv -f install_sire.sh $HOME/sire.app/share/Sire/build/

# Set environment variables for installation directory and run file.
ENV SIRE_RUN_FILE $HOME/biosimspace_devel_latest_linux.run

# Create a self-extracting binary file.
RUN $HOME/sire.app/bin/package_sire > package.log 2> package.err
RUN chmod a+x $HOME/biosimspace_devel_latest_linux.run

ENTRYPOINT ["bash"]
