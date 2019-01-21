# This image is used to package the latest version of the devel
# branch of BioSimSpace.

FROM biosimspace/biosimspace-devel:latest

WORKDIR $HOME

# Add the modified Sire installation script.
ADD install_sire.sh .

# Move modified installation script into Sire application.
RUN mv -f install_sire.sh $HOME/sire.app/pkgs/sire-*/share/Sire/build/

# Set environment variables for installation directory and run file.
ENV SIRE_RUN_FILE biosimspace_devel_latest_linux.run

# Create a self-extracting binary file.
RUN $HOME/sire.app/bin/package_sire
RUN chmod a+x biosimspace_devel_latest_linux.run

ENTRYPOINT ["scl", "enable", "devtoolset-6", "bash"]
