# Workshop image that builds on the biosimspace docker image
# This starts from biosimspace-base and adds on any late
# fixes or other things. This prevents us from having to
# rebuild everything built on top of biosimspace-base
# if we just need to update biosimspace or make small
# changes or fixes

# start from the BioSimSpace image
FROM chryswoods/biosimspace-base:latest

USER $NB_USER
WORKDIR /home/$NB_USER

# Add any missing environment variables
ENV OPENMM_PLUGIN_DIR=$SIRE/lib/plugins

# Add a file into the path that allows the user to update
# BioSimSpace in place easily
COPY update_biosimspace /opt/conda/bin/

# Finally, run 'update_workshops' so that everything is 
# fully updated (or we want to update without rebuilding
# everything)
RUN /opt/conda/bin/update_biosimspace

# Ensure that the container runs as $NB_USER
USER $NB_USER
