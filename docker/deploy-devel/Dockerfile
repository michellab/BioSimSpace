# This image is used to deploy the latest BioSimSpace binary.

FROM biosimspace/package-devel:latest

WORKDIR $HOME

USER $FN_USER

ARG par_url
ENV PAR_URL=$par_url

RUN $HOME/sire.app/bin/conda install -y pycurl

ADD deploy.py .
RUN $HOME/sire.app/bin/python deploy.py

ADD deploy_release.sh .
RUN ./deploy_release.sh

ENTRYPOINT ["bash"]
