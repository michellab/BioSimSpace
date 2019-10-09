# This image is used to build the documentation and deploy the website for the
# latest version of the devel branch of BioSimSpace.

FROM biosimspace/biosimspace-devel:latest

WORKDIR $HOME

ARG github_token
ARG github_email
ENV GITHUB_TOKEN=$github_token
ENV GITHUB_EMAIL=$github_email

# Install additional Python dependencies.
RUN $HOME/sire.app/bin/pip install sphinx sphinx_issues sphinx_rtd_theme sphinxcontrib.youtube

# Patch the sphinxcontrib.youtube package.
RUN cp $HOME/BioSimSpace/docker/document-devel/youtube.py $HOME/sire.app/lib/python3.7/site-packages/sphinxcontrib/youtube

# Build the html docs.
RUN cd $HOME/BioSimSpace/doc && SPHINXBUILD=$HOME/sire.app/bin/sphinx-build make html

# Clone the BioSimSpaceWebsite repository, copy across latest docs, commit and push.
RUN git clone https://github.com/michellab/BioSimSpaceWebsite.git && \
    cd $HOME/BioSimSpaceWebsite && \
    cp -a $HOME/BioSimSpace/doc/build/html/* docs && \
    git config user.name "BioSimSpaceBot" && \
    git config user.email "$GITHUB_EMAIL" && \
    git add docs && \
    git commit -a -m "Updating website." && \
    git push --repo https://biosimspacebot:$GITHUB_TOKEN@github.com/michellab/BioSimSpaceWebsite.git > /dev/null 2>&1

ENTRYPOINT ["bash"]
