# This image is used to build the documentation and deploy the website for the
# latest version of the devel branch of BioSimSpace.

FROM biosimspace/biosimspace-devel:latest

WORKDIR $HOME

ARG github_token
ARG github_email
ENV GITHUB_TOKEN=$github_token
ENV GITHUB_EMAIL=$github_email

# Install additional Python dependencies.
RUN $HOME/sire.app/bin/pip install docutils==0.17.1 sphinx==2.2.2 sphinx_issues sphinx_rtd_theme sphinxcontrib-youtube==0.1.2

# Patch the sphinxcontrib.youtube package.
RUN cp $HOME/BioSimSpace/docker/document-devel/youtube.py $HOME/sire.app/lib/python3.7/site-packages/sphinxcontrib/youtube

# Inject the dynamically create force field function names into the Sphinx
# documentation in the source repository.
ADD inject_force_fields.sh .
RUN $HOME/inject_force_fields.sh

# Build the html docs against the patched source code.
RUN cd $HOME/BioSimSpace/doc && PYTHONPATH=$HOME/BioSimSpace/python SPHINXBUILD=$HOME/sire.app/bin/sphinx-build make html

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
