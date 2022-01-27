#!/usr/bin/env bash

# Get the website repository user token.
WEBSITE_TOKEN=${WEBSITE_TOKEN:-$1}

git push --repo https://biosimspacebot:$WEBSITE_TOKEN@github.com/michellab/BioSimSpaceWebsite.git > /dev/null 2>&1
