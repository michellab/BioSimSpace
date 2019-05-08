This directory contains the Dockerfile needed to build the 
BioSimSpace notebook image.

It also contains the yaml file used to configure a jupyterhub
instance to use this image.

# Build the docker file

```
$ docker build . -t chryswoods/biosimspace-notebook:v5
```

# Test the docker file locally

```
$ docker run -p 8888:8888 -it chryswoods/biosimspace-notebook:latest
```

(and then go to the link specified)

# Push the docker file to dockerhub

```
$ docker push chryswoods/biosimspace-notebook:v5
```

# Installing into k8s

```
$ helm install jupyterhub/jupyterhub --version=0.8.0 --timeout=36000 --debug --install=notebook --namespace=notebook --values notebook.yaml

# Upgrading jupyterhub

```
$ helm upgrade --install=notebook --namespace=notebook jupyterhub/jupyterhub --timeout=36000 --debug --version=0.8.0 --values notebook.yaml
```
