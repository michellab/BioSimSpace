# Creating a new JupyterHub from scratch

This is often necessary when kubernetes, helm or JupyterHub are upgraded.

In these cases, it is better to start again and build everything from scratch.

Here is the process I went through the last time I did this.

(note that I did all of this in an Oracle Cloud Shell, having already created
a k8s cluster and followed Oracle's instructions to activate kubectl
and connect it to the cluster)

(note also that I followed the instructions 
`from here <https://zero-to-jupyterhub.readthedocs.io/en/latest/jupyterhub/installation.html>`__)

## Default JupyterHub

Start from a default JupyterHub. Use the default `workshop.yaml`

::

   # This file can update the JupyterHub Helm chart's default configuration values.
   #
   # For reference see the configuration reference and default values, but make
   # sure to refer to the Helm chart version of interest to you!
   #
   # Introduction to YAML:     https://www.youtube.com/watch?v=cdLNKUoMc6c
   # Chart config reference:   https://zero-to-jupyterhub.readthedocs.io/en/stable/resources/reference.html
   # Chart default values:     https://github.com/jupyterhub/zero-to-jupyterhub-k8s/blob/HEAD/jupyterhub/values.yaml
   # Available chart versions: https://jupyterhub.github.io/helm-chart/
   #

Use the links to find the current chart versions. For today, the latest JupyterHub
version is 1.2.0

## Install JupyterHub Repo

.. code-block:: bash

    helm repo add jupyterhub https://jupyterhub.github.io/helm-chart/
    helm repo update

## Install JupyterHub

.. code-block:: bash

    helm upgrade --cleanup-on-fail --install workshop jupyterhub/jupyterhub --namespace workshop --create-namespace --version=1.2.0 --values workshop.yaml

## Access the hub

Get the IP address of the service using

.. code-block:: bash

    kubectl --namespace=workshop get services

e.g. I see

::

    NAME           TYPE           CLUSTER-IP      EXTERNAL-IP   PORT(S)        AGE
    hub            ClusterIP      10.96.244.150   <none>        8081/TCP       25s
    proxy-api      ClusterIP      10.96.110.222   <none>        8001/TCP       25s
    proxy-public   LoadBalancer   10.96.212.29    <pending>     80:32237/TCP   25s

You need to wait until the `EXTERNAL-IP` of `proxy-public` goes from `<pending>` 
to an actual IP address.

Once you have the IP address, put it into a browser to navigate to the site.

If everything worked, then you are ready to go onto the next stage.

## Building the docker image

The next step is to build the docker image.

I am starting from the `Dockerfile` we have used for other workshops.

To build, run

.. code-block:: bash

    docker build . -t biosimspace/biosimspace-workshop:v1

## Updating the config

.. code-block:: bash

    helm upgrade --debug workshop jupyterhub/jupyterhub --namespace workshop --version=1.2.0 --values workshop.yaml
