This directory contains the Dockerfile needed to build the
BioSimSpace workshop image.

It also contains the yaml file used to configure a jupyterhub
instance to use this image.

Build the docker file
---------------------

.. code-block:: bash

    docker build . -t biosimspace/biosimspace-workshop:v1

(We tag workshop images with a version so that we can revert to a working state.)

Test the docker file locally
----------------------------

(Note that you can't test locally if the ``jupyter_notebook_config.py``
file is included - you need to comment out.)

.. code-block:: bash

    docker run -p 8888:8888 -it biosimspace/biosimspace-workshop:latest

(and then go to the link specified)

Push the docker file to dockerhub
---------------------------------

.. code-block:: bash

    docker push biosimspace/biosimspace-workshop:v1

Installing into k8s
-------------------

Install the jupyterhub helm chart

.. code-block:: bash

    helm repo add jupyterhub https://jupyterhub.github.io/helm-chart/
    helm repo update

Remember to update `workshop.yaml` to set the TAG for the workshop image you want,
plus to set a SECRET_KEY. You can generate the secret key using the command

.. code-block:: bash

    openssl rand -hex 32

Once you are ready, you can install jupyterhub.

.. code-block:: bash

    helm install workshop jupyterhub/jupyterhub --version=1.0.0 --timeout=36000s --debug --values workshop.yaml

Once it has installed, you need to get the IP address of the server, e.g. via

.. code-block:: bash

    kubectl get services

Then, assign this IP address to the domain name that you've set in `workshop.yaml`.

Note that it may take some time for letsencrypt to create a SSL certificate for you. Be patient.

Upgrading jupyterhub
--------------------

.. code-block:: bash

    helm upgrade workshop jupyterhub/jupyterhub --timeout=36000 --debug --version=1.0.0 --values workshop.yaml

Deleting user pods
------------------

Assuming your kubectl works, then you can delete all active user workshop pods
using

.. code-block:: bash

    python killjupyter.py

This will kill all logged in pods as well - so only use it if you need to
reset everything because the server is overloaded.
