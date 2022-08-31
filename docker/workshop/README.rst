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

.. code-block:: bash

    helm install jupyterhub/jupyterhub --version=0.8.0 --timeout=36000 --debug --install=workshop --namespace=workshop --values workshop.yaml

Upgrading jupyterhub
--------------------

.. code-block:: bash

    helm upgrade --install workshop --namespace workshop jupyterhub/jupyterhub --timeout=36000 --debug --version=0.8.0 --values workshop.yaml

Deleting user pods
------------------

Assuming your kubectl works, then you can delete all active user workshop pods
using

.. code-block:: bash

    python killjupyter.py

This will kill all logged in pods as well - so only use it if you need to
reset everything because the server is overloaded.
