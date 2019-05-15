This directory contains the Dockerfile needed to build the
BioSimSpace notebook image.

It also contains the yaml file used to configure a jupyterhub
instance to use this image.

Build the docker file
---------------------

.. code-block:: bash

    docker build . -t biosimspace/biosimspace-notebook:v5

(We tag notebook images with a version so that we can revert to a working state.)

Test the docker file locally
----------------------------

(Note that you can't test locally if the ``jupyter_notebook_config.py``
file is included - you need to comment out.)

.. code-block:: bash

    docker run -p 8888:8888 -it biosimspace/biosimspace-notebook:latest

(and then go to the link specified)

Push the docker file to dockerhub
---------------------------------

.. code-block:: bash

    docker push biosimspace/biosimspace-notebook:v5

Installing into k8s
-------------------

.. code-block:: bash

    helm install jupyterhub/jupyterhub --version=0.8.0 --timeout=36000 --debug --install=notebook --namespace=notebook --values notebook.yaml

Upgrading jupyterhub
--------------------

.. code-block:: bash

    helm upgrade --install notebook --namespace notebook jupyterhub/jupyterhub --timeout=36000 --debug --version=0.8.0 --values notebook.yaml

Deleting user pods
------------------

Assuming your kubectl works, then you can delete all active user notebook pods
using

.. code-block:: bash

    python killjupyter.py

This will kill all logged in pods as well - so only use it if you need to
reset everything because the server is overloaded.
