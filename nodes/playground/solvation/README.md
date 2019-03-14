Running BioSimSpace:
--------------------

All the documentation can be found [here](www.biosimspace.org).

1. Get BioSimSpace   
   Easiest is to download the development binaries from [here](https://objectstorage.eu-frankfurt-1.oraclecloud.com/p/ZH4wscDHe59T28yVJtrMH8uqifI_ih0NL5IyqxXQjSo/n/chryswoods/b/biosimspace_releases/o/biosimspace_devel_latest_linux.run).
   Then run:

   ```
   chmod +x biosimspace_devel_latest_linux.run
   ./biosimspace_devel_latest_linux.run
   ```
   Unless put somewhere else you should now be able to find a directory called biosimspace.app in your home directory. The python version that is available with it, should be used to run any BioSimSpace script. 

   There are future plans to fully incorporate BioSimSpace into conda, but they have not yet been finalised. 

2. Run the examples   
   BioSimSpace has some dependencies in order to run smoothly. For full functionality you need to have AmberTools18 and Gromacs installed as well as BioSimSpace. This will allow for solvation and dynamics to be run smoothly. 

   If you perfer working with docker containers, why not checkout the BioSimSpace docker image [here](https://cloud.docker.com/u/biosimspace/repository/docker/biosimspace/biosimspace-devel).

   The example that solvates and minimises benzene from a pdb file can be found in the `solvation` directory and readily be run as:

   ```
   ~/biosimspace.app/bin/python setup_molecule.py  --molecule benzene.pdb --boxsize 40 --parametrisation_base params --minimisation_base mini --verbose --forcefield gaff 
   ```

   To get more detailed help on setup_molecule.py you can simply run:
   ```
   ~/biosimspace.app/bin/python setup_molecule.py --help
   ```
