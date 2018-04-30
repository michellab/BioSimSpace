# Make sure we're using the Sire python interpreter.
try:
    import Sire.Base
    bin_dir = Sire.Base.getBinDir()
    lib_dir = Sire.Base.getLibDir()
except ModuleNotFoundError:
    raise ModuleNotFoundError("BioSimSpace currently requires the Sire Python interpreter: www.siremol.org")

from setuptools import setup, find_packages

import subprocess
import versioneer

# A list of authors and their email addresses.
authors=("Lester Hedges <lester.hedges@gmail.com, "
         "Christopher Woods <chryswoods@gmail.com>, "
         "Antonia Mey <antonia.mey@gmail.com")

# Run the setup.
try:
    setup(name='BioSimSpace',
          version=versioneer.get_version(),
          cmdclass=versioneer.get_cmdclass(),
          description='BioSimSpace: Making biomolecular simulation a breeze.',
          author=authors,
          url='https://github.com/michellab/BioSimSpace',
          license='GPLv2',
          packages=find_packages(),
          zip_safe=True
        )

# Post setup configuration.
finally:
    print("\nSetting up python environment...")

    # Install Python dependencies and enable Jupyter widget extensions.

    print("Updating conda")
    command = "%s/conda update -y -q -n base conda" % bin_dir
    subprocess.run(command, shell=True, stdout=subprocess.PIPE)

    print("Adding conda-forge channel")
    command = "%s/conda config --system --prepend channels conda-forge" % bin_dir
    subprocess.run(command, shell=True, stdout=subprocess.PIPE)

    print("Disabling conda auto update")
    command = "%s/conda config --system --set auto_update_conda false" % bin_dir
    subprocess.run(command, shell=True, stdout=subprocess.PIPE)

    print("Installing package: 'mdtraj'")
    command = "%s/conda install -y -q -c omnia mdtraj" % bin_dir
    subprocess.run(command, shell=True, stdout=subprocess.PIPE)

    print("Upgrading pip")
    command = "%s/pip install --upgrade pip" % bin_dir
    subprocess.run(command, shell=True, stdout=subprocess.PIPE)

    print("Installing package: 'watchdog'")
    command = "%s/pip install watchdog" % bin_dir
    subprocess.run(command, shell=True, stdout=subprocess.PIPE)

    print("Installing package: 'mdanalysis'")
    command = "%s/pip install mdanalysis" % bin_dir
    subprocess.run(command, shell=True, stdout=subprocess.PIPE)

    print("Installing package: 'jupyter'")
    command = "%s/pip install jupyter" % bin_dir
    subprocess.run(command, shell=True, stdout=subprocess.PIPE)

    print("Installing package: 'duecredit'")
    command = "%s/pip install duecredit" % bin_dir
    subprocess.run(command, shell=True, stdout=subprocess.PIPE)

    print("Installing package: 'pygtail'")
    command = "%s/pip install pygtail" % bin_dir
    subprocess.run(command, shell=True, stdout=subprocess.PIPE)

    print("Installing package: 'matplotlib'")
    command = "%s/pip install matplotlib" % bin_dir
    subprocess.run(command, shell=True, stdout=subprocess.PIPE)

    print("Installing package: 'fileupload'")
    command = "%s/pip install fileupload" % bin_dir
    subprocess.run(command, shell=True, stdout=subprocess.PIPE)

    print("Activating notebook extension: 'fileupload'")
    command = "%s/jupyter-nbextension install fileupload --py --sys-prefix" % bin_dir
    subprocess.run(command, shell=True, stdout=subprocess.PIPE)

    command = "%s/jupyter-nbextension enable fileupload --py --sys-prefix" % bin_dir
    subprocess.run(command, shell=True, stdout=subprocess.PIPE)

    print("Installing package: 'nglview'")
    command = "%s/pip --no-cache-dir install nglview" % bin_dir
    subprocess.run(command, shell=True, stdout=subprocess.PIPE)

    print("Activating notebook extension: 'nglview'")
    command = "%s/jupyter-nbextension install nglview --py --sys-prefix" % bin_dir
    subprocess.run(command, shell=True, stdout=subprocess.PIPE)

    command = "%s/jupyter-nbextension enable nglview --py --sys-prefix" % bin_dir
    subprocess.run(command, shell=True, stdout=subprocess.PIPE)

    print("Patching 'mdtraj' NetCDF writer")
    command = ("%s/curl -k -s https://raw.githubusercontent.com/scipy/scipy/master/scipy/io/netcdf.py"
               "> %s/python3.5/site-packages/scipy/io/netcdf.py"
              ) % (bin_dir, lib_dir)
    subprocess.run(command, shell=True, stdout=subprocess.PIPE)

    print("\nDone!")

    print("\nIf you have problems with Jupyter permissions, try removing '$HOME/.jupyter' or '$HOME/.local/share/jupyter'")
