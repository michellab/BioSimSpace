# Set the minimum allowed Sire version.
min_ver = "2019.1.0"
min_ver_int = int(min_ver.replace(".", ""))

# Make sure we're using the Sire python interpreter.
try:
    import Sire.Base
    bin_dir = Sire.Base.getBinDir()
    lib_dir = Sire.Base.getLibDir()
except ModuleNotFoundError:
    raise ModuleNotFoundError("BioSimSpace currently requires the Sire Python interpreter: www.siremol.org")

# Check the Sire version.
if int(Sire.__version__.replace(".", "")) < min_ver_int:
    raise ImportError("BioSimSpace requires Sire version '%s' or above." % min_ver)

from setuptools import setup, find_packages

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
    import os
    import sys

    if "install" in sys.argv and not os.getenv("BSS_SKIP_INSTALL"):
        import subprocess

        # Install Python dependencies and enable Jupyter widget extensions.
        print("\nSetting up python environment...")

        # Open files for stdout/stderr.
        stdout = open("setup.out", "w")
        stderr = open("setup.err", "w")

        print("Updating conda")
        command = "%s/conda update -y -q -n base conda" % bin_dir
        subprocess.run(command, shell=True, stdout=stdout, stderr=stderr)

        print("Adding conda-forge channel")
        command = "%s/conda config --system --prepend channels conda-forge" % bin_dir
        subprocess.run(command, shell=True, stdout=stdout, stderr=stderr)

        print("Disabling conda auto update")
        command = "%s/conda config --system --set auto_update_conda false" % bin_dir
        subprocess.run(command, shell=True, stdout=stdout, stderr=stderr)

        print("Installing package: mdtraj")
        command = "%s/conda install -y -q -c omnia mdtraj" % bin_dir
        subprocess.run(command, shell=True, stdout=stdout, stderr=stderr)

        print("Installing package: mdanalysis")
        command = "%s/conda install -y -q mdanalysis" % bin_dir
        subprocess.run(command, shell=True, stdout=stdout, stderr=stderr)

        print("Upgrading pip")
        command = "%s/pip install --upgrade pip" % bin_dir
        subprocess.run(command, shell=True, stdout=stdout, stderr=stderr)

        print("Installing package: configargparse")
        command = "%s/pip install configargparse" % bin_dir
        subprocess.run(command, shell=True, stdout=stdout, stderr=stderr)

        print("Installing package: pyyaml")
        command = "%s/pip install pyaml" % bin_dir
        subprocess.run(command, shell=True, stdout=stdout, stderr=stderr)

        print("Installing package: pypdb")
        command = "%s/pip install pypdb" % bin_dir
        subprocess.run(command, shell=True, stdout=stdout, stderr=stderr)

        print("Installing package: watchdog")
        command = "%s/pip install watchdog" % bin_dir
        subprocess.run(command, shell=True, stdout=stdout, stderr=stderr)

        print("Installing package: jupyter")
        command = "%s/pip install jupyter" % bin_dir
        subprocess.run(command, shell=True, stdout=stdout, stderr=stderr)

        print("Installing package: duecredit")
        command = "%s/pip install duecredit" % bin_dir
        subprocess.run(command, shell=True, stdout=stdout, stderr=stderr)

        print("Installing package: mock")
        command = "%s/pip install mock" % bin_dir
        subprocess.run(command, shell=True, stdout=stdout, stderr=stderr)

        print("Installing package: pygtail")
        command = "%s/pip install pygtail" % bin_dir
        subprocess.run(command, shell=True, stdout=stdout, stderr=stderr)

        print("Installing package: pymbar")
        command = "%s/pip install pymbar" % bin_dir
        subprocess.run(command, shell=True, stdout=stdout, stderr=stderr)

        print("Installing package: fileupload")
        command = "%s/pip install fileupload" % bin_dir
        subprocess.run(command, shell=True, stdout=stdout, stderr=stderr)

        print("Activating notebook extension: fileupload")
        command = "%s/jupyter-nbextension install fileupload --py --sys-prefix --log-level=0" % bin_dir
        subprocess.run(command, shell=True, stdout=stdout, stderr=stderr)
        command = "%s/jupyter-nbextension enable fileupload --py --sys-prefix" % bin_dir
        subprocess.run(command, shell=True, stdout=stdout, stderr=stderr)

        print("Installing package: nglview")
        command = "%s/pip --no-cache-dir install nglview" % bin_dir
        subprocess.run(command, shell=True, stdout=stdout, stderr=stderr)

        print("Activating notebook extension: nglview")
        command = "%s/jupyter-nbextension install nglview --py --sys-prefix --log-level=0" % bin_dir
        subprocess.run(command, shell=True, stdout=stdout, stderr=stderr)
        command = "%s/jupyter-nbextension enable nglview --py --sys-prefix" % bin_dir
        subprocess.run(command, shell=True, stdout=stdout, stderr=stderr)

        print("Cleaning conda environment")
        command = "%s/conda clean -all -y -q" % bin_dir
        subprocess.run(command, shell=True, stdout=stdout, stderr=stderr)

        # Close the file handles.
        stdout.close()
        stderr.close()

        try:
            import BioSimSpace

            # Installation worked. Remove the stdout/stderr files.
            os.remove("setup.out")
            os.remove("setup.err")
        except:
            print("\nPossible installation issues. Please check output in 'setup.out' and 'setup.err'")
            sys.exit()

        print("\nDone!")

        print("\nIf you have problems with Jupyter permissions, try removing '$HOME/.jupyter' or '$HOME/.local/share/jupyter'")

        print("\nFor optional package support...")
        print("AMBER:   http://ambermd.org")
        print("GROMACS: http://www.gromacs.org")
        print("NAMD:    http://www.ks.uiuc.edu/Research/namd")
