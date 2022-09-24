import os
import platform

if not os.getenv("BSS_CONDA_INSTALL"):
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

# Function to check if a conda dependency has been installed
def is_installed(dep: str, conda: str):
    p = subprocess.Popen([conda, "list", dep], stdout=subprocess.PIPE)
    lines = str(p.stdout.read())

    return lines.find(dep) != -1


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
          include_package_data=True,
          zip_safe=False
        )

# Post setup configuration.
finally:
    import sys

    if "install" in sys.argv and not (os.getenv("BSS_CONDA_INSTALL") or os.getenv("BSS_SKIP_DEPENDENCIES")):
        import shlex
        import subprocess

        # Install Python dependencies and enable Jupyter widget extensions.
        print("\nSetting up python environment...")

        # Open files for stdout/stderr.
        stdout = sys.stdout
        stderr = sys.stderr

        # Create a list of the conda dependencies.
        conda_deps = ["configargparse",
                      "pygtail",
                      "pytest",
                      "pyyaml",
                      "watchdog",
                      "pydot",
                      "networkx",
                      "nglview",
                      "ipywidgets<8",
                      "py3dmol",
                      "pypdb",
                      "rdkit",
                      "parmed",
                      "lomap2",
                      "mdtraj",             # known not available on aarch64
                      "mdanalysis",         # known not available on aarch64
                      "openff-toolkit"      # known not available on aarch64
                     ]

        # Don't try to install things that are already installed...
        to_install_deps = []

        print("Checking for dependencies that are already installed...")

        for dep in conda_deps:
            if not is_installed(dep, conda="%s/conda" % bin_dir):
                to_install_deps.append(dep)
            else:
                print("Already installed %s" % dep)

        conda_deps = to_install_deps

        print("Adding conda-forge channel")
        command = "%s/conda config --system --prepend channels conda-forge" % bin_dir
        subprocess.run(shlex.split(command), shell=False, stdout=stdout, stderr=stderr)

        print("Disabling conda auto update")
        command = "%s/conda config --system --set auto_update_conda false" % bin_dir
        subprocess.run(shlex.split(command), shell=False, stdout=stdout, stderr=stderr)

        print("Installing conda dependencies: %s" % ", ".join(conda_deps))
        command = "%s/conda install -y -q %s" % (bin_dir, " ".join(conda_deps))

        all_installed_ok = True

        try:
            subprocess.run(shlex.split(command), shell=False,
                           stdout=stdout, stderr=stderr, check=True)
        except Exception:
            all_installed_ok = False

        if not all_installed_ok:
            print("There were errors installing some of the dependencies.")
            print("We will now try to install them one-by-one. This may take some time...")

            failures = []

            for dep in conda_deps:
                if not is_installed(dep, conda="%s/conda" % bin_dir):
                    print("Trying again to install '%s'" % dep)
                    command = "%s/conda install -y -q %s" % (bin_dir, dep)

                    try:
                        subprocess.run(shlex.split(command), shell=False,
                                       stdout=stdout, stderr=stderr, check=True)
                    except Exception:
                        failures.append(dep)

            if len(failures) == 0:
                print("All dependencies installed successfully!")
            else:
                print("\n** Failed to install these dependencies: %s" % ", ".join(failures))
                print("** BioSimSpace will still install and run, but some functionality may not be available.\n")
        else:
            print("All dependencies install successfully first time!")


        print("Activating notebook extension: nglview")
        command = "%s/jupyter-nbextension install nglview --py --sys-prefix --log-level=0" % bin_dir
        subprocess.run(shlex.split(command), shell=False, stdout=stdout, stderr=stderr)
        command = "%s/jupyter-nbextension enable nglview --py --sys-prefix" % bin_dir
        subprocess.run(shlex.split(command), shell=False, stdout=stdout, stderr=stderr)

        print("Cleaning conda environment")
        command = "%s/conda clean --all --yes --quiet" % bin_dir
        subprocess.run(shlex.split(command), shell=False, stdout=stdout, stderr=stderr)

        try:
            import BioSimSpace
        except:
            print("\nPossible installation issues.")
            sys.exit()

        print("\nDone!")

        print("\nIf you have problems with Jupyter permissions, try removing '$HOME/.jupyter' or '$HOME/.local/share/jupyter'")

        print("\nFor optional package support...")
        print("AMBER:   http://ambermd.org")
        print("GROMACS: http://www.gromacs.org")
        print("NAMD:    http://www.ks.uiuc.edu/Research/namd")
        print("FKCOMBU: https://pdbj.org/kcombu")
