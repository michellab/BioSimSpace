import os
import sys
import glob

script = os.path.abspath(sys.argv[0])

# go up one directories to get the source directory
# (this script is in Sire/actions/)
srcdir = os.path.dirname(os.path.dirname(script))

print(f"BioSimSpace source is in {srcdir}\n")

# Get the anaconda token to authorise uploads
if "ANACONDA_TOKEN" in os.environ:
    conda_token = os.environ["ANACONDA_TOKEN"]
else:
    conda_token = "TEST"

# Get the anaconda channel labels.
if "ANACONDA_LABEL" in os.environ:
    conda_label = os.environ["ANACONDA_LABEL"]
else:
    conda_label = "dev"

# get the root conda directory
conda = os.environ["CONDA"]

# Set the path to the conda-bld directory.
conda_bld = os.path.join(conda, "envs", "bss_build", "conda-bld")

print(f"conda_bld = {conda_bld}")

# Find the packages to upload
bss_pkg = glob.glob(os.path.join(conda_bld, "*-*", "biosimspace-*.tar.bz2"))

if len(bss_pkg) == 0:
    print("No BioSimSpace packages to upload?")
    sys.exit(-1)

packages = bss_pkg

print(f"Uploading packages:")
print(" * ", "\n *  ".join(packages))

packages = " ".join(packages)


def run_cmd(cmd):
    import subprocess

    p = subprocess.Popen(cmd.split(), stdout=subprocess.PIPE)
    return str(p.stdout.read().decode("utf-8")).lstrip().rstrip()


gitdir = os.path.join(srcdir, ".git")

tag = run_cmd(f"git --git-dir={gitdir} --work-tree={srcdir} tag --contains")

# Upload the packages to the openbiosim channel on Anaconda Cloud.
cmd = f"anaconda --token {conda_token} upload --user openbiosim --label {conda_label} --force {packages}"

print(f"\nUpload command:\n\n{cmd}\n")

# Label release packages with main and dev so that dev is at least as new as
# main. Only need to uncomment the libcpuid and fkcombu package uploads when
# there new versions are released.
if conda_token == "TEST":
    print("Not uploading as the ANACONDA_TOKEN is not set!")
    sys.exit(-1)

output = run_cmd(cmd)

print(output)

print("Package uploaded!")
