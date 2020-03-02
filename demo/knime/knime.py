# Example script showing how to interface BioSimSpace with KNIME.

import pandas
import subprocess
import yaml

# Set the BioSimSpace Python interpreter.
bss_python = "/home/lester/sire.app/bin/python"

# Set the script name.
script = "minimisation.py"

# Set the name of the BioSimSpace output YAML file. (This is the default.)
output = "output.yaml"

# The KNIME Python extension assumes Pandas DataFrame objects as input/output
# of nodes. While this works well for tabular data, it is a poor format for
# the unstructured data required by BioSimSpace. As a hack, we assume the input
# consists of a DataFrame containing a single Series called "input". This
# is simply a dictionary containing key/value pairs for all of the required
# input to the script. BioSimSpace is able to convert strings to any of its
# built in types, so strings can be used for the value of all input variables.

# (Note that the input DataFrame, "input_table", should be generated elsewhere
# and passed as input to this node.)

# For example, create a dictionary to hold the input arguments.
input_dict = {"steps" : 1000, "files" : ["amber/ala/ala.crd", "amber/ala/ala.top"]}

# Insert the dictionary into the input DataFrame with label "input".
input_table = pandas.DataFrame(data={"input" : input_dict})

# Write a YAML configuration file for the BioSimSpace script.
with open("input.yaml", "w") as file:
    yaml.dump(input_table["input"].to_dict(), file, default_flow_style=False)

# Generate the shell command.
command = "%s %s -c input.yaml" % (bss_python, script)

# Run the BioSimSpace script as a subprocess and raise any error if it fails.
try:
    subprocess.run(command, shell=True, check=True)
except:
    raise

# Read the output YAML file into a dictionary.
with open(output, "r") as file:
    output_dict = yaml.safe_load(file)

# Insert the dictionary into the output DataFrame with label "output".
output_table = pandas.DataFrame(data={"output" : output_dict})
