import BioSimSpace as BSS

# A script that is run as a subprocess by test_node to allow testing
# of command-line argument parsing and validation.

# Create a node object with a description of its use.
node = BSS.Gateway.Node("A node to test command-line argument parsing.")

# Set some input requirements.
node.addInput("bool", BSS.Gateway.Boolean(help="A boolean requirement", default=False))
node.addInput("int", BSS.Gateway.Integer(help="An integer requirement."))
node.addInput("float", BSS.Gateway.Float(help="A float requirement."))
node.addInput("string", BSS.Gateway.String(help="A string requirement."))
node.addInput("file", BSS.Gateway.File(help="A file requirement."))
node.addInput("fileset", BSS.Gateway.FileSet(help="A file set requirement."))
node.addInput("temperature", BSS.Gateway.Temperature(help="A temperature requirement.", unit="kelvin"))
node.addInput("time", BSS.Gateway.Time(help="A time requirement.", unit="minute"))
node.addInput("length", BSS.Gateway.Length(help="A length requirement.", unit="angstrom"))
node.addInput("area", BSS.Gateway.Area(help="An area requirement.", unit="angstrom squared"))
node.addInput("volume", BSS.Gateway.Volume(help="A volume requirement.", unit="angstrom cubed"))
node.addInput("angle", BSS.Gateway.Angle(help="An angle requirement.", unit="radian"))
node.addInput("charge", BSS.Gateway.Charge(help="A charge requirement.", unit="electron charge"))
node.addInput("energy", BSS.Gateway.Energy(help="An energy requirement.", unit="kcal/mol"))
node.addInput("pressure", BSS.Gateway.Pressure(help="A pressure requirement.", unit="atmosphere"))

# Set some output requirements.
node.addOutput("out1", BSS.Gateway.Float(help="a float requirement"))
node.addOutput("out2", BSS.Gateway.String(help="a string requirement"))

# Get the values of the parsed inputs.
my_bool = node.getInput("bool")
my_int = node.getInput("int")
my_float = node.getInput("float")
my_string = node.getInput("string")
my_file = node.getInput("file")
my_fileset = node.getInput("fileset")
my_temperature = node.getInput("temperature")
my_time = node.getInput("time")
my_length = node.getInput("length")
my_area = node.getInput("area")
my_volume = node.getInput("volume")
my_charge = node.getInput("charge")
my_energy = node.getInput("energy")
my_pressure = node.getInput("pressure")

# Create some variables based on the input.
out1 = my_int * my_float
out2 = "%s %.1f" % (my_string, my_float)

# Set the output.
node.setOutput("out1", out1)
node.setOutput("out2", out2)

# Validate the output of the node.
node.validate()
