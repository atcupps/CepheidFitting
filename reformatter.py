'''
This program reformats a SIMBAD output file as a new file compatible with
IRSA multi-object search. Some manual checking of the output file may be
necessary after running the program to eliminate occasional errors due to
SIMBAD formatting inconsistencies.
Written by Andrew Cupps.

In order to run this program: Replace the `input_file_name` and
`output_file_name` variables; then run.
'''
input_file_name = "input.txt"
output_file_name = "output.txt"

# Open SIMBAD file and output file
input = open(input_file_name, "r")
output = open(output_file_name, "w")

# Write headers for ouput file
output.write("\ EQUINOX = 'J2000.0'\n")
output.write("\ Query\n")
output.write("|   ra            |   dec           |\n")
output.write("|   double        |   double        |\n")
output.write("|   deg           |   deg           |\n")

# Read through and discard first two line (header with column names)
_ = input.readline()
_ = input.readline()

# Read through all data
lines = list(input.readlines())
for line in lines:
    coords = line[46:81].split(' ')
    ra = coords[0]
    dec = coords[1]
    if dec[-1] == '|':
        dec = dec[0:-1]
    output.write(f" {ra:17} {dec:17}\n")

input.close()
output.close()