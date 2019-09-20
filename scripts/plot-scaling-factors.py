import argparse
import math
import numpy
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(prog='plot-scaling-factors.py', description=__doc__)
parser.add_argument('histo', metavar='HISTO', help='.histo file output of PGGTyper with scaling factors.')
parser.add_argument('--max-value', default='100', metavar='MAX_VALUE', help='max value to plot on x-axis (default: 100).')
args = parser.parse_args()

# read the histogram from input file
x_bins = []
abundances = []

x_range = 0.0
x_unit = 0.0

print('Reading file ...')

for line in open(args.histo, 'r'):
	parsed_line = line.split();
	if parsed_line[0] == "xvalues:":
		x_range = float(parsed_line[1])
		x_unit = float(parsed_line[2])
		break
	x_bin = int(parsed_line[0])
	y_value = int(parsed_line[1])
	x_bins.append(x_bin);
	abundances.append(y_value);
print('Plotting ..')
# plot distributions
x_values = [x_range + x_unit*i for i in x_bins]
plt.plot(x_values, abundances, marker='o', linestyle=':')
print('1')
plt.xlim(-int(args.max_value), int(args.max_value))
print(2)
#plt.yticks(range(0, max(abundances)+1, 2))
print(3)
plt.xticks(numpy.arange(-int(args.max_value), int(args.max_value)+1, 5*x_unit))

#plt.yscale('log')
plt.title('histogram of kmer scaling factors')
plt.xlabel('scaling factor')
plt.ylabel('count')
plt.savefig(args.histo + '.pdf')
plt.close()
