__author__ = 'Sarah, heavily borowed from Jessica'

# This program will take a template file and fit a sigmoidal curve to each 
# set of data and their average. The curve will be shown and the user asked
# whether the fit is acceptable. If yes, the averages for each concentration,
# standard deviations and curve parameters will be appended to an existing
# csv file compatible with the 'plot_from_csv.py' program. Otherwise, program
# closes to give the user a chance to modify the input data or parameters.
#
# Script should be run in IDLE directly.
#
# NOTE: The file utilitiesJ.py and DRC_TEMPLATE.txt must be in the same folder.
#
# Template file 'DRC_TEMPLATE.txt' should be an eleven-line plain text file:
# 1.  absolute path to directory containing the .fcs files
# 2.  channel name (choose from Cy3-A, DAPI-A, PE-Texas Red-A, CFP-A, Alexa Fluor
#       488-A. The channel names are listed in the FCS file header when it is
#       opened as a text file. No quotes. 
# 3.  well rows (capital letters, separated by spaces) of interest
# 4.  well columns (with leading 0 if < 10, separated by spaces) of interest
# 5.  concentrations (in nM), use 1 instead of 0 for the control
# 6.  orientation (byrow or bycolumn, quotes not necessary)
# 7.  the slope of the linear correction for plates being stored at 4C prior to
#       being read. (float) (to leave data uncorrected, specify slope of 0.0,
#       intercept of 1.0 below)
# 8.  the intercept of the correction line (float)
# 9.  the average 0nM fluorescence of the reference strain
# 10. the name of the strain
# 11. the name of the output file to be appended (not in quotes. must be in the
#       same folder. Do NOT include a \n character on this line.)  

import utilitiesJ as util
import numpy as np
import matplotlib.pylab as plt
import csv
from FlowCytometryTools import FCPlate
from sys import argv
from scipy.optimize import curve_fit

def calculate_adjusted_mean(sample):
	data = sample.get_data()
	FP_data = data[channel]
	return np.mean(FP_data)

# Parse template file.
with open('DRC_TEMPLATE.txt','r') as f:
	data_dir = util.rm_newline(f.readline())
	channel = f.readline()
	channel = str(channel.strip())
	rows = f.readline()
	cols = f.readline()
	conc = f.readline()
	conc = [float(x) for x in conc.split()]
	direction = util.rm_newline(f.readline())
	slope = float(f.readline())  
	b = float(f.readline())
	WT0 = float(f.readline())
	strainname = f.readline()
	strainname = strainname.strip()
	outfile = f.readline()

# Parse all data files into plate.
print '\n>> Getting data for each well'
plate, first_well = util.parse_plate(data_dir, rows, cols)

print '>> Calculating arithmetic mean for normalized GFP'
means = plate.apply(calculate_adjusted_mean)

# Extract values from means array.
m = means.values;
m = m[~np.all(np.isnan(m), axis = 1)]

# Correct data to ySS008 controls within the plate (added Sarah 10/18/2017)
for counter in range(0,len(conc)):
        if direction == 'bycolumn':
                conccorr = util.correction(np.log10(conc[counter]),slope,b)
                for counter2 in range(0,len(cols.split())):
                    m[counter,counter2]=m[counter,counter2]/conccorr
        elif direction == 'byrow':
                print 'Program is cannot handle data in by row format.'

# Normalize data to the wild-type 0nM fluorescence (added by Sarah 11/4/2017)

m[:,:]= [x/WT0 for x in m]
for i in range(0,len(conc)):
        print m[i,0], '\t', m[i,1],'\t', m[i,2]

# Calculate average and standard deviation for means across replicates.
m_avg = [x for x in np.average(m, axis = 1) if str(x) != 'nan']
m_std = [x for x in np.std(m, axis = 1) if str(x) != 'nan']

# Fit data to curve for average.
print '\n','>> Fitting curve(s)'

n = len(cols.split())
params = np.empty([len(conc), n + 1])
for i in range(0,  n + 1):
	if i == 0:              
		y_data = m_avg
        # this makes the first column be averages, while following columns are
        # calculated for individual columns.
	else:
		y_data = m[:,i - 1]
	x_data = np.log(conc)/np.log(10)
	popt, pcov = curve_fit(util.func, x_data, y_data, maxfev = 100000)
	params[0,i] = popt[0]
	params[1,i] = popt[1]
	params[2,i] = 10**popt[2]
	params[3,i] = popt[3]
	for counter in range(4,len(conc)): params[counter,i]=0

# Print parameters.
param_names = ['MIN','MAX','EC50','HILL']
for i in range(0,4):
	print '\n',param_names[i],'\t',params[i,0],'\t',np.std(params[i,1::]),'\t',
	for j in range (1, n + 1):
		print params[i,j] , '\t' ,
print '\n'

# Plot data points and fitted curve.
x = np.linspace(-1, 5, 100)
y = util.func(x, params[0,0], params[1,0], np.log10(params[2,0]), params[3,0])
plt.plot(x, y, marker = 'None', lw = 2, color = 'k')
plt.errorbar(x_data, m_avg, yerr = m_std, mfc = 'b', ls = 'None', marker = 'o')
plt.xlabel('[ligand]', size = 15)
plt.ylabel('Normalized Response', size = 15)
plt.show()


# If user accepts the fit as shown, appends data to csv file 
accept=input('Accept this fit? y or n')
if str(accept) == 'y':                          
        print 'writing to ', str(outfile)
        with open(outfile,'ab') as ofile:
                firstrow = strainname,0,0
                csv.writer(ofile).writerow(firstrow)
                for counter in range(0,len(conc)):
                        newrow = m_avg[counter],m_std[counter],params[counter,0]
                        csv.writer(ofile).writerow(newrow)
  
if str(accept) == 'n':
        print 'closing program.'                       
