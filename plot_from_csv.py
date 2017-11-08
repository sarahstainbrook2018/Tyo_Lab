__author__ = 'Sarah'
__update__ = '11/8/2017'

# This program will take an instruction file and a csv data file and output
# graphs of the dose response curves with both the points and errorbars, and
# the smoothed fit.
#
# Script should be run in IDLE directly.
#
# NOTE: The file plotdirections.txt must be in the same folder.
#
# Template file 'plotdirections.txt' should be a plain text file consisting of:
#   the input file, with path if it's not located in the same folder (no quotes)
#   one line for each strain to be plotted:
#       the name of the strain in the csv data file,
#       the line color (no quote needed), the line width (int), the line style
#       of the smoothed fit (- for solid, . for no line), the line style of the
#       non-smoothed raw data plot line. No quotes on anything in this line.
#       Ex: ySS008_RFP,k,1,-,.
#   the title of the plot (no quotes)

import csv
import numpy as np
import matplotlib.pylab as plt


def func(x, bot, top, ec50, hill):
	return bot + (top - bot)/(1 + 10**((np.log10(ec50) - x)*hill))


params=[]
x_data = [0.0, 1.69897000434, 2.0, 2.69897000434, 3.0, 3.69897000434, 4.0]
# If you use alpha-factor concentrations other than standard, modify x_data.
x=np.linspace(0,4.0,1000)

# reads in the instructions file and its arguments
with open('plotdirections.txt','r') as f:
    linecount = sum(1 for line in f)
    f.seek(0)        #returns reader to the beginning of the file
    datafile=str(f.readline())
    datafile = datafile.strip()
    for counter in range(0,linecount-2):
        args=f.readline()
        args = args[0:-1]
        argssplit=[str(item) for item in args.split(',',5)]
        params.append(argssplit)
    title=f.readline()

# reads in the datafile
with open(datafile,'rU') as infile:
    reader=csv.reader(infile,dialect='excel',delimiter=',')
    newrow = []
    sampledict={}
    for row in reader:
        newrow.append(row)
    for index in range(0,len(newrow),8):
        sampledict[str(newrow[index][0])]=[[],[],[]]
        for index2 in range(1,(len(x_data)+1)):
            sampledict[str(newrow[index][0])][0].append(float
                                                        (newrow[((index)+index2)][0]))
            sampledict[str(newrow[index][0])][1].append(float
                                                        (newrow[((index)+index2)][1]))
            sampledict[str(newrow[index][0])][2].append(float
                                                        (newrow[((index)+index2)][2]))

# plots the strains specified in the instructions file
for counter in range(0,linecount-2):
    strainID = str(params[counter][0])
    fit = func(x,sampledict[strainID][2][0],sampledict[strainID][2][1],
               sampledict[strainID][2][2],sampledict[strainID][2][3])
    plt.plot(x,fit,linewidth=int(params[counter][2]),linestyle=
             params[counter][3],color=params[counter][1],label=strainID)
    plt.errorbar(x_data,sampledict[strainID][0],sampledict[strainID][1],
                 linewidth=int(params[counter][2]),color=params[counter][1],
                 marker='o',linestyle=params[counter][4])#, label=strainID)

plt.xlabel('[alpha factor], log nM', size=20)
plt.ylabel('Intensity, arbitrary units', size=20)
plt.suptitle(str(title), fontsize=25)
plt.xlim(-0.01,4.01)
plt.legend(loc=2)
plt.show()
