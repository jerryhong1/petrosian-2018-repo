# used to create csv of data to visualize in d3

from __future__ import division
import csv
import numpy as np
import quasars as quas
import math
#plt.style.use('mystyle.mplstyle')
#plt.style.use('poster.mplstyle')

# In[1]: import file and put into QuasarData anLarad band objects.
file = open("SDSSFIRSTComplete.dat","r") 
lines = file.readlines()
data = np.empty((len(lines) - 5, 3))
for i in range(5, len(lines)):
    linearr = np.array(lines[i].strip().split())
    for j in [0, 5]:
        data[i - 5, 0] = float(linearr[0])
        data[i - 5, 1] = float(linearr[6])

# z = 0, L_r = 1, Lrmin = 2
header = ['z', 'L', 'Lmin']

fmin_rad = 1.0e-26 #1 mJy        
def lmin(z):
    k = quas.kcorrect_rad(z);
    return 4 * math.pi * (quas.d_lum(z)**2) * fmin_rad / k

Lmin = [lmin(z) for z in data[:, 0]]

data[:, 2] = Lmin

with open("sampledata.csv", "w") as csvfile:
    writer = csv.writer(csvfile, delimiter=',',
                            quotechar='|', quoting=csv.QUOTE_MINIMAL)
    writer.writerow(['Redshift', 'L', 'Lmin'])
    for i in range(len(data[:,0]))[::50]:
        writer.writerow([data[i, 0], data[i, 1], data[i, 2]])