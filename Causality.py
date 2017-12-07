# -*- coding: utf-8 -*-
"""
Created on Mon Dec  4 14:53:19 2017

@author: Max
"""


## Import all the librairies

import csv;
import math as m;
import numpy as np;
import matplotlib.pylab as plt;
from statsmodels.tsa.stattools import grangercausalitytests;
from scipy.stats.stats import pearsonr;
from scipy.stats.stats import pointbiserialr;


## Prepare the data
with open('C:/Users/Max/Documents/UPM/Data_Analysis/Mini_pro/Allianz.csv', newline='') as csvfile:
    assur1_x = [];
    assur1_y = [];
    
    reader = csv.DictReader(csvfile);
    for row in reader:      
            assur1_x.append(row['Date']);
            if (row['AdjClose']!="null"):
                assur1_y.append([float(row['AdjClose'])]);
            else:
                assur1_y.append([float(0)]);
    newAssur1_y = [];
    newAssur1_y.append([float(assur1_y[0][0])]);
    for i in range(1,len(assur1_y)):
        newAssur1_y.append([m.log10(float(assur1_y[i][0]))/float(assur1_y[i-1][0])]);
        
with open('C:/Users/Max/Documents/UPM/Data_Analysis/Mini_pro/AssicurazioniGenerali.csv', newline='') as csvfile:
    assur2_x = [];
    assur2_y = [];
    reader = csv.DictReader(csvfile)
    for row in reader:      
        if (row['AdjClose']!="null"):
            assur2_x.append(row['Date']);
            assur2_y.append([float(row['AdjClose'])]);
    newAssur2_y = [];
    newAssur2_y.append([float(assur2_y[0][0])]);
    for i in range(1,len(assur2_y)):
        newAssur2_y.append([m.log10(float(assur2_y[i][0]))/float(assur2_y[i-1][0])]);
        
with open('C:/Users/Max/Documents/UPM/Data_Analysis/Mini_pro/Bayer.csv', newline='') as csvfile:
    pharma1_x = [];
    pharma1_y = [];
    reader = csv.DictReader(csvfile)
    for row in reader:      
        if (row['AdjClose']!="null"):
            pharma1_x.append(row['Date']);
            pharma1_y.append([float(row['AdjClose'])]);
    newPharma1_y = [];
    newPharma1_y.append([float(pharma1_y[0][0])]);
    for i in range(1,len(pharma1_y)):
        newPharma1_y.append([m.log10(float(pharma1_y[i][0]))/float(pharma1_y[i-1][0])]);
        
with open('C:/Users/Max/Documents/UPM/Data_Analysis/Mini_pro/LOreal.csv', newline='') as csvfile:
    pharma2_x = [];
    pharma2_y = [];
    reader = csv.DictReader(csvfile)
    for row in reader:      
        if (row['AdjClose']!="null"):
            pharma2_x.append(row['Date']);
            pharma2_y.append([float(row['AdjClose'])]);
    newPharma2_y = [];
    newPharma2_y.append([float(pharma2_y[0][0])]);
    for i in range(1,len(pharma2_y)):
        newPharma2_y.append([m.log10(float(pharma2_y[i][0]))/float(pharma2_y[i-1][0])]);
        
with open('C:/Users/Max/Documents/UPM/Data_Analysis/Mini_pro/BNP.csv', newline='') as csvfile:
    bank1_x = [];
    bank1_y = [];
    reader = csv.DictReader(csvfile)
    for row in reader:      
        if (row['AdjClose']!="null"):
            bank1_x.append(row['Date']);
            bank1_y.append([float(row['AdjClose'])]);
    newBank1_y = [];
    newBank1_y.append([float(bank1_y[0][0])]);
    for i in range(1,len(bank1_y)):
        newBank1_y.append([m.log10(float(bank1_y[i][0]))/float(bank1_y[i-1][0])]);
    
with open('C:/Users/Max/Documents/UPM/Data_Analysis/Mini_pro/BBVA.csv', newline='') as csvfile:
    bank2_x = [];
    bank2_y = [];
    reader = csv.DictReader(csvfile)
    for row in reader:      
        if (row['AdjClose']!="null"):
            bank2_x.append(row['Date']);
            bank2_y.append([float(row['AdjClose'])]);
    newBank2_y = [];
    newBank2_y.append([float(bank2_y[0][0])]);
    for i in range(1,len(bank2_y)):
        newBank2_y.append([m.log10(float(bank2_y[i][0]))/float(bank2_y[i-1][0])]);
        
a1 = set(pharma1_x);
b1 = set(pharma2_x);
c1 = a1.intersection(b1);
c1n = list(c1);

a2 = set(bank1_x);
b2 = set(bank2_x);
c2 = a2.intersection(b2);
c2n = list(c2);

a3 = set(assur1_x);
b3 = set(assur2_x);
c3 = a3.intersection(b3);
c3n = list(c3);

for i in range(0,len(c1n)-1):
    if(pharma2_x[i] not in c1n):
        del pharma2_x[i]
        del newPharma2_y[i]
    if(pharma1_x[i] not in c1n):
        del pharma1_x[i]
        del newPharma1_y[i]

for i in range(0,len(c2n)-1):
    if(bank2_x[i] not in c2n):
        del bank2_x[i]
        del newBank2_y[i]
    if(bank1_x[i] not in c2n):
        del bank1_x[i]
        del newBank1_y[i]
        
for i in range(0,len(c3n)-1):
    if(assur2_x[i] not in c3n):
        del assur2_x[i]
        del newAssur2_y[i]
    if(assur1_x[i] not in c3n):
        del assur1_x[i]
        del newAssur1_y[i]


## Plot the data insurance
plt.figure(1)

x1 = np.linspace(0, 1, len(newAssur1_y)-1);

mean11 = np.mean(newAssur1_y[1:len(newAssur1_y)]);
std11 = np.std(newAssur1_y[1:len(newAssur1_y)]);

assur1_plot = [];

for e in newAssur1_y[1:len(newAssur1_y)]:
    assur1_plot.append((e[0]-mean11)/std11);

mean12 = np.mean(newAssur2_y[1:len(newAssur1_y)]);
std12 = np.std(newAssur2_y[1:len(newAssur1_y)]);

assur2_plot = [];

for e in newAssur2_y[1:len(newAssur1_y)]:
    assur2_plot.append((e[0]-mean12)/std12);
    
plt.plot(x1, assur1_plot, x1, assur2_plot)

plt.savefig("C:/Users/Max/Documents/UPM/Data_Analysis/Mini_pro/insurance.png")


## Plot the data pharmaceutic industry
plt.figure(2)

x2 = np.linspace(0, 1, len(newPharma1_y)-1);

mean21 = np.mean(newPharma1_y[1:len(newPharma1_y)]);
std21 = np.std(newPharma1_y[1:len(newPharma1_y)]);

pharma1_plot = [];

for e in newPharma1_y[1:len(newPharma1_y)]:
    pharma1_plot.append((e[0]-mean21)/std21);

mean22 = np.mean(newPharma2_y[1:len(newPharma1_y)]);
std22 = np.std(newPharma2_y[1:len(newPharma1_y)]);

pharma2_plot = [];

for e in newPharma2_y[1:len(newPharma1_y)]:
    pharma2_plot.append((e[0]-mean22)/std22);
    
plt.plot(x2, pharma1_plot, x2, pharma2_plot)

plt.savefig("C:/Users/Max/Documents/UPM/Data_Analysis/Mini_pro/pharma.png")


## Plot the data bank
plt.figure(3)

x3 = np.linspace(0, 1, len(newBank1_y)-1);

mean31 = np.mean(newBank1_y[1:len(newBank1_y)]);
std31 = np.std(newBank1_y[1:len(newBank1_y)]);

bank1_plot = [];

for e in newBank1_y[1:len(newBank1_y)]:
    bank1_plot.append((e[0]-mean31)/std31);

mean32 = np.mean(newBank2_y[1:len(newBank1_y)]);
std32 = np.std(newBank2_y[1:len(newBank1_y)]);

bank2_plot = [];

for e in newBank2_y[1:len(newBank1_y)]:
    bank2_plot.append((e[0]-mean32)/std32);
    
plt.plot(x3, bank1_plot, x3, bank2_plot)

plt.savefig("C:/Users/Max/Documents/UPM/Data_Analysis/Mini_pro/bank.png")


## Granger causality

hstack1 = np.hstack((newAssur1_y[1:len(newAssur1_y)], newAssur2_y[1:len(newAssur1_y)]));
hstack2 = np.hstack((newPharma1_y[1:len(newPharma1_y)], newPharma2_y[1:len(newPharma1_y)]));
hstack3 = np.hstack((newBank1_y[1:len(newBank1_y)], newBank2_y[1:len(newBank1_y)]));

g1 = grangercausalitytests(hstack1, maxlag = 2, verbose = True);
g2 = grangercausalitytests(hstack2, maxlag = 2, verbose = True);
g3 = grangercausalitytests(hstack3, maxlag = 2, verbose = True);

## Pearson correlation
gp1 = pearsonr(newAssur1_y, newAssur2_y);
gp2 = pearsonr(newPharma1_y, newPharma2_y);
gp3 = pearsonr(newBank1_y, newBank2_y);

## Cross-correlation
gb1 = pointbiserialr(newAssur1_y, newAssur2_y);
gb2 = pointbiserialr(newPharma1_y, newPharma2_y);
gb3 = pointbiserialr(newBank1_y, newBank2_y);

## Comparison with simple causality
np.correlate()