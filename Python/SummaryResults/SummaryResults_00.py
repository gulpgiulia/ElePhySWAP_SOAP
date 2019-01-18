#------------------------------------------------------
#  Copyright 2017 Giulia De Bonis @ INFN, Rome - Italy
#  Version: 1.0 - April 24, 2017
#-----------------------------------------------------

# IMPORT

import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
import scipy.io as sio 
from operator import itemgetter, attrgetter, methodcaller
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import math

#------------------------------------------------------------------------------------------
# PATH

BaseDir = '/apotto/home1/homedirs/gdebonis/mouseExpData/'
InputDir = BaseDir + 'ISS-DataAnalysis/Results/Nov16/Spontaneous/161101_rec07_Spontaneous_RH/'
OutputDir = BaseDir + 'Python-DataAnalysis/Results/Spontaneous/Nov16/161101_rec07_Spontaneous_RH/'

#------------------------------------------------------------------------------------------
debug=False
#------------------------------------------------------------------------------------------
if debug:
# LOAD  Results.mat
    data = sio.loadmat(InputDir + 'Results.mat')
    datastruct = data['results']
    datastruct.shape
    nCh=datastruct.shape[1]

    ### (datastruct[0,i]['UpState']['numel_duration']) --> number of transitions
    ### (datastruct[0,i]['DownState']['mean_duration']) --> mean suration of the DownState

    # DownState -- MEAN Duration
    DownState = []
    for i in range(0,nCh):
        DownState.append(datastruct[0,i]['DownState']['mean_duration'])

    DownState = np.concatenate(DownState).astype(None)
    Ch = np.arange(1, nCh+1, 1)

    plt.plot(Ch,DownState,'bo')
    #plt.errorbar(Ch, DownState, yerr=.1, fmt=None, color='b')
    plt.xlabel('Channel')
    plt.ylabel('DownState [s]')
    plt.show()

    # BoxPlot
    plt.figure()
    plt.boxplot(DownState)
    plt.show()

#------------------------------------------------------------------------------------------
# LOAD Duration.mat

data = sio.loadmat(InputDir + 'Duration.mat')
datastruct = data['duration']
datastruct.shape
nCh=datastruct.shape[1]

# Cortical Areas
#area = {'M':[0,1,2,5,8],'S':[3,4,6,7,9,10,11],'P':[13,14,15],'R':[12,16,20,25],
#        'V':[17,18,19,21,22,23,24,26,27,28,29,30,31]} # (dictionary)

area = [['M',[0,1,2,5,8]],['S',[3,4,6,7,9,10,11]],['P',[13,14,15]],['R',[12,16,20,25]],
        ['V',[17,18,19,21,22,23,24,26,27,28,29,30,31]]] # (list)
Ch = np.arange(1, nCh+1, 1)

if debug: # (Channels are NOT grouped by Cortical Area)
#--- DownStateLen ---
    DownStateLen = []
    for i in range(0,nCh):
        DownStateLen.append([datastruct[0,i]['DownStateLen']])
   
    # BoxPlot_Channels
    plt.figure(1,figsize=(12, 6))
    plt.boxplot(DownStateLen,sym='+')
    plt.xlabel('Channel')
    plt.ylabel('Duration of the Down State [s]')
    plt.show(0)
    plt.savefig(OutputDir + 'DownState_BoxPlot_Channels.pdf')

# *** BoxPlot_Channels *** (Channels ARE grouped by Cortical Area)

DownStateLen = []
meanCh = np.empty((0,nCh))
ChList = np.empty((0,nCh))
Separator=np.array([0.])
for i in area:
    Separator=np.append(Separator,Separator[-1]+len(i[1]))
    for j in i[1]:
        DownStateLen.append([datastruct[0,j]['DownStateLen']])
        meanCh=np.append(meanCh,np.mean([datastruct[0,j]['DownStateLen']]))
        ChList=np.append(ChList,str(j+1))  # The labels for the bottom x axis

separator=Separator[1:]

fig, ax1 = plt.subplots(figsize=(12, 6))
plt.boxplot(DownStateLen,sym='+')
ax1.set_xticklabels(ChList)
ax1.set_xlabel('Channel') 

plt.plot(ax1.get_xticks(),meanCh,linestyle='none',color='b',marker='*')

#width=[]
height=ax1.get_ylim()[1]-ax1.get_ylim()[0]
TickPos = []
offset=0.5
for i in range(0,len(area)):
    TickPos.append(Separator[i] + (Separator[i+1]-Separator[i])/2)
#    width.append(Separator[i+1]-Separator[i]-0.5)
    width=(Separator[i+1]-Separator[i])
    if i%2 == 0:
        Alpha=0.2
    else:
        Alpha=0.4
    ax1.add_patch(
        patches.Rectangle(
            (0.+offset,ax1.get_ylim()[0]),   # (x,y)
            width,                    # width
            height,                   # height
            alpha=Alpha
        )
    )
    offset=offset+width

plt.ylabel('Duration of the Down State [s]')


labels = ['Motor','Somatosensory','Parietal','Retrosplenial','Visual'] # The labels for the top x axis

TickPosNorm=[i/nCh for i in TickPos]

ax2 = ax1.twiny()
ax2.set_xticks(TickPosNorm)
ax2.set_xticklabels(labels)
#ax2.set_xlabel('Cortical Area')


# Legend
plt.figtext(0.80, 0.92, '*', color='b', weight='roman', size='medium')
plt.figtext(0.815, 0.92, 'Mean Value', color='black', weight='roman', size='medium',
            bbox=dict(facecolor='none', edgecolor='k', boxstyle='round'))

plt.show(0)
plt.savefig(OutputDir + 'DownState_BoxPlot_Areas-and-Channels.pdf')



# *** BoxPlot_CorticalAreas ***

data_to_plot = []
meanA=[]
stdA=[]
stderrA=[]
#for i in sorted(area.keys(),key=itemgetter(0)): #(areas sorted from front to back)
    #print(i)
for i in area:
    A=[]
    for j in (i[1]):
        A=np.append(A,[datastruct[0,j]['DownStateLen']])
        #print(len(A))
    data_to_plot.append(A)
    meanA.append(np.mean(A))
    stdA.append(np.std(A))
    stderrA.append(np.std(A)/math.sqrt(len(A)))
    #print(len(A))
#print(meanA)
#print(stdA)
#print(stderrA)

fig, ax1 = plt.subplots(figsize=(10, 6))
fig.canvas.set_window_title('DownStateDuration - BoxPlot - CorticalAreas')
plt.subplots_adjust(left=0.075, right=0.95, top=0.9, bottom=0.1)

plt.boxplot(data_to_plot,sym='+', vert=1, whis=1.5)
plt.xlabel('Cortical Area')
plt.xticks([1,2,3,4,5], ['M','S','P','R','V'])
plt.xticks([1,2,3,4,5], ['Motor','Somatosensory','Parietal','Retrosplenial','Visual'])
xtickNames = plt.setp(ax1, xticklabels=['Motor','Somatosensory','Parietal','Retrosplenial','Visual'])
#plt.setp(xtickNames,rotation=45,fontsize=8)
plt.setp(xtickNames,fontsize=8)

plt.ylabel('Duration of the Down State [s]')

ml = MultipleLocator(0.2)
ax1.yaxis.set_minor_locator(ml)
ax1.yaxis.grid(True, linestyle='-', which='minor', color='lightgrey', alpha=0.5)
plt.plot([1,2,3,4,5],meanA,linestyle='none',color='b',marker='*')
#plt.errorbar([1,2,3,4,5],meanA,yerr=stderrA,fmt='o',color='b',marker='*') #stderr SMALL (large sample)
#plt.errorbar([1,2,3,4,5],meanA,yerr=stdA,fmt='o',color='b',marker='*') 
# overplot the sample averages (markeredgecolor='k')

# Legend
plt.figtext(0.80, 0.92, '*', color='b', weight='roman', size='medium')
plt.figtext(0.815, 0.92, 'Mean Value', color='black', weight='roman', size='medium',
            bbox=dict(facecolor='none', edgecolor='k', boxstyle='round'))

plt.show(0)
plt.savefig(OutputDir + 'DownState_BoxPlot_Areas.pdf')
