#!/usr/bin/python
#--------1---------2---------3---------4---------5---------6---------7--------X
#  Copyright 2017 Giulia De Bonis @ INFN, Rome - Italy
#  Version: 1.0 - June 6, 2017
#--------1---------2---------3---------4---------5---------6---------7--------X

#*** DEBUG Flux Control ***
#    if debug:
#        continue # new item in the for loop
#    else:
#        pass #continue
#
#if debug:
#     sys.exit()
#**************************

#-----------------------------------------------------
# IMPORT
#-----------------------------------------------------
execfile("importPackages.py")
print("...importPackages...")
#-----------------------------------------------------
# FUNCTIONS
#-----------------------------------------------------
execfile("defFunctions.py") 
#import test_defFunctions as my
print("...defineFunctions...")
#[NB] currently the module is *executed* and *not imported* in the main
#-----------------------------------------------------
#--- Channels
nCh = 32
Ch = np.arange(1, nCh+1, 1)

#--- N.B. Constants should be moved in importPackage or in defFunctions
#-----------------------------------------------------
# DEBUG
#-----------------------------------------------------
debug = True
ddebug = False #(deep-debug)
#-----------------------------------------------------
BaseDir = '/apotto/home1/homedirs/gdebonis/mouseExpData/'
OutputDir = BaseDir + 'SWAP/RESULTS/Nov2016/Spontaneous/comboRes/'
#-----------------------------------------------------
# load WT data
#-----------------------------------------------------
PythonDir = 'SWAP/RESULTS/Nov2016/Spontaneous/WT/Python/'
Dir = BaseDir + PythonDir
inputF = Dir + 'SD/SD.npy'
F = open(inputF,'r')
F.seek(0)
npzfile = np.load(F)
wtSD=npzfile['SD']
wtSDr=npzfile['SDr']
wtSD0=npzfile['SD0']
wtOUT=npzfile['OUT']
wtOUTsd=npzfile['OUTsd']
F.close()
#-----------------------------------------------------
# load KO data
#-----------------------------------------------------
PythonDir = 'SWAP/RESULTS/Nov2016/Spontaneous/KO/Python/'
Dir = BaseDir + PythonDir
inputF = Dir + 'SD/SD.npy'
F = open(inputF,'r')
F.seek(0)
npzfile = np.load(F)
koSD=npzfile['SD']
koSDr=npzfile['SDr']
koSD0=npzfile['SD0']
koOUT=npzfile['OUT']
koOUTsd=npzfile['OUTsd']
F.close()

#-----------------------------------------------------
# Figure (1)
#-----------------------------------------------------
fig1, ax1 = plt.subplots(121,figsize=(10, 6))  
plt.suptitle('[log(MUA)] Standard Deviation (SD) of the Down State')
#-----------------------------------------------------
plt.subplot(121)
plt.hist(MaskNaN(wtSD0),color='orange')
plt.hist(MaskNaN(koSD0),color='green')
plt.title('WT and KO mice (outliers exlcuded)')
plt.xlabel('SD')
plt.ylabel('#entries')
plt.figtext(0.4,0.85,'WT mice',color='orange',weight='roman',size='medium',
            bbox=dict(facecolor='white', edgecolor='black', boxstyle='round'))
plt.figtext(0.4, 0.815,'KO mice',color='green',weight='roman',size='medium',
            bbox=dict(facecolor='white', edgecolor='black', boxstyle='round'))
#-----------------------------------------------------
plt.subplot(122)
SDtot=np.append(wtSD0,koSD0)
N=len(MaskNaN(SDtot));
m=np.mean(MaskNaN(SDtot));
std=np.std(MaskNaN(SDtot));
M=np.median(MaskNaN(SDtot));
Q1=np.percentile(MaskNaN(SDtot),25)
Q3=np.percentile(MaskNaN(SDtot),75)
IQR=Q3-Q1
sigma3n=m-3*std
sigma2n=m-2*std
sigma2p=m+2*std
sigma3p=m+3*std

plt.hist(MaskNaN(SDtot))
plt.axvline(sigma3n,color='#8B0000')
plt.axvline(sigma2n,color='r')
plt.axvline(sigma2p,color='r')
plt.axvline(sigma3p,color='#8B0000')

plt.title('Total sample')
plt.xlabel('SD')
plt.ylabel('#entries')

plt.figtext(0.85, 0.85, 
                'N=%d\nmean=%f\nstd=%f\nmedian=%f\nQ1=%f\nQ3=%f'%(N,m,std,M,Q1,Q3),
                fontsize=8, 
                bbox=dict(facecolor='white', edgecolor='black', boxstyle='round'))
plt.figtext(0.95,0.855,'3*std',color='#8B0000',weight='roman',size='medium',
            bbox=dict(facecolor='white', edgecolor='black', boxstyle='round'))
plt.figtext(0.95, 0.82,'2*std',color='r',weight='roman',size='medium',
            bbox=dict(facecolor='white', edgecolor='black', boxstyle='round'))
#-----------------------------------------------------
plt.show(0)
plt.savefig(OutputDir + 'SDhisto.pdf')

#-----------------------------------------------------
# OUTPUT file
#-----------------------------------------------------
outfile = OutputDir + 'SDtotal.npy'
out = open(outfile,'w')
np.save(out,SDtot)
out.close()

################################################################
#-----------------------------------------------------
# How to load/read the output file
#-----------------------------------------------------
###execfile("importPackages.py")
#import numpy as np
#BaseDir = '/apotto/home1/homedirs/gdebonis/mouseExpData/'
#Dir = BaseDir + 'SWAP/RESULTS/Nov2016/Spontaneous/comboRes/'
#SDtot= np.load(Dir + 'SDtotal.npy')
################################################################

