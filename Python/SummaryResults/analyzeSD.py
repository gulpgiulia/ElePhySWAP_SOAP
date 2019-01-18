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
#--- Channels
nCh = 32
Ch = np.arange(1, nCh+1, 1)

#--- N.B. Constants should be moved in importPackage or in defFunctions

#-----------------------------------------------------
# PATH
#-----------------------------------------------------
###[ikkio]BaseDir   = '/apotto/home1/homedirs/gdebonis/mouseExpData/'
BaseDir   = '/Users/giuliadebonis/Desktop/'
MatlabDir = 'SWAP/RESULTS/Nov2016/Spontaneous/WT/Matlab/v2/'
PythonDir = 'SWAP/RESULTS/Nov2016/Spontaneous/WT/Python/v2/'

#-----------------------------------------------------
# FUNCTIONS
#-----------------------------------------------------
execfile("defFunctions.py") 
#import test_defFunctions as my
print("...defineFunctions...")
#[NB] currently the module is *executed* and *not imported* in the main

#-----------------------------------------------------
# DEBUG
#-----------------------------------------------------
debug = True
ddebug = False #(deep-debug)
testFile = False

#-----------------------------------------------------
# Select the FileList to be analysed
#-----------------------------------------------------
FileList = []
if testFile:
#--- subset of experiments with Frequency > 0.5 Hz (Apr2018)
    ExpList=['161111_rec10_Spontaneous_LH']
else:
    ExpList=os.listdir(BaseDir+MatlabDir) 
# N.B. listdir lists elements in arbitrary order
    unwanted = {'.DS_Store','Icon\r','logbook'}
    ExpList = [e for e in ExpList if e not in unwanted]
    
print(ExpList)
nFile = len(ExpList)

nExp=nFile # set the maximum number of experiments to process
if nFile<nExp:
     nExp=nFile
print 'nExp =',nExp

#-----------------------------------------------------
# INITIALIZATION (stacking of Experiments)
#-----------------------------------------------------
SD=np.empty((nExp,nCh)); SD.fill(np.nan); OUT=[None]*nExp;

#-----------------------------------------------------
# MAIN LOOP over the Experiments 
#-----------------------------------------------------
for n in range(0,nExp):
    FileName=ExpList[n]
    print(FileName)
    FileList=np.append(FileList,FileName)
    InputDir = BaseDir + MatlabDir + FileName +'/'
    OutputDir = BaseDir + PythonDir + FileName +'/'
    if not os.path.isdir(OutputDir):
        os.makedirs(OutputDir)
#-----------------------------------------------------
    # *** Results ***
    data = sio.loadmat(InputDir + 'Results.mat') #(scipy.io)
    print data.keys()
    datastruct = data['results'] # (each SummaryFile.mat has its own structure)
    print datastruct.dtype #(dtype = data type)
    
    if(ddebug):
        print '***'
        print 'type(data):',type(data)
        print 'type(datastruct):',type(datastruct)
        print '***'

    if datastruct.shape[1] != nCh:
        print 'WARNING --- Check your data! some channels are missing: nCh =', \
                                                             datastruct.shape[1] 
    outliers = data['OUT']
    if outliers.size != 0:
    #[OR]if outliers.any() != 0: #(for Matlab v1_Dec2017)
        print 'OUTLIER CHANNELS: ', outliers
#-----------------------------------------------------    
    OUT[n]=outliers
        
    SD[n]=np.squeeze(datastruct['Delta'])
    for i in range(nCh):
        SD[n,i]=np.asscalar(SD[n,i])
        if i in outliers-1:
            SD[n,i]=np.nan
    SD[n]=SD[n].astype(np.float)

#-----------------------------------------------------
# Figure (1)
#-----------------------------------------------------
fig1, ax1 = plt.subplots(figsize=(10, 6))  
for n in range(nExp):
    plt.stem(Ch, SD[n],markerfmt='o') 
    
plt.title('[log(MUA)] Standard Deviation (SD) of the Down State \
\n as $\sigma$ of the first peak (Gaussian fit)')
plt.figtext(0.025,0.025,'(for each experiment, outlier channels are already excluded)',\
            weight='roman',size='small')
plt.xlabel('Channel')
plt.ylabel('SD')
ax1.set_xticks(Ch)
plt.show(0)
plt.savefig(BaseDir+PythonDir + 'SD/SD_stem.pdf')

#-----------------------------------------------------
# Figure (2)
#-----------------------------------------------------
SDr=SD.reshape(SD.shape[0]*SD.shape[1])
N=len(MaskNaN(SDr));
m=np.mean(MaskNaN(SDr));
std=np.std(MaskNaN(SDr));
M=np.median(MaskNaN(SDr));
Q1=np.percentile(MaskNaN(SDr),25)
Q3=np.percentile(MaskNaN(SDr),75)
IQR=Q3-Q1
sigma3n=m-3*std
sigma2n=m-2*std
sigma2p=m+2*std
sigma3p=m+3*std
#-----------------------------------------------------
fig2, ax2 = plt.subplots(121,figsize=(10, 6))  
plt.suptitle('[log(MUA)] Standard Deviation (SD) of the Down State (full data sample)',y=0.9)
#plt.figtext(0.81, 0.67, 
#                'N=%d\nmean=%.3f\nstd=%.3f\nmedian=%.3f\nQ1=%.3f\nQ3=%.3f'%(N,m,std,M,Q1,Q3),
#                fontsize=8, 
#                bbox=dict(facecolor='white', edgecolor='black', boxstyle='round'))
plt.figtext(0.768, 0.615, 
                'N=%d\nmean=%.3f\nstd=%.3f\nmedian=%.3f\nQ1=%.3f\nQ3=%.3f'%(N,m,std,M,Q1,Q3),
                fontsize=12, 
                bbox=dict(facecolor='white', edgecolor='black', boxstyle='round'))
#plt.figtext(0.43,0.765,'3*std',color='#8B0000',weight='roman',size='medium',
#            bbox=dict(facecolor='white', edgecolor='black', boxstyle='round'))
#plt.figtext(0.43, 0.73,'2*std',color='r',weight='roman',size='medium',
#            bbox=dict(facecolor='white', edgecolor='black', boxstyle='round'))

plt.figtext(0.4,0.765,'3*std',color='#8B0000',weight='roman',size='large',
            bbox=dict(facecolor='white', edgecolor='black', boxstyle='round'))
plt.figtext(0.4, 0.725,'2*std',color='r',weight='roman',size='large',
            bbox=dict(facecolor='white', edgecolor='black', boxstyle='round'))
#-----------------------------------------------------
plt.subplot(121)
plt.hist(MaskNaN(SDr))
plt.axvline(sigma3n,color='#8B0000')
plt.axvline(sigma2n,color='r')
plt.axvline(sigma2p,color='r')
plt.axvline(sigma3p,color='#8B0000')
plt.title('Histogram')
plt.xlabel('SD')
plt.ylabel('#entries')
#-----------------------------------------------------
plt.subplot(122)
meanpoint = dict(marker='*',markerfacecolor='blue',markeredgecolor='blue')
plt.boxplot(MaskNaN(SDr),sym='+',showmeans=True,meanprops=meanpoint)
###[ScatterPlot]
###y=MaskNaN(SDr)
###x=np.random.normal(1, 0.04, size=len(y))
###plt.plot(x,y,'k.',alpha=0.2)
plt.title('BoxPlot')
plt.ylabel('SD')
#-----------------------------------------------------
plt.subplots_adjust( wspace=0.2, top=0.8)
plt.show(0)
plt.savefig(BaseDir+PythonDir + 'SD/SD_fullsample.pdf')

#-----------------------------------------------------
# Figure (3)
#-----------------------------------------------------
fig3=plt.subplots(figsize=(10, 6)) 
plt.suptitle('SD distribution after the exclusion of SD outlier channels')
gs = gridspec.GridSpec(1,2,width_ratios=[3,1])  
#-----------------------------------------------------
ax31=plt.subplot(gs[0])

TickLabels = [0 for x in range(nExp)]
for i in range(0,nExp):
    TickLabels[i] = FileList[i][0:7]
    TickLabels[i] = TickLabels[i] + FileList[i][25:27]

for n in range(nExp):
    if OUT[n].size !=0:
        ax31.plot(n,OUT[n],linestyle='None',marker='o',color='red',markersize=10)

moreOUT=np.array([],dtype=np.uint8)
OUTsd=[[] for n in range(nExp)]
for i in range(len(SDr)):
    if ~np.isnan(SDr[i]) and (SDr[i]>Q3+1.5*IQR or SDr[i]<Q1-1.5*IQR):
        exp=i/nCh
        chID=i%nCh
        OUTsd[exp]=np.append(OUTsd[exp],chID+1)
        print i,exp,chID
        ax31.plot(exp,chID+1,linestyle='None',marker='o',color='b',markersize=10)
        moreOUT=np.append(moreOUT,i)

ax31.set_yticks(Ch)
ax31.set_xticks(range(nExp))
plt.grid(True)
xtickNames = plt.setp(ax31,xticklabels=TickLabels)
plt.setp(xtickNames,fontsize=6)
plt.title('Outlier Channels')
plt.xlabel('Experiment')
plt.ylabel('Outlier Channels')

for n in range(nExp):
    OUTsd[n]=np.array(OUTsd[n], dtype=np.uint8)
    OUT[n]=np.squeeze(OUT[n])

plt.figtext(0.025,0.95,'Experiment outliers',color='r',weight='roman',size='medium')
plt.figtext(0.025, 0.92,'Full data sample SD outliers',color='blue',weight='roman',size='medium')
#-----------------------------------------------------
ax32=plt.subplot(gs[1])

SD0=SDr
SD0[moreOUT[:]]=np.nan

N0=len(MaskNaN(SD0));
m0=np.mean(MaskNaN(SD0));
std0=np.std(MaskNaN(SD0));
M0=np.median(MaskNaN(SD0));
Q10=np.percentile(MaskNaN(SD0),25)
Q30=np.percentile(MaskNaN(SD0),75)
IQR0=Q3-Q1

sigma3n0=m0-3*std0
sigma2n0=m0-2*std0
sigma2p0=m0+2*std0
sigma3p0=m0+3*std0

ax32.hist(MaskNaN(SD0))
plt.axvline(sigma3n0,color='#8B0000')
plt.axvline(sigma2n0,color='r')
plt.axvline(sigma2p0,color='r')
plt.axvline(sigma3p0,color='#8B0000')

plt.title('Histogram')
plt.xlabel('SD')
plt.ylabel('#entries')
plt.figtext(0.84, 0.65, 
                'N=%d\nmean=%.3f\nstd=%.3f\nmedian=%.3f\nQ1=%.3f\nQ3=%.3f'%(N0,m0,std0,M0,Q10,Q30),
                fontsize=8, 
                bbox=dict(facecolor='white', edgecolor='black', boxstyle='round'))

plt.figtext(0.84,0.85,'3*std',color='#8B0000',weight='roman',size='medium',
            bbox=dict(facecolor='white', edgecolor='black', boxstyle='round'))
plt.figtext(0.84,0.81,'2*std',color='r',weight='roman',size='medium',
            bbox=dict(facecolor='white', edgecolor='black', boxstyle='round'))
#-----------------------------------------------------
plt.show(0)
plt.savefig(BaseDir+PythonDir + 'SD/outliers.pdf')

#-----------------------------------------------------
# OUTPUT file
#-----------------------------------------------------
OutputDir = BaseDir + PythonDir
outfile = OutputDir + 'SD/SD.npz'
out = open(outfile,'w')
OUTm=np.asarray(OUT,object)
np.savez(out, SD=SD, SDr=SDr, SD0=SD0, OUT=OUTm, OUTsd=OUTsd)
out.close()

outfile = OutputDir + 'SD/FileList.txt'
F = open(outfile,'w')
for item in FileList:
  F.write("%s\n" % item)
F.close()

################################################################
#-----------------------------------------------------
# How to load/read the output file
#-----------------------------------------------------
###execfile("importPackages.py")
#import numpy as np
#BaseDir   = '/apotto/home1/homedirs/gdebonis/mouseExpData/'
#PythonDir = 'SWAP/RESULTS/Nov2016/Spontaneous/WT/Python/'
#Dir = BaseDir + PythonDir
#inputF = Dir + 'SD/SD.npy'
#F = open(inputF,'r')
#F.seek(0)
#npzfile = np.load(F)
#SD=npzfile['SD']
#SDr=npzfile['SDr']
#SD0=npzfile['SD0']
#OUT=npzfile['OUT']
#OUTsd=npzfile['OUTsd']
#F.close()
################################################################
