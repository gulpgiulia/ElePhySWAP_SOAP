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
SummaryFile = ['Duration', 'Results', 'DownwardTrans', 'UpwardTrans']

#--- Cortical Areas (list)
area = [['M',[0,1,2,5,8]],['S',[3,4,6,7,9,10,11]],['P',[13,14,15]],
        ['R',[12,16,20,25]],['V',[17,18,19,21,22,23,24,26,27,28,29,30,31]]]
#--- coordinates in the grid
ArrayMap = [(0,0),(1,0),(0,1),(0,2),(0,3),(1,1),(2,1),(1,2),(2,2),(1,3),(2,3),
            (3,3),(1,4),(2,4),(3,4),(0,4),(0,5),(0,6),(0,7),(1,5),(2,5),(3,5),
            (1,6),(2,6),(3,6),(4,6),(1,7),(2,7),(3,7),(1,8),(2,8),(3,8)] 

labels = ['Motor','Somatosensory','Parietal','Retrosplenial','Visual']

#--- Channels
nCh = 32
Ch = np.arange(1, nCh+1, 1)

#--- N.B. Constants should be moved in importPackage or in defFunctions

#-----------------------------------------------------
# PATH
#-----------------------------------------------------
BaseDir   = '/apotto/home1/homedirs/gdebonis/mouseExpData/'
MatlabDir = 'SWAP/RESULTS/Nov2016/Spontaneous/WT/Matlab/'
PythonDir = 'SWAP/RESULTS/Nov2016/Spontaneous/WT/Python/'

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
#--- Produce PLOTs for each Variable and Experiment 

FileList = []
if testFile:
#--- subset of experiments with Frequency > 0.5 Hz (Apr2018)
    ExpList=['161111_rec10_Spontaneous_LH']
else:
    ExpList=os.listdir(BaseDir+MatlabDir) 
# N.B. listdir lists elements in arbitrary order
    print(ExpList)

nFile = len(ExpList)

nExp=nFile # set the maximum number of experiments to process
if nFile<nExp:
     nExp=nFile
print 'nExp =',nExp

#-----------------------------------------------------
# INITIALIZATION (stacking of Experiments)
#-----------------------------------------------------
MuTot = np.empty((nExp,nCh)); MuTot.fill(np.nan)
SkewTot = np.empty((nExp,nCh)); SkewTot.fill(np.nan)
MuTail = np.empty((nExp,nCh)); MuTail.fill(np.nan)
SkewTail = np.empty((nExp,nCh)); SkewTail.fill(np.nan)
SD=np.empty((nExp,nCh)); SD.fill(np.nan)

#-----------------------------------------------------
# MAIN LOOP over the Experiments 
#-----------------------------------------------------
#--- Produce PLOTs for each Variable and Experiment 

for n in range(0,nExp):
    FileName=ExpList[n]
    print(FileName)
    FileList=np.append(FileList,FileName)
    InputDir = BaseDir + MatlabDir + FileName +'/'
    OutputDir = BaseDir + PythonDir + FileName +'/'
    if not os.path.isdir(OutputDir):
        os.makedirs(OutputDir)

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
                
#--- 1) ModeParams.SkewnessTot [Mu and Skew]
    for i in range(nCh):
        if i not in outliers-1:
            MuTot[n,i]=np.squeeze(datastruct['ModeParams'][0,i][0,0][5][0][0][0])
            SkewTot[n,i]=np.squeeze(datastruct['ModeParams'][0,i][0,0][5][0][0][7])
            MuTail[n,i]=np.squeeze(datastruct['ModeParams'][0,i][0,0][6][0][0][0])
            SkewTail[n,i]=np.squeeze(datastruct['ModeParams'][0,i][0,0][6][0][0][7])
    
    SD[n]=np.squeeze(datastruct['Delta'])
    for i in range(nCh):
        SD[n,i]=np.asscalar(SD[n,i])
        #if i in outliers-1:
        #    SD[n,i]=np.nan
    SD[n]=SD[n].astype(np.float)


SDr=SD.reshape(SD.shape[0]*SD.shape[1])
mmm=np.mean(MaskNaN(SDr))
sss=np.std(MaskNaN(SDr))

for i in range(len(SDr)):
    if ~np.isnan(SDr[i]):
        scartoQ=math.sqrt((SDr[i]-mmm)*(SDr[i]-mmm))
        if scartoQ>3*sss:
            print(i)




m=0
for n in range(nExp):
    m=m+np.mean(MaskNaN(SD[n]))
media=m/nExp


    SDmean=np.mean(MaskNaN(SD))
    SDmed=np.median(MaskNaN(SD))



    for i in range(nCh):
        if i not in outliers-1:
            scartoQ=math.sqrt((SD[i]-SDmean)*(SD[i]-SDmean))
            if scartoQ>3*STD:
                print(i)


    ###SD=np.delete(SD,outliers-1)
        
    [X,Y,Z]=Contour(MuTot[n],'MuTot')
    Contour(SkewTot[n],'SwekTot')
    Contour(MuTail[n],'MuTail')
    Contour(SkewTail[n],'SkewTail')
    plt.close('all')

FileName='messageA'
OutputDir = BaseDir + PythonDir
mMuTot = np.empty((nCh))
mSkewTot = np.empty((nCh))
mMuTail = np.empty((nCh))
mSkewTail = np.empty((nCh))
for i in range(0,nCh):
    mMuTot[i]=np.mean(MuTot[:,i])
    mSkewTot[i]=np.mean(SkewTot[:,i])
    mMuTail[i]=np.mean(MuTail[:,i])
    mSkewTail[i]=np.mean(SkewTail[:,i])
Contour(mMuTot,'meanMuTot')
Contour(mSkewTot,'meanSkewTot')
Contour(mMuTail,'meanMuTail')
Contour(mSkewTail,'meanSkewTail')
plt.close('all')
