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
#-----------------------------------------------------

#-----------------------------------------------------
# Geometry of the Electrode Array
#-----------------------------------------------------
#--- Cortical Areas (list)
area = [['M',[0,1,2,5,8]],['S',[3,4,6,7,9,10,11]],['P',[13,14,15]],
        ['R',[12,16,20,25]],['V',[17,18,19,21,22,23,24,26,27,28,29,30,31]]]

Separator=np.array([0.]) # (AreaSeparator, for visualization)
ChList = np.array([])    # (string) --- channel numbering
for i in area:
    Separator=np.append(Separator,Separator[-1]+len(i[1]))
    ### (for data representation in BoxPlotChannels)
    for j in i[1]:
        ChList=np.append(ChList,str(j+1))
        ### (labels for channels grouped by cortical areas)
    
labels = ['Motor','Somatosensory','Parietal','Retrosplenial','Visual']


#--- coordinates in the grid
ArrayMap = [(0,0),(1,0),(0,1),(0,2),(0,3),(1,1),(2,1),(1,2),(2,2),(1,3),(2,3),
            (3,3),(1,4),(2,4),(3,4),(0,4),(0,5),(0,6),(0,7),(1,5),(2,5),(3,5),
            (1,6),(2,6),(3,6),(4,6),(1,7),(2,7),(3,7),(1,8),(2,8),(3,8)] 

#--- Channels
nCh = 32
Ch = np.arange(1, nCh+1, 1)
#-----------------------------------------------------

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
# INITIALIZATION (stacking of Experiments)
#-----------------------------------------------------
###SD=np.empty((nExp,nCh)); SD.fill(np.nan); OUT=[None]*nExp;
#-------------
# per CHANNELS
#-------------
#--- median per channel, state duration [s]
Ch_Mdown=np.empty((0,32)) 
#[NB]Ch_Mdown=np.array([]) cannot be used, because the np.vstack requires
#--- that all the input array dimensions must match exactly
#[DEBUG]print Ch_Mdown,Ch_Mdown.shape,len(Ch_Mdown)
Ch_Mup=np.empty((0,32)) 
Ch_Mud=np.empty((0,32))

#--- mean per channel, frequency [Hz] (=1/<LenUD>)
Freq = np.empty((0,32))
errFreq = np.empty((0,32))

#---  slope (of the mean transition), per channel (dynamics of transition)
Ch_slopeD=np.empty((0,32)) # downward
Ch_slopeU=np.empty((0,32)) # upward

#--- peak  in the UP state, per channel
Ch_peak=np.empty((0,32))

#---------
# per AREA
#---------
#--- median values
M_down = np.empty((0,5)) # (5 is the number of Cortical Areas)
M_up = np.empty((0,5))
M_ud = np.empty((0,5))
M_Sd = np.empty((0,5))
M_Su = np.empty((0,5))
M_peak = np.empty((0,5))
m_down = np.empty((0,5))
#--- mean values + SEM + SD (...not used)
m_up = np.empty((0,5))
m_ud = np.empty((0,5))
m_Sd = np.empty((0,5))
m_Su = np.empty((0,5))
m_peak = np.empty((0,5))
SEM_down = np.empty((0,5))
SEM_up = np.empty((0,5))
SEM_ud = np.empty((0,5))
SEM_Sd = np.empty((0,5))
SEM_Su = np.empty((0,5))
SEM_peak = np.empty((0,5))
SD_down = np.empty((0,5))
SD_up = np.empty((0,5))
SD_ud = np.empty((0,5))
SD_Sd = np.empty((0,5))
SD_Su = np.empty((0,5))
SD_peak = np.empty((0,5))

# mean frequency per area
FreqA = np.empty((0,5))
errFreqA = np.empty((0,5))

#-----------------------------------------------------
# DEBUG
#-----------------------------------------------------
debug = True
ddebug = False #(deep-debug)
testFile = True

#-----------------------------------------------------
# ExpList
#-----------------------------------------------------
if testFile:
#--- subset of experiments with Frequency > 0.5 Hz (Apr2018)
    ExpList=['161101_rec07_Spontaneous_RH','161103_rec08_Spontaneous_LH',
             '161107_rec07_Spontaneous_LH','161114_rec10_Spontaneous_RH',
             '161115_rec10_Spontaneous_RH','161116_rec08_Spontaneous_RH',
             '161117_rec07_Spontaneous_RH','161120_rec07_Spontaneous_RH']
    ExpList=['161101_rec07_Spontaneous_RH']
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
# OULIERS
#-----------------------------------------------------
OUTall=True
outTOT=[None]*nExp;
infoOUT=BaseDir+PythonDir+'SD/SD.npz'
fff=np.load(infoOUT); out1=fff['OUT']; out2=fff['OUTsd']
for i in range(nExp):
    outTOT[i]=sorted(np.concatenate((out1[i], out2[i]), axis=None))
    outTOT[i]=np.asarray(outTOT[i])
# WARNING ! ! ! if testFile option is used

#-----------------------------------------------------
# MAIN LOOP over the Experiments (in ExpList)
#-----------------------------------------------------
#--- Produce PLOTs for each Variable and Experiment
#--- Store information in SummaryData, to be used for SummaryPlots

FileList=[]
n=0
for FileName in ExpList:
    n=n+1
    if n <= nExp:

        print(FileName)
        FileList=np.append(FileList,FileName)
        InputDir = BaseDir + MatlabDir + FileName +'/'
        OutputDir = BaseDir + PythonDir + FileName +'/'
        if not os.path.isdir(OutputDir):
            os.makedirs(OutputDir) # create OutputDir if not existing
# ----------------
# | Duration.mat | (DownStateLen, UpStateLen, UDcycleLen)
# ----------------
        data = sio.loadmat(InputDir + 'Duration.mat') #(scipy.io)
        print data.keys()
 
#>>> OUTLIERS
        outliers = data['OUT']
        if OUTall:
            outliers = outTOT[n-1]
        if outliers.size != 0:
        #[OR]if outliers.any() != 0: #(for Matlab v1_Dec2017)
            print 'OUTLIER CHANNELS: ', outliers

#>>> DATA
        datastruct = data['duration'] # (each SummaryFile.mat has its own structure)
        print datastruct.dtype #(dtype = data type)
        if datastruct.shape[1] != nCh:
            print 'WARNING --- Check your data! some channels are missing: nCh =', \
                                                                 datastruct.shape[1]

            
###################################################################################
# OPTIMIZE the CODE           
###################################################################################
#        DownState=[None]*nCh; UpState=[None]*nCh; UDcycle=[None]*nCh         
#        for i in range(nCh):    
#            DownState[i]=np.squeeze(datastruct['DownStateLen'][0,i])
#            UpState[i]=np.squeeze(datastruct['UpStateLen'][0,i])
#            UDcycle[i]=np.squeeze(datastruct['UDcycleLen'][0,i])
#        ### Number of entries for each channel: var[i].size OR len(var[i])
#
#        variable='DownStateLen'
#        yLabel='Duration of the Down State [s]'
#        (A,B,C)=ArrangeChannelsByArea(DownState,outliers,area)
#        BoxPlotChByArea(A,yLabel,ChList,Separator,FileName)
#        plt.savefig(OutputDir + variable +'_BoxPlot_ChannelsByAreas.pdf')
#        Ch_Mdown=np.vstack([Ch_Mdown,B]) # channels arranged by cortical areas
###################################################################################

        
#--- 1) DOWNSTATE
#[OR]outcome = my.BoxPlotCh('DownStateLen','Duration of the Down State [s]',debug) 
        #(as module)
        outcome = BoxPlotCh('DownStateLen','Duration of the Down State [s]',debug)
        Ch_Mdown = np.vstack([Ch_Mdown,outcome[0]]) # outcome[0] is the median for each channel
        ###Contour(outcome[0],'DownStateLen')
        (M,m,SD,SEM) = BoxPlotArea('DownStateLen','Duration of the Down State [s]')
        M_down = np.vstack([M_down, M])
        m_down = np.vstack([m_down, m])
        SD_down = np.vstack([SD_down, SD])
        SEM_down = np.vstack([SEM_down, SEM])
 
#--- 2) UPSTATE
        outcome = BoxPlotCh('UpStateLen','Duration of the Up State [s]',debug)
        Ch_Mup = np.vstack([Ch_Mup,outcome[0]]) # outcome[0] is the median for each channel
        ###Contour(outcome[0],'UpStateLen')
        (M,m,SD,SEM) = BoxPlotArea('UpStateLen','Duration of the Up State [s]')
        M_up = np.vstack([M_up, M])
        m_up = np.vstack([m_up, m])
        SD_up = np.vstack([SD_up, SD])
        SEM_up = np.vstack([SEM_up, SEM])

#--- 3) UD-CYCLE 
        outcome = BoxPlotCh('UDcycleLen','Duration of the Up-Down cycle [s]',debug) 
        Ch_Mud = np.vstack([Ch_Mud,outcome[0]]) # outcome[0] is the median for each channel
        ###Contour(outcome[0],'UDStateLen')
        (M,m,SD,SEM) = BoxPlotArea('UDcycleLen','Duration of the Up-Down cycle [s]')
        M_ud = np.vstack([M_ud, M])
        m_ud = np.vstack([m_ud, m])
        SD_ud = np.vstack([SD_ud, SD])
        SEM_ud = np.vstack([SEM_ud, SEM])
       
#--- 4) FREQUENCY
        F = 1./outcome[2] #(meanCh, array of nCh elements)
        errF = F*F*outcome[3] #(semCh)
        # --- PLOT mean frequency per channel
        StemPlotCh(F,'F','Frequency [Hz]')
        PlotCh(F,errF,'F','Frequency [Hz]') #(maybe not so useful..)
        ###Contour(F,'Frequency')
        # --- PLOT mean frequency per area
        f=1./m
        errf = f*f*SEM
        StemPlotArea(f,'F','Frequency [Hz]')
        PlotArea(f,errf,'F','Frequency [Hz]')
        # --- store data (stacking)
        Freq = np.vstack([Freq,F]) #(array of nExp el, each el is an array of nCh elements)
        errFreq = np.vstack([errFreq,errF]) 
        FreqA = np.vstack([FreqA, f]) # SummaryData (experiments, areas)
        errFreqA = np.vstack([errFreqA,errf])

        # *** Results ***

        # *** DownwardTrans *** --> SLOPE
        data = sio.loadmat(InputDir + 'DownwardTrans.mat')
        print data.keys()
        datastruct = data['downtrans'] # (each SummaryFile.mat has its own structure)
        print datastruct.dtype
        if datastruct.shape[1] != nCh:
            print 'WARNING --- Check your data! some channels are missing: nCh =', \
                                                                 datastruct.shape[1] 
### DO NOT REDEFINE OUTLIERS!
##        outliers = data['OUT']
##        if outliers.size != 0:
##            print 'OUTLIER CHANNELS: ', outliers

#--- 1) DOWNWARD transition SLOPE
        (slope,peak)=GetSlope([-0.025,0.010]) #(peak not used for downward transitions)
        StemPlotCh(slope,'DownwardSlope','Up to Down state transition [s-1]')
        ###Contour(slope,'DownwardSlope') 
        Ch_slopeD=np.vstack([Ch_slopeD,slope])

        (M,m,SD,SEM)=GetAreaInfo(slope)
        M_Sd = np.vstack([M_Sd, M])
        m_Sd = np.vstack([m_Sd, m])
        SD_Sd = np.vstack([SD_Sd, SD])
        SEM_Sd = np.vstack([SEM_Sd, SEM])

        # *** UpwardTrans *** --> SLOPE
        data = sio.loadmat(InputDir + 'UpwardTrans.mat')
        print data.keys()
        datastruct = data['uptrans'] # (each SummaryFile.mat has its own structure)
        print datastruct.dtype
        if datastruct.shape[1] != nCh:
            print 'WARNING --- Check your data! some channels are missing: nCh =', \
                                                                 datastruct.shape[1] 


### DO NOT REDEFINE OUTLIERS!
##        outliers = data['OUT']
##        if outliers.size != 0:
##            print 'OUTLIER CHANNELS: ', outliers

#--- 1)2) UPWARD transitionSLOPE and Relative Firing Rate Max
        (slope,peak)=GetSlope([-0.010,0.025])
        StemPlotCh(slope,'UpwardSlope','Down to Up state transition [s-1]')
        ###Contour(slope,'UpwardSlope') 
        Ch_slopeU=np.vstack([Ch_slopeU,slope])
        StemPlotCh(peak,'PeakY','Relative Firing Rate Max')                  
        ###Contour(peak,'PeakY') 
        Ch_peak=np.vstack([Ch_peak,peak])

        (M,m,SD,SEM)=GetAreaInfo(slope)
        M_Su = np.vstack([M_Su, M])
        m_Su = np.vstack([m_Su, m])
        SD_Su = np.vstack([SD_Su, SD])
        SEM_Su = np.vstack([SEM_Su, SEM])

        (M,m,SD,SEM)=GetAreaInfo(peak)
        M_peak = np.vstack([M_peak, M])
        m_peak = np.vstack([m_peak, m])
        SD_peak = np.vstack([SD_peak, SD])
        SEM_peak = np.vstack([SEM_peak, SEM])

        plt.close('all')

#------------- END of LOOP over the Exps -------------

#-----------------------------------------------------
# STORE DATA
#-----------------------------------------------------
A=np.array([M_down,m_down,SD_down,SEM_down])
B=np.array([M_up,m_up,SD_up,SEM_up])
C=np.array([M_ud,m_ud,SD_ud,SEM_ud])
D=np.array([M_Sd,m_Sd,SD_Sd,SEM_Sd])
E=np.array([M_Su,m_Su,SD_Su,SEM_Su])
F=np.array([M_peak,m_peak,SD_peak,SEM_peak])
G=np.array([FreqA,errFreqA])

# SummaryData: store the information of the nExp in a unique array    
SummaryDataArea=np.array([A,B,C,D,E,F,G])
SummaryDataCh=np.array([Ch_Mdown,Ch_Mup,Ch_Mud,Freq,Ch_slopeD,Ch_slopeU,Ch_peak])

#-----------------------------------------------------
# --- SAVE DATA to files
#-----------------------------------------------------
#outfile = BaseDir + PythonDir + 'SummaryData.npy'
#np.save(outfile, SummaryData)

OutputDir = BaseDir + PythonDir

outfile = OutputDir + 'SummaryDataArea.npy'
F = open(outfile,'w')
np.save(F, SummaryDataArea)

outfile = OutputDir + 'SummaryDataCh.npy'
F = open(outfile,'w')
np.save(F, SummaryDataCh)

outfile = OutputDir + 'Frequency.npy'
np.save(outfile,Freq)
np.save(outfile,errFreq) # save twice and load twice

outfile = OutputDir + 'Frequency.npz'
F = open(outfile,'w')
np.savez(F, Freq=Freq, errFreq=errFreq)
F.close()

outfile = OutputDir + 'FileList.txt'
F = open(outfile,'w')
for item in FileList:
  F.write("%s\n" % item)
F.close()

#np.load(outfile)

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# f = [0] * 30
# a = np.ones((3,2))        # a 2D array with 3 rows, 2 columns, filled with ones
# b = np.array([1,2,3])     # a 1D array initialised using a list [1,2,3]
#c = np.linspace(2,3,100)  # an array with 100 points beteen (and including) 2 and 3
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


