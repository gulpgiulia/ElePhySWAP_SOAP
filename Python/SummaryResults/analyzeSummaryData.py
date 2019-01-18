#!/usr/bin/python
#--------1---------2---------3---------4---------5---------6---------7--------X
#  Copyright 2017 Giulia De Bonis @ INFN, Rome - Italy
#  Version: 1.0 - July 10, 2017
#--------1---------2---------3---------4---------5---------6---------7--------X

# ------------------------------------------------------
# IMPORT
# ------------------------------------------------------
execfile("importPackages.py")
print("...importPackages...")

#-----------------------------------------------------
# PATH
# ------------------------------------------------------
###[ikkio]BaseDir   = '/apotto/home1/homedirs/gdebonis/mouseExpData/'
BaseDir   = '/Users/giuliadebonis/Desktop/'
MatlabDir = 'SWAP/RESULTS/Nov2016/Spontaneous/WT/Matlab/v2/'
PythonDir = 'SWAP/RESULTS/Nov2016/Spontaneous/WT/Python/v2/'

#-----------------------------------------------------
# FUNCTIONS
#-----------------------------------------------------
execfile("defFunctions.py")
print("...defineFunctions...")

#-----------------------------------------------------
nCh=32

#--- Cortical Areas (list)
area = [['M',[0,1,2,5,8]],['S',[3,4,6,7,9,10,11]],['P',[13,14,15]],
        ['R',[12,16,20,25]],['V',[17,18,19,21,22,23,24,26,27,28,29,30,31]]]
#--- coordinates in the grid
ArrayMap = [(0,0),(1,0),(0,1),(0,2),(0,3),(1,1),(2,1),(1,2),(2,2),(1,3),(2,3),
            (3,3),(1,4),(2,4),(3,4),(0,4),(0,5),(0,6),(0,7),(1,5),(2,5),(3,5),
            (1,6),(2,6),(3,6),(4,6),(1,7),(2,7),(3,7),(1,8),(2,8),(3,8)] 

ChList = [1,2,3,6,9,4,5,7,8,10,11,12,14,15,16,13,
          17,21,26,18,19,20,22,23,24,25,27,28,29,30,31,32]
ChListStr = [str(i) for i in ChList]

#-----------------------------------------------------
# Step-by-Step analysis
#-----------------------------------------------------
FreqPlot = False
BoxPlot = False
ContourPlot = True
StatTest = False #(Hypothesis Tests are done in Matlab)

debug = False
# ------------------------------------------------------
# *** LOAD DATA ***
# ------------------------------------------------------
Dir = BaseDir + PythonDir

#--- SummaryData
SummaryDataArea = np.load(Dir + 'SummaryDataArea.npy')
SummaryDataCh = np.load(Dir + 'SummaryDataCh.npy') #(used for ContourPlot)

#--- FileList
infile = Dir + 'FileList.txt'
with open(infile) as F:
    FileList = F.read().splitlines() 
    print FileList
F.close()
nExp=len(FileList)
print 'nExp =',nExp

#--- FileList
infile = Dir + 'Frequency.npz'
F = open(infile,'r')
F.seek(0)
npzfile = np.load(F)
Freq=npzfile['Freq']
errFreq=npzfile['errFreq']
F.close()

#--- SD and OUTLIERS
infile = Dir + 'SD/SD.npz'
F = open(infile,'r')
F.seek(0)
npzfile = np.load(F)
wtSD=npzfile['SD']
wtSDr=npzfile['SDr']
wtSD0=npzfile['SD0']
wtOUT=npzfile['OUT']
wtOUTsd=npzfile['OUTsd']
F.close()

# ------------------------------------------------------
# *** MeanFrequency ***
# ------------------------------------------------------

if(FreqPlot):

    #--- MEAN FREQUENCY for each Experiment

    meanF=[np.nan]*nExp
    stdF=[np.nan]*nExp
    semF=[np.nan]*nExp
    for j in range(0,nExp):
        S=0;STD=0;n=0
        for i in range(0,nCh): # meanF = mean of the meanCh
            if ~np.isnan(Freq[j,i]):
                S=S+Freq[j,i]; n=n+1
                STD=STD + Freq[j,i]* Freq[j,i];
        meanF[j]=S/n
        stdF[j]=math.sqrt(STD/n)-meanF[j]*meanF[j]
        semF[j]=stdF[j]/math.sqrt(n)

    mm=min(meanF)
    MM=max(meanF)
    m=np.mean(meanF)
    SD=np.std(meanF)
    SEM=stats.sem(meanF,ddof=0)
    M=np.median(meanF)
    Q1=np.percentile(meanF,25)
    Q3=np.percentile(meanF,75)

    TickLabels = [0 for x in range(nExp)]
    for i in range(0,nExp):
        TickLabels[i] = FileList[i][0:7]
        TickLabels[i] = TickLabels[i] + FileList[i][25:27]

    #--- Frequency PLOTS

    #--- 1) STEM PLOT of the Mean Frequency
    fig, ax1 = plt.subplots(figsize=(10, 6))
    plt.stem(meanF)
    ax1.set_xticks(range(0,nExp))
    xtickNames = plt.setp(ax1,xticklabels=TickLabels)
    plt.setp(xtickNames,fontsize=6)
    plt.xlabel('Experiment')
    plt.ylabel('Mean Frequency [Hz]')
    plt.figtext(0.025, 0.9, 'For each exp., Frequency is computed for each channel \
(as 1/<UDcycle>) and MeanFrequency is the mean of Frequencies across channels',
                fontsize=9)
    plt.figtext(0.822, 0.765, 
                'mean=%.2f\nSD=%.2f\nSEM=%.2f\nmedian=%.2f\nQ1=%.2f\nQ3=%.2f'%(m,SD,SEM,M,Q1,Q3),
                fontsize=8, 
                bbox=dict(facecolor='white', edgecolor='black', boxstyle='round'))
    plt.show(0)
    plt.savefig(BaseDir + PythonDir +'MeanFrequency_StemPlot.pdf')

    #--- 2) PLOT of the Mean Frequency with Error Bars
    fig, ax1 = plt.subplots(figsize=(10, 6))
    x = range(nExp)
    plt.errorbar(x,meanF,xerr=0,yerr=semF,fmt='o',capsize=5) #fmt='.'
    ax1.set_xticks(range(0,nExp))
    xtickNames = plt.setp(ax1,xticklabels=TickLabels)
    plt.setp(xtickNames,fontsize=6)
    plt.xlabel('Experiment')
    plt.ylabel('Mean Frequency [Hz]')
    plt.figtext(0.025, 0.95, 'For each exp., Frequency is computed for each channel \
(as 1/<UDcycle>) and MeanFrequency is the mean of Frequencies across channels.',
                fontsize=9)
    plt.figtext(0.025, 0.9, 'For each value, the error is the s.e.m \
(standard error of the mean) computed across channels.',fontsize=9)
    plt.figtext(0.822, 0.76, 
                'mean=%.2f\nSD=%.2f\nSEM=%.2f\nmedian=%.2f\nQ1=%.2f\nQ3=%.2f'%(m,SD,SEM,M,Q1,Q3),
                fontsize=8, 
                bbox=dict(facecolor='white', edgecolor='black', boxstyle='round'))
    #plt.figtext(0.825, 0.75, 
    #            'mean=%s\nSD=%s\nSEM=%s\nmedian=%s\nQ1=%s\nQ3=%s'%(m,SD,SEM,M,Q1,Q3),
    #            fontsize=8,
    #            bbox=dict(facecolor='white', edgecolor='black', boxstyle='round'))
    plt.show(0)
    plt.savefig(BaseDir + PythonDir +'MeanFrequency_Plot.pdf')

    #--- 3) HISTOGRAM of the Mean Frequency 
    #    --> What is the probability distribution of the Mean Frequency?
    fig, ax1 = plt.subplots(figsize=(10, 6))
    #plt.hist(meanF) # (default: 10 bins, from Fmin to Fmax)
    #plt.hist(meanF,bins=np.arange(0.2, 1.2, 0.2),edgecolor = "k")
    # Histogram Range
    hstep=0.1
    hmax=(math.modf(max(meanF)/hstep)[1]+1)*hstep
    hmin=(math.modf(min(meanF)/hstep)[1])*hstep
    plt.hist(meanF,bins=np.arange(hmin,hmax+hstep,hstep),edgecolor = "k")
    plt.title("Mean Frequency per Experiment")
    plt.xlabel("Value [Hz]")
    plt.ylabel("N")
    plt.show(0)
    plt.savefig(BaseDir + PythonDir +'MeanFrequency_Histo.pdf')

    plt.close('all')


# ------------------------------------------------------
#--- SUMMARY DATA Areas (BoxPlot)
# ------------------------------------------------------

if(BoxPlot):
    #NormData = np.array(SummaryData) # First, clone. Then, normalize
    #NormData = np.copy(SummaryData)
    #NormData[:] = SummaryData
    NormDataArea = copy.deepcopy(SummaryDataArea)

    # *** SummaryData[variable][0=median, 1=mean][experiment] ***
    variable=['_down','_up','_UD','_slopeDown','_slopeUp','_peak','frequency']
    Variable = ['DownStateLen', 'UpStateLen', 'UDcycleLen', 'Frequency',
                'SlopeDown','SlopeUp','Peak']
    #what=['M','m'] # SD and SEM do not need to be normalized
    what=['M'] # MEDIAN only (Apr2018)

    for k in range(len(what)):
        name=what[k]*len(variable[0:6])
        varname = ([str(a) + b for a,b in zip(name,variable[0:6])]) 
        for j in range(len(variable[0:6])):
            for i in range(nExp):
                # normalization of data
                NormDataArea[j][k][i]=(SummaryDataArea[j][k][i])/np.mean(SummaryDataArea[j][k][i])

            v=varname[j]  
            print ">>> ",v
# 6 plots for each variable: (3 types of plot) x 2 (SummaryData & NormData) 
# BoxPlot, BoxPlotNotches, ErrorPlot, and each plot is for SummaryData and NormData
            SummaryBoxPlotArea(NormDataArea[j][k],v + '_norm','normalized ' + v,False)
#            SummaryBoxPlotArea(NormDataArea[j][k],v + '_norm','normalized ' + v,True)
            SummaryBoxPlotArea(SummaryDataArea[j][k],v,v,False)
#            SummaryBoxPlotArea(SummaryDataArea[j][k],v,v,True) #(notches)
#            SummaryErrorPlotArea(NormDataArea[j][k],v + '_norm','normalized ' + v)
#            SummaryErrorPlotArea(SummaryDataArea[j][k],v,v)

            plt.close('all') #---close plots (after each variable)

    #--- Plot for the Frequency (FreqA)
    for i in range(nExp):
        NormDataArea[6][0][i]=(SummaryDataArea[6][0][i])/np.mean(SummaryDataArea[6][0][i])
        v=variable[6]  
    print ">>> ",v  

    SummaryBoxPlotArea(NormDataArea[6][0],v + '_norm','normalized ' + v,False)
    #SummaryBoxPlotArea(NormDataArea[6][0],v + '_norm','normalized ' + v,True)
    SummaryBoxPlotArea(SummaryDataArea[6][0],v,v,False)
    #SummaryBoxPlotArea(SummaryDataArea[6][0],v,v,True) #(notches)
    #SummaryErrorPlotArea(NormDataArea[6][0],v + '_norm','normalized ' + v)
    #SummaryErrorPlotArea(SummaryDataArea[6][0],v,v)

    plt.close('all')

    sio.savemat(Dir+'NormDataArea.mat', {'NormDataArea':NormDataArea})
    sio.savemat(Dir+'SummaryDataArea.mat', {'SummaryDataArea':SummaryDataArea})

    nArea=len(area)
# --------------------------
# *** Hypothesis Testing ***
# --------------------------
    alpha0 = 0.05
    Bonferroni = nArea*(nArea-1)/2
    alphaB = alpha0/Bonferroni
    alpha = alpha0
    pvalMatrixArea=np.empty((nVar,nArea,nArea)); pvalMatrixArea.fill(np.nan);

    for i in range(nVar):
        for a1 in range(nArea):
            for a2 in range(a1+1,nArea):
                test=stats.ranksums(NormDataArea[i][0][:,a1],
                                    NormDataArea[0][0][:,a2])
                pvalMatrixArea[i,a1,a2]=test[1]

# ---------------------------------------------------------------
# ***   FDR (False Discovery Rate) ***
# --- Benjamini-Hockberg procedure ---
    pvalArrayArea=[]
    CpvalArrayArea=[] # CORRECTED p-values
    BoolArrayArea=[]
    for i in range(nVar):
        pvalArrayArea.append(pvalMatrixArea[i].ravel())

    for i in range(nVar):
        pvalArrayArea[i]=MaskNaN(pvalArrayArea[i])
        if len(pvalArrayArea[i]) != (nArea*(nArea-1))/2:
            print('Hey! Check the number of multiple hypotheses!)')
        T=multitest.multipletests(pvalArrayArea[i],alpha=alpha,method='fdr_bh')
#        T=multitest.multipletests(pvalArray[i],alpha=alpha,method='bonferroni')
        CpvalArrayArea.append(T[1])
        BoolArrayArea.append(T[0])

# ---------------------------------------------------------------    
#--- Build the BoolMatrix and the CorrectedPvalMatrix
    BoolMatrixArea=np.empty((nVar,nArea,nArea)); BoolMatrixArea.fill(np.nan);
    CpvalMatrixArea=np.empty((nVar,nArea,nArea)); CpvalMatrixArea.fill(np.nan);
    for i in range(nVar):
        j=0
        for a1 in range(nArea):
            for a2 in range(a1+1,nArea):
                BoolMatrixArea[i][a1][a2]=BoolArrayArea[i][j]
                CpvalMatrixArea[i][a1][a2]=CpvalArrayArea[i][j]  
                j=j+1
#--- Symmetrize the BoolMatrix and the CorrectedPvalMatrix
    for i in range(nVar):
        for a2 in range(nArea):
            for a1 in range(a2+1,nArea):
                BoolMatrixArea[i][a1][a2]=BoolMatrixArea[i][a2][a1]
                CpvalMatrixArea[i][a1][a2]=CpvalMatrixArea[i][a2][a1]  

# ------------------------------------------------------    
#--- Plot the BoolMatrix
    for i in range(nVar):
        fig,ax = plt.subplots(figsize=(6,6))
        plt.figtext(0.025,0.025,
                    'pvalue correction: Benjamini-Hochberg (alpha='+str(alpha)+')',
                    bbox=dict(facecolor='none', edgecolor='k', boxstyle='round'))
        ax.grid(True,which='minor',axis='both',linestyle='-',
                color='lightgrey',alpha=0.5)
        ax.set_xticks(range(nArea), minor=True)
        ax.set_yticks(range(nArea), minor=True)
        labels = ['Motor','Somatosensory','Parietal','Retrosplenial','Visual']
        plt.xticks([0.5,1.5,2.5,3.5,4.5],labels)
        plt.yticks([0.5,1.5,2.5,3.5,4.5],labels)
        plt.setp(ax.xaxis.get_majorticklabels(),fontsize=8)
        plt.setp(ax.yaxis.get_majorticklabels(),rotation=45,fontsize=8)
        ax.pcolor(BoolMatrixArea[i], cmap='Greys')
        ax.xaxis.set_ticks_position('none')
        ax.yaxis.set_ticks_position('none') 
        plt.title(Variable[i])
#    plt.colorbar()
#    plt.clim(0,alpha)
        plt.show(0)
        plt.savefig(OutputDir+Variable[i]+'_NORM_MatrixArea_'+str(alpha)+'.pdf')
    plt.close('all')
    
# ------------------------------------------------------    
#--- Plot the BoolMatrix and the BoxPlot
    for i in range(nVar):
# ------------------------------------------------------    
        fig,ax=plt.subplots(figsize=(12, 6))
        plt.subplots_adjust(left=0.075, right=0.95, top=0.9, bottom=0.1)
        gs = gridspec.GridSpec(1,2,width_ratios=[1,1])
        plt.figtext(0.025,0.025,
                    'pvalue correction: Benjamini-Hochberg (alpha='+str(alpha)+')',
                    bbox=dict(facecolor='none', edgecolor='k', boxstyle='round'))
# ------------------------------------------------------            
        ax=plt.subplot(gs[0])
        ###fig,ax = plt.subplots(1,1)
        ax.grid(True,which='minor',axis='both',linestyle='-',
                color='lightgrey',alpha=0.5)
        ax.set_xticks(range(nArea), minor=True)
        ax.set_yticks(range(nArea), minor=True)
        labels = ['Motor','Somatosensory','Parietal','Retrosplenial','Visual']
        plt.xticks([0.5,1.5,2.5,3.5,4.5],labels)
        plt.yticks([0.5,1.5,2.5,3.5,4.5],labels)
        plt.setp(ax.xaxis.get_majorticklabels(),fontsize=8)
        plt.setp(ax.yaxis.get_majorticklabels(),rotation=45,fontsize=8)
        ax.pcolor(BoolMatrixArea[i], cmap='Greys')
        ax.xaxis.set_ticks_position('none')
        ax.yaxis.set_ticks_position('none') 
        plt.title(Variable[i])
#    plt.colorbar()
#    plt.clim(0,alpha)
# ------------------------------------------------------            
        ax2=plt.subplot(gs[1])
        ax2=SummaryBoxPlotAreaNEW(NormDataArea[i][0],Variable[i] + '_norm',
                                  'normalized '+Variable[i],ax2)
# ------------------------------------------------------
        plt.show(0)
        plt.savefig(OutputDir+Variable[i]+
                    '_NORM_MatrixBoxPlotArea_'+str(alpha)+'.pdf')
    plt.close('all')

# ------------------------------------------------------    
#--- Plot the CpvalMatrix and the BoxPlot
    for i in range(nVar):
# ------------------------------------------------------    
        fig,ax=plt.subplots(figsize=(12, 6))
        plt.subplots_adjust(left=0.075, right=0.95, top=0.9, bottom=0.1)
        gs = gridspec.GridSpec(1,2,width_ratios=[1,1])
        plt.figtext(0.025,0.025,
                    'pvalue correction: Benjamini-Hochberg',
                    bbox=dict(facecolor='none', edgecolor='k', boxstyle='round'))
# ------------------------------------------------------            
        ax=plt.subplot(gs[0])
        ###fig,ax = plt.subplots(1,1)
        ax.grid(True,which='minor',axis='both',linestyle='-',
                color='lightgrey',alpha=0.5)
        ax.set_xticks(range(nArea), minor=True)
        ax.set_yticks(range(nArea), minor=True)
        labels = ['Motor','Somatosensory','Parietal','Retrosplenial','Visual']
        plt.xticks([0.5,1.5,2.5,3.5,4.5],labels)
        plt.yticks([0.5,1.5,2.5,3.5,4.5],labels)
        plt.setp(ax.xaxis.get_majorticklabels(),fontsize=8)
        plt.setp(ax.yaxis.get_majorticklabels(),rotation=45,fontsize=8)
        plt.pcolor(CpvalMatrixArea[i], cmap='hot')
        ax.xaxis.set_ticks_position('none')
        ax.yaxis.set_ticks_position('none') 
        plt.title(Variable[i])
        plt.colorbar()
        plt.clim(0,alpha)
# ------------------------------------------------------            
        ax2=plt.subplot(gs[1])
        ax2=SummaryBoxPlotAreaNEW(NormDataArea[i][0],Variable[i] + '_norm',
                                  'normalized '+Variable[i],ax2)
# ------------------------------------------------------
        plt.show(0)
        plt.savefig(OutputDir+Variable[i]+
                    '_NORM_pvalMatrixBoxPlotArea.pdf')
    plt.close('all')

    
# ------------------------------------------------------
#--- SUMMARY DATA Channels (ContourPlot)
# ------------------------------------------------------

if(ContourPlot):
    NormDataCh = copy.deepcopy(SummaryDataCh)
    nVar = len(SummaryDataCh) # number of variables
    variable = ['DownStateLen', 'UpStateLen', 'UDcycleLen', 'Frequency',
                'SlopeDown','SlopeUp','Peak']
    unit = [' [s]', ' [s]', ' [s]', ' [Hz]',
                ' [1/s]',' [1/s]',' [a.u.]']

#--- SummaryDataCh[variable][experiment] --> NormDataCh
    for i in range(nVar):                   
        for j in range(nExp):
            #--- exclude SD outliers (for comparison of summary data)
            k=[ChList.index(out) for out in wtOUTsd[j]]
            SummaryDataCh[i][j][k]=np.nan
            ### N.B. already done in analyzeRecordingSetResults if OUTall=True
#--- clean with the Mask that identify NaN (i.e. missing) channels
            Clean = MaskNaN(SummaryDataCh[i][j])
#--- compute the mean considering only non-Nan channels
            Mean = Clean.mean() 
            NormDataCh[i][j] = SummaryDataCh[i][j]/Mean

    storeA=np.empty((nVar,nCh,nExp))
    MapDataA = np.empty((nVar,nCh)) # data for the ContourPlot (average map)
    #--- for each variable and each channel, mean of SummaryDataCh across exps   
    StdDataA = np.empty((nVar,nCh))
    SemDataA = np.empty((nVar,nCh))
    RelErrA = np.empty((nVar,nCh))

    storeB=np.empty((nVar,nCh,nExp))
    MapDataB = np.empty((nVar,nCh)) # data for the ContourPlot
    #--- for each variable and each channel, mean of NormDataCh across experiments
    StdDataB = np.empty((nVar,nCh))
    SemDataB = np.empty((nVar,nCh))
    RelErrB = np.empty((nVar,nCh))

    numel = np.empty((nCh),dtype=np.uint8) 

    for i in range(nVar):
        for k in range(nCh):
            
            A=np.array([])
            for j in range(nExp):
                A=np.append(A,SummaryDataCh[i][j][k])
            storeA[i][k]=A

            B=np.array([])
            for j in range(nExp):
                B=np.append(B,NormDataCh[i][j][k])
            storeB[i][k]=B
                
            CleanA=MaskNaN(A)
            CleanB=MaskNaN(B)

            numel[k]=len(CleanA)

            MapDataA[i][k]=CleanA.mean()
            StdDataA[i][k]=CleanA.std()
            SemDataA[i][k]=CleanA.std()/math.sqrt(numel[k])
            
            MapDataB[i][k]=CleanB.mean()
            StdDataB[i][k]=CleanB.std()
            SemDataB[i][k]=CleanB.std()/math.sqrt(numel[k])

    RelErrA=abs(StdDataA/MapDataA)
    RelErrB=abs(StdDataB/MapDataB)

    OutputDir = BaseDir + PythonDir

# --------------------------
# *** Hypothesis Testing ***
# --------------------------
    alpha0 = 0.05
    Bonferroni = nCh*(nCh-1)/2
    alphaB = alpha0/Bonferroni
    alpha = alpha0
    pvalMatrix=np.empty((nVar,nCh,nCh)); pvalMatrix.fill(np.nan);

    for i in range(nVar):
        print ">>> ", variable[i]

# 1) SummaryData
        fig,ax1=plt.subplots(figsize=(5, 6)) 
        FileName='messageA'
        ax1=Contour(MapDataA[i],RelErrA[i],variable[i],unit[i],ax1)
### Significance [WilcoxonRankSum]
        for k1 in range(nCh):
            for k2 in range(k1+1,nCh):
                test=stats.ranksums(MaskNaN(storeA[i][k1]),MaskNaN(storeA[i][k2]))
        plt.savefig(OutputDir + variable[i] + '_ContourPlot.pdf')
        
# 2) NormData
        FileName='messageB'
        fig,ax2=plt.subplots(figsize=(5, 6)) 
        ax2=Contour(MapDataB[i],RelErrB[i],variable[i],unit[i],ax2)
### Significance [WilcoxonRankSum]
        for k1 in range(nCh):
            for k2 in range(k1+1,nCh):
                test=stats.ranksums(MaskNaN(storeB[i][k1]),MaskNaN(storeB[i][k2]))
                pvalMatrix[i,k1,k2]=test[1]
        plt.savefig(OutputDir + variable[i] + '_NORM_ContourPlot.pdf')

# ---------------------------------------------------------------
# ***   FDR (False Discovery Rate) ***
# --- Benjamini-Hockberg procedure ---
    pvalArray=[]
    CpvalArray=[] # CORRECTED p-values
    BoolArray=[]
    for i in range(nVar):
        pvalArray.append(pvalMatrix[i].ravel())

    for i in range(nVar):
        pvalArray[i]=MaskNaN(pvalArray[i])
        if len(pvalArray[i]) != (nCh*(nCh-1))/2:
            print('Hey! Check the number of multiple hypotheses!)')
        T=multitest.multipletests(pvalArray[i],alpha=alpha,method='fdr_bh')
#        T=multitest.multipletests(pvalArray[i],alpha=alpha,method='bonferroni')
        CpvalArray.append(T[1])
        BoolArray.append(T[0])

# ---------------------------------------------------------------    
#--- Build the BoolMatrix
    BoolMatrix=np.empty((nVar,nCh,nCh)); BoolMatrix.fill(np.nan);
    CpvalMatrix=np.empty((nVar,nCh,nCh)); CpvalMatrix.fill(np.nan);
    for i in range(nVar):
        j=0
        for k1 in range(nCh):
            for k2 in range(k1+1,nCh):
                BoolMatrix[i][k1][k2]=BoolArray[i][j]
                CpvalMatrix[i][k1][k2]=CpvalArray[i][j]  
                j=j+1
#--- Symmetrize the BoolMatrix
    for i in range(nVar):
        for k2 in range(nCh):
            for k1 in range(k2+1,nCh):
                BoolMatrix[i][k1][k2]=BoolMatrix[i][k2][k1]
                CpvalMatrix[i][k1][k2]=CpvalMatrix[i][k2][k1]

# -----------------------------------------------------------------                   
#--- Plot the ContourPlot and the BoolMatrix
    for i in range(nVar):
# ------------------------------------------------------    
        fig,ax=plt.subplots(figsize=(10, 6)) 
        gs = gridspec.GridSpec(1,2,width_ratios=[1.5,1])
        plt.figtext(0.025,0.025,
                    'pvalue correction: Benjamini-Hochberg (alpha='+str(alpha)+')',
                    bbox=dict(facecolor='none', edgecolor='k', boxstyle='round'))
# ------------------------------------------------------    
        ax1=plt.subplot(gs[0]) 
#        fig,ax = plt.subplots(1,1,figsize=(6,6))
        ax1.grid(True,which='minor',axis='both',\
                   linestyle='-',color='lightgrey',alpha=0.5)
        ax1.set_xticks(range(nCh), minor=True)
        ax1.set_yticks(range(nCh), minor=True)
        ax1.pcolor(BoolMatrix[i], cmap='Greys')

        ax1.tick_params(axis=u'both',direction="in",length=0)

    ##... red lines to separate electrodes of different cortical areas
        ax1.plot((5,5),(0,32),'r')
        ax1.plot((0,32),(5,5),'r')
        ax1.plot((12,12),(0,32),'r')
        ax1.plot((0,32),(12,12),'r')
        ax1.plot((15,15),(0,32),'r')
        ax1.plot((0,32),(15,15),'r')
        ax1.plot((19,19),(0,32),'r')
        ax1.plot((0,32),(19,19),'r')

#        ax1.set_yticklabels([])
#        ax1.set_xticklabels([])
#        ax1.set_yticklabels(ChListStr,minor=True,fontsize=6)
#        ax1.set_xticklabels(ChListStr,minor=True,fontsize=6)

        for axis in [ax1.xaxis, ax1.yaxis]:
            axis.set(ticks=np.arange(0.5, len(ChListStr)), ticklabels=ChListStr)
            axis.set_tick_params(labelsize=6)

#        labels = ['Motor','Somatosensory','Parietal','Retrosplenial','Visual']
#        ax1.text(0.175,0.06,'M'); ax1.text(0.07,0.16,'M')
#        ax1.text(0.32,0.06,'S');  ax1.text(0.07,0.30,'S')
#        ax1.text(0.445,0.06,'P'); ax1.text(0.07,0.425,'P')
#        ax1.text(0.53,0.06,'R');  ax1.text(0.07,0.51,'R')
#        ax1.text(0.735,0.06,'V'); ax1.text(0.07,0.717,'V')

        plt.title(variable[i])
# ------------------------------------------------------            
        ax2=plt.subplot(gs[1])
        FileName='messageB'
        ax2=Contour(MapDataB[i],RelErrB[i],variable[i],unit[i],ax2)

        #--- Find the most significant channels                
        SigCh=np.zeros(nCh)            
        for k in range(nCh):
            SigCh[k]=sum(MaskNaN(BoolMatrix[i][k]))
        idx = (-SigCh).argsort()[:3]


        myColor=['r','b','k']
        for k1 in idx:
            for k2 in range(nCh):
                if BoolMatrix[i,k1,k2]:
                    ax2.plot((ArrayMap[k1][0]*550e-3, ArrayMap[k2][0]*550e-3), 
                             (ArrayMap[k1][1]*550e-3, ArrayMap[k2][1]*550e-3),
                             color=myColor[np.where(idx==k1)[0][0]])
                             ##color=myColor[np.squeeze(np.where(idx==k1)[0][0]]))
# ------------------------------------------------------    
        plt.show(0)
        plt.savefig(OutputDir+variable[i]+'_NORM_Matrix-and-Map_'+str(alpha)+'.pdf')
        
    plt.close('all')
 

# -----------------------------------------------------------------                   
#--- Plot the ContourPlot and the CpvalMatrix
    for i in range(nVar):
# ------------------------------------------------------    
        fig,ax=plt.subplots(figsize=(10, 6))
        plt.subplots_adjust(left=0.05, right=0.95)
        gs = gridspec.GridSpec(1,2,width_ratios=[1.5,1])
        plt.figtext(0.025,0.025,
                    'pvalue correction: Benjamini-Hochberg',
                    bbox=dict(facecolor='none', edgecolor='k', boxstyle='round'))
# ------------------------------------------------------    
        ax1=plt.subplot(gs[0]) 
#        fig,ax = plt.subplots(1,1,figsize=(6,6))
        ax1.grid(True,which='minor',axis='both',\
                   linestyle='-',color='lightgrey',alpha=0.5)
        ax1.set_xticks(range(nCh), minor=True)
        ax1.set_yticks(range(nCh), minor=True)

        ax1.tick_params(axis=u'both',direction="in",length=0)
        
        plt.pcolor(CpvalMatrix[i], cmap='hot')
        plt.colorbar()
        plt.clim(0,alpha)
        
    ##... red lines to separate electrodes of different cortical areas
        ax1.plot((5,5),(0,32),'r')
        ax1.plot((0,32),(5,5),'r')
        ax1.plot((12,12),(0,32),'r')
        ax1.plot((0,32),(12,12),'r')
        ax1.plot((15,15),(0,32),'r')
        ax1.plot((0,32),(15,15),'r')
        ax1.plot((19,19),(0,32),'r')
        ax1.plot((0,32),(19,19),'r')

        #ax1.set_yticklabels([])
        #ax1.set_xticklabels([])
        #ax1.set_yticklabels(ChListStr,minor=True,fontsize=6)
        #ax1.set_xticklabels(ChListStr,minor=True,fontsize=6)

        for axis in [ax1.xaxis, ax1.yaxis]:
            axis.set(ticks=np.arange(0.5, len(ChListStr)), ticklabels=ChListStr)
            axis.set_tick_params(labelsize=6)

#        labels = ['Motor','Somatosensory','Parietal','Retrosplenial','Visual']
#        ax1.text(0.175,0.06,'M'); ax1.text(0.07,0.16,'M')
#        ax1.text(0.32,0.06,'S');  ax1.text(0.07,0.30,'S')
#        ax1.text(0.445,0.06,'P'); ax1.text(0.07,0.425,'P')
#        ax1.text(0.53,0.06,'R');  ax1.text(0.07,0.51,'R')
#        ax1.text(0.735,0.06,'V'); ax1.text(0.07,0.717,'V')

        plt.title(variable[i])
# ------------------------------------------------------            
        ax2=plt.subplot(gs[1])
        FileName='messageB'
        ax2=Contour(MapDataB[i],RelErrB[i],variable[i],unit[i],ax2)

        #--- Find the most significant channels                
        SigCh=np.zeros(nCh)            
        for k in range(nCh):
            SigCh[k]=sum(MaskNaN(BoolMatrix[i][k]))
        idx = (-SigCh).argsort()[:3]


        myColor=['r','b','k']
        for k1 in idx:
            for k2 in range(nCh):
                if BoolMatrix[i,k1,k2]:
                    ax2.plot((ArrayMap[k1][0]*550e-3, ArrayMap[k2][0]*550e-3), 
                             (ArrayMap[k1][1]*550e-3, ArrayMap[k2][1]*550e-3),
                             color=myColor[np.where(idx==k1)[0][0]])
                             ##color=myColor[np.squeeze(np.where(idx==k1)[0][0]]))
# ------------------------------------------------------    
        plt.show(0)
        plt.savefig(OutputDir+variable[i]+'_NORM_pvalMatrix-and-Map.pdf')
        
    plt.close('all')

    
### Boolean Matrix    
#    pvalMatrixBool=copy.deepcopy(pvalMatrix)
#    for i in range(nVar):
#        for k1 in range(nCh):
#            for k2 in range(k1+1,nCh):
#                if pvalMatrix[i][k1][k2]<alpha:
#                    pvalMatrixBool[i][k1][k2]=1
#                else:
#                    pvalMatrixBool[i][k1][k2]=0

#    for i in range(nVar):
#        for k2 in range(nCh):
#            for k1 in range(k2+1,nCh):
#                pvalMatrixBool[i][k1][k2]=pvalMatrixBool[i][k2][k1]

#    fig,ax = plt.subplots(1,1,figsize=(6,6))
#    ax.grid(True,which='minor',axis='both',\
#               linestyle='-',color='lightgrey',alpha=0.5)
#    ax.set_xticks(range(nCh), minor=True)
#    ax.set_yticks(range(nCh), minor=True)
#    ax.pcolor(pvalMatrixBool[0], cmap='Greys')

##... red lines to separate electrodes of different cortical areas
#    ax.plot((5,5),(0,32),'r')
#    ax.plot((0,32),(5,5),'r')
#    ax.plot((12,12),(0,32),'r')
#    ax.plot((0,32),(12,12),'r')
#    ax.plot((15,15),(0,32),'r')
#    ax.plot((0,32),(15,15),'r')
#    ax.plot((19,19),(0,32),'r')
#    ax.plot((0,32),(19,19),'r')
    
#    ax.set_yticklabels([])
#    ax.set_xticklabels([])
    
#    plt.show(0)                  

    
##    for i in range(nVar):
##        #fig,ax = plt.subplots(figsize=(10,10))
##        plt.figure()
##        plt.pcolor(pvalMatrix[i], cmap='Greys_r')
##        plt.colorbar()# colorbar = add a color bar
##        plt.clim(0,alpha)
##        plt.title(variable[i])
##        plt.show(0)
##        plt.savefig(OutputDir + variable[i]+'_pvalMatrix.pdf')
##    plt.close('all') 
    
    
#    plt.figure()
#    plt.imshow(pvalMatrix[0])
#    plt.colorbar()
#    plt.show(0)

##    plt.figure()
#    fig,ax = plt.subplots(1,1)
#    ax.grid(True,which='minor',axis='both',linestyle='-',color='lightgrey',alpha=0.5)
#    ax.set_xticks(range(nCh), minor=True)
#    ax.set_yticks(range(nCh), minor=True)
#    plt.pcolor(pvalMatrix[0], cmap='greys_r')
#    plt.colorbar()
#    plt.clim(0,alpha)
#    plt.show(0)

#    fig,axs=plt.subplots(2,4,figsize=(15,6))
#    fig.subplots_adjust(hspace = .5, wspace=.001)
#    axs=axs.ravel()
#    for i in range(nVar):
#        axs[i].pcolor(pvalMatrix[i], cmap='Greys_r')
#        #plt.colorbar()# colorbar = add a color bar
#        #plt.clim(0,0.05)
#        axs[i].set_title(variable[i])
#    plt.show(0)
#    plt.savefig(OutputDir + 'pvalMatrix.pdf')

#    fig,ax = plt.subplots(1,2)
#    ax=ax.ravel()
#    ax[0].grid(True,which='minor',axis='both',\
#               linestyle='-',color='lightgrey',alpha=0.5)
#    ax[0].set_xticks(range(nCh), minor=True)
#    ax[0].set_yticks(range(nCh), minor=True)
#    ax[0].pcolor(pvalMatrixBool[0], cmap='Greys')
#    ax[1].plot(ArrayMap)
    


## 'Greys' = from white to Black (or: cmap='binary')  (cmap=cm.Greys)       
## 'grey' = from Black to White

#plt.title('Title')
#fig.suptitle('Title', fontsize=14, fontweight='bold')
#fig.text(0.21,0.9,'text')

# ------------------------------------------------------
### Significance ###
## Ho: the value is the same at each electrode position 
## (and equal to the mean across electrodes)

#    for i in range(nVar):
#        ax1=Contour(MapDataA[i],RelErrA[i],variable[i],unit[i])
#        H0=MapDataA[i].mean()
#        sigma0=MapDataA[i].std()
#        for k in range(nCh):
#            if abs(MapDataA[i][k]-H0)>2*sigma0:
#                print ChList[k]
#                ax1.add_patch(patches.Circle(
#                    (ArrayMap[k][0]*550e-3,ArrayMap[k][1]*550e-3),
#                    radius=0.2,edgecolor='w',facecolor='None'))   
# ------------------------------------------------------
### Significance [ttest] ###

#    for i in range(nVar):
#        ax1=Contour(MapDataA[i],RelErrA[i],variable[i],unit[i])
#        for k1 in range(nCh):
#            for k2 in range(k1+1,nCh):
#                test=stats.ttest_ind(MaskNaN(storeA[i][k1]),MaskNaN(storeA[i][k2]))

### SignificanceLines on the Map (without pval correction)
#                 if test[1]< alpha:
#                    ax1.plot((ArrayMap[k1][0]*550e-3, ArrayMap[k2][0]*550e-3), 
#                             (ArrayMap[k1][1]*550e-3, ArrayMap[k2][1]*550e-3),'r') 
# ------------------------------------------------------
### [for check] ###
#    sort_index = np.argsort(ChList)
#    SampleSize=numel[sort_index]
#    nOut=nExp-SampleSize
# ------------------------------------------------------


if StatTest:

    # *** ANOVA *** (STATISTICAL SIGNIFCICANCE of the RESULTS)

    # Focus on MEDIAN per area. Two approaches: 
    # 1) T-Student Test (TTest)
    # 2) Wilcoxon Rank Sum Test (WTest)--> 
    #        can be used as an alternative to the paired TTest, 
    #        when the population cannot be assumed to be normally distributed


    # 1) --- TTtest ---
    # d.o.f. = n1+n2-2 = 11+11-2 = 20 (n=11 experiments)
    dof=20
    alpha = 0.001 
    halpha = alpha/2
    par=1-halpha
    TC = stats.t.ppf(par,dof)

    outfile1 = BaseDir + PythonDir + 'TTest.txt'
    F1 = open(outfile1,'w')
    F1.write("alpha=%s (2-sided)\n" %alpha)
    F1.write("T_CriticalValue=%s\n" %TC)

    # 2) --- Wilcoxon RankSum Test ---
    sig = 3 # threshold significance: 3sigma 
    alpha2=0.001; # Level of Significance
    outfile2 = BaseDir + PythonDir + 'WilcoxonTest.txt'
    F2 = open(outfile2,'w')
    F2.write("alpha=%s (2-sided)\n" %alpha2)
    #F2.write("significance=%s sigma\n" %sig)

    for j in range(len(variable)):
        print(variable[j])
        F1.write("*********************\n%s\n" %variable[j])
        F2.write("*********************\n%s\n" %variable[j])
        for k in range(len(area)):
            for i in range(k+1,len(area)):
                pair=[str(area[k][0])+'-'+str(area[i][0])]
                #print(pair)
                F1.write("%s\n" %pair)
                F2.write("%s\n" %pair)

                # 1) TTest
                TTest=(stats.ttest_ind(NormData[j][0][:,k],NormData[j][0][:,i])) 
                # [][0][]--> consider always the median
                #print(TTest)
                F1.write(str(TTest)+"\n")
                if abs(TTest[0])>TC:
                    significance=stats.norm.ppf(1-TTest[1]/2)
                    print(str(pair)+" >>> H0 rejected")
                    print("significance (in units of sigma): %s" %significance)
                    F1.write(">>> H0 rejected\n")
                    F1.write("significance (in units of sigma): %s\n" %significance)
                    F1.write("p-value (2-sided): %s\n" %TTest[1])

                # 2) Wilcoxon RankSum Test
                Wilco=(stats.ranksums(NormData[j][0][:,k],NormData[j][0][:,i])) 
                F2.write(str(Wilco)+"\n")
                #if abs(Wilco[0])>3:
                if Wilco[1]<=alpha:
                    significance=abs(Wilco[0])
                    print(str(pair)+" >>> H0 rejected")
                    print("significance (in units of sigma): %s" %significance)
                    F2.write(">>> H0 rejected\n")
                    F2.write("significance (in units of sigma): %s\n" %significance)
                    F2.write("p-value (2-sided): %s\n" %Wilco[1])


                # 3) Mann-Whitney Test
                MW=(stats.mannwhitneyu(NormData[j][0][:,k],NormData[j][0][:,i],alternative='two-sided')) 
                F3.write(str(MW)+"\n")
                #if abs(Wilco[0])>3:
                if MW[1]<=alpha:
                    significance=abs(MW[0])
                    print(str(pair)+" >>> H0 rejected")
                    print("significance (in units of sigma): %s" %significance)
                    F3.write(">>> H0 rejected\n")
                    F3.write("significance (in units of sigma): %s\n" %significance)
                    F3.write("p-value (2-sided): %s\n" %MW[1])

    F1.close()
    F2.close()
    F3.close()
