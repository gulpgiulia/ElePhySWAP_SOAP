#!/usr/bin/python
#--------1---------2---------3---------4---------5---------6---------7--------X
#  Copyright 2017 Giulia De Bonis @ INFN, Rome - Italy
#  Version: 1.0 - June 8, 2017
#--------1---------2---------3---------4---------5---------6---------7--------X
"Module with user-defined functions (and some constants)"

#[NB] currently the module is *executed* and *not imported* in the main
#--- For a 'proper' import: 
#--- 1) 'importPackages' should include the constants, and should be executed 
#--- by the module 'defFunction'
#--- 2) datastruct (aut similia) should be passed as parameters

#execfile("importPackages.py") 
#print("...importPackages...")

#*** DEBUG Flux Control ***
#    if debug:
#        return
#        #pass
#    else:
#        pass
#**************************

#-----------------------------------------------------
# CONSTANTS
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

#-----------------------------------------------------
# FUNCTIONS
#-----------------------------------------------------
def MaskNaN(variable):
    "Apply a Mask to exclude NaN entries \
(... useful for computing mean and median and percentiles...)"

#>>> RETURN: variable without NaN entries    

    mask = np.ones(len(variable), dtype=bool)
    
    for i in range(len(variable)):
        if np.isnan(variable[i]):
            mask[[i]] = False
    var=variable[mask]
    return(var)

#--------1---------2---------3---------4---------5---------6---------7--------X
def ArrangeChannelsByArea(variable,outliers,area):
    "Compute median and mean of the variable for each channel.\
Outlier channels are excluded. Channels are numbered and grouped by cortical areas"
    
    Variable= [] # LIST (data re-arranged and channels grouped by cortical areas)

    #[OR]M = np.empty((0,nCh))
    M = np.array([])      # median (of the channel)
    meanCh = np.array([]) # mean (of the channel)
    stdCh = np.array([])  # standard deviation (of the channel)
    semCh = np.array([])  # standard error of the mean
    ChList = np.array([]) # (string) --- channel numbering
    Separator=np.array([0.])   # (AreaSeparator, for visualization)
  

    for i in area:
        for j in i[1]:
            if j not in outliers-1:
                Variable.append(variable[j])
                M=np.append(M,np.median(variable[j]))
                meanCh=np.append(meanCh,np.mean(variable[j]))   
                stdCh=np.append(stdCh,np.std(variable[j]))
                semCh=np.append(semCh,np.std(variable[j])/math.sqrt(len(variable[j])))
            else:
                Variable.append(np.array([])) #(list-of-arrays)               
                M=np.append(M,np.nan)
                meanCh=np.append(meanCh,np.nan)
                stdCh=np.append(stdCh,np.nan)
                semCh=np.append(semCh,np.nan)

    return(Variable,M,meanCh)

#--------1---------2---------3---------4---------5---------6---------7--------X
def BoxPlot(variable):
    "BoxPlot of a variable, with mean values.\
Returns 'ax' (axes) for setting graphical properties"
                
# --- BOXPLOT (with mean values)
    fig, ax = plt.subplots(figsize=(10,6))
    plt.subplots_adjust(left=0.075, right=0.95, top=0.9, bottom=0.1)
    meanpoint = dict(marker='*',markerfacecolor='blue',markeredgecolor='blue')
    plt.boxplot(variable,sym='+',showmeans=True,meanprops=meanpoint) 
    plt.show(0)

# --- legend
    plt.figtext(0.85,0.92,'*',color='b',weight='roman',size='medium')
    plt.figtext(0.865, 0.92,'Mean Value',color='black',weight='roman',size='medium',
                bbox=dict(facecolor='none',edgecolor='k',boxstyle='round'))
    
    return(ax)
    
#--------1---------2---------3---------4---------5---------6---------7--------X
def BoxPlotChByArea(variable,yLabel,ChList,Separator,Filename):
    "Produce a BoxPlot with the distribution of the variable for each channel. \
Channels are numbered and grouped by Cortical Area"

    ax1=BoxPlot(variable)
 
# --- x-axis
    ax1.set_xticklabels(ChList)
    ax1.set_xlabel('Channel') # *** x-axis --- Label

# --- separators
    height=ax1.get_ylim()[1]-ax1.get_ylim()[0]
    TickPos = []
    offset=0.5
    for i in range(len(Separator)-1):
        TickPos.append(Separator[i] + (Separator[i+1]-Separator[i])/2)
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

# --- y-axis
    ax1.set_ylabel(yLabel) # *** y-axis --- Label

# --- secondary x-axis
    labels = ['Motor','Somatosensory','Parietal','Retrosplenial','Visual']
    TickPosNorm=[i/len(ChList) for i in TickPos]
    ax2 = ax1.twiny()
    ax2.set_xticks(TickPosNorm)
    ax2.set_xticklabels(labels)

# --- FileName in the canvas
    plt.figtext(0.05,0.02,FileName,color='blue',weight='roman',size='medium',
                bbox=dict(facecolor='none', edgecolor='b', boxstyle='round'))

#--------1---------2---------3---------4---------5---------6---------7--------X
def BoxPlotCh(variable,yLabel,debug):
    "Produce a BoxPlot with the distribution of the variable for each channel. \
Channels are numbered and grouped by Cortical Area"

#>>> RETURN:
#[0] array of medians at each channels (channels numbered by cortical areas);
#    nan for outlier channels
#[1] ChList (...to be removed... it's a constant... non-sense to be computed here)
#[2] array of means at each channels (channels numbered by cortical areas);
#    nan for outlier channels
#[3] array of sem at each channels (channels numbered by cortical areas);
#    nan for outlier channels

    testVar = False
    ddebug = False #(deep-debug)
    
    #[WRONG]Variable= np.empty([nCh]) 
# --- Variable must be a LIST and not an array: otherwise...
# --- "ValueError: setting an array element with a sequence."
    Variable= [] # LIST
    if(ddebug):
        print 'type(Variable):',type(Variable)
        #print len(Variable)
    #[OR]M = np.empty((0,nCh))
    M = np.array([])      # median (of the channel)
    meanCh = np.array([]) # mean (of the channel)
    stdCh = np.array([])  # standard deviation (of the channel)
    semCh = np.array([])  # standard error of the mean
    ChList = np.array([]) # (string) --- channel numbering
    Separator=np.array([0.])   # (AreaSeparator, for visualization)
  
    outliers_reShape = outliers-1 #(Matlab vs Python array indexing)

    for i in area:
        Separator=np.append(Separator,Separator[-1]+len(i[1]))
        for j in i[1]:
            if j not in outliers_reShape:
                #[WRONG]Variable[j]=datastruct[0,j][variable] #(array)
                #[WRONG]Variable.append([datastruct[0,j][variable]]) #(list-of-lists)
                #Variable.append(np.squeeze(datastruct[0,j][variable])) #(list-of-arrays)
# --- NOTE the np.squeeze. Alternatively, extract the dimension:                
                Variable.append(datastruct[variable][0,j][0])

                M=np.append(M,np.median(datastruct[variable][0,j][0]))                
                meanCh=np.append(meanCh,np.mean(datastruct[variable][0,j][0]))
                stdCh=np.append(stdCh,np.std(datastruct[variable][0,j][0]))
                semCh=np.append(semCh,np.std(datastruct[variable][0,j][0])
                                /math.sqrt(len(datastruct['DownStateLen'][0,j][0])))
                #a = np.array(np.array([datastruct[0,j][variable]][0])[0])
                #semCh=np.append(semCh,stats.sem(a))
            else:
                #[WRONG]Variable[j]=np.nan #(array --- but produces ValueError)
                #[WRONG]Variable.append([]) #(list-of-lists)
                Variable.append(np.array([])) #(list-of-arrays)
                
                M=np.append(M,np.nan)
                meanCh=np.append(meanCh,np.nan)
                stdCh=np.append(stdCh,np.nan)
                semCh=np.append(semCh,np.nan)
            ChList=np.append(ChList,str(j+1))  # The labels for the bottom x axis
            # print len(Variable)

    if(testVar):
        A=[]
        for i in range(5):
            if i !=3:            
                A.append(np.random.rand(int(np.random.rand(1)*100)))
            else:
                A.append(np.empty(1))

        for el in A:
            print(type(el))
            print el.shape
            
        print('***',type(A))
        A=np.asarray(A)
        print('***',type(A))

    if(ddebug):
        print 'len(Variable)',len(Variable)
        i=0
        for el in Variable:
            i=i+1
            print '--- element:',i,type(el),el.shape
            #print np.squeeze(el).shape
        
        print 'type(Variable)',type(Variable)
        try:
            Variable=np.asarray(Variable)
            print('np.asarray WORKS')
            print 'type(Variable)',type(Variable)
        except:
            print('np.asarray FAILS')
            print 'type(Variable)',type(Variable)


# --- very rough BUG FIXING --- (before the squeeze solution)
#    if outliers.any() == 0:
#    if outliers.size == 0:
#        Variable.append([])
    

# --- BOXPLOT (with mean values)
    fig, ax1 = plt.subplots(figsize=(10, 6))
    plt.subplots_adjust(left=0.075, right=0.95, top=0.9, bottom=0.1)
    meanpoint = dict(marker='*',markerfacecolor='blue',markeredgecolor='blue')
    plt.boxplot(Variable,sym='+',showmeans=True,meanprops=meanpoint) 
 
# --- x-axis
    ax1.set_xticklabels(ChList)
    ax1.set_xlabel('Channel') # *** x-axis --- Label

# --- separators
    height=ax1.get_ylim()[1]-ax1.get_ylim()[0]
    TickPos = []
    offset=0.5
    for i in range(0,len(area)):
        TickPos.append(Separator[i] + (Separator[i+1]-Separator[i])/2)
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

# --- y-axis
    plt.ylabel(yLabel) # *** y-axis --- Label

# --- secondary x-axis
    labels = ['Motor','Somatosensory','Parietal','Retrosplenial','Visual']
    TickPosNorm=[i/nCh for i in TickPos]
    ax2 = ax1.twiny()
    ax2.set_xticks(TickPosNorm)
    ax2.set_xticklabels(labels)
    #ax2.set_xlabel('Cortical Area')

# --- legend
    plt.figtext(0.85,0.92,'*',color='b',weight='roman',size='medium')
    plt.figtext(0.865, 0.92,'Mean Value',color='black',weight='roman',size='medium',
                bbox=dict(facecolor='none',edgecolor='k',boxstyle='round'))

# --- FileName in the canvas
    plt.figtext(0.05,0.02,FileName,color='blue',weight='roman',size='medium',
                bbox=dict(facecolor='none', edgecolor='b', boxstyle='round'))
                
    plt.show(0)
    if ~debug:
        plt.savefig(OutputDir + variable +'_BoxPlot_Areas-and-Channels.pdf')

    return(M,ChList,meanCh,semCh)

#--------1---------2---------3---------4---------5---------6---------7--------X
def Contour(variable,error,name,unit,ax1):
    "Produce a CountourPlot of the variable on the Map of the ElectrodeArray"

    ddebug = False #(deep-debug)

#--- make data (Electrode Coordinates in the Map)
    X = np.empty((0,nCh))
    Y = np.empty((0,nCh))

    for i in range(nCh):
        X=np.append(X,ArrayMap[i][0]*550e-3) # unit: mm
        Y=np.append(Y,ArrayMap[i][1]*550e-3)
    
    if(ddebug):
        print 'X:', X
        print 'Y:', Y
        print 'Variable:', variable

#--- "clean" data = exclude NaN entries (and the corresponding electrode position XY)
    mask = np.ones(len(variable), dtype=bool)
    for i in range(len(variable)):
        if np.isnan(variable[i]):
            mask[[i]] = False
            #print 'NaN at channel ',ChList[i]
    z=variable[mask]
    x=X[mask]
    y=Y[mask]

    if(ddebug):
        print 'x:', x
        print 'y:', y
        print 'z:', z

# ------------------------------------------------------    
###    fig,ax=plt.subplots(figsize=(5, 8)) 
###    gs = gridspec.GridSpec(2,1,height_ratios=[3,1])  
# ------------------------------------------------------    

#--- Plot of the ArrayMap   
###    fig, ax1 = plt.subplots(figsize=(5, 6))
###    ax1=plt.subplot(gs[0])
###    plt.plot(X,Y,"k.",markersize=10,markeredgecolor='w') # missing channels
###    plt.plot(x,y,"k.",markersize=10) # positions of the electrodes in the map 

###[move afterwards to keep the marker in front]
###    plt.scatter(X,Y,s=error*100, c='k',edgecolor='w')
###    plt.scatter(x,y,s=error*100, c='k')
#(superimposed)

    ArrayStep = 0.55
    step = 0.55/2
    GridStep = 0.55/10
    e = GridStep/2 #(epsilon)

# ------------------------------------------------------
#--- Plot the borders of CorticalAreas

#*** MOTOR *** (polygon)
    coord = [[0.-step+e,0.-step+e], [0.55+step-e,0.-step+e], [0.55+step-e,0.+step-1.5*e], 
             [0.+step-1.5*e,0.+step-1.5*e], [0.+step-1.5*e,1.65+step-e], [0.-step+e,1.65+step-e]]
    coord.append(coord[0]) #repeat the first point to create a 'closed loop'
    xs, ys = zip(*coord) #create lists of x and y values
    plt.plot(xs,ys,color='k',linewidth=2)

#*** SOMATOSENSORY *** (polygon)
    coord = [[0.55-step+0.5*e,0.55-step+e], [1.1+step-e,0.55-step+e],
             [1.1+step-e,1.1+step+e], [1.65+step,1.65-step+e],
             [1.65+step,1.65+step-e], [0.55-step+0.5*e,1.65+step-e]]
    coord.append(coord[0]) #repeat the first point to create a 'closed loop'
    xs, ys = zip(*coord) #create lists of x and y values
    plt.plot(xs,ys,color='k',linewidth=2)

#*** PARIETAL *** (rectangle)
    ax1.add_patch(
        patches.Rectangle(
            (0.55-step+0.5*e, 2.2-step+e),   # (x,y)
            3*0.55-0.5*e,          # width
            1*0.55,          # height
            fill=False,
            linewidth=2
        ))

#*** RETROSPLENIAL *** (rectangle)
    ax1.add_patch(
        patches.Rectangle(
            (0.-step+e, 2.2-step+e),   # (x,y)
            1*(0.55-2.5*e),          # width
            4*0.55,          # height
            fill=False,
            linewidth=2
        ))

#*** VISUAL *** (polygon)
    coord = [[0.55-step+0.5*e,2.75-step+3*e], [1.65+step,2.75-step+3*e], 
             [1.65+step,2.75+step+3*e], [2.2+step-e,2.75+step+3*e], [2.2+step-e,3.3+step+e], 
             [1.65+step,3.3+step+e], [1.65+step,4.4+step+e], [0.55-step+0.5*e,4.4+step+e]]
    coord.append(coord[0]) #repeat the first point to create a 'closed loop'
    xs, ys = zip(*coord) #create lists of x and y values
    plt.plot(xs,ys,color='k',linewidth=2) 

# ------------------------------------------------------
#--- Labels
    ax1.set_xlabel('X [mm]')
    ax1.set_ylabel('Y [mm]')
    ax1.set_title('Electrode Array')
# ------------------------------------------------------    
 
#--- INTERPOLATOR
#    interp = scipy.interpolate.Rbf(x,y,z, function='thin_plate')
#[THIS]    interp = scipy.interpolate.Rbf(x,y,z, function='linear')
    interp = scipy.interpolate.Rbf(x,y,z, function='linear')

#--- MESH-GRID (denser than the electrode array grid)
#[NB] np.mgrid = construct a meshgrid
#--- if the step length is a **complex number** (e.g. 5j), then the 
#--- integer part of its magnitude is interpreted as specifying the
#--- number of points to create between the start and stop values, 
#--- where the stop value **is inclusive**.

#    yi, xi = np.mgrid[0:1:100j, 0:1:100j] # (the step length is a complex number)
#    yi, xi = np.mgrid[0:1.5:100j, 0:3:200j] # (the step length is a complex number)
    xi, yi = np.mgrid[min(X)-275e-3:max(X)+275e-3:50j, min(Y)-275e-3:max(Y)+275e-3:90j] 
#[DEBUG]print min(X),max(X),min(Y),max(Y)
#[NB] [50,90][25,45][15,27][10,18][5,9]

#--- INTERPOLATED VALUES (at the points of the MESH-GRID)
#[THIS]    yi, xi = np.mgrid[min(X)-275e-3:max(X)+275e-3:50j, min(Y)-275e-3:max(Y)+275e-3:90j] 
    zi = interp(xi,yi) # function evaluation
# --- interpolate at points that are on a regular grid. 

    for i in range(zi.shape[0]):
        for j in range(zi.shape[1]):
            if(xi[i][j]>(0.55+step) and yi[i][j]<(0.55-step)):
                #print i,j
                zi[i][j]=np.nan
            if(xi[i][j]>(1.1+step) and yi[i][j]<(1.65-step)):
                zi[i][j]=np.nan
            if(xi[i][j]>(1.65+step) and yi[i][j]<(3.3-step)):
                zi[i][j]=np.nan
            if(xi[i][j]>(1.65+step) and yi[i][j]>(3.3+step)):
                zi[i][j]=np.nan
            if(xi[i][j]<(0.55-step) and yi[i][j]>(3.85+step)):
                zi[i][j]=np.nan


#--- CONTOUR PLOT
    plt.contourf(zi.T, extent=[-0.275, 2.475, -0.275, 4.75], cmap='gist_earth')

#    plt.contourf(zi, extent=[-0.275, 2.5, -0.275, 4.75], cmap='magma')
#    plt.imshow(zi, extent=[-0.5, 2.5, -0.5, 4.75], cmap='magma')
#    plt.imshow(zi, extent=[-0.275, 2., -0.275, 4.75], cmap='magma')

#[THIS]plt.contourf(zi, extent=[-0.275, 2.5, -0.275, 4.75], cmap='magma')

    cbar=plt.colorbar()
    ###cbar.ax.set_yticklabels(['0','1','2','>3'])
    ###cbar.set_label(name+unit, rotation=270)    
    cbar.ax.set_title(unit)
    
    plt.scatter(X,Y,s=error*100, c='k',edgecolor='w')
    plt.scatter(x,y,s=error*100, c='k')


#--- To exclude non-physical regions in the mouse cortex --> 
#--- GRID = concatenation of rectangular grids

#x1, y1 = np.mgrid[min(X)-275e-3:0.55+275e-3:21j, min(Y)-275e-3:3.85+275e-3:81j]
#x2, y2 = np.mgrid[0.55-275e-3:1.65+275e-3:31j, 4.4-275e-3:4.4+275e-3:11j]
#x3, y3 = np.mgrid[0.55-275e-3:1.65+275e-3:31j, 4.4-275e-3:4.4+275e-3:11j]
#z1 = interp(x1,y1) # function evaluation
#z2 = interp(x2,y2)

# --- FileName and VariableName in the canvas

    if FileName=='messageA':
        plt.figtext(0.05,0.97,'SummaryData',color='blue',weight='roman',size='medium',
                bbox=dict(facecolor='none', edgecolor='b', boxstyle='round'))
        plt.figtext(0.05,0.93,'(averaged across experiments)',
                    color='blue',weight='roman',size='small')
        note='';

    elif FileName=='messageB':
        plt.figtext(0.05,0.97,'NormData',color='blue',weight='roman',size='medium',
                bbox=dict(facecolor='none', edgecolor='b', boxstyle='round'))
        plt.figtext(0.05,0.93,'(normalized across channels and averaged across experiments)',
                    color='blue',weight='roman',size='small')
        note='_norm';

    else:
        plt.figtext(0.05,0.97,FileName,color='blue',weight='roman',size='medium',
                bbox=dict(facecolor='none', edgecolor='b', boxstyle='round'))
        note='';

    plt.figtext(0.75,0.97,name,color='red',weight='roman',size='medium',
                bbox=dict(facecolor='none', edgecolor='r', boxstyle='round'))

# ------------------------------------------------------    
###    ax2=plt.subplot(gs[1])
###    plt.hist(variable)
#--- Labels
###    ax2.set_xlabel(name+unit)
###    ax2.set_ylabel('N')
    #ax2.set_title('Histogram of mean values at the electrode positions',fontsize=10)

# ------------------------------------------------------                    
    ###fig.tight_layout()
    plt.show(0)

##    if ~debug:
##        plt.savefig(OutputDir + name + note + '_ContourPlot.pdf')
 
    return(ax1)
    ##return(fig,ax1)


# GRID
#    x=np.arange(5) # OR x = np.linspace(0, 4, 5)
#    y=np.arange(9) # OR y=np.linspace(0, 8, 9)
#    x,y=np.meshgrid(x,y)

#-----------------------------------------------------
def BoxPlotArea(variable,yLabel):
    "Produce a BoxPlot with the distribution of the variable for each Cortical Area"
    
    data_to_plot = []
    meanA=[]
    stdA=[]
    stderrA=[]
    nValidCh=[]

    outliers_reShape = outliers-1 #(array indexing different in Matlab and in Python)

    for i in area:
        A=[]
        nChA=0  # Number of channels per area
        for j in (i[1]):
            if j not in outliers_reShape:
                A=np.append(A,[datastruct[0,j][variable]])
                nChA=nChA+1
                #print(len(A))
            else:
                A=np.append(A,[])
        data_to_plot.append(A)
        nValidCh.append(nChA)
        meanA.append(np.mean(A))
        stdA.append(np.std(A))
        stderrA.append(np.std(A)/math.sqrt(len(A)))

    #print(meanA)

    M = np.empty((0,len(data_to_plot)))
    m = np.empty((0,len(data_to_plot)))
    #SEM = np.empty((0,len(data_to_plot)))
    for i in range(0,len(data_to_plot)):
        M = np.append(M, np.median(data_to_plot[i]))
        m = np.append(m, np.mean(data_to_plot[i])) # (m and meanA are the same)
        #SEM=np.append(SEM,stats.sem((data_to_plot[i])) # standar error of the mean

    #print M
    #print len(M)
    #print nValidCh

    fig, ax1 = plt.subplots(figsize=(10, 6))
    #fig.canvas.set_window_title(variable+' - BoxPlot - CorticalAreas')
    plt.subplots_adjust(left=0.075, right=0.95, top=0.9, bottom=0.1)

    meanpointprops = dict(marker='*',markerfacecolor='blue',markeredgecolor='blue')
    plt.boxplot(data_to_plot,notch=True,sym='+',showmeans=True,meanprops=meanpointprops) # BOXPLOT (with mean values)
    #plt.plot([1,2,3,4,5],meanA,linestyle='none',color='b',marker='*')
    #plt.errorbar([1,2,3,4,5],meanA,yerr=stderrA,fmt='o',color='b',marker='*') #stderr SMALL (large sample)
    #plt.errorbar([1,2,3,4,5],meanA,yerr=stdA,fmt='o',color='b',marker='*') 
    # overplot the sample averages (markeredgecolor='k')
    plt.xlabel('Cortical Area')

    #plt.xticks([1,2,3,4,5], ['M','S','P','R','V'])
    #plt.xticks([1,2,3,4,5], ['Motor','Somatosensory','Parietal','Retrosplenial','Visual'])
    #xtickNames = plt.setp(ax1, xticklabels=['Motor','Somatosensory','Parietal','Retrosplenial','Visual'])
    TickLabels = ['Motor','Somatosensory','Parietal','Retrosplenial','Visual']
    for i in range(0,5):
        TickLabels[i] = TickLabels[i]+' (nCh='+str(nValidCh[i])+')'
    xtickNames = plt.setp(ax1,xticklabels=TickLabels)

    #plt.setp(xtickNames,rotation=45,fontsize=8)
    plt.setp(xtickNames,fontsize=8)

    plt.ylabel(yLabel)

    ml = MultipleLocator(0.2)
    ax1.yaxis.set_minor_locator(ml)
    ax1.yaxis.grid(True, linestyle='-', which='minor', color='lightgrey', alpha=0.5)
    
    # Legend
    plt.figtext(0.85, 0.92, '*', color='b', weight='roman', size='medium')
    plt.figtext(0.865, 0.92, 'Mean Value', color='black', weight='roman', size='medium',
                bbox=dict(facecolor='none', edgecolor='k', boxstyle='round'))

    # FileName in the canvas
    plt.figtext(0.05, 0.92, FileName, color='blue', weight='roman', size='medium',
                bbox=dict(facecolor='none', edgecolor='b', boxstyle='round'))

    plt.show(0)
    plt.savefig(OutputDir + variable +'_BoxPlot_Areas.pdf')

    return(M,m,stdA,stderrA)
                      
#-----------------------------------------------------
def SummaryBoxPlotArea(variable,varname,yLabel,option):
    "Produce a BoxPlot with the distribution of values of the Variable across Experiments for each Cortical Area"

    fig, ax1 = plt.subplots(figsize=(10, 6))
    #fig.canvas.set_window_title(variable+' - BoxPlot - CorticalAreas')
    plt.subplots_adjust(left=0.075, right=0.95, top=0.9, bottom=0.1)

    meanpointprops = dict(marker='*',markerfacecolor='blue',markeredgecolor='blue')
    #plt.boxplot(variable,notch=True,sym='+',showmeans=True,meanprops=meanpointprops) # BOXPLOT (with mean values)
    #plt.boxplot(variable,notch=False,sym='+',showmeans=True,meanprops=meanpointprops) # BOXPLOT (with mean values)
    plt.boxplot(variable,notch=option,sym='+',showmeans=True,meanprops=meanpointprops) # BOXPLOT (with mean values)
    #plt.plot([1,2,3,4,5],meanA,linestyle='none',color='b',marker='*')
    #plt.errorbar([1,2,3,4,5],meanA,yerr=stderrA,fmt='o',color='b',marker='*') #stderr SMALL (large sample)
    #plt.errorbar([1,2,3,4,5],meanA,yerr=stdA,fmt='o',color='b',marker='*') 
    # overplot the sample averages (markeredgecolor='k')
    plt.xlabel('Cortical Area')

    # --- SCATTER PLOT superimposed
    #for i in range(5):
    #    y = variable[:,i]
    #    x = np.random.normal(1+i, 0.04, size=len(y))
    #    plt.plot(x, y, 'r.', alpha=0.6)

    #plt.xticks([1,2,3,4,5], ['M','S','P','R','V'])
    #plt.xticks([1,2,3,4,5], ['Motor','Somatosensory','Parietal','Retrosplenial','Visual'])
    #xtickNames = plt.setp(ax1, xticklabels=['Motor','Somatosensory','Parietal','Retrosplenial','Visual'])
    TickLabels = ['Motor','Somatosensory','Parietal','Retrosplenial','Visual']
    xtickNames = plt.setp(ax1,xticklabels=TickLabels)

    #plt.setp(xtickNames,rotation=45,fontsize=8)
    plt.setp(xtickNames,fontsize=8)

    plt.ylabel(yLabel)

    ml = MultipleLocator(0.2)
    ax1.yaxis.set_minor_locator(ml)
    ax1.yaxis.grid(True, linestyle='-', which='minor', color='lightgrey', alpha=0.5)

    # Legend
    plt.figtext(0.85, 0.92, '*', color='b', weight='roman', size='medium')
    plt.figtext(0.865, 0.92, 'Mean Value', color='black', weight='roman', size='medium',
                bbox=dict(facecolor='none', edgecolor='k', boxstyle='round'))

    # FileName in the canvas
    #plt.figtext(0.05, 0.92, 'nExp = 3', color='blue', weight='roman', size='medium',
    #            bbox=dict(facecolor='none', edgecolor='b', boxstyle='round'))

    plt.show(0)
    #plt.savefig(PythonDir + variable +'_BoxPlot_Areas.pdf')
    if option:
        plt.savefig(BaseDir + PythonDir + varname +'_BoxPlot_Areas_notches.pdf')
    else:
        plt.savefig(BaseDir + PythonDir + varname +'_BoxPlot_Areas.pdf')

#-----------------------------------------------------
def SummaryBoxPlotAreaNEW(variable,varname,yLabel,ax):
    "Produce a BoxPlot with the distribution of values of the Variable across Experiments for each Cortical Area"

    ##fig, ax1 = plt.subplots(figsize=(10, 6))
    ##fig.canvas.set_window_title(variable+' - BoxPlot - CorticalAreas')
    ##plt.subplots_adjust(left=0.075, right=0.95, top=0.9, bottom=0.1)

    meanpointprops = dict(marker='*',markerfacecolor='blue',markeredgecolor='blue')
    # BOXPLOT (with mean values)
    plt.boxplot(variable,sym='+',showmeans=True,meanprops=meanpointprops) 
    plt.xlabel('Cortical Area')

    # --- SCATTER PLOT superimposed
    #for i in range(5):
    #    y = variable[:,i]
    #    x = np.random.normal(1+i, 0.04, size=len(y))
    #    plt.plot(x, y, 'r.', alpha=0.6)
    
    TickLabels = ['Motor','Somatosensory','Parietal','Retrosplenial','Visual']
    xtickNames = plt.setp(ax,xticklabels=TickLabels)
    plt.setp(xtickNames,fontsize=8)

    plt.ylabel(yLabel)

    ml = MultipleLocator(0.2)
    ax.yaxis.set_minor_locator(ml)
    ax.yaxis.grid(True, linestyle='-', which='minor', color='lightgrey', alpha=0.5)

    # Legend
    plt.figtext(0.85, 0.92, '*', color='b', weight='roman', size='medium')
    plt.figtext(0.865, 0.92, 'Mean Value',
                color='black', weight='roman', size='medium',
                bbox=dict(facecolor='none', edgecolor='k', boxstyle='round'))

    return(ax)

    ##plt.show(0)
    ##plt.savefig(PythonDir + variable +'_BoxPlot_Areas.pdf')
    ##if option:
    ##    plt.savefig(BaseDir + PythonDir + varname +'_BoxPlot_Areas_notches.pdf')
    ##else:
    ##    plt.savefig(BaseDir + PythonDir + varname +'_BoxPlot_Areas.pdf')

#-----------------------------------------------------


def SummaryErrorPlotArea(variable,varname,yLabel):
    "Produce a Error Plot of the meanValue of the Variable across Experiments for each Cortical Area"
    x=range(variable.shape[1])
    y=[]
    err=[]
    for j in range(variable.shape[1]):
        #print('***************')
        A=[]
        for i in range(variable.shape[0]):
        #print(variable[i][j])
            A.append(variable[i][j])
        #print(np.mean(A))
        y.append(np.mean(A))
        err.append(stats.sem(A,ddof=0))
        #err.append(np.std(A)/math.sqrt(variable.shape[0]))
    #print(y)
    #print(err)

    fig, ax1 = plt.subplots(figsize=(10, 6))
    #fig.canvas.set_window_title(variable+' - BoxPlot - CorticalAreas')
    plt.subplots_adjust(left=0.075, right=0.95, top=0.9, bottom=0.1)
    plt.errorbar(x,y,xerr=0,yerr=err,fmt='o',capsize=5) #fmt='.'
    plt.xlabel('Cortical Area')
    ax1.set_xticks(range(0,5))
    TickLabels = ['Motor','Somatosensory','Parietal','Retrosplenial','Visual']
    xtickNames = plt.setp(ax1,xticklabels=TickLabels)
    plt.setp(xtickNames,fontsize=8)
    plt.ylabel(yLabel)
   
    #plt.ylim((25,250))

    #ml = MultipleLocator(0.1*abs(max(y)-min(y)))
    #ax1.yaxis.set_minor_locator(ml)
    #ax1.yaxis.grid(True, linestyle='-', which='minor', color='lightgrey', alpha=0.5)

    plt.figtext(0.025, 0.95, 'Mean Value of the variable across the 11 experiments, for each area. The error is the s.e.m across experiments.',fontsize=9)

    plt.show(0)
    plt.savefig(BaseDir + PythonDir + varname +'_ErrorPlot_Areas.pdf')

#-----------------------------------------------------
def StemPlotCh(variable,varname,yLabel):
    "Produce a StemPlot of values per channel"
    ChList = np.empty((0,nCh))
    Separator=np.array([0.])

    for i in area:
        Separator=np.append(Separator,Separator[-1]+len(i[1]))
        for j in i[1]:
            ChList=np.append(ChList,str(j+1))  # The labels for the bottom x axis
    #print(ChList)
    #print(Separator)

    fig, ax1 = plt.subplots(figsize=(10, 6))
    #fig.canvas.set_window_title(variable+' - BoxPlot - CorticalAreas')
    plt.subplots_adjust(left=0.075, right=0.95, top=0.9, bottom=0.1)

    plt.stem(variable)
    ax1.set_xlabel('Channel') # *** x-axis --- Label
    plt.ylabel(yLabel)

    #separator=Separator[1:]
    #width=[]
    height=ax1.get_ylim()[1]-ax1.get_ylim()[0]
    TickPos = []
    offset=-0.5
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


    # *** secondary x-axis
    labels = ['Motor','Somatosensory','Parietal','Retrosplenial','Visual'] # The labels for the top x axis
    TickPosNorm=[i/nCh for i in TickPos]
    ax2 = ax1.twiny()
    ax2.set_xticks(TickPosNorm)
    ax2.set_xticklabels(labels)
    #ax2.set_xlabel('Cortical Area')

    ax1.set_xlim(-0.5,31.5)
    ax1.set_xticks(range(0,nCh))
    ax1.set_xticklabels(ChList)

    # FileName in the canvas
    plt.figtext(0.05, 0.02, FileName, color='blue', weight='roman', size='medium',
                bbox=dict(facecolor='none', edgecolor='b', boxstyle='round'))
     
    # *** mean, Median and percentiles in the canvas (exclude entries == nan)
    mask = np.ones(len(variable), dtype=bool)
    
    for i in range(len(variable)):
        if np.isnan(variable[i]):
            mask[[i]] = False
    var=variable[mask]

    N=len(var)
    m=np.mean(var)
    SD=np.std(var)
    SEM=stats.sem(var,ddof=0)
    M=np.median(var)
    Q1=np.percentile(var,25)
    Q3=np.percentile(var,75)
    plt.figtext(0.85, 0.85, 'N=%s\nmean=%s\nSD=%s\nSEM=%s\nmedian=%s\nQ1=%s\nQ3=%s'%(N,m,SD,SEM,M,Q1,Q3),fontsize=8,
                bbox=dict(facecolor='white', edgecolor='black', boxstyle='round'))

#pylab.title('Minimal Energy Configuration of %s Charges on Disc W = %s'%(N, W))

    plt.show(0)

    plt.savefig(OutputDir + varname +'_StemPlotCh.pdf')

#-----------------------------------------------------
def PlotCh(variable,err,varname,yLabel):
    "Produce a Plot of values per channels and with errors (mean values + SEM)"

    ChList = np.empty((0,nCh))
    Separator=np.array([0.])

    for i in area:
        Separator=np.append(Separator,Separator[-1]+len(i[1]))
        for j in i[1]:
            ChList=np.append(ChList,str(j+1))  # The labels for the bottom x axis
    #print(ChList)
    #print(Separator)

    fig, ax1 = plt.subplots(figsize=(10, 6))
    #fig.canvas.set_window_title(variable+' - BoxPlot - CorticalAreas')
    plt.subplots_adjust(left=0.075, right=0.95, top=0.9, bottom=0.1)

    x = range(nCh)
    plt.errorbar(x,variable,xerr=0,yerr=err,fmt='o',capsize=5) #fmt='.'
    ax1.set_xlabel('Channel') # *** x-axis --- Label
    plt.ylabel(yLabel)

    #separator=Separator[1:]
    #width=[]
    height=ax1.get_ylim()[1]-ax1.get_ylim()[0]
    TickPos = []
    offset=-0.5
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


    # *** secondary x-axis
    labels = ['Motor','Somatosensory','Parietal','Retrosplenial','Visual'] # The labels for the top x axis
    TickPosNorm=[i/nCh for i in TickPos]
    ax2 = ax1.twiny()
    ax2.set_xticks(TickPosNorm)
    ax2.set_xticklabels(labels)
    #ax2.set_xlabel('Cortical Area')

    ax1.set_xlim(-0.5,31.5)
    ax1.set_xticks(range(0,nCh))
    ax1.set_xticklabels(ChList)

    # FileName in the canvas
    plt.figtext(0.05, 0.02, FileName, color='blue', weight='roman', size='medium',
                bbox=dict(facecolor='none', edgecolor='b', boxstyle='round'))
     
    # *** mean, Median and percentiles in the canvas (exclude entries == nan)
    mask = np.ones(len(variable), dtype=bool)
    for i in range(len(variable)):
        if np.isnan(variable[i]):
            mask[[i]] = False
    var=variable[mask]

    N=len(var)
    m=np.mean(var)
    SD=np.std(var)
    SEM=stats.sem(var,ddof=0)
    M=np.median(var)
    Q1=np.percentile(var,25)
    Q3=np.percentile(var,75)
    plt.figtext(0.825, 0.85, 'N=%s\nmean=%s\nSD=%s\nSEM=%s\nmedian=%s\nQ1=%s\nQ3=%s'%(N,m,SD,SEM,M,Q1,Q3),fontsize=8,
                bbox=dict(facecolor='white', edgecolor='black', boxstyle='round'))

#pylab.title('Minimal Energy Configuration of %s Charges on Disc W = %s'%(N, W))

    plt.show(0)

    plt.savefig(OutputDir + varname +'_PlotCh.pdf')

#-----------------------------------------------------
def StemPlotArea(variable,varname,yLabel):
    "Produce a StemPlot of values per area"

    fig, ax1 = plt.subplots(figsize=(10, 6))
    #fig.canvas.set_window_title(variable+' - BoxPlot - CorticalAreas')
    plt.subplots_adjust(left=0.075, right=0.95, top=0.9, bottom=0.1)

    plt.stem(variable)
    ax1.set_xlabel('Cortical Area') # *** x-axis --- Label
    plt.ylabel(yLabel)

    ax1.set_xticks(range(0,5))
    TickLabels = ['Motor','Somatosensory','Parietal','Retrosplenial','Visual']
#    for i in range(0,5):
#        TickLabels[i] = TickLabels[i]+' (nCh='+str(nValidCh[i])+')' 
    xtickNames = plt.setp(ax1,xticklabels=TickLabels)
    plt.setp(xtickNames,fontsize=8)

    # FileName in the canvas
    plt.figtext(0.05, 0.92, FileName, color='blue', weight='roman', size='medium',
                bbox=dict(facecolor='none', edgecolor='b', boxstyle='round'))

    plt.show(0)
    plt.savefig(OutputDir + varname +'_StemPlotAreas.pdf')

#-----------------------------------------------------
def PlotArea(variable,err,varname,yLabel):
    "Produce a Plot of values per area and with errors (mean values + SEM)"

    fig, ax1 = plt.subplots(figsize=(10, 6))
    #fig.canvas.set_window_title(variable+' - BoxPlot - CorticalAreas')
    plt.subplots_adjust(left=0.075, right=0.95, top=0.9, bottom=0.1)

    x = range(5)
    plt.errorbar(x,variable,xerr=0,yerr=err,fmt='o',capsize=5) #fmt='.'
    ax1.set_xlabel('ChannelCortical Area') # *** x-axis --- Label
    plt.ylabel(yLabel)

    ax1.set_xticks(range(0,5))
    TickLabels = ['Motor','Somatosensory','Parietal','Retrosplenial','Visual']
#    for i in range(0,5):
#        TickLabels[i] = TickLabels[i]+' (nCh='+str(nValidCh[i])+')' 
    xtickNames = plt.setp(ax1,xticklabels=TickLabels)
    plt.setp(xtickNames,fontsize=8)

    # FileName in the canvas
    plt.figtext(0.05, 0.92, FileName, color='blue', weight='roman', size='medium',
                bbox=dict(facecolor='none', edgecolor='b', boxstyle='round'))

    plt.show(0)
    plt.savefig(OutputDir + varname +'_PlotAreas.pdf')

#-----------------------------------------------------
def GetSlope(interval):#,varname,yLabel):
    "Compute the slope of the average transition in the specified interval"

    X=[];Y=[]
    for i in area:
        for j in i[1]:
            #print j
            if j not in outliers-1: #(array indexing different in Matlab and in Python)
                Y.append(np.squeeze(datastruct[0,j]['MeanY']))
                X.append(np.squeeze(datastruct[0,j]['t']))
            else:
                X.append([])
                Y.append([])
    
    peak=[]; z=[]
    for j in range(nCh):
        if len(X[j]) is not 0:
            xx=[]; yy=[]; zz=[]; yyy=[]
### xx and yy subset for the computation of the slope
### yyy subset for the computation of the peak
            for i in range(len(X[j])):
                if (X[j][i]>=0): # select time range to find the Peak
                    yyy.append(Y[j][i])
                if (X[j][i]>=interval[0]) and (X[j][i]<=interval[1]): # select time interval and samples
                    xx.append(X[j][i])
                    yy.append(Y[j][i])
            peak.append(max(yyy)) # find the Peak
            zz = np.polyfit(xx, yy, 3) # fit the transition with a cubic
            z.append(zz)
        else:
            peak.append(np.nan)
            z.append([])


    slope=np.empty((nCh))
    for i in range(len(z)):
        if len(z[i])==4:
            slope[i]=(z[i][2]) # SLOPE = derivative of the fitted transition computed at x = 0
        else:
            slope[i]=np.nan

    peak=np.array(peak)

    return(slope,peak)

#-----------------------
#    X=[];Y=[]

#    for i in area:
#        for j in i[1]:
#            #print j
#            if j not in outliers-1: #(array indexing different in Matlab and in Python)
#                Y.append([datastruct[0,j]['MeanY']])
#                X.append([datastruct[0,j]['t']])
#            else:
#                X.append([])
#                Y.append([])
    
#    x=[];y=[];z=[];peak=[]
#    for j in range(nCh):
#        if j not in outliers-1:
#            x.append(np.array(np.array(X[j][0][0])))
#            y.append(np.array(np.array(Y[j][0][0])))
#            xx=[];yy=[];yyy=[];zz=[]
#            for i in range(len(np.squeeze(X[j]))):
#                if (x[j][i]>=0): # select time range to find the Peak
#                    #print(i,x[j][i])
#                    yyy.append(y[j][i])
#                if (x[j][i]>=interval[0]) and (x[j][i]<=interval[1]): # select time interval and samples
#                    xx.append(x[j][i])
#                    yy.append(y[j][i])
#            peak.append(max(yyy)) # find the Peak
#            zz = np.polyfit(xx, yy, 3) # fit the transition with a cubic
#            z.append(zz)
#        else:
#            x.append([])
#            y.append([])
#            z.append([])
#            peak.append(np.nan)

#    slope=np.empty((nCh))
#    for i in range(len(z)):
#        if len(z[i])==4:
#            slope[i]=(z[i][2]) # SLOPE = derivative of the fitted transition computed at x = 0
#        else:
#            slope[i]=np.nan

#    peak=np.array(peak)

#    return(slope,peak)

#-----------------------------------------------------
def GetAreaInfo(variable):#,varname,yLabel):
    "Get the median, mean, SD and SED per area"

    M = np.empty(len(area))
    m = np.empty(len(area))
    std =  np.empty(len(area))
    SEM = np.empty(len(area))
    
    for i in range(len(area)):
        if i==0:
            start=0
        else:
            start=offset
        end=start+len(area[i][1])
        #print(start,end)
        
        subset = variable[start:end]
        mask = np.ones(len(subset), dtype=bool)
        for j in range(len(subset)):
            if np.isnan(subset[j]):
                mask[[j]] = False
        var=subset[mask]

        M[i]=np.median(var)
        m[i]=np.mean(var)
        std[i]=np.std(var)
        SEM[i]=stats.sem(var,ddof=0)
        offset=end

    return(M,m,std,SEM)
