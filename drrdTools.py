def med2tec(fNAME,flag='A'):
    # written by Marcelo Bussotti Reyes - Universidade Federal do ABC
    # Based on the Matlab code writen by Marcelo Caetano.
    # 2019-01-10

    TIME_RESOLUTION = 2e-3 # (in seconds) Time resolutuion in the med associates box = 2ms

    try:
        fid = open(fNAME, 'r')
    except IOError:
        print("Could not read file:", fNAME)
        return([])

    ##cans the entire line in as a string, with spaces (to preserve number info)
    fileString = fid.readlines() 

    # looking for the line positions that contain the string
    # notice that the range starts in 1 (not in zero) avoiding the first line 
    # (the file name is always the first line)
    indexes = [i for i in range(1,len(fileString)) if flag in fileString[i]]

    #checks if letter is at the beginning of the string and posceded by a :
    for i in indexes:
        if fileString[i].index(flag) == 0 and fileString[i][1] == ':':
            startParse = i+1

    # now let's look for the next letter of end of file to check how many lines we should read
    stopParse = -1                                   # just a code in case it finds nothing
    for i in range(startParse+1,len(fileString)):    # for from startParse until end of file
        if str.isalpha(fileString[i][0]) and fileString[i][1]==':': # looks if the first element of the line is alphabetic and the next is :
            stopParse = i                            # saves the last position that did not start with a letter
            break 
    # checking if did not find any letters, then just record the last line of the file
    if stopParse == -1:
        stopParse = len(fileString)


    M = []
    for i in range(startParse,stopParse):
        for j in str.split(fileString[i])[1:]:
            data = str.split(j,'.')
            M.append([round(int(data[0])*TIME_RESOLUTION,3), int(data[1])])
    
    fid.close()

    return(M)

#_____________________________________________________________________

def plotDrrd(D,title_label):
    # ___________________________________________________________________________________________
    # File:             plotDrrd.py
    # File type:        Function
    # Created on:       January, 14, 2019
    # Created by:       Marcelo Bussotti Reyes
    # Purpose:          A function that analyzes the performance of rats in the drrd
    #                   procedure by plotting the trials 
    #
    # Input:            A matrix with 6 columns of data according to 
    #                   the one produced by the function drrd.
    #                   Each line of the matrix D is a trial
    #                   Column 0 is the duration of the lever press
    #                   Column 1 is the the time between the lever release and the next lever press (ITI)
    #                   Column 2 is 1 for the reinforced trials
    #                   Column 3 is 1 for trials where the light was on (valid trials)
    #                   Column 4 shows the criterion (prime time) for each trial
    #                   Column 5 is the session number
    #
    # Output:           A plot of the data contained in D and some statistics
    #
    #
    # Format:           plotDrrd(D,'tittle label')
    # Example:          d = plotDrrd(D,'191014-Session1'); 

    #if nargin == 1
    #  title_label = [];
    #end

    import numpy as np
    import matplotlib.pyplot as plt


    #D = np.load('D.npy')

    primed  = 2
    valid   = 3
    primeT  = 4
    session = 5 	# column with the session number

    N = D.shape[0]

    # --- looking for the specific trials ---
    #validPrimed    =  find(D(:,primed)==1 	& D(:,valid)==1);
    validPrimed    = [i for i in range(N) if (D[i,primed]==1 and D[i,valid]==1)]
    validNonPrimed = [i for i in range(N) if (D[i,primed]==0 and D[i,valid]==1)]
    invalid        = [i for i in range(N) if D[i,valid]==0] 

    # --- plotting the prime times ---
    plt.plot(D[:,primeT],range(N))#,'r','linewidth', 1.5);

    # --- alternative: patch ---
    #patch([ D(:,5); D(end,5); 0.00; 0.00], [1:N N+5 N+5 0], [.7 .8 .7] ,'EdgeColor' ,'none');% % [.7 .8 .7]
    #patch([ D(:,5); D(end,5); 0.00; 0.00], [1:N N+5 N+5 0], [0 110 144]/255 ,'EdgeColor' ,'none');% % [.7 .8 .7]
    #patch([ D(:,5); D(end,5); 0.00; 0.00], [1:N N+5 N+5 0], [0.8 0.8 0.8] ,'EdgeColor' ,'none');% % [.7 .8 .7]

    # --- Plotting each trial in a different style ---
    #mycolor = [0,.8,.9]
    mysize  = 30
    lw      = 0.5
    plt.scatter(D[validPrimed,0]   ,validPrimed,   s=mysize,  linewidths=lw, marker='o',c='w', edgecolors='k' )
    plt.scatter(D[validNonPrimed,0],validNonPrimed,s=mysize,  linewidths=lw, marker='o',c='w', edgecolors='k' )
    plt.scatter(D[invalid,0]       ,       invalid,s=mysize,  linewidths=lw, marker='o',c='w', edgecolors='r' )


    # --- printing the lines dividing the sessions ---
    div = np.where(np.diff(D[:,session]))[0]
    for i in range(len(div)):
        plt.plot(plt.xlim,[div[i],div[i]])
        
    # --- Plotting the moving average of the lever press durations ---
    #plot(movingAverage(D(:,1),20),1:N,'linewidth',2);


    # --- setting up the scale and title ---
    plt.xlim((0,5))
    plt.ylim((0,N+5))
    plt.xlabel('time (s)',fontsize=20)
    plt.ylabel('trial',fontsize=20)
    plt.title(title_label,fontsize=22)

    plt.show()
    plt.ion()
    # --- mounting return variable ---
    return ([len(validPrimed)/N, len(validNonPrimed)/N, len(invalid)/N] *100)

#___________________________________________________________________________
def drrd(prefix='AB1',animalID = 64,session  = 1, plotFlag = True):
    # Created on:       January, 14, 2019
    # Created by:       Marcelo Bussotti Reyes
    # Purpose:          A function that analyzes the performance of rats in the drrd
    #                   procedure. The animals are supposed to press a lever for a
    #					duration longer than the prime time in order to receive food.
    #
    # Input:            File name prefix as a string, animal id as a number and 
    #					session as a number. Filene name are in the format:
    #					[prefix00A.00S] where A is the animal and S is the sesson 
    #
    # Output:           D, a matrix with 6 columns containing data from a trial in 
    #					each line.
    #                   In case saveMatFlag is parsed as true, the program will
    #                   save a matlab file (.mat) with the same name as the
    #                   original data file.
    # Coments:          Uses functions: med2tec.m, and plotDrrd.m
    #
    # Format:           drrd(prefix,animalID,session))
    # Example:          D = drrd('AB1',1,2'); 
    #					this will analyze the file AB1001.002, animal 1 and session 2 
    # Previous Modifications:
    # Modification:     Included the option of saving a matlab file with the
    #                   matrix D. This helps to speed up the analysis when
    #                   several sessions and animals are to be analyzed.
    #                   Included session column (6th) and the input format
    #                   (Marcelo B. Reyes 2013)

    #from drrdTools import med2tec
    import numpy as np

    # builds file name from parsed parameters - uses zfill to pad with zeros
    filename = prefix + str(animalID).zfill(3) + '.' + str(session).zfill(3)

    # --- Indexes for each of the output variable columns
    dtCol      = 0
    itiCol     = 1
    primedCol  = 2
    validCol   = 3
    phaseCol   = 4
    sessionCol = 5              # session variable column index

    dPh        = 0.1            # initiates the variable just in case it cannot be 
                                # obtained from the data, as for example when no 
                                # trials were primed (reinforced) in one session

    data = med2tec(filename)    # reads data from medpc format to time-event code
    data = np.array(data)


    if len(data) == 0:
        print('Empty file or file not found: no data analyzed\n')
        #return([])        

    # small correction for a bug in the med-pc file ---
    # if the animal presses the lever at the same cycle of the Start command in the
    # box, the first time can be registered wrong, so this sets it to zero as it 
    # should be
    if data[0,0] > data[1,0]:
        data[1,1] = 0
        
    # --- look for indexes of temporal events ---
    #startIndex      = np.find(data[:,2]== 1)  
    startIndex    = np.array([i for i in range(len(data)) if data[i,1] == 1])
    endIndex      = np.array([i for i in range(len(data)) if data[i,1] == 3])
    primeIndex    = np.array([i for i in range(len(data)) if data[i,1] == 18]) # trials that lasted longer than criterion (primed trials)
    lightOnIndex  = np.array([i for i in range(len(data)) if data[i,1] == 11]) # when light was turned on 
    lightOffIndex = np.array([i for i in range(len(data)) if data[i,1] == 21]) # when light was turned off
    phaseAdvIndex = np.array([i for i in range(len(data)) if data[i,1] == 17]) # indexes of trials where phase was advanced
    phaseBckIndex = np.array([i for i in range(len(data)) if data[i,1] == 17]) # indexes of trials where phase was retreated

    startIndex    = startIndex[range(len(endIndex))] #eliminates the last trial in case it was incomplete

    # --- searching for trials in which the animals received food. We call these "primed" ----
    #primedTrials = findTrial(startIndex,primeIndex);
    primedTrials = np.array([startIndex[startIndex<primeIndex[i]].size-1 for i in range(len(primeIndex))])

    # --- searching for trials in which animals progressed or retreated phase
    phaseAdvTrials = np.array([startIndex[startIndex<phaseAdvIndex[i]].size-1 for i in range(len(phaseAdvIndex))])
    phaseBckTrials = np.array([startIndex[startIndex<phaseBckIndex[i]].size-1 for i in range(len(phaseBckIndex))]) 

    # --- searching for trials in which the animals responded with the light on 
    # (not in timeout). We'll call these trials "valid". 
    #validTrials     = findValidTrial(startIndex,lightOnIndex,lightOffIndex);

    if len(lightOnIndex)!= len(lightOffIndex):
        if len(lightOnIndex) == len(lightOffIndex)+1:
            #print(startIndex[-1])
            lightOffIndex = np.append(lightOffIndex,startIndex[-1]+1)
        else: 
            print('Incompatible number of events')
            #exit(-2)


    validTrials = np.array([], dtype=np.int64).reshape(0,)
    for i in range(len(lightOnIndex)):
        Nu = np.array(len(startIndex[startIndex<lightOnIndex[i] ]))
        Nv = np.array(len(startIndex[startIndex<lightOffIndex[i]]))
        #print(validTrials.size,Nu.size,Nv.size)
        validTrials = np.hstack([validTrials,range(Nu,Nv)])

    # --- search for valid trials in which animals were and were not reinforded ---
    validPrimed    = np.intersect1d(validTrials,primedTrials)
    validNonPrimed = np.setdiff1d  (validTrials,primedTrials)
    invalid        = np.setdiff1d  (range(len(startIndex)), validTrials)

    # --- gets the initial prime time ---
    if len(primeIndex) >= 1:
        iniPh = data[primeIndex[0],0] - data[primeIndex[0]-1,0]
    else:
        iniPh = 1.2
        print('Warning: initial criterion cound not be retrieved from file, assumed 1.2 s')


    D = np.zeros(len(startIndex)*6)    # initiates the vector for includind the data     

    # --- Organizing data in one single matriz: D --- 
    D = np.zeros((len(startIndex),6))             # Initiates the vector for speed

    # --- Calculating the duration of the lever presses ---
    D[:  , dtCol]               = data[endIndex,0] - data[startIndex,0]
    D[:-1,itiCol]               = data[startIndex[1:],0] - data[endIndex[:-1],0]
    D[-1 ,itiCol]               = np.nan

    if len(primedTrials)>0: D[primedTrials,primedCol]   = 1  # sets to 1 all the trials that were primed
    if len(validTrials) >0: D[validTrials , validCol]   = 1  # sets to 1 all the trials that were primed
    if len(phaseAdvTrials)>0: D[phaseAdvTrials,phaseCol]  =  dPh
    if len(phaseBckTrials)>0: D[phaseBckTrials,phaseCol]  = -dPh

    D[:,phaseCol]               = np.cumsum(D[:,phaseCol])+iniPh
    D[:,sessionCol]             = session        # adds the session number to data (same for all lines)

    # --- printing output with the summary of the results ---
    print(f'Rat{D.shape[0]}  Trials:{len(startIndex)}  Reinforced:{len(validPrimed)}  Non-Reinforced:{len(validNonPrimed)}  Invalid:{len(invalid)}\n')


    # --- graphical part ---
    if plotFlag:
        plotDrrd(D,filename)


    #np.save('D',D) # remove this
    #print(D)
    return(D)

