from sklearn import mixture
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import leastsq
from scipy.optimize import curve_fit
import drrdTools as dr
import os

# --- General Settings ----
homePath = os.path.realpath('../data/raw/pharmac_mPFC/') + '/'

## parameters for the graphical output
font = {'family': 'sans',
        'weight': 'medium',
        'size': 12}
# plt.rc('font', **font)#
# plt.rc('text', usetex=True)

# parameters for histograms and for plotting the distributions
(xmi, xma, dt) = 0, 6, 0.1
dtFine = dt / 10

# determining the bin edges for the plots
x = np.arange(xmi, xma, dt)
xfine = np.arange(xmi, xma, dtFine)


# ---------------------------------------------------


def med2tec(fNAME, flag='A'):
    # written by Marcelo Bussotti Reyes - Universidade Federal do ABC
    # Based on the Matlab code writen by Marcelo Caetano.
    # 2019-01-10

    TIME_RESOLUTION = 2e-3  # (in seconds) Time resolutuion in the med associates box = 2ms

    try:
        fid = open(fNAME, 'r')
    except IOError:
        print("Could not read file:", fNAME)
        return ([])

    ##reads the entire line in as a string, with spaces (to preserve number info)
    fileString = fid.readlines()

    # looking for the line positions that contain the string
    # notice that the range starts in 1 (not in zero) avoiding the first line 
    # (the file name is always the first line)
    indexes = [i for i in range(1, len(fileString)) if flag in fileString[i]]

    # checks if letter is at the beginning of the string and posceded by a :
    for i in indexes:
        if fileString[i].index(flag) == 0 and fileString[i][1] == ':':
            startParse = i + 1

    # now let's look for the next letter of end of file to check how many lines we should read
    stopParse = -1  # just a code in case it finds nothing
    for i in range(startParse + 1, len(fileString)):  # for from startParse until end of file
        if str.isalpha(fileString[i][0]) and fileString[i][1] == ':':  # looks if the first element of the line is alphabetic and the next is :
            stopParse = i  # saves the last position that did not start with a letter
            break
            # checking if did not find any letters, then just record the last line of the file
    if stopParse == -1:
        stopParse = len(fileString)

    M = []
    for i in range(startParse, stopParse):
        for j in str.split(fileString[i])[1:]:
            data = str.split(j, '.')
            M.append([round(int(data[0]) * TIME_RESOLUTION, 3), int(data[1])])

    fid.close()

    return (M)


def plotDrrd(D, title_label='Title_Label'):
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

    # --- checking if data variable contains any data. Otherwise exit
    if len(D) == 0:
        print('No data found, exiting\n')
        return ([])

    primed = 2
    valid = 3
    primeT = 4
    session = 5  # column with the session number
    N = D.shape[0]

    # --- looking for the specific trials ---
    validPrimed = [i for i in range(N) if (D[i, primed] == 1 and D[i, valid] == 1)]
    validNonPrimed = [i for i in range(N) if (D[i, primed] == 0 and D[i, valid] == 1)]
    invalid = [i for i in range(N) if D[i, valid] == 0]

    # activating iteractive mode
    plt.ion()

    # --- plotting the prime times ---
    plt.plot(D[:, primeT], range(N), 'grey', linewidth=2, alpha=0.5)

    # --- setting up the scale and title ---
    plt.xlim((0,4))
    plt.ylim((0, N + 5))
    plt.xlabel('time (s)', fontsize=20)
    plt.ylabel('trial', fontsize=20)
    plt.title(title_label, fontsize=22)

    # plt.plot(D[:,primeT],range(N),'r','linewidth', 5.5)
    # plt.plot([1.5]*N,[i for i in range(1,N+1)],'b','linewidth', 1.5)

    # --- alternative: patch ---
    # patch([ D(:,5); D(end,5); 0.00; 0.00], [1:N N+5 N+5 0], [.7 .8 .7] ,'EdgeColor' ,'none');% % [.7 .8 .7]
    # patch([ D(:,5); D(end,5); 0.00; 0.00], [1:N N+5 N+5 0], [0 110 144]/255 ,'EdgeColor' ,'none');% % [.7 .8 .7]
    # patch([ D(:,5); D(end,5); 0.00; 0.00], [1:N N+5 N+5 0], [0.8 0.8 0.8] ,'EdgeColor' ,'none');% % [.7 .8 .7]

    # --- printing the lines dividing the sessions ---
    div = np.where(np.diff(D[:, session]))[0]
    for i in div:
        plt.plot(plt.xlim(), [div, div], lw=0.5, alpha=0.5)

    # --- Plotting the response distributions ---
    xmin = 0
    xmax = 10
    # xma = np.round(np.mean(D[:,0])+4*np.std(D[:,0]))
    Nx_grid = 100
    x_grid = np.linspace(xmin, xmax, Nx_grid)
    split = True
    thisKDE = KDE(D[:, 0], x_grid, split=split, bw=0.5)

    # --- plotting lines splitting the sessions ---
    if not split:
        # plt.plot(x_grid,thisKDE,)
        plt.plot(x_grid, thisKDE, color='blue', alpha=0.5, lw=4)
    else:
        for i in range(0, 2):
            plt.plot(x_grid, thisKDE[i], alpha=0.5, lw=4)

    # plt.fill(,thisKDE,color='blue',alpha=0.5,lw=3,baseline=0)
    # --- Plotting the moving average of the lever press durations ---
    # plot(movingAverage(D(:,1),20),1:N,'linewidth',2);

    # --- Plotting each trial in a different style ---
    # mycolor = [0,.8,.9]
    mysize = 5
    lw = 0.8
    face_color = 'black'
    plt.scatter(D[validPrimed, 0], validPrimed, s=mysize, linewidths=lw, marker='o', c=face_color, edgecolors='k')
    plt.scatter(D[validNonPrimed, 0], validNonPrimed, s=mysize, linewidths=lw, marker='o', c=face_color, edgecolors='k')
    plt.scatter(D[invalid, 0], invalid, s=mysize, linewidths=lw, marker='o', c=face_color, edgecolors='r')

    plt.show()
    # --- mounting return variable ---
    return ([len(validPrimed) / N, len(validNonPrimed) / N, len(invalid) / N] * 100)


def individual_drrd(prefix='AB1', animalID=64, session=1, plotFlag=True, dataPath='', events_to_eliminate= None):
    # ________________________________________________________
    # Created on:       January, 14, 2019
    # Created by:       Marcelo Bussotti Reyes
    # Purpose:          A function that analyzes the performance of rats in the drrd
    #                   procedure. The animals are supposed to press a lever for a
    #					duration longer than the prime time in order to receive food.
    #
    # Input:            File name prefix as a string, animal id as a number and
    #					session as a number. Filene name are in the format:
    #					[prefix00A.00S] where A is the animal and S is the sesson
    #     prefix:       usually two letters and a digit, e.g. AB0
    #     animalID:     the number of the rat whose data you want to analyze
    #     session:      is the number of the session to be analyzed (starts in 1)
    #     plotFlag:     True if one desires to plot the results False otherwise
    #     dataPath:     path to where the data is. Default is the current directory
    #     events_to_eliminate: gives the option to eliminate events that are not analyzed, 
    #                   for example, responses to a non-rewarded lever, nosepokes, etc.
    #                   Needs to be a vector with all events e.g. [5,9] or [5], or None (default)
    #
    # Output:           D, a matrix with 6 columns containing data from a trial in
    #					each line.
    #                   In case saveMatFlag is parsed as true, the program will
    #                   save a matlab file (.mat) with the same name as the
    #                   original data file.
    # Coments:          Uses functions: med2tec.m, and plotDrrd.m
    #
    # Format:           drrd(prefix,animalID,session))
    # Example:          D = drrd('AB1',1,2);
    #					this will analyze the file AB1001.002, animal 1 and session 2
    # Previous Modifications:
    # Modification:     Included the option of saving a matlab file with the
    #                   matrix D. This helps to speed up the analysis when
    #                   several sessions and animals are to be analyzed.
    #                   Included session column (6th) and the input format
    #                   (Marcelo B. Reyes 2013)

    
    # defining a function meant to eliminate events that are not analyzed, for 
    # example, nose pokes, responses to a lever that is not rewarded, etc.
    # we suggest the experimentalist
    def eliminate_events(M, events_to_eliminate = None):
        
        # tests if there are events to be eliminated
        if events_to_eliminate != None:
        
            # loop over all events to be eliminated
            for event in events_to_eliminate:
            
                # keeps all events but the one to be eliminated
                # stores the results to the same variable
                M = [M[k] for k in range(len(M)) if M[k][1] != event]
        
        # returns the updated variable (the same parsed if events_to_elimnate == None)
        return(M)



    # builds file name from parsed parameters - uses zfill to pad with zeros
    filename = prefix + str(animalID).zfill(3) + '.' + str(session).zfill(3)

    # --- Indexes for each of the output variable columns
    dtCol = 0
    itiCol = 1
    primedCol = 2
    validCol = 3
    phaseCol = 4
    sessionCol = 5  # session variable column index
    Ncols = 6

    data = med2tec(dataPath + filename)  # reads data from medpc format to time-event code
    
    data = eliminate_events(data, events_to_eliminate = events_to_eliminate)
        
    data = np.array(data)                # transforms data to numpy array

    # Checking if data array is empty
    if len(data) == 0:
        print('Empty file or file not found: no data analyzed\n')
        return ([])
    
    # small correction for a bug in the med-pc file ---
    data = fix_clock_reset(data)


    # identifying in which index the first trial begins.
    index_first_start = np.where(data[:,1]==1)[0][0]
    index_first_end   = np.where(data[:,1]==3)[0][0]

#---------------------------------------------------------------------------------
    # ******** Attention **********
    # this line is a potential problem for data with valid an invalid trials
    # Still lacks tests

    
    # if the first trial begins after an end of trial, we need to eliminate 
    # the "fragments" of the first trial. Everything before (including) the trial
    # end must be eliminated
    if index_first_start > index_first_end:
        data = data[index_first_end+1:,:]
#---------------------------------------------------------------------------------


    # --- look for indexes of temporal events ---
    startIndex    = np.array([i for i in range(len(data)) if data[i, 1] ==  1])  # code that indicates that the lever was pressed
    endIndex      = np.array([i for i in range(len(data)) if data[i, 1] ==  3])  # instant when the lever was released
    primeIndex    = np.array([i for i in range(len(data)) if data[i, 1] == 18])  # moment when the trial reached the criterion
    lightOnIndex  = np.array([i for i in range(len(data)) if data[i, 1] == 11])  # when light was turned on (beginning of session)
    lightOffIndex = np.array([i for i in range(len(data)) if data[i, 1] == 21])  # when light was turned off (end of session)
    phaseAdvIndex = np.array([i for i in range(len(data)) if data[i, 1] == 17])  # moment of phase advance (criterion increased)
    phaseBckIndex = np.array([i for i in range(len(data)) if data[i, 1] == 27])  # moment of phase retreat (criterion decreased)

    # eliminating the last trial in case it was incomplete
    startIndex = startIndex[0:len(endIndex)]  
    
    print('_' * 80)
    print(len(endIndex),' trials found')
    
    # --- checking if there was at least one trial, 
    # --- otherwise stops the routine end return empty vector
    if len(startIndex) == 0:
        print('No trials recorded')
        return ([])

    # --- searching for trials in which the animals received food. We call these "primed" ----
    primedTrials = np.array([startIndex[startIndex < i].size - 1 for i in primeIndex])

    # --- searching for trials in which animals progressed or retreated phase
    phaseAdvTrials = np.array([startIndex[startIndex < i].size - 1 for i in phaseAdvIndex])
    phaseBckTrials = np.array([startIndex[startIndex < i].size - 1 for i in phaseBckIndex])

    # completing data in case the lightOff event wasn't found
    # includes a lightOff event in the last trial
    if len(lightOnIndex) != len(lightOffIndex):
        if len(lightOnIndex) == len(lightOffIndex) + 1:
            lightOffIndex = np.append(lightOffIndex, startIndex[-1] + 1)
        else:
            print('Incompatible number of events')
            return (-2)

    validTrials = np.array([], dtype=np.int64).reshape(0, )   # initiates empty vector to store all valid trials
    
    for i in range(len(lightOnIndex)):
        Nu = np.array(len(startIndex[startIndex < lightOnIndex[i]]))
        Nv = np.array(len(startIndex[startIndex < lightOffIndex[i]]))
        validTrials = np.hstack([validTrials, range(Nu, Nv)])

    validTrials = np.array(validTrials)  # transforms to numpy array

    # --- search for valid trials in which animals were and were not reinforded ---
    validPrimed    = np.intersect1d(validTrials, primedTrials)
    validNonPrimed = np.setdiff1d  (validTrials, primedTrials)
    invalid        = np.setdiff1d  (range(len(startIndex)), validTrials)

    # --- Organizing data in one single matriz: D --- 
    D = np.zeros((len(startIndex), Ncols))  # Initiates the vector for speed

    # --- Calculating the duration of the lever presses ---
    D[:, dtCol]    = data[endIndex, 0]       - data[startIndex, 0]     # time elapsed from press until relesase
    D[:-1, itiCol] = data[startIndex[1:], 0] - data[endIndex[:-1], 0]  # time between release and press ( beginning of next trial)
    D[-1, itiCol]  = np.nan                                            # the last trial has no ITI, stores NaN

    # --- saving each data in respective column
    if len(primedTrials) > 0: D[primedTrials, primedCol] = 1  # sets to 1 all the trials that were primed

    if len(validTrials) > 0:
        for k in validTrials: D[int(k), validCol] = 1

    if len(phaseAdvTrials) > 0: D[phaseAdvTrials, phaseCol] = 1
    if len(phaseBckTrials) > 0: D[phaseBckTrials, phaseCol] = -1

    # D[:,phaseCol]               = np.cumsum(D[:,phaseCol])+iniPh
    D[:, sessionCol] = session  # adds the session number to data (same for all lines)

    # Getting the time where the rats exceeded the criterion time
    # these are called primed trials, and occurr at the prime time.
    # we can get the primed times by subtracting the data where the prime happened 
    # from the event before, which necessarily is the lever press
    primeTimes = [0] * D.shape[0]

    for primed in primedTrials:
        thisIndex = startIndex[primed]
        primeTimes[primed] = np.round(data[thisIndex + 1, 0] - data[thisIndex, 0], decimals=5)

    print(np.unique(primeTimes))
    crit = extractCriterion(phAdv=D[:, phaseCol], primed=D[:, primedCol], primeTimes=primeTimes)
    D[:, phaseCol] = crit

    # --- printing output with the summary of the results ---
    print(
        f'Rat{animalID}  Trials:{len(startIndex)}  Reinforced:{len(validPrimed)}  Non-Reinforced:{len(validNonPrimed)}  Invalid:{len(invalid)}\n')

    # --- graphical part ---
    if plotFlag:
        plotDrrd(D, filename)

    # --- returning results ---    
    return (D)


def extractCriterion(phAdv, primed, primeTimes):
    def fixPrimeTime(ptvalues):  # for each prime time value
        print('Trying to fix inconsistencies in prime times')

        # first let'see if there are values close to zero
        t = [t for t in ptvalues if t > 0.1]

        # second, values that are too close
        t = [round(k, 1) for k in t]  # narrows down to round values
        t = list(set(t))              # and eliminate the values that are equal
        
        if len(t) == 1:
            print('Successfully fixed')
            print(t)
            return (t[0])
        else:
            print('Unable to fix')
            return ([])

    def getPrimeTimes(i, j, x):
        subset = list(set(x[i:j]))  # find the unique values in the list x
        subset.sort()  # organizes starting from zero

        if len(subset) > 2:  # there can be only zero and the prime time.
            print("Warning: More than one prime time found:", subset[1:])
            return (fixPrimeTime(subset[1:]))  # If there are more - something is wrong - return empty vector
        elif len(subset) == 2:  # if there are two, values are probably correct
            # TODO: here there is a loose end: there can be two values but one of them not be zero!!!
            return subset[1]  # return the last value of the list (the positive value)
        elif int(subset[0]) == 0:  # if there are only zeros in the list
            return (-1)  # no way to determine the prime time: returns -1
        else:  # if  there is one non-zero value, all trials were primed
            return (subset[0])  # return the prime time

    newPrimeTimes = primeTimes.copy()  # make a copy of prime time to return later

    if len(phAdv) == len(primed) == len(primeTimes):  # checks if all vectors are of the same size

        ind = [i + 1 for i, j in enumerate(phAdv) if j == 1]  # gets trials indices to find the trials of phase advance
        ind = [0] + ind + [len(phAdv)]  # add the beginning and the end indices
        ind = list(set(ind))  # eliminates redundancies - incase the animal advanced in the very last trial, for example
        ind.sort()

        for i, j in zip(ind[0:-1], ind[1:]):  # loop for initial and final indices for each phase

            pt = getPrimeTimes(i, j, primeTimes)  # get the prime time for the trials and checks for inconsistencies

            if not pt:  # if the prime time is empty, smthg is wrong
                print("inconsistency found")
                print("i, j = ", i, j)  # alert user!
                newPrimeTimes = [0] * len(phAdv)  # make all values equal to zero
                break
            elif pt == -1:  # if there were no correct trials
                if i > 0:  # and if it is not first phase
                    pt = newPrimeTimes[i - 1]  # gets the same criterion as last phase
                else:  # if it is first phase, just make it zero
                    pt = 0
            for k in range(i, j):  # finally uniformizes all values in the phase
                newPrimeTimes[k] = pt

    else:
        print("vectors are not of the same size")
        print(len(phAdv))
        newPrimeTimes = [0] * len(phAdv)  # make all values equal to zero
    return (newPrimeTimes)


def drrd(prefix='AB1', animalID=7, sessions=[7], plotFlag=True, dataPath='', events_to_eliminate = None):
    # function D = drrd(prefix,animalID,sessions,plotFlag)
    # prefix = the code for the experiment
    # example: D = drrd('AB1',1,1:9,True)
    # runs the sessions 1 through 9 for animal 1 of the AB1 experiment.

    D = np.array([])

    if type(sessions) == int: sessions = [sessions]  # transforms the int into a list

    for session in sessions:
        if len(D) > 0:
            D = np.vstack((D, individual_drrd(prefix, animalID, session, plotFlag=False, dataPath=dataPath, events_to_eliminate = events_to_eliminate)))
        else:
            D =               individual_drrd(prefix, animalID, session, plotFlag=False, dataPath=dataPath, events_to_eliminate = events_to_eliminate)

    if plotFlag:
        plotDrrd(D, title_label='Rat:' + str(animalID) + ' Sess:' + str(sessions))

    return (D)


def evaluateKDE(x, x_grid, bw=0.2):
    from statsmodels.nonparametric.kde import KDEUnivariate

    kde = KDEUnivariate(x)
    kde.fit(bw=bw)
    return (kde.evaluate(x_grid))


def scaleKDE(x, y, scale=''):
    # determines the amplitude of the KDE to fit to the graph range (y-axis)
    # usually it will scale to the number of trials

    if scale == '': scale = len(x)

    factor = scale / np.max(y)
    return (y * factor)


def KDE(x, x_grid=range(0, 8), bw=0.25, split=True, NSplit=100, normalize=False):
    # function that calculates the kernel density estimate
    # parameters: 
    #   xmi,xma = minimum (maximum) value of the range for x_grid
    #   Nx_grid   = number of x_grid for kernel estimation
    #   bw      = width of the kernel
    #   split   = True if the routine splits in two phases, beginning and end of the session
    #   plotFlag= True to also print output to a graph otherwise just return the calculations
    #   Nsplit  = number of trials to use at the beginning and end of the session    

    if 2 * NSplit > len(x):
        NSplit = int(np.floor(len(x) / 2))
        print('Number of trials to split too large, trials split in half')

    if not split:
        y = evaluateKDE(x, x_grid, bw=bw)
        return (scaleKDE(x, y))
    else:
        xSplit = x[0:NSplit]
        y = evaluateKDE(xSplit, x_grid, bw=bw)
        yscaled1 = scaleKDE(xSplit, y, scale=len(x))

        xSplit = x[-NSplit - 1:-1]
        y = evaluateKDE(xSplit, x_grid, bw=bw)
        yscaled2 = scaleKDE(xSplit, y, scale=len(x))

        return ([yscaled1, yscaled2])


def fix_clock_reset(data):
    # sometimes the clock reset of the hardware (medPC) happens after
    # the beginning of the experiment. Typically it happens in the second
    # trial, but we detected that is happened later. This small routine 
    # fix this bug by eliminating all the trials before the reset. 

    # --- detecting if the event times are allways increasing
    # --- using the np.diff and then the where
    ind = np.where(np.diff(data[:, 0]) < 0)  # looks for moments where the time decreased (diff<0)
    ind = ind[0]  # gets only the ndarray (peculiar to numpy.where and its way to return the variables)

    # detects if there were one or more events
    # keeping only the last one
    if len(ind) >= 1:
        print('Warning: eliminated first ', ind, ' trials due to clock reset issue')
        ind = int(ind[-1])
        if data[ind + 1, 1] == 3 and len(data[:, 0] > 2):  # if it starts with a 3 (which is a lever release), eliminate this trial
            return (data[(ind + 2):, :])
        else:
            return (data[(ind + 1):, :])

            # if it is always increasing, just returns the same data
    return (data)


def separate_gaussians(x, params, number):
    (gamma, mu1, sigma1, mu2, sigma2) = params
    dt = np.unique(np.round(np.diff(x), 4))
    res = (1 - gamma) * np.exp(- (x - mu1) ** 2.0 / (2.0 * sigma1 ** 2.0)) \
          + gamma * np.exp(- (x - mu2) ** 2.0 / (2.0 * sigma2 ** 2.0))
    norm = np.sum(res) * dt

    if number == 0:
        res = (1 - gamma) * np.exp(- (x - mu1) ** 2.0 / (2.0 * sigma1 ** 2.0)) / norm
    else:
        res = gamma * np.exp(- (x - mu2) ** 2.0 / (2.0 * sigma2 ** 2.0)) / norm
    return res


def double_gaussian(x, gamma, mu1, sigma1, mu2, sigma2):
    dt = np.unique(np.round(np.diff(x), 4))  # getting the value of dt as difference between x steps
    res = (1 - gamma) * np.exp(- (x - np.abs(mu1)) ** 2.0 / (2.0 * sigma1 ** 2.0)) \
          + gamma * np.exp(- (x - np.abs(mu2)) ** 2.0 / (2.0 * sigma2 ** 2.0))
    return (res / np.sum(res) / dt)


def single_gaussian(x, mu, sigma):
    dt = np.unique(np.round(np.diff(x), 4))
    res = np.exp(- (x - mu) ** 2.0 / (2.0 * sigma ** 2.0)) / sigma / np.sqrt(2 * np.pi)
    return (res / np.sum(res) / dt)


def get_data_from_rat(animal, session, x, indexes=(0, None), prefix='AL0', dataPath=homePath):
    # indexes are the trials that should be considered 
    # example: indexes=(0,100) only the first 100 trials will be analyezed
    # indexes = (-100, None) the last 100 trials will be analyzed

    # get the indexes values to slice the vector
    i, j = indexes

    # reading data from raw file
    D = dr.drrd(prefix=prefix, animalID=animal, sessions=session, plotFlag=False, \
                dataPath=dataPath)

    # Trying to slice the data matrix according to indexes
    try:
        D = D[i:j, :]
    except:
        print('Could not reduce the data matrix to indexes')
        return ([], [])

    y = calc_histogram(D[:, 0], x)

    return (D, y)


def calc_histogram(T, x):
    # calculating the histogram counts for the data
    import seaborn as sns

    pal = sns.cubehelix_palette(6, rot=-.45, light=.7)
    data = plt.hist(T, bins=x - dt / 2, histtype='stepfilled', density=True, align='mid', color=pal[0], alpha=1)

    # selecting only the counts and ignoring the bins edges
    data = data[0]

    # adding zero to the end to generate a vector with the same size as x
    y = np.append(data, 0)

    # normalizing the histogram
    y = y / np.sum(y) / dt

    return (y)


def plot_all_curves(x, fit, fit_single, dt, fontDict=font, plotSingle=False):
    # plotting the separate gaussians

    plt.plot(x, separate_gaussians(x, fit, number=0), '-.', c='k', lw=1.0, label='Gauss1')
    plt.plot(x, separate_gaussians(x, fit, number=1), '--', c='k', lw=1.0, label='Gauss2')
    plt.plot(x, double_gaussian(x, *fit), '-', c='k', lw=5.0, label='dGauss', alpha=0.3)
    # plt.xlabel('Duration (s)')

    if plotSingle:
        # plotting the single gaussian fitted
        plt.plot(x, single_gaussian(x, *fit_single), '-.', c='k', lw=1.5, label='single gauss.')

    # xlim and adding axes labels
    plt.xlim([-dt / 2, 3])
    plt.xlabel('duration (s)')
    plt.ylabel('probability density')

    # returning the axes for further plotting
    return (plt.gca)


def calculate_distance(x, y, fit, model=0):
    if model == 'double':  # double gaussian
        ycalc = double_gaussian(x, *fit)
        sqDist = np.sum((y - ycalc) ** 2)

    else:
        ycalc = single_gaussian(x, *fit)
        sqDist = np.sum((y - ycalc) ** 2)

    return (sqDist)


def fix_parameters_order(p):
    # p is a vector with 5 entries: gamma, mu1, sigma1, mu2, sigma2
    # if mu1>mu2 (the first gaussian has a greater mean), we invert the order of the parameter
    if p[1] > p[3]:
        print('Parameters inverted')
        p[0] = 1 - p[0]  # we invert the order, and gamma needs to be replaced by 1-gamma
        return (p[[0, 3, 4, 1, 2]])
    else:
        return (p)


def fit_single_animal(animal, session, plotFlag=True, indexes=(0, None), ax=None, prefix='AL0', dataPath=homePath):
    # get data from file: 
    # y is the probabilty density with bins defined by x
    # D is the data from the session obtained from drrdTools.drrd() method

    # getting data from file
    D, y = get_data_from_rat(animal, session, x, indexes, prefix=prefix, dataPath=dataPath)

    # Checking if D is not empty    
    if len(D) == 0:
        print('No data found, quiting.')
        return ([np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan])

    # fit single gaussian
    popt_sngl, pcov_sngl = curve_fit(single_gaussian, x, y, bounds=(0, [10, 5]), p0=[1, 0.5])
    d_sngl = calculate_distance(x, y, popt_sngl, model='single')

    # fit the double gaussian
    initParsDouble = (0.5, 0.2, 0.1, 1, 0.5)
    popt, pcov = curve_fit(double_gaussian, x, y, bounds=(0, [1, 5, 5, 5, 5]), p0=initParsDouble)

    popt = fix_parameters_order(popt)

    d = calculate_distance(x, y, popt, model='double')

    # bundling all values to return them together
    r = np.concatenate((popt, [d], popt_sngl, [d_sngl]))

    # drawing all curves into the graph
    if plotFlag:
        plot_all_curves(xfine, popt, popt_sngl, dtFine)
    #       add_info_to_graph(animal,session,r)
    #

    # That'it. Let's return the values and get out of here
    return (r)


def add_info_to_graph(rat, session, rr, addSingleInfo=False):
    # ax = axis handle to plot
    # rat = animal id
    # session = session to plot
    # rr the parameters of the fit: gamma mu1 sig1 mu2 sig2 chi mu sig chi
    # parameters for the graphical output

    font = {'family': 'sans',
            'weight': 'medium',
            'size': 12}
    plt.rc('font', **font)
    plt.rc('text', usetex=True)

    # adding a title for the graph
    plt.title('animal:' + str(rat) + '\tsession: ' + str(session))

    # some mambo jambo with the coordinates to find the right coordinates for text
    yi, yf = plt.ylim()
    xi, xf = plt.xlim()
    xt = (xf - xi) * 0.70
    yt = (yi + yf) * 0.55
    dy = (yf - yi) * 0.06

    # now let's finally add the information to the graph
    rr = np.round(rr, 2)
    plt.text(xt, yt, '\mu_1=' + str(rr[1]) + ' \sigma_1=' + str(rr[2]))
    plt.text(xt, yt - dy, '\mu_2=' + str(rr[3]) + ' \sigma_2=' + str(rr[4]))
    plt.text(xt, yt - 2.1 * dy, '\gamma=' + str(rr[0]) + ' sd_d=' + str(rr[5]))
    plt.text(xt, yt - 2.1 * dy, '\gamma$=' + str(rr[0]))

    if addSingleInfo:
        plt.text(xt, yt - 4 * dy, '$\mu$= ' + str(rr[6]) + ' $\sigma$=' + str(rr[7]))
        plt.text(xt, yt - 5.3 * dy, '$sd_s$= ' + str(rr[8]))

    plt.legend()
    plt.show()


if __name__ == "__main__":
    D = drrd(prefix='AB1', animalID=64, sessions=[1], plotFlag=True, dataPath='../sampleData/')
    
