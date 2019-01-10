def med2tec(fNAME,flag='A'):
    # written by Marcelo Bussotti Reyes - Universidade Federal do ABC
    # Based on the Matlab code writen by Marcelo Caetano.
    # 2019-01-10

    TIME_RESOLUTION = 2e-3 # (in seconds) Time resolutuion in the med associates box = 2ms

    try:
        fid = open(fNAME, 'r')
    except IOError:
        print("Could not read file:", fName)

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
    count = 0
    for i in range(startParse,stopParse):
        for j in str.split(fileString[i])[1:]:
            data = str.split(j,'.')
            M.append([round(int(data[0])*TIME_RESOLUTION,3), int(data[1])])
    
    fid.close()

    return(M)

