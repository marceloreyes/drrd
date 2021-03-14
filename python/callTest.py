import matplotlib.pyplot as plt
import drrdTools as dt

prefix = 'AB1'
id   = 46
#sessions = 1#,2,3,4,5]
sessions = [1,2,3,4,5]
D = dt.drrd(prefix=prefix,animalID=id,sessions=sessions,dataPath='/home/mbreyes/ufabc/dados/AB/data/raw/')
