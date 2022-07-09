from pickle import FALSE
import pandas as pd
import statistics as st
import itertools
import scipy as sp
from collections import defaultdict
import math     #for sqrt

###READING IN THE DATA####
hmat = pd.read_csv(snakemake.input[hmat], delimiter=",")
#hmat = pd.read_csv("/storage/mathelierarea/processed/petear/analysis/test/hmat.csv", delimiter=",")
ClusterDF = pd.read_csv(snakemake.input[cluster], delimiter=",")
#ClusterDF = pd.read_csv("/storage/mathelierarea/processed/petear/analysis/test/ClusterDF.csv", delimiter=",")

hmat = hmat.set_index("Unnamed: 0")
hmat.index.names = ["TopicNo"]

####SORTING THE DATA####
#Sort hmat by clusterDF
ClusterDF["Cluster"] = ClusterDF["Cluster"].str.replace(r"Cluster", '')

sorter = list(ClusterDF["Patient"])
cluCol = list(ClusterDF["Cluster"]) 

hmat = hmat[sorter]


thmat = hmat.T

thmat["Cluster"] = cluCol   

simpThmat = thmat[[1,"Cluster"]]


#####CREATING DICTIONARY OF CLUSTERS####
uniques = simpThmat.Cluster.unique()

clusterDict = {}  
varianceDict={}

"""
Create a dictionary of clusters with the following structure: Key = Cluster, Value = List of patients

"""
for value in uniques:
    n=int(value)
    clusterList = []
    varList = []
    ta = 0
    for index, row in simpThmat.iterrows():
        if ta < len(simpThmat):
            if int(simpThmat["Cluster"].iat[ta]) == n:
                var = simpThmat[1].iat[ta]    
                clusterList.append(index)
                varList.append(var)
            ta += 1
            clusterDict[n] = clusterList                               
    varianceDict[n] = st.variance(varList)         
      

"""
Function creates chunks over a list, to be used later for selecting low variance values. 
WARNING had to change from tuples to list because downstream function was not working - Gave error: "TypeError: 'tuple' object is not callable"
"""   
def makeChunks(iters, chunks):
    return ([
        list([x for x in y if x])
        for y in list(itertools.zip_longest(*[iter(iters)] * chunks))
    ])

"""
Make dictionary with keys being columns where the first digit is the column number and the second digit is the patients numbered so patients can be found again. Values are the variance of the patients: 7_71 is column 7, patients in the 71st chunk. 
"""
tupVarDict = {}
nint = 3
for clusKey, clusVal in clusterDict.items():
    df = thmat.loc[clusVal]
    for column in df.columns[:-1]:
        vector = df[column]
        chunkedVector = makeChunks(vector, nint)
        CVlen = 1
        col = column
        for element in chunkedVector:
                ident = (f"{col}_{CVlen}")
                tupVarDict[ident] = element
                CVlen += 1



    
for valKey, valVal in tupVarDict.items():
    for element in valVal:
        var1 = st.var(element))
        




"""
Create a dictionary of columns and the respective values from the hmat-dataframe. Also the clusters are preserved in the keys name so that the first digit is the cluster number and the number after the hyphen is the column number: 1-8 is cluster 1 for column 8. As value for the dictionary are a list of the values in the column for that cluster.
"""
valDict = {}

for col in thmat.columns[:-1]:
    for key, value in clusterDict.items():
        #print(f"this is key:{key}, this is value: {value}, len is {len(value)}, and this is column: {col}")
        valval = value
        valList = []
        for element in valval:
            vals = thmat.loc[[element], [col]]
            vals = vals[col].values[0]
            valList.append(vals)
            keyID = (f"{key}-{col}")
        valDict[keyID] = valList   
  



########################TESTING AREA############################
"""
For this funcion add a test by feks having a dictionary with n1n2 for the variance of n1 and n2. Feks n1n2 as key and the variance as value. Then check n2n3 and the variance. Continue until end of column. The lowest variance is the starting point for finding clusterblocks. 
"""
nnDictClusKey = {}
nnDictFull = {}

for clusKey, clusVal in clusterDict.items():
    df = thmat.loc[clusVal] 
    m = len(df)
    nnList = []
    nnDictClusKey = {}
    for column in df.columns[:-1]:
        nnDictColCluster = {}
        varList = []
        for n in range(m):
            if n == m-1:
                break
            else:
                val_1 = df.iloc[n][column]
                o = n+1
                val_2 = df.iloc[o][column]
                keyIndex1 = df.index[n]
                keyIndex2 = df.index[o]
                key = (f"{keyIndex1}")
                var = st.variance([val_1, val_2])
                nnDictFull[key] = var
                nnDictColCluster[key] = var
        minVarValCluster = min(nnDictColCluster, key=nnDictFull.get)
        nnList.append(minVarValCluster)  
        nnDictClusKey[clusKey] = nnList
    #the part where I crawl over the values:
    for donoKey, donoVal in nnDictClusKey.items():
        for element in donoVal:
            var1 = st.variance(df.iloc[n][column]), df.iloc[n+1][column]) 
            var2 = var1 
            while var1 <= var2
                
"""
Because variation will continously decrease as long as the new value is under sqrt(1+(1/n)) * stdev ,  from the mean. So negative if: sqrt(1+1/n) ⋅ stdev  > x > mean  
"""
def varLenTest(newval,lenlist):
    newVal = newval
    nlen = len(lenlist)
    meanmean = st.mean(lenlist)
    standDev = st.stdev(lenlist)
    trhsHold = math.sqrt(1+(1/nlen)) * standDev
    if trhsHold > newVal and newVal > meanmean:
        return True


#Trying to crawl by checking variance by every couple 
for clusKey, clusVal in clusterDict.items():
    df = thmat.loc[clusVal] 
    for colu in df.columns[:-1]:
        nnDictColCluster = {}
        valList = []
        nnDictPatients = {}
        patList = []
        m = len(df)
        n = 0
        o = 0
        while n in range(m):
            if n == m-2:
                break
            else:
                val_1 = df.iloc[n][column]
                pat_1 = df.index[n]
                valList.append(val_1)
                patList.append(pat_1)
                o = n+1
                val_2 = df.iloc[o][column]
                pat_2 = df.index[o]
                valList.append(val_2)
                patList.append(pat_2)
                #keyIndex1 = df.index[n] not needed thus commented out
                #keyIndex2 = df.index[o] not needed thus commented out
                key = (f"{keyIndex1}")
                var1 = st.variance(valList)
                val_3 = df.iloc[o+1][column]
                valList.append(val_3)
                var2 = st.variance(valList)
            
                while var1 >= var2 and o < m-1 and varLenTest(var2, valList) == True:
                    o += 1
                    val_2 = df.iloc[o][column]
                    pat2 = df.index[o]
                    valList.append(val_2)
                    patList.append(pat2)
                    var2 = st.variance(valList)
            
            ##next n!:::
            keyName = (f"{clusKey}-{colu}")
            nnDictColCluster[keyName] = valList
            n=o+1    
            ###    









for varKey, varVal in nnDictFull.items():           
     




















            minVarVal  = min(nnDictFull, key=nnDictFull.get)


            if n < 2:
                
                val_1 = df.iloc[n][column]
                varList.append(val_1)
                
            elif n > 3:
                
                var1 = st.variance(varList)
                print(f"is{n}")



########same as over but copied for n1n2
for clusKey, clusVal in clusterDict.items():
    df = thmat.loc[clusVal]
    
    m = len(df)
    for column in df.columns[:-1]:
        varList = []
        for n in range(m):
            if n < 2:
                
                val_1 = df.iloc[n][column]
                varList.append(val_1)
                
            elif n > 3:
                
                var1 = st.variance(varList)
                print(f"is{n}")

####################


        n -= 1
        val_1 = df.iloc[n][column]
        while n > 0:
            n -= 1
            val_2 = df.iloc[n][column]
            varList = []
            varList.append(val_1)
            varList.append(val_2)
            var1 = st.variance(varList)











for clusKey, clusVal in clusterDict.items():
    df = thmat.loc[clusVal]
    for column in df.columns[:-1]:
        n = len(df)
        n -= 1
        val_1 = df.iloc[n][column]
        while n > 0:
            n -= 1
            val_2 = df.iloc[n][column]
            varList = []
            varList.append(val_1)
            varList.append(val_2)
            var1 = st.variance(varList)

#Make dictionary with values being values in the dataframe

keyN = 1

#work

        
        if key == thmat["Cluster"].iat[ta]:
            print(key)
            colValList = []
            for i in value:
                cellVal = thmat.loc[i][col]
                colValList.append(cellVal)
                print(colValList)
            colVar = st.variance(colValList)    
                

        n += 1



thmat.loc[thmat['Clusters'] == some_value]

def crawler(list): 
    x = 3
    for 





###################4


uniques = simpThmat.Cluster.unique()

clusterDict = {}  
clusterList = []                        #R1-R4

for column in simpThmat.columns:
    FullColVariance = simpThmat.var()[column]
    
for value in uniques:
    n=int(valu
    clusterList = []
    ta = 0
    for index, row in simpThmat.iterrows():
        if ta < len(simpThmat):
            if int(simpThmat["Cluster"].iat[ta]) == n:
                clusterList.append(index)
            ta += 1
            clusterDict[n] = clusterList                               

def getCluster(cluster):
    return clusterDict[cluster]





######################


colList1 = list(hmat.columns)    

colList2 = CluDFSorted.Cluster

hmat.columns = [colList1, colList2]



n=1
for ind in thmat.index.unique("Cluster"):
	print(ind)	
























maxVal = hmat.max().max()
minVal = hmat.min().min() 
"""
function for 
"""
n = 5
for column, row in hmat.iterrows():
    

tests to do in the row: 
for row in hmat.iterrows(): 
        
    o = 4
    var = st.var(row.iloc[n:o])
    if low variance in n first columns:
        while variance not increasing by more than 5 % adding +1 to nn, is TRUE:
            test if the variance doesnt increase noteworthy if adding another
        while variance not increasing by more than 5 % addding +2 to n is TRUE

    elif :  
     
                


    else (low variance of first n columns is FALSE)
        skip 1 further for new n



"""
Test area :
"""        

#itertuples becasuse it iterates over DataFrame rows as namedtuples of the values.
#starts at 1 because itertuples starts at rowname as 0. 
for topic in hmat.itertuples():

topic = hmat.loc["Topic6"]
n=0
o=9
fullVar = st.variance(topic)
var = st.variance(topic[n:o])


c = 1
for row in hmat.iterrows():
    for hmat.columns.get_level_values(1) = 1:
      print()