"""
Make heatmap based on the probabilites of CpGs belonging to a certain topic. 
"""
"s"

"""
The methylation level of CpGs in the topics, you could plot the complete distribution of values for a given set of patients associated with the corresponding topic and compare to other CpGs in the other topics. For instance, using violin plots with strips
"""
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import lzma
import scipy
import statistics
import numpy as np

####Read tables#######################
"""
methTable = pd.read_csv("/storage/mathelierarea/processed/petear/analysis/test/methtable.csv.xz", delimiter=",", compression="gzip")
regAssigUnormal = pd.read_csv("/storage/mathelierarea/processed/petear/analysis/test/RegAssigUnormal.csv", delimiter=",")
regScrPrtopic = pd.read_csv("/storage/mathelierarea/processed/petear/analysis/test/RegScrPrtopic.csv", delimiter=",")
topicAssigToPatient = pd.read_csv("/storage/mathelierarea/processed/petear/analysis/test/topicAssigToPatient.csv", delimiter=",")
bina = pd.read_csv("/storage/mathelierarea/processed/petear/analysis//test/bina.csv", delimiter=",")
meta = pd.read_csv("/storage/mathelierarea/processed/petear/analysis/test/meta.csv", delimiter=",")
hmat = pd.read_csv("/storage/mathelierarea/processed/petear/analysis/test/hmat.csv", delimiter=",")
"""
methTable = pd.read_csv(snakemake.input[methTab], delimiter=",", compression="gzip")
regScrPrtopic = pd.read_csv(snakemake.input[regScrNorm], delimiter=",")
regAssigUnormal = pd.read_csv(snakemake.input[regScrUnrm], delimiter=",")
topicAssigToPatient = pd.read_csv(snakemake.input[topicAssig], delimiter=",")
bina = pd.read_csv(snakemake.input[bina] delimiter=",")
meta = pd.read_csv(snakemake.input[methTab], delimiter=",")
hmat = pd.read_csv(snakemake.input[hmat], delimiter=",")
####################################


#assign patients to topics, and topics to CpGs, so you can go from patients to CpGs
methTabWithRegAsRowname = methTable
methTabWithRegAsRowname.index = methTabWithRegAsRowname["Unnamed: 0"] 
methTabWithRegAsRowname.index.names = [None]

#regTable = methTabWithRegAsRowname["Unnamed: 0"]

methTabWithRegAsRowname = methTabWithRegAsRowname.drop(["Unnamed: 0"], axis=1)




"""
Create a dict with topics as key and donors as values. Which topic the donors belong to is found by subtracting each value in the row of topicAssigToPatient by the mean of that row. Then the row with the highest value after the subtraction is found and this is the topic for that column/donor 
Using raw values will just select the same topic for all donors as some topics have a higher value for all donors compared to other topics.  
Gå over rad og finn donorer som har størst contribution. 
"""

#Add Topic to int in Unnamed: 0 column
topicAssigToPatient["Unnamed: 0"] = "Topic" + topicAssigToPatient["Unnamed: 0"].astype(str)

#Remove unnamed column with topic numbers and set those as rownames
topicAsPatRownames = topicAssigToPatient.set_index("Unnamed: 0")
topicAsPatRownames.index.names = [None]
#topicAsPatRownames["mean"] = topicAsPatRownames.mean(axis=1)

topAssigPatT = topicAsPatRownames.T

#Create a dict with with topics as key and donors' probabilites for belonging to that topic as values, but only values after a cutoff of median absolute deviation + median of that c
topicDict={}
for column in topAssigPatT:
    colMad = topAssigPatT[column].mad(axis=0)
    colMedian = topAssigPatT[column].median(axis=0)
    cutOff = colMad+colMedian
    maxMin = topAssigPatT[column].max(axis=0) - topAssigPatT[column].min(axis=0)
    colMean = topAssigPatT[column].mean(axis=0)
    highInsurance = maxMin - 1.5*colMean 
    Q1 = topAssigPatT[column].quantile(0.25)
    Q3 = topAssigPatT[column].quantile(0.75)
    IQR = Q3 - Q1
    print(f"{column}: {IQR}")
    if highInsurance < 0:
        topicDict[column] = topAssigPatT[column][topAssigPatT[column]>cutOff]

#function to find CpGs with similar values in each topic in topAssigPatTx


  

#Virker men må se om dette stemmer med det som sees i heatmap





for key, value in topicDict.items():
    print(key)

   

#topicAsPatRownames["MaxMin"] = np.ptp(topicAsPatRownames.values, axis=1)

highRows = topicAsPatRownames.eq(topicAsPatRownames.max()).stack()


#subract mean of each 
topicPatMean = topicAsPatRownames.sub(topicAsPatRownames.mean(axis=1), axis=0)

highRows = topicPatMean.eq(topicPatMean.max()).stack()
topicPatDict = highRows[highRows].reset_index(level=1).groupby(level=0)['level_1'].agg(list).to_dict()




"""
Use binarized to find high rank CpGs: Patients -> Topics -> CpGs :: Patients -> CpG
"""

bina = bina.drop("Unnamed: 0", axis=1)
bina = bina.set_index("chrPos")
bina.dropna(axis = 0, how = 'all', inplace = True)

    binaDict = {col : bina[bina[col].notnull()].index.tolist() for col in bina.columns}

#find number of values per key
nr_vals_binaDict = {key : len(value) for key, value in binaDict.items()}

for key, value in topicDict.items():
    print(key)


"""
Select rows and columns from methTabWithRegAsRowname that have high contribution to topics (no differentiating between topics here)
"""			
col_names = []
for key, col_list in topicPatDict.items():
    col_names += col_list

row_names = []
for key, row_list in binaDict.items():
    row_names += row_list

df = methTabWithRegAsRowname[col_names].filter(items=row_names, axis='index')
df = df.sort_index()

#df = df.T


"""
The goal is to find metadata for the heatmap. To do this assign metadata to every donor in topicPatDict
"""

meta = meta.fillna('NA', inplace=False)
meta = meta.rename(columns={"Unnamed: 0" : "Donors"})

meta.index = meta["Donors"] 
meta.index.names = [None]
meta = meta[meta.index.isin(df.index)]

"""
Because the dataframe 'df' is too big for seaborn to make a clustermap from I have to write the dataframes to csv to use in R with complexheatmap 
"""


df.to_csv(snakemake.output[py_DF]) 
meta.to_csv(snakemake.output[py_meta])

"""
df.to_csv("py_df.csv") 
meta.to_csv("py_meta.csv")
df.to_csv("py_df.tsv", sep="\t") 
"""



"""
Because 'df' is very large I will attempt to reduce the size here by removing rows with low difference in values
"""


#subtract min value from max value 
#df["mean"] = df.mean(axis=1)
#DFMean = df.sub(df.mean(axis=1), axis=0)

#Add column with Max - Min value in that row

df["MaxMin"] = np.ptp(df.values, axis=1)

#Remove rows under Min and over Max
Min = 0.3
Max = 1 
df = df[df["MaxMin"].between(Min, Max)]

"""
hmat is the matrix, topAss
"""     
hmat["Unnamed: 0"] = "Topic" + hmat["Unnamed: 0"].astype(str)


hmat = hmat.set_index("Unnamed: 0")
hmat.index.names = [None]

for column, row in hmat.iterrows():
    print(row)
    hmgrouped = hmat.groupby(np.ceil(hmat[column] * 100) // 10).mean()    
#group each row by similar values


maxVal = hmat.max().max()
minVal = hmat.min().min() 



thamat = hmat.T