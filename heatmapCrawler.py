from pickle import FALSE
import pandas as pd
import statistics as st
import itertools as it 
import scipy as sp

hmat = pd.read_csv(snakemake.input[hmat], delimiter=",")
#hmat = pd.read_csv("/storage/mathelierarea/processed/petear/analysis/test/hmat.csv", delimiter=",")
ClusterDF = pd.read_csv(snakemake.input[cluster], delimiter=",")
#ClusterDF = pd.read_csv("/storage/mathelierarea/processed/petear/analysis/test/ClusterDF.csv", delimiter=",")

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

