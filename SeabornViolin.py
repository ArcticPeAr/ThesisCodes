"""
The methylation level of CpGs in the topics, you could plot the complete distribution of values for a given set of patients associated with the corresponding topic and compare to other CpGs in the other topics. For instance, using violin plots with strips
"""
    from matplotlib.backends.backend_pdf import PdfPages
    import matplotlib.pyplot as plt
    import seaborn as sns
    import pandas as pd
    import pyreadr
#####

####Read tables#######################
methTable = snakemake.input["methTable2rds"]
regScrNorm = snakemake.input["RegScrPrtopicOut"]
regScrUnrm = snakemake.input["RegAssigUnormalOut"]
topicAssig = snakemake.input["topicAssigToPatientOut"]

#read RDS in python 
RDSreadMT = pyreadr.read_r(methTable) 
methTable = list(RDSreadMT.values())[0]
methTable.index.names = [None]

#read RDS in python
#remove rownames 
RDSreadRSPT = pyreadr.read_r(regScrNorm)
regScrPrtopic = list(RDSreadRSPT.values())[0]

RDSreadRSU = pyreadr.read_r(regScrUnrm)
regAssigUnormal = list(RDSreadRSU.values())[0]

RDSreadTA = pyreadr.read_r(topicAssig)
topicAssigToPatient = list(RDSreadTA.values())[0]

####################################


"""
Following to find the donors that have the highest topic assignment:
"""
#set row as index
topicAssigToPatient["index"] = range(1, len(topicAssigToPatient) + 1)
topicAssigToPatient.set_index('index', inplace=True)

#subract mean of each 
topicPatMean = topicAsPatRownames.sub(topicAsPatRownames.mean(axis=1), axis=0)

#topicAbsMean = topicAbsMean.abs()


#Make dictionary-class .
class topDict(dict):
	def __init__(self):
		self = dict()
	def add(self, col, row):
		self[col] = row	





# Make a list of column names
colList = list(topicPatMean)



topicPatDict = topDict()

#Add key (patientID) with n topics (as values), after subtracting mean of that topic - this because it can find the topics contributing most for that patient compared to other topics.
for item in colList:
	topicPatDict.col = item
	patient = topicPatMean[item]
	topicPatDict.row = patient.nlargest(4)
	topicPatDict.add(topicPatDict.col, list(topicPatDict.row.index))





#Find number of topics 
tpcNo = len(topicAssigToPatient)

#Make lists for each topic, named topic_ + topic number 
nr = 0

while nr < tpcNo: 
	nr = nr + 1
	locals()[f"Topic_{str(nr)}"] = []
	

#number of topics and add 1 because this is used for setting columns to be worked on on following for loop. (First column is region: "Unnamed: 0")
tcpNoPlus1 = tpcNo+1

#Make list with Topic_X as name of list (x = topic number) and patients as items where the patients had that topic as highest contributing topic. This will cluster patients that has the same "highest contributing topic" together
for key, value in topicPatDict.items():
    for n in range(1,tcpNoPlus1,1):
        if n in value:
            locals()[f"Topic_{str(n)}"].append(key)



"""
Following to find regions most contributing to a topic
"""
#Set region as rownumber methTable and drop the column. Also make regTable as DF with only region as column so that columns from topics in next while-loop can be added
methTabWithRegAsRowname = methTable
methTabWithRegAsRowname.index = methTable["Unnamed: 0"] 
methTabWithRegAsRowname.index.names = [None]

regTable = methTabWithRegAsRowname["Unnamed: 0"]

methTabWithRegAsRowname = methTabWithRegAsRowname.drop(["Unnamed: 0"], axis=1)


###################################    RegScrPrTopic    #######################################################
#Remove unneeded columns
regScrPrTopicDrpd = regScrPrtopic.drop(["seqnames", "start", "end", "width", "nCounts", "nCells"], axis=1)

#Make dict where keys are columns(topicScores) and values are columns of "Unnamed: 0 and the values, respectively. For the values, those below 0.5, are discarded "
LookupDict = {}
for col in regScrPrTopicDrpd.columns[1:]:
    newDF = regScrPrTopicDrpd[["Unnamed: 0", col]]
    newDF = newDF[newDF[col] > 0.5]
    LookupDict[col] = list(newDF["Unnamed: 0"])
"""
{'Scores_Topic1':                       Unnamed: 0
								205     chrX:134569361-134569362
								227       chrX:89177497-89177498
								314       chrX:75249221-75249222
								387     chrX:111874417-111874418
								389     chrX:145077290-145077291
"""
#####################################################################################
###


#Make dataframe named TopFrame_X with X_patients as columns to ease further downstream analysis (X = topic number). Make new DF of 30 % of columns in each Topic
nr = 0
while nr < tpcNo:
	nr = nr + 1
	locals()[f"Scores_Topic{str(nr)}"] = methTabWithRegAsRowname[methTabWithRegAsRowname.columns & locals()[f"Topic_{str(nr)}"]]
	locals()[f"Scores_Topic{str(nr)}"] = locals()[f"Scores_Topic{str(nr)}"].add_prefix(f"{str(nr)}_") 
	locals()[f"Scores_Topic{str(nr)}"] = pd.concat([methTable["Unnamed: 0"], (locals()[f"Scores_Topic{str(nr)}"].sample(frac=1, axis=1))], axis=1)


#Remove rows from Scores_TopicX that are not present in the LookupDict
for key in LookupDict.keys():
    locals()[f"{str(key)}"] = (locals()[key].loc[LookupDict[key]])






nr = 0
while nr < tpcNo:
    nr = nr + 1
    df = locals()[f"Scores_Topic{str(nr)}"].round(decimals=2)
    for column in df.columns[1:]:
        df1 = pd.concat([df["Unnamed: 0"], df[column]], axis=1)
        df1 = df1.drop_duplicates(subset=[column])
        g = sns.violinplot(y=column, data=df1)
        fig = g.get_figure()
        name = f"{column}.png"
        fig.savefig(name)
        plt.clf()







#less decimals in the table to easily drop duplicates in the next for-loop
deciRegTab = regTable.round(decimals=3)	 
dRTSansUn = deciRegTab.drop(["Unnamed: 0"], axis=1)

#drop rows that have duplicate values in dataframe. Then make new dataframe of each colum with that column and column named "Unnamed: 0"
for column in dRTSansUn.columns:
	locals()[f"DF_{column}"] = pd.concat([deciRegTab["Unnamed: 0"], dRTSansUn[column]], axis=1)
	locals()[f"DF_{column}"] = locals()[f"DF_{column}"].drop_duplicates(subset=[column])
	g = sns.violinplot(y=column, data=locals()[f"DF_{column}"])
	fig = g.get_figure()
	name=column+".png"
	fig.savefig(name)
	plt.clf()

