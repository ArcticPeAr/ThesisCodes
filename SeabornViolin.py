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
meth_table = snakemake.input["methTable2rds"]
reg_scr_norm = snakemake.input["RegScrPrtopicOut"]
reg_scr_unrm = snakemake.input["RegAssigUnormalOut"]
topic_assig = snakemake.input["topicAssigToPatientOut"]

# read RDS from R into Python
rds_read_mt = pyreadr.read_r(meth_table) 
meth_table = list(rds_read_mt.values())[0]
meth_table.index.names = [None]

# read RDS in python
# remove rownames 
rds_read_rspt = pyreadr.read_r(reg_scr_norm)
reg_scr_prtopic = list(rds_read_rspt.values())[0]

rds_read_rsu = pyreadr.read_r(reg_scr_unrm)
reg_assig_unormal = list(rds_read_rsu.values())[0]

rds_read_ta = pyreadr.read_r(topic_assig)
topic_assig_to_patient = list(rds_read_ta.values())[0]

####################################

"""
Following to find the donors that have the highest topic assignment:
"""
# set row as index
topic_assig_to_patient["index"] = range(1, len(topic_assig_to_patient) + 1)
topic_assig_to_patient.set_index('index', inplace=True)

# Find number of topics 
tpc_no = len(topic_assig_to_patient)

nr = 0
topics_dict = {}

for nr in range(1, tpc_no + 1):
    topics_dict[f"Topic_{nr}"] = []
    
# Make list with Topic_X as name of list (x = topic number) and patients as items where the patients had that topic as highest contributing topic. This will cluster patients that has the same "highest contributing topic" together
transposed_df = topic_assig_to_patient.transpose()

# Iterate through the transposed DataFrame
for patient, row in transposed_df.iterrows():
    # Identify the topic with the maximum value
    highest_contrib_topic = row.idxmax()
    topic_number = highest_contrib_topic + 1  # Adding 1 because the index starts from 0
    # Append the patient to the corresponding list in the topics_dict
    topics_dict[f"Topic_{topic_number}"].append(patient)


    
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



#####################################################################################
#SNS
#####################################################################################


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

#Topics on x-axis
import pandas as pd

melted_df = topicAssigToPatient.reset_index().melt(id_vars="index", 
                                                   value_vars=topicAssigToPatient.columns, 
                                                   var_name="Patient", 
                                                   value_name="Value")
melted_df.rename(columns={"index": "Topic"}, inplace=True)

import seaborn as sns
import matplotlib.pyplot as plt

plt.figure(figsize=(15,10))
sns.violinplot(x="Topic", y="Value", data=melted_df)
plt.title("Distribution of Topic Assignments Across Patients")
plt.xticks(rotation=90)
plt.tight_layout()
plt.show()

#Patients on x-axis


melted_df = topicAssigToPatient.reset_index().melt(id_vars="index", 
                                                   value_vars=topicAssigToPatient.columns, 
                                                   var_name="Patient", 
                                                   value_name="Value")
melted_df.rename(columns={"index": "Topic"}, inplace=True)

plt.figure(figsize=(20,10))
sns.violinplot(x="Patient", y="Value", hue="Topic", data=melted_df, inner="quartile")
plt.title("Distribution of Topic Assignments for Each Patient")
plt.xticks(rotation=90)
plt.tight_layout()
plt.legend(title='Topic', bbox_to_anchor=(1.05, 1), loc='upper left')
plt.show()

plt.figure(figsize=(20,10))
sns.violinplot(x="Patient", y="Value", data=melted_df, inner="quartile", scale="width")
plt.title("Distribution of Topic Assignments for Each Patient")
plt.xticks(rotation=90)
plt.tight_layout()
plt.show()

