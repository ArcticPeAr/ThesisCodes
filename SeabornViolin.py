"""
The methylation level of CpGs in the topics, you could plot the complete distribution of values for a given set of patients associated with the corresponding topic and compare to other CpGs in the other topics. For instance, using violin plots with strips
"""
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import pyreadr

######################################################################################################
#### Define functions
######################################################################################################
def read_rds(file_path):
    """
    Reads an RDS file and returns it as a DataFrame.
    
    Parameters:
        file_path (str): The path to the RDS file.
    
    Returns:
        pd.DataFrame: The DataFrame containing the data.
    """
    rds_data = pyreadr.read_r(file_path)
    return list(rds_data.values())[0]



def get_topic_patient_dict(df):
    """
    Creates a dictionary where keys are the topics and the values are lists of patients 
    that have that corresponding topic as their highest contributing topic.
    
    Parameters:
        df (pd.DataFrame): DataFrame containing topic assignment for patients.
        
    Returns:
        dict: Dictionary with topics as keys and list of patients as values.
    """
    
    topics_dict = {f"Topic_{i+1}": [] for i in range(len(df))}
    
    for col in df:
        # Get the topic number that the patient (column) has the highest score for
        topic_number = df[col].idxmax()
        topic_label = f"Topic_{topic_number + 1}"
        
        # Append the patient to the list of patients for the topic they score highest on
        topics_dict[topic_label].append(col)
        
    return topics_dict



def remove_empty_lists_from_dict(input_dict):
    """
    Removes keys that have empty lists as values from the dictionary.
    
    Parameters:
        input_dict (dict): The input dictionary.
        
    Returns:
        dict: The dictionary with keys having empty lists removed.
    """
    return {k: v for k, v in input_dict.items() if v}



def collect_values_from_dict_values(dict_values, df):
    """
    Collects values from a DataFrame based on the values in a dictionary.
    
    Parameters:
        dict_values (dict): The dictionary containing column names as values in lists.
        df (pd.DataFrame): The pandas DataFrame to collect values from.
    
    Returns:
        dict: A dictionary with the same keys as dict_values, but with lists of values as values. 
              Each list contains the values from the corresponding columns in df specified in dict_values, 
              with the values from different columns concatenated into a single list.
    """
    collected_values = {}
    
    for key, values in dict_values.items():
        collected_values[key] = [val for col in values if col in df.columns for val in df[col].tolist()]
    
    return collected_values



def create_violin_plots(data_dict):
    """
    Create a series of violin plots from a dictionary of data.
    
    The function takes a dictionary where keys are topics and values are lists 
    of servations for that topic. The dictionary 
    is converted to a pandas DataFrame, then melted to create a long-form 
    DataFrame, which is used to create the violin plots using Seaborn's 
    violinplot function.
    
    Parameters:
    data_dict (dict): A dictionary where the keys are strings representing 
                      topic names and the values are lists of numeric values 
                      representing observations for that topic.
    
    Returns:
    None: The function displays a violin plot using matplotlib but does not 
          return any values.
    
    Example:
    data_dict = {"Topic1": [1, 2, 3, 4, 5], "Topic2": [5, 6, 7, 8, 9]}
    create_violin_plots(data_dict)
    """
    # Create a new DataFrame from the dictionary
    df = pd.DataFrame.from_dict(data_dict, orient='index').transpose()
    
    # Melt the DataFrame to create a long-form DataFrame
    df_melted = df.melt(var_name='Topics', value_name='Values')
    
    # Create a violin plot using the long-form DataFrame
    plt.figure()
    sns.violinplot(x='Topics', y='Values', data=df_melted)
    
    # Set plot title and labels
    plt.title('Violin plots for all topics')
    plt.xlabel('Topics')
    plt.ylabel('Values')
    
    # Rotate x labels for better visibility
    plt.xticks(rotation=90)
    
    # Display the plot
    plt.show()



def save_current_plot():
    """Save the current plot to the PDF and close the figure."""
    pdf.savefig()
    plt.close()


def plot_topic_distribution_violin(df, scale=None):
    """
    Plot a violin plot showing topic distributions across documents.
    
    Parameters:
    - df (pd.DataFrame): DataFrame containing topic distributions.
    - scale (str): If "log", will plot the logarithm of the values.
    """
    melted_df = df.reset_index().melt(id_vars='index', var_name='Document', value_name='Value')
    ylabel = 'Value'
    
    if scale == 'log':
        melted_df['Value'] = np.log(melted_df['Value'])
        ylabel = 'ln(Value)'
        
    plt.figure(figsize=(12, 8))
    sns.violinplot(x='index', y='Value', data=melted_df)
    plt.xlabel('Topics', labelpad=10, fontsize=10, fontweight='bold')
    plt.ylabel(ylabel, labelpad=10, fontsize=10, fontweight='bold')
    plt.title('Topic Distribution with Violin Plot')
    
    save_current_plot()


def plot_topic_distribution_boxplot(df, scale=None):
    """
    Plot a single boxplot showing topic distributions across multiple documents for all topics.
    
    Parameters:
    - df (pd.DataFrame): DataFrame where each row corresponds to a topic and each column to a document. 
                         Cells contain numeric values representing the relevance score of a topic in a document.
    - scale (str): If "log", will plot the logarithm of the values.
    """
    melted_df = df.T.melt(var_name='Topic', value_name='Value')
    ylabel = 'Value'
    
    if scale == 'log':
        melted_df['Value'] = np.log(melted_df['Value'])
        ylabel = 'ln(Value)'
        
    plt.figure(figsize=(12, 8))
    sns.boxplot(x='Topic', y='Value', data=melted_df)
    plt.xlabel('Topics', labelpad=10, fontsize=10, fontweight='bold')
    plt.ylabel(ylabel, labelpad=10, fontsize=10, fontweight='bold')
    plt.title('Topic Distributions Across Documents (Boxplot)')
    
    save_current_plot()

def plot_topic_distribution_histogram(df):
    """
    Plot histograms showing topic distributions across multiple documents for each topic.
    
    Parameters:
    - df (pd.DataFrame): DataFrame where each row corresponds to a topic and each column to a document. 
                         Cells contain numeric values representing the relevance score of a topic in a document.
    """
    topics = df.index.tolist()
    n = len(topics)
    max_val = df.values.max()
    
    fig, axes = plt.subplots(n, 1, figsize=(10, 6 * n))
    
    for i, topic in enumerate(topics):
        sns.histplot(df.loc[topic], kde=False, ax=axes[i], bins=20)
        axes[i].set_xlim(0, max_val)
        axes[i].set_title(f'Topic: {topic}', loc='right')
        axes[i].set_ylabel('Frequency')
        axes[i].set_xlabel('Value')
    
    fig.suptitle('Topic Distributions Across Documents', fontsize=16)
    plt.subplots_adjust(left=0.15, hspace=0.3)
    
    save_current_plot()

######################################################################################################
#### Load data
######################################################################################################

meth_table = read_rds(snakemake.input["methTable2rds"])
meth_table.index.names = [None]

reg_scr_prtopic = read_rds(snakemake.input["RegScrPrtopicOut"])
reg_scr_prtopic.index.names = [None]

reg_assig_unormal = read_rds(snakemake.input["RegAssigUnormalOut"])
reg_assig_unormal.index.names = [None]

cell_model_mat <- read_rds(snakemake.input["cellModelMatOut"]) 

region_model_mat <- read_rds(snakemake.input["regionModelMatOut"])

######################################################################################################
#### Create plots
######################################################################################################
pdf = PdfPages('all_plots.pdf')

plot_topic_distribution_violin(region_model_mat, "log")
plot_topic_distribution_violin(cell_model_mat)

plot_topic_distribution_boxplot(region_model_mat, "log")
plot_topic_distribution_boxplot(cell_model_mat)

plot_topic_distribution_histogram(region_model_mat)
plot_topic_distribution_histogram(cell_model_mat)





pdf.close()


# topics_dict = get_topic_patient_dict(topic_assig_to_patient)

# topics_dict = remove_empty_lists_from_dict(topics_dict)

# meth_table.index.names = [None]
# topics_dict = collect_values_from_dict_values(topics_dict, meth_table)

# #create violin plots
# create_violin_plots(topics_dict)

save_plot_as_pdf("violin_plots.pdf")




######################################################################################################
#### Create plots
######################################################################################################

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

