import pandas as pd
import re
import sys

metaphlan_output = sys.argv[1]
df = pd.read_table(metaphlan_output, sep = '\t')
#df = pd.read_table("/Users/brandonbuenohernandez/Downloads/subsamples_merge.txt", sep = '\t')

df_taxa = df['ID'].str.split( '|', expand=True)
taxa_cols = ["Kingdom","Phylum","Class","Order","Family","Genus","Species","Strain"]
taxa_dict = {'Kingdom':1, 'Phylum':2, 'Class':3, 'Order':4, 'Family':5, 'Genus':6, 'Species':7, 'Strain':8}
rank = "Strain"
value = taxa_dict.get(rank)
taxa_cols = taxa_cols[0:value]
df_taxa.columns = taxa_cols

def trim_taxa_names(x):
    match = re.sub(r'^[kpcofgs]__',"",str(x))
    return match

for col in df_taxa.columns:
    df_taxa[col]=df_taxa[col].apply(trim_taxa_names)

df_taxa.to_csv("taxa_table.txt", index=False, sep = '\t')
print(df_taxa)
