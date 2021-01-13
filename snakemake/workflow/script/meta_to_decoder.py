import pandas as pd  
df=pd.read_excel("meta/meta.xlsx",  engine='openpyxl')
df.columns = ["sample.ID", "group.ID", "batch.ID"] 
df.to_csv("meta/decoder.txt", sep="\t", index=False)  
