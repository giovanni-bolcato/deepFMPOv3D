import pandas as pd
from rdkit.Chem import PandasTools
import ast
from ast import literal_eval

df=pd.read_csv('df.csv')

PandasTools.AddMoleculeColumnToFrame(df, smilesCol='smile', molCol='ROMol', includeFingerprints=False)

PandasTools.WriteSDF(df, 'sim.sdf', molColName='ROMol', idName=None, properties=None, allNumeric=False)

with open("History/decodings.txt", "r") as data:
    dictionary = ast.literal_eval(data.read())

l=[]
for i in range(len(dictionary)):
    l.append(i)

list_element_dictionary=[]
list_key=[]
for value in dictionary.items():
    print(value[1])
    list_element_dictionary.append(value[1])
    list_key.append(value[0])

df_3 = pd.DataFrame(list_element_dictionary,columns=['smile'])
df_4 = pd.DataFrame(list_key,columns=['decodings'])

df_merge=pd.concat([df_3, df_4], axis=1)

PandasTools.AddMoleculeColumnToFrame(df_merge, smilesCol='smile', molCol='ROMol', includeFingerprints=False)

PandasTools.WriteSDF(df_merge, 'sim_ordered_decoding.sdf', molColName='ROMol', idName=None, properties=list(df_merge.columns), allNumeric=False)


