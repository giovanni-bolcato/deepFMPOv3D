import rdkit.Chem as Chem

from rdkit.Chem import PandasTools
import pandas as pd

from espsim import EmbedAlignScore



#insert reference as smile string, example 'c1ccccc1'
##############
smi_reference='Brc1n2c(C(NCc3cncnc3)=CC(c3ccccc3)=C2)nc1'
##############

#insert substructures to delete, example 'c1ccccc1'. If this filter is not desired leave "None"
#############
substructure=None
#############

#n12c(ncc1)C=CC=C2

df=PandasTools.LoadSDF('results_uniq.sdf', idName='Name', molColName='mol')

ref=Chem.rdmolops.AddHs(Chem.MolFromSmiles(smi_reference))

l=[]
for i in range(len(df)):
    l.append(df['mol'].iloc[[i]].item())


if substructure is not None:
    substructure=Chem.MolFromSmiles(substructure)
    matches=[x for x in l if not x.HasSubstructMatch(substructure)]
    df = pd.DataFrame(matches, columns=['mol'])

def calc_espsim(mol):
    simShape,simEsp=EmbedAlignScore(ref,mol)
    print(simEsp)
    return simEsp[0]

df['ESP_sim_reference']=df.mol.apply(calc_espsim)

Chem.PandasTools.WriteSDF(df, 'query_comparison.sdf', molColName='mol', properties=['ESP_sim_reference'])
