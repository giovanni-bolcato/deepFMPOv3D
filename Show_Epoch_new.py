import sys
sys.path.insert(0, './Modules/')

from global_parameters import EPOCHS
from build_encoding import read_decodings, decode
from rewards import evaluate_chem_mol
import rdkit.Chem as Chem

from rdkit.Chem import PandasTools
import pandas as pd
from rdkit.Chem import Descriptors

import numpy as np
import argparse
from rdkit import rdBase
rdBase.DisableLog('rdApp.error')


def safe_decode(x, decodings):
    try:
        m = decode(x,decodings)
        Chem.Kekulize(m)
        return m
    except:
        return None



def main(epoch):
    decodings2 = read_decodings()


    in_mols = np.load("History/in-{}.npy".format(epoch)) 
    out_mols = np.load("History/out-{}.npy".format(epoch)) 
    in_mols = [decode(m, decodings2) for m in in_mols] 
    out_mols = [safe_decode(m, decodings2) for m in out_mols] 
    use = [(not out_mols[i] is None) and \
                  Chem.MolToSmiles(out_mols[i]) != Chem.MolToSmiles(in_mols[i])
           for i in range(len(out_mols))]
    
    plot_mols = [[m1,m2] for m1,m2,u in zip(in_mols,out_mols,use) if u] 
    
    df = pd.DataFrame(columns=['Name'])
    l_in=[]
    l_out=[]
    for i in plot_mols:
        l_in.append(Chem.MolToSmiles(i[0]))
        l_out.append(Chem.MolToSmiles(i[1]))

    df['in_smiles']=l_in
    df['out_smiles']=l_out

    Name=[]
    for i in range(len(l_in)):
       Name.append(epoch)
    df['Name']=Name

    df.drop_duplicates(subset=['out_smiles'],inplace=True)    
    return df


def calc_molWT(mol):
    molwt=Descriptors.MolWt(mol)
    return molwt

def calc_logp(mol):
    logP=Descriptors.MolLogP(mol)
    return logP

def calc_tpsa(mol):
    tpsa=Descriptors.TPSA(mol)
    return tpsa

def check_smile(smi):
    m = Chem.MolFromSmiles(smi)
    if m is None:
        check='invalid'
    else:
        check ='ok'
    return check

if __name__ == "__main__":

    df_list=[]
    for i in range(1,EPOCHS):
        df=main(i)
        df_list.append(df)
    df = pd.concat(df_list)
    PandasTools.AddMoleculeColumnToFrame(df, smilesCol='out_smiles', molCol='Mol_out', includeFingerprints=False)

    df['check'] = df.out_smiles.apply(check_smile)
    df = df[df.check != 'invalid']
    del df['check']

    df['logP'] = df.Mol_out.apply(calc_logp)
    df['TPSA'] = df.Mol_out.apply(calc_tpsa)
    df['MolWT'] = df.Mol_out.apply(calc_molWT)

    PandasTools.WriteSDF(df, 'tot_duplicati.sdf',molColName='Mol_out' ,idName=None, allNumeric=False, properties=list(df.columns))
    header = ["Name","in_smiles", "out_smiles", "MolWT", "logP","TPSA"]
    del df['Mol_out']
    df.to_csv('tot.csv')
