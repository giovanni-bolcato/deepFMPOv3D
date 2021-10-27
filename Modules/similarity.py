from global_parameters import ETA
import Levenshtein
from rdkit.Chem import rdFMCS
from global_parameters import partialCharges
from global_parameters import basisPsi4
from global_parameters import methodPsi4
from mol_utils import neutralize_atoms

import os

import shutil

from rdkit import Chem
from rdkit.Chem import AllChem


from espsim import ConstrainedEmbedMultipleConfs,GetEspSim

import psi4
import resp

import pandas as pd

import numpy as np

import dask.dataframe as dd


def calculateMCStanimoto(ref_mol, target_mol):

    numAtomsRefCpd = float(ref_mol.GetNumAtoms())
    numAtomsTargetCpd = float(target_mol.GetNumAtoms())

    if numAtomsRefCpd < numAtomsTargetCpd:
        leastNumAtms = int(numAtomsRefCpd)
    else:
        leastNumAtms = int(numAtomsTargetCpd)

    pair_of_molecules = [ref_mol, target_mol]
    numCommonAtoms = rdFMCS.FindMCS(pair_of_molecules, 
                                    atomCompare=rdFMCS.AtomCompare.CompareElements,
                                    bondCompare=rdFMCS.BondCompare.CompareOrderExact, matchValences=True).numAtoms
    mcsTanimoto = numCommonAtoms/((numAtomsTargetCpd+numAtomsRefCpd)-numCommonAtoms)

    return mcsTanimoto, leastNumAtms


def EmbedAlignConstrainedScore(prbMol,
                               refMols,
                               core,
                               prbNumConfs = 10,
                               refNumConfs = 10,
                               prbCharge = [],
                               refCharges = [],
                               metric = "carbo",
                               integrate = "gauss",
                               partialCharges = "gasteiger",
                               renormalize = False,
                               customrange = None,
                               marginMC = 10,
                               nMC = 1,
                               basisPsi4 = '3-21G',
                               methodPsi4 = 'scf',
                               gridPsi4 = 1,):
    
    if type(refMols) != list:
        refMols=[refMols]

    if refCharges == []:
        refCharges=[[]]*len(refMols)
        
    prbMol=ConstrainedEmbedMultipleConfs(prbMol, core, numConfs=prbNumConfs)
    for refMol in refMols:
        refMol=ConstrainedEmbedMultipleConfs(refMol, core, numConfs=refNumConfs)
        
    prbMatch = prbMol.GetSubstructMatch(core)
    allShapeDist = []
    allEspSim = []
    
    for idx,refMol in enumerate(refMols):
        shapeDist=1
        prbBestConf=0
        refBestConf=0
        refMatch = refMol.GetSubstructMatch(core)
        for i in range(refNumConfs):
            for j in range(prbNumConfs):
                AllChem.AlignMol(prbMol,refMol,atomMap=list(zip(prbMatch,refMatch)),prbCid=j,refCid=i)
                shape = AllChem.ShapeTanimotoDist(prbMol,refMol,confId1=j,confId2=i)
                if shape<shapeDist:
                    shapeDist=shape
                    prbBestConf=j
                    refBestConf=i
                AllChem.AlignMol(prbMol,refMol,atomMap=list(zip(prbMatch,refMatch)),prbCid=prbBestConf,refCid=refBestConf) 
        allShapeDist.append(1-shapeDist)
        best_conf_1_0=Chem.MolToMolBlock(refMol,confId=refBestConf)
        best_conf_2_0=Chem.MolToMolBlock(prbMol,confId=prbBestConf)
        best_conf_1=Chem.MolFromMolBlock(best_conf_1_0)
        best_conf_2=Chem.MolFromMolBlock(best_conf_2_0)
    return allShapeDist,best_conf_1,best_conf_2

def replace_dummies(smi1):
    frag1_smi=smi1
    frag1_not_cleaned=Chem.AddHs(Chem.MolFromSmiles(frag1_smi))
    chain = Chem.MolFromSmiles('C[n+]1nnnnn1')
    methyl = Chem.MolFromSmiles('C')
    frag1_A=Chem.rdmolops.ReplaceSubstructs(frag1_not_cleaned,Chem.MolFromSmarts('[Yb]'),chain)
    frag1_B=Chem.rdmolops.ReplaceSubstructs(frag1_A[0],Chem.MolFromSmarts('[Lu]'),methyl)
    frag1=Chem.rdmolops.ReplaceSubstructs(frag1_B[0],Chem.MolFromSmarts('[Ta]'),methyl)
    Chem.SanitizeMol(frag1[0])
    frag=Chem.AddHs(frag1[0])
    Chem.SanitizeMol(frag)
    return frag

def core():
    patt=Chem.MolFromSmiles("C[n+]1nnnnn1")
    helper=Chem.AddHs(Chem.MolFromSmiles("C[n+]1nnnnn1"))
    AllChem.EmbedMolecule(helper,AllChem.ETKDG()) 
    AllChem.UFFOptimizeMolecule(helper) 
    core = AllChem.DeleteSubstructs(AllChem.ReplaceSidechains(helper,patt),Chem.MolFromSmiles('*')) 
    core.UpdatePropertyCache()
    return core

def delete_anchor(mol1):
    sub=Chem.MolFromSmarts('[n+]1nnnnn1')
    frag1=AllChem.DeleteSubstructs(mol1,sub)
    Chem.SanitizeMol(frag1)
    return frag1

def delete_dummies(smi):
    mol_a=AllChem.DeleteSubstructs(Chem.MolFromSmiles(smi),Chem.MolFromSmarts('[Yb]'))
    mol_b=AllChem.DeleteSubstructs(mol_a,Chem.MolFromSmarts('[Lu]'))
    mol_c=AllChem.DeleteSubstructs(mol_b,Chem.MolFromSmarts('[Ta]'))
    smi_out=Chem.rdmolfiles.MolToSmiles(mol_c)
    return smi_out

def list_dummies(smi):
    list_dummies=['Yb']
    if 'Lu' in smi:
        list_dummies.append('Lu')
    if 'Ta' in smi:
        list_dummies.append('Ta')
    return list_dummies

def similarity(smi1,smi2,partialCharges,methodPsi4,basisPsi4):
    number_dummies_1=len(list_dummies(smi1))
    number_dummies_2=len(list_dummies(smi2))

    smi1_no_dummies=delete_dummies(smi1)
    smi2_no_dummies=delete_dummies(smi2)


    if number_dummies_1 != number_dummies_2:
        tot_sim=0

    
    elif smi1_no_dummies == smi2_no_dummies:
        tot_sim=0

    else:

        frag1_with_anchor=replace_dummies(smi1)
        frag2_with_anchor=replace_dummies(smi2)
    
        Chem.SanitizeMol(frag1_with_anchor)
        Chem.SanitizeMol(frag2_with_anchor)
    
        frag1_with_anchor_hs_charges=Chem.AddHs(frag1_with_anchor)
        frag2_with_anchor_hs_charges=Chem.AddHs(frag2_with_anchor)

        frag1_with_anchor_hs=neutralize_atoms(frag1_with_anchor_hs_charges)
        frag2_with_anchor_hs=neutralize_atoms(frag2_with_anchor_hs_charges)
    
        simShape,frag1_with_anchor_aligned,frag2_with_anchor_aligned=EmbedAlignConstrainedScore(frag1_with_anchor_hs,frag2_with_anchor,core())
    
        frag1_aligned=delete_anchor(frag1_with_anchor_aligned)
        frag2_aligned=delete_anchor(frag2_with_anchor_aligned)
    
        frag1_aligned_san=Chem.AddHs(frag1_aligned,addCoords=True)
        frag2_aligned_san=Chem.AddHs(frag2_aligned,addCoords=True)
    
        sim_esp=round(GetEspSim(frag1_aligned_san,frag2_aligned_san,partialCharges = partialCharges,metric="tanimoto",methodPsi4=methodPsi4,basisPsi4=basisPsi4,renormalize=True),3)
        simShape_value=simShape[0]
        tot_sim=(sim_esp+simShape_value)/2
    print('tot_sim is',tot_sim, 'between',smi1,'and',smi2)
    return tot_sim
