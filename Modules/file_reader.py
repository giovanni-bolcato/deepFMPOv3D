from mol_utils import drop_salt
from rdkit import Chem

# Read a file containing SMILES
# The file should be a .smi or a .csv where the first column should contain a SMILES string
def read_file(file_name, drop_first=True):
    mol_list=[]
    molObjects = []
    with open(file_name) as f:
        for l in f:
            if drop_first:
                drop_first = False
                continue
            l = l.strip().split(",")[0]
            smi = drop_salt(l.strip())
            print('reading', smi)
            if "/" in smi:
                print(smi,'discarded')
            else:
                molObjects.append(Chem.MolFromSmiles(smi))
                mol_list.append(smi)

    with open("readed.txt", "w") as output:
        output.write(str(mol_list))

    return molObjects

