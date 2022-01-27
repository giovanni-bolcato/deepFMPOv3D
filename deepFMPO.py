import sys
sys.path.insert(0, './Modules/')

from argsPassing import read_args
from Show_Epoch_new import write_results


def main(fragment_file, lead_file, args):
    import numpy as np
    from file_reader import read_file
    from mol_utils import get_fragments
    from build_encoding import get_encodings, encode_molecule, decode_molecule, encode_list, save_decodings
    from models import build_models
    from training import train
    from rewards import clean_good
    from rdkit import rdBase
    rdBase.DisableLog('rdApp.error')






    print ("Reading files")
    fragment_mols = read_file(fragment_file)
    lead_mols = read_file(lead_file)
    fragment_mols += lead_mols

    print ("Fragmenting molecules")
    fragments, used_mols = get_fragments(fragment_mols)

    print ("Encoding molecules")
    encodings, decodings = get_encodings(fragments)
    save_decodings(decodings)

    lead_mols = np.asarray(fragment_mols[-len(lead_mols):])[used_mols[-len(lead_mols):]]


    X = encode_list(lead_mols, encodings)


    print ("Building models and training")
    actor, critic = build_models(X.shape[1:])

    X = clean_good(X, decodings)

    history = train(X, actor, critic, decodings)

    np.save("History/history.npy", history)




if __name__ == "__main__":


    args = read_args()


    fragment_file = args.fragment_file
    lead_file = args.lead_file

    main(fragment_file, lead_file, args)

    write_results(args.out_file)
