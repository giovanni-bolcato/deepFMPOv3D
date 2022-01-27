# DeepFMPO v3D

Code accompanying the paper "On the value of using 3D-shape and electrostatic similarities in deep generative methods". 
The paper can be found at https://doi.org/10.33774/chemrxiv-2021-sqvv9

## Instructions
To run the code on the default dataset use the command:
```sh
python Main.py -f ./Data/molecules.smi -l ./Data/lead.smi -o results.sdf
```

#####
To use other sets of molecules (one of molecules that will be optimized and one of molecules that will be fragmented), use the command:
```sh
python Main.py -f your_fragment.smi -l your_lead.smi -o your_results.sdf
```

#####
In both cases several parameters can be set in the global_parameters.py file.
The most important parameter is the type of charges to be used: gasteiger, mmff or psi4. If psi4 charges are used, methodPsi4 and basisPsi4 must also be set.
Please note that using psi4 charges the calculation requires a lot of time, so this option should be used only with small datasets. 
In The same file the target values can also be setted.

#####
To save the generated fragments in an sdf file use the script decoding_to_sdf.py, with the command:

```sh
python decoding_to_sdf.py
```

The fragments are reported in the same way in which they are present in the similarity tree.

## Requirements

The program is written in Python 3.7 using the following Python libraries:
- rdkit
- scipy
- resp
- psi4
- dask
- ESP-Sim
- numpy
- sklearn
- keras
- pandas
- bisect
- Levenshtein
- A backend to keras, such as theano, tensorflow or CNTK

Please note that the ESP-sim can be obtained from the "Comparison of electrostatic potential and shape" GitHub repo:
https://github.com/hesther/espsim
