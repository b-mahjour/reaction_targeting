This repository contains all the data and scripts needed generate the figures shown in the manuscript "One-Step Retrosynthesis of Drug Molecules Leveraging C–H Coupling Reactions with Commercially Available Building Blocks". 

Each folder named 'figure' contains the data and script needed to generate the respective figure, as well as the illustrator file used to generate the graphic for the paper.

The manuscript folder contains the paper, the supporting information, and the references file.

The SI folder contains folders corresponding to the figures shown in the SI, with associated data and scripts, if relevant.

The example_notebook folder contains the reaction targeting algorithm used to generate the data used to write the manuscript. The algorithm is self contained in the file "reaction_targeter.py", which can be run from a notebook or via the command line. An example of running the script from the command line is shown in "script_writer.ipynb". The algorithm requires a set of user defined information, which is detailed in the notebook. A catalog or database of molecules is required to serve as the building block pool for synthons to be matched against (not provided). 

```
input_params = {
    "target_molecule": "O=C(C1=C(OC(C)=O)C=CC=C1)O",
    "disconnection_map": [1,2],
    "target_substructure": "[c:1]-[C:2]",
    "substructure_atom_map": [1,2],
    "building_blocks_a": ["[O:3]", "[N:3]", "[C:3](=O)[O]", "[Cl:3]", "[Br:3]", "[I:3]", "[B:3](O)[O]"],
    "building_blocks_b": ["[O:4]", "[N:4]", "[C:4](=O)[O]", "[Cl:4]", "[Br:4]", "[I:4]", "[B:4](O)[O]"],
    "building_blocks_labels_a": ["hydrogen", "alcohol", "amine", "acid", "chloride", "bromide", "iodide", "boronate"],
    "building_blocks_labels_b": ["hydrogen", "alcohol", "amine", "acid", "chloride", "bromide", "iodide", "boronate"],
    "output_path":f"out.xlsx",
    "commercial_data_path": "desalted_sigma_lt_500_mw.xlsx",
    "debug_flag": False
}
```

target_molecule: the SMILES of the molecule that is being analyzed for single step syntheses
disconnection_map: the atom map labels of the atoms associated with the bond being disconnected
target_substructure: the SMARTS of the bond being targeted (sp2 sp3 C-C in this example)
substructure_atom_map: the atom map labels to be used to label the each atom at the disconnected bond (mainly for bookkeeping, keep same as disconnection_map if unsure)
building_blocks_a and building_blocks_b: the allowed building blocks to be enumerated at each disconnection point, encoded as SMARTS
building_blocks_labels_a and building_blocks_labels_b: the associated label for each building block (same length as corresponding building_block list)
output_path: the file to save all output data to
commercial_data_path: the location of the commercial catalog
debug_flag: whether to output debuggin information
