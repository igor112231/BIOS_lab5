import Bio.PDB
import numpy
from Bio.PDB import Selection
from Bio import PDB

def main():
    pdb_file_ref = "R1107_reference"
    pdb_filename_ref = pdb_file_ref + ".pdb"
    ref_structure = Bio.PDB.PDBParser().get_structure(pdb_file_ref, pdb_filename_ref)
    ref_atom_list = Selection.unfold_entities(ref_structure, "A")

    pdb_file = "R1107TS081_4"
    pdb_filename = pdb_file + ".pdb"
    structure = Bio.PDB.PDBParser().get_structure(pdb_file, pdb_filename)
    atom_list = Selection.unfold_entities(structure, "A")

    superimposer = PDB.Superimposer()
    superimposer.set_atoms(ref_atom_list, atom_list)
    superimposer.apply(structure.get_atoms())

    RMSD=0.0
    for x in range(len(ref_atom_list)):
        diff_vector = atom_list[x].coord - ref_atom_list[x].coord
        diff_vector = numpy.sqrt(numpy.sum(diff_vector * diff_vector))
        diff_vector = diff_vector * diff_vector
        RMSD=RMSD+diff_vector
    RMSD=RMSD/len(ref_atom_list)
    RMSD=numpy.sqrt(RMSD)

    print(RMSD)

if __name__ == '__main__':
    main()