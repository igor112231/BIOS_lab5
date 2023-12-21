import Bio.PDB
import numpy
from Bio.PDB import Selection


def main():

    pdb_file = input("Provide filename of RNA structure: ")
    pdb_filename = pdb_file + ".pdb"
    structure = Bio.PDB.PDBParser().get_structure(pdb_file, pdb_filename)
    atom_list = Selection.unfold_entities(structure, "A")
    Clash_score = Clash_Score_Calculation(atom_list)

    print(f"Clash Score for model {pdb_file}: {Clash_score}")


def Clash_Score_Calculation(atom_list):
    bad_overlaps = 0.0
    for x in range(len(atom_list)):
        for y in range(x+1, len(atom_list)):
            diff_vector = atom_list[x].coord - atom_list[y].coord
            diff_vector = numpy.sqrt(numpy.sum(diff_vector * diff_vector))
            if (diff_vector < 0.4):
                bad_overlaps += 1

    clash_score = 1000*bad_overlaps/len(atom_list)
    return clash_score


if __name__ == '__main__':
    main()
