from Bio.PDB import PDBParser, Superimposer


def main():
    # Load the reference structure
    parser = PDBParser(QUIET=True)
    reference_structure = parser.get_structure("reference", "R1107_reference.pdb")

    # Load the file with multiple models
    # all_models_structure = parser.get_structure("all_models", "R1107TS081.pdb")

    # Initialize Superimposer
    superimposer = Superimposer()

    # Get coordinates of the reference structure
    atoms_reference = list(reference_structure.get_atoms())
    atoms_ids_reference = [atom.get_id() for atom in atoms_reference]

    # iterate all models
    for itr in range(1, 5):
        iterator = "_" + str(itr)  # custom filenames
        print(iterator)
        model_structure = load_model(iterator, parser)
        model_atoms_copy = list(model_structure.get_atoms())
        sub_atoms = []
        for ref in atoms_reference:
            for sub in model_atoms_copy:
                if ref.get_parent() == sub.get_parent() and ref.get_id() == sub.get_id():
                    sub_atoms.append(sub)

        # Sort atoms in the model
        sub_atoms.sort(key=lambda x: x.get_id())
        atoms_reference.sort(key=lambda x: x.get_id())

        # Align structures using Superimposer
        superimposer.set_atoms(atoms_reference, sub_atoms)
        superimposer.apply(model_structure.get_atoms())

        # Calculate and print RMSD
        rmsd_value = calculate_rmsd([atom.get_coord() for atom in sub_atoms],
                                    [atom.get_coord() for atom in atoms_reference])
        print(f"RMSD for Model {itr}: {rmsd_value:.4f} Angstrom")
        del sub_atoms

        # Write each model to a separate PDB file
        # output_file = f"model_{}_aligned.pdb"
        # PDBIO().set_structure(model)
        # PDBIO().save(output_file)


def load_model(iterator, parser):
    file_name = f"R1107TS081{iterator}.pdb"
    model_structure = parser.get_structure(iterator, file_name)
    return model_structure


def calculate_rmsd(coords1, coords2):
    # Calculate the RMSD between two sets of coordinates
    n = len(coords1)
    rmsd = 0.0
    for i in range(n):
        rmsd += sum((coords1[i][j] - coords2[i][j]) ** 2 for j in range(3))
    rmsd = (rmsd / n) ** 0.5
    return rmsd


if __name__ == "__main__":
    main()
