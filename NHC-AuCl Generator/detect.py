from rdkit import Chem


def find_carbanion_carbon(mol):
    pattern = Chem.MolFromSmarts("[C-]")
    matches = mol.GetSubstructMatches(pattern)
    return [m[0] for m in matches]