from rdkit import Chem
from rdkit.Chem import AllChem
from nhc_generator.config import settings

def optimize_nhc(smiles, carbene_idx=None):
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)

    if carbene_idx is not None:
        rw = Chem.RWMol(mol)
        for nbr in mol.GetAtomWithIdx(carbene_idx).GetNeighbors():
            if nbr.GetSymbol() == 'H':
                rw.RemoveAtom(nbr.GetIdx())
        mol = rw.GetMol()

    try:
        return generate_lowest_energy_conformer(mol)
    except:
        AllChem.EmbedMolecule(mol, AllChem.ETKDG())
        AllChem.MMFFOptimizeMolecule(mol)
        return mol

def generate_lowest_energy_conformer(mol, num_confs=10, maxIters=200):
    params = AllChem.ETKDGv3()
    ids = AllChem.EmbedMultipleConfs(mol, numConfs=num_confs, params=params)
    energies = []
    for cid in ids:
        try:
            ff = AllChem.UFFGetMoleculeForceField(mol, confId=cid)
            ff.Minimize(maxIts=maxIters)
            energies.append((ff.CalcEnergy(), cid))
        except:
            continue
    best_id = sorted(energies)[0][1]
    mol.RemoveAllConformers()
    mol.AddConformer(mol.GetConformer(best_id), assignId=True)
    return mol