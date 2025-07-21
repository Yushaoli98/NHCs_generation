import os
from nhc_generator.detect import find_carbanion_carbon
from nhc_generator.conformer import optimize_nhc
from nhc_generator.geometry import mol_to_coords, select_best_direction
from nhc_generator.aucl_attach import generate_aucl
from nhc_generator.io_utils import write_xyz
from nhc_generator.xtb_wrapper import optimize_with_xtb
from nhc_generator.config import settings
from rdkit import Chem

def process_row(row):
    mol_id, smi = row['ID'], row['Processed_SMILES']
    mol = Chem.AddHs(Chem.MolFromSmiles(smi))
    carbene_idx = find_carbanion_carbon(mol)[0]
    mol_nhc = optimize_nhc(smi, carbene_idx)
    coords = mol_to_coords(mol_nhc)
    elems = [a.GetSymbol() for a in mol_nhc.GetAtoms()]
    nhc_xyz = os.path.join(settings.OUTPUT_DIR, f"{mol_id}_nhc.xyz")
    write_xyz(nhc_xyz, elems, coords)
    direction = select_best_direction(coords[carbene_idx], coords, elems)
    new_atoms, new_coords = generate_aucl(coords[carbene_idx], direction)
    all_coords = np.vstack([coords, new_coords])
    all_elems = elems + new_atoms
    raw_xyz = os.path.join(settings.OUTPUT_DIR, f"{mol_id}.xyz")
    write_xyz(raw_xyz, all_elems, all_coords)
    final_xyz = os.path.join(settings.OUTPUT_DIR, f"{mol_id}_opt.xyz")
    optimize_with_xtb(raw_xyz, final_xyz)