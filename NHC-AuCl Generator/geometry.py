import numpy as np
from scipy.spatial import KDTree
from nhc_generator.config import settings

def mol_to_coords(mol):
    return np.array([list(mol.GetConformer().GetAtomPosition(i)) for i in range(mol.GetNumAtoms())])

def get_vdw_radius(elem):
    return {"H":1.2,"C":1.7,"N":1.55,"O":1.52,"Cl":1.75,"Au":1.66}.get(elem, 1.7)

def select_best_direction(carb_xyz, coords, elems):
    dirs = uniform_sphere_directions(settings.NUM_DIRS)
    nbrs = np.argsort(np.linalg.norm(coords - carb_xyz, axis=1))[1:3]
    normal = np.cross(coords[nbrs[0]] - carb_xyz, coords[nbrs[1]] - carb_xyz)
    normal /= np.linalg.norm(normal)
    kdt = KDTree(coords)
    best_score, best_dir = -1, normal
    for u in dirs:
        if np.abs(np.degrees(np.arccos(np.dot(u, normal)))) > 60:
            continue
        au_pos = carb_xyz + u * settings.AU_C_DIST
        cl_pos = au_pos + u * settings.CL_DIST
        min_dist = min(kdt.query(au_pos)[0], kdt.query(cl_pos)[0])
        if min_dist > best_score:
            best_score, best_dir = min_dist, u
    return best_dir

def uniform_sphere_directions(n):
    phi = np.arccos(1 - 2 * (np.arange(n) + 0.5) / n)
    theta = np.pi * (1 + 5**0.5) * (np.arange(n) + 0.5)
    return np.vstack([np.cos(theta)*np.sin(phi), np.sin(theta)*np.sin(phi), np.cos(phi)]).T