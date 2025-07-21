import numpy as np
from nhc_generator.config import settings

def generate_aucl(carb_xyz, direction):
    u = direction / np.linalg.norm(direction)
    au = carb_xyz + u * settings.AU_C_DIST
    cl = au + u * settings.CL_DIST
    return ['Au', 'Cl'], np.vstack([au, cl])