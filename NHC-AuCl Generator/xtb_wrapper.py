import os, shutil, subprocess
from nhc_generator.config import settings

def optimize_with_xtb(in_xyz, out_xyz):
    workdir = os.path.dirname(out_xyz)
    base = os.path.splitext(os.path.basename(out_xyz))[0]
    tmp_in = os.path.join(workdir, base + "_xtb_in.xyz")
    shutil.copy(in_xyz, tmp_in)
    cmd = [settings.XTB_EXECUTABLE, tmp_in, "--opt", "--chrg", "0", "--gfn", "2"]
    result = subprocess.run(cmd, cwd=workdir, capture_output=True, text=True)
    out_file = os.path.join(workdir, "xtbopt.xyz")
    if os.path.exists(out_file):
        shutil.move(out_file, out_xyz)
    else:
        raise RuntimeError("XTB optimization failed")
