def write_xyz(path, elems, coords):
    with open(path, 'w') as f:
        f.write(f"{len(elems)}\nGenerated\n")
        for e, (x, y, z) in zip(elems, coords):
            f.write(f"{e} {x:.6f} {y:.6f} {z:.6f}\n")