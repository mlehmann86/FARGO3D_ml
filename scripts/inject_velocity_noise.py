import numpy as np
import os


def read_parameters(filename, IDEFIX=False):
    """
    Read and parse simulation parameters.
    Supports FARGO3D (default) and IDEFIX (if IDEFIX=True).
    """
    print('READING SIMULATION PARAMETERS')
    parameters = {}
    with open(filename, 'r', encoding='utf-8') as f:
        lines = f.readlines()

    if IDEFIX:
        in_params = False
        for line in lines:
            line = line.strip()
            if "Input Parameters using input file" in line:
                in_params = True
                continue
            if in_params:
                if line.startswith("---") or line == "":
                    continue
                if line.startswith("["):  # Section header
                    current_section = line.strip("[]")
                    continue

                # Special handling for smoothing line with plummer
                if line.startswith("smoothing") and "plummer" in line:
                    parts = line.split()
                    try:
                        smoothing_val = float(parts[2])  # This is the 0.02 value
                        parameters['smoothing'] = smoothing_val
                    except (IndexError, ValueError):
                        parameters['smoothing'] = None
                    continue

                # Split into key-value, special case for grid lines
                if 'grid' in line.lower():
                    parts = line.split()
                    key = parts[0]  # e.g., X1-grid
                    parameters[key] = parts[1:]  # store the full list of values
                elif "=" in line:
                    key, val = map(str.strip, line.split("=", 1))
                    try:
                        parameters[key] = float(val)
                    except ValueError:
                        parameters[key.upper()] = val
                else:
                    parts = line.split()
                    if len(parts) >= 2:
                        key = parts[0]
                        val = parts[1]
                        try:
                            parameters[key] = float(val)
                        except ValueError:
                            parameters[key.upper()] = val

    else:
        # FARGO3D parser logic (left unchanged)
        in_parameters_section = False
        for line in lines:
            line = line.strip()
            if "PARAMETERS SECTION:" in line:
                in_parameters_section = True
            elif in_parameters_section and line:
                if line.startswith('=') or line == "":
                    continue
                parts = line.split(maxsplit=1)
                if len(parts) == 2:
                    key, value = parts
                    value = value.strip().split()[0]  # Extract the numeric value, strip any comment
                    try:
                        value = float(value)
                    except ValueError:
                        pass
                    if key.upper() not in parameters:
                        parameters[key.upper()] = value

    return parameters

########################################################################################################################################

def reconstruct_grid(parameters, IDEFIX=False):
    """Reconstruct grid based on parameters. Supports both FARGO3D and IDEFIX modes."""
    print('RECONSTRUCTING GRIDS')

    if IDEFIX:
        # Extract from X1-grid, X2-grid, etc.
        if 'X1-grid' not in parameters:
            raise ValueError("X1-grid not found in parameters.")
        if 'X2-grid' not in parameters:
            parameters['X2-grid'] = ['1', '0.0', '1', 'u', '1.0']
        if 'X3-grid' not in parameters:
            parameters['X3-grid'] = ['1', '0.0', '1', 'u', '1.0']

        x1grid = parameters['X1-grid']
        x2grid = parameters['X2-grid']
        x3grid = parameters['X3-grid']

        nx = int(x1grid[2])
        ny = int(x2grid[2])
        nz = int(x3grid[2])

        r_min = float(x1grid[1])
        r_max = float(x1grid[4])
        y_min = float(x2grid[1])
        y_max = float(x2grid[4])
        z_min = float(x3grid[1])
        z_max = float(x3grid[4])
    else:
        ny = int(parameters['NX'])
        nx = int(parameters['NY'])
        nz = int(parameters['NZ'])

        y_min = parameters['XMIN']
        y_max = parameters['XMAX']
        r_min = parameters['YMIN']
        r_max = parameters['YMAX']
        z_min = parameters['ZMIN']
        z_max = parameters['ZMAX']

    xgrid = np.linspace(r_min, r_max, nx)
    ygrid = np.linspace(y_min, y_max, ny)
    zgrid = np.linspace(z_min, z_max, nz)

    print(f"xgrid: [{r_min:.3f}, {r_max:.3f}], nx={nx}")
    print(f"ygrid: [{y_min:.3f}, {y_max:.3f}], ny={ny}")
    print(f"zgrid: [{z_min:.3f}, {z_max:.3f}], nz={nz}")

    return xgrid, ygrid, zgrid, ny, nx, nz




def inject_velocity_noise(sim_path, strength=0.1, verbose=True):
    logfile = os.path.join(sim_path, "inject_cos_noise.debug.log")
    sys.stdout = open(logfile, "a")
    sys.stderr = sys.stdout
    print(f"\n--- New injection run ---")
    print(f"[inject_velocity_noise] sim_path = {sim_path}")
    print(f"[inject_velocity_noise] strength = {strength}")
    """
    Injects white noise (strength * c_s(r)) into the initial velocity fields of a Fargo3D simulation.
    Requires data_reader.py utilities for param and grid reading.
    """

    print(f"[inject_velocity_noise] sim_path = {sim_path}")
    print(f"[inject_velocity_noise] strength = {strength}")

    # === Read parameters and grid ===
    summaryfile = os.path.join(sim_path, "summary0.dat")
    if not os.path.isfile(summaryfile):
        raise FileNotFoundError(f"Parameter file not found: {summaryfile}")
    parameters = read_parameters(summaryfile)
    xgrid, ygrid, zgrid, ny, nx, nz = reconstruct_grid(parameters)

    shape_disk = (nz, nx, ny)  # FARGO3D binary order

    # === File paths ===
    vxfile = os.path.join(sim_path, "gasvx0.dat")
    vyfile = os.path.join(sim_path, "gasvy0.dat")
    vzfile = os.path.join(sim_path, "gasvz0.dat")
    densfile = os.path.join(sim_path, "gasdens0.dat")
    enefile = os.path.join(sim_path, "gasenergy0.dat")
    for f in [vxfile, vyfile, vzfile, densfile, enefile]:
        if not os.path.isfile(f):
            raise FileNotFoundError(f"{f} not found!")

    # === Read arrays and transpose to (ny, nx, nz) ===
    dens = np.fromfile(densfile, dtype=np.float64).reshape(shape_disk).transpose(2, 1, 0)
    ene  = np.fromfile(enefile, dtype=np.float64).reshape(shape_disk).transpose(2, 1, 0)
    vx   = np.fromfile(vxfile, dtype=np.float64).reshape(shape_disk).transpose(2, 1, 0)
    vy   = np.fromfile(vyfile, dtype=np.float64).reshape(shape_disk).transpose(2, 1, 0)
    vz   = np.fromfile(vzfile, dtype=np.float64).reshape(shape_disk).transpose(2, 1, 0)

    # === Compute c_s(r) profile ===
    cs2 = ene / dens
    cs2_profile = np.mean(cs2, axis=(0, 2))  # (nx,)
    cs_profile = np.sqrt(cs2_profile)

    if verbose:
        print(f"Sound speed profile shape: {cs_profile.shape} (nx={nx})")
        print(f"Example c_s values: min={cs_profile.min():.3e}, max={cs_profile.max():.3e}")

    # === Inject noise ===
    rng = np.random.default_rng()
    for vfield, fname in zip([vx, vy, vz], ["gasvx0.dat", "gasvy0.dat", "gasvz0.dat"]):
        noise = (2.0 * rng.random(size=vfield.shape) - 1.0)  # shape (ny, nx, nz)
        for ix in range(nx):
            noise[:, ix, :] *= strength * cs_profile[ix]
        vfield += noise
        # Write back as (nz, nx, ny)
        outarr = vfield.transpose(2, 1, 0)
        outpath = os.path.join(sim_path, fname)
        outarr.astype(np.float64).tofile(outpath)
        if verbose:
            print(f"Injected noise and overwrote {fname}")

    print("All velocity fields updated with random noise.")

if __name__ == "__main__":
    import sys
    print(f"[DEBUG] sys.argv = {sys.argv}")
    if len(sys.argv) < 2:
        print("Usage: python3 inject_cos_noise.py <sim_path>")
        sys.exit(1)
    sim_path = sys.argv[1]
    inject_velocity_noise(sim_path, strength=0.005)
