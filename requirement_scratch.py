import numpy as np

KG_FAC = 1.67e-27


def calc_mass(atomicity, reader):
    """
    Python version of calc_mass from mod_1_requirement_scratch.f90.
    atomicity == 0 is treated as monatomic.
    reader must provide next_data_line().
    """
    num_atom = int(reader.next_data_line()[0])
    if atomicity == 0:
        # The original Fortran discards the number-of-atoms line and then reads one mass.
        masses = [float(reader.next_data_line()[0])]
    else:
        masses = [float(reader.next_data_line()[0]) for _ in range(num_atom)]

    mass_tot = sum(masses)
    mass_tot_kg = mass_tot * KG_FAC
    return num_atom, masses, mass_tot, mass_tot_kg


def calc_iner(atomicity, num_atom, mass_tot, mass_atom, reader):
    """
    Python version of calc_iner.
    Returns COM, principal moments I_A, I_B, I_C, and coordinate data.
    I_A/I_B/I_C are in amu Angstrom^2 before unit conversion.
    """
    if atomicity == 0:
        return np.zeros(3), 0.0, 0.0, 0.0, []

    symbols = []
    coords = []

    for _ in range(num_atom):
        parts = reader.next_data_line()
        symbols.append(parts[0])
        coords.append([float(parts[1]), float(parts[2]), float(parts[3])])

    coords = np.array(coords, dtype=float)
    masses = np.array(mass_atom, dtype=float)

    com = np.sum(coords * masses[:, None], axis=0) / mass_tot
    shifted = coords - com

    x = shifted[:, 0]
    y = shifted[:, 1]
    z = shifted[:, 2]

    i_xx = np.sum(masses * (y**2 + z**2))
    i_yy = np.sum(masses * (x**2 + z**2))
    i_zz = np.sum(masses * (x**2 + y**2))
    i_xy = -np.sum(masses * x * y)
    i_xz = -np.sum(masses * x * z)
    i_yz = -np.sum(masses * y * z)

    inertia_tensor = np.array(
        [[i_xx, i_xy, i_xz],
         [i_xy, i_yy, i_yz],
         [i_xz, i_yz, i_zz]],
        dtype=float,
    )

    if i_xy == 0.0 and i_xz == 0.0 and i_yz == 0.0:
        moments = np.array([i_xx, i_yy, i_zz], dtype=float)
    else:
        # LAPACK dsyev in Fortran returns sorted eigenvalues for symmetric matrices.
        moments = np.linalg.eigvalsh(inertia_tensor)

    atoms = [
        {"symbol": sym, "x": float(xyz[0]), "y": float(xyz[1]), "z": float(xyz[2])}
        for sym, xyz in zip(symbols, coords)
    ]

    return com, float(moments[0]), float(moments[1]), float(moments[2]), atoms
