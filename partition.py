import math

PI = math.pi
K_B = 1.380649e-23
H = 6.626e-34
C = 2.99e10
R = 8.314
KJ_2_J_FAC = 1000.0
ANG_2_M = 1.0e-10
KG_FAC = 1.67e-27
TOLERANCE = 0.3


def trans_parti(mass_tot_kg, vol_fac, temp):
    return (((2 * PI * mass_tot_kg * K_B * temp) ** 1.5) / (H**3)) / vol_fac


def rot_parti(atomicity, sym_fac, temp, i_a, i_b, i_c):
    i_a_unit = i_a * KG_FAC * (ANG_2_M**2)
    i_b_unit = i_b * KG_FAC * (ANG_2_M**2)
    i_c_unit = i_c * KG_FAC * (ANG_2_M**2)

    if atomicity == 0:
        rot_val = 1.0
    elif abs(i_a) < TOLERANCE and abs(i_b - i_c) < TOLERANCE:
        # Linear rotor
        rot_val = (8 * PI**2 * i_c_unit * K_B * temp) / (sym_fac * H**2)
    elif abs(i_a - i_b) < TOLERANCE and abs(i_b - i_c) < TOLERANCE:
        # Spherical top
        rot_val = ((8 * PI**2 * i_c_unit * K_B * temp) / H**2) ** 1.5
        rot_val *= math.sqrt(PI) / sym_fac
    elif i_a < i_b and abs(i_b - i_c) < TOLERANCE:
        # Prolate symmetric top
        rot_val = ((8 * PI**2 * K_B * temp) / H**2) ** 1.5
        rot_val *= math.sqrt(PI * i_a_unit * i_c_unit**2) / sym_fac
    elif i_a < i_c and abs(i_a - i_b) < TOLERANCE:
        # Oblate symmetric top
        rot_val = ((8 * PI**2 * K_B * temp) / H**2) ** 1.5
        rot_val *= math.sqrt(PI * i_c_unit * i_b_unit**2) / sym_fac
    elif i_a < i_b < i_c and (i_b - i_a) > 0.4 and (i_c - i_a) > 0.4:
        # Asymmetric top
        rot_val = ((8 * PI**2 * K_B * temp) / H**2) ** 1.5
        rot_val *= math.sqrt(PI * i_c_unit * i_b_unit * i_a_unit) / sym_fac
    else:
        # Fallback for cases not covered by the original if/else tree.
        rot_val = ((8 * PI**2 * K_B * temp) / H**2) ** 1.5
        rot_val *= math.sqrt(max(0.0, PI * i_c_unit * i_b_unit * i_a_unit)) / sym_fac

    return rot_val, i_a_unit, i_b_unit, i_c_unit


def vib_parti(atomicity, temp, reader):
    if atomicity == 0:
        return 1.0, []

    freq_num = int(reader.next_data_line()[0])
    freq = [float(reader.next_data_line()[0]) for _ in range(freq_num)]

    factors = [(H * (f * C)) / (K_B * temp) for f in freq]

    if freq_num == 1:
        vib_val = 1.0 / (1.0 - math.exp(-factors[0]))
    elif freq_num == 2:
        # This follows the original Fortran: second frequency is squared.
        vib_val = 1.0 / ((1.0 - math.exp(-factors[0])) * (1.0 - math.exp(-factors[1])) ** 2)
    else:
        # Generalized product for more than two frequencies.
        vib_val = 1.0
        for fac in factors:
            vib_val *= 1.0 / (1.0 - math.exp(-fac))

    return vib_val, freq


def parti(trans_val, rot_val, vib_val):
    return trans_val * rot_val * vib_val


def rate(parti_fac, e_0, temp):
    return ((K_B * temp) / H) * parti_fac * math.exp(-((e_0 * KJ_2_J_FAC) / (R * temp)))
