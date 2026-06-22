from pathlib import Path

from requirement_scratch import calc_mass, calc_iner
from partition import trans_parti, rot_parti, vib_parti, parti, rate, TOLERANCE


class InputReader:
    def __init__(self, filename):
        self.lines = Path(filename).read_text().splitlines()
        self.index = 0

    def next_data_line(self):
        while self.index < len(self.lines):
            raw = self.lines[self.index]
            self.index += 1
            clean = raw.split("!", 1)[0].strip()
            if clean:
                return clean.split()
        raise EOFError("Unexpected end of input file.")


def top_type(atomicity, i_a, i_b, i_c):
    if atomicity == 0 and i_a == 0.0 and i_b == 0.0 and i_c == 0.0:
        return "single atom"
    if atomicity > 0 and abs(i_a) < TOLERANCE and abs(i_b - i_c) < TOLERANCE:
        return "linear"
    if atomicity > 0 and abs(i_a - i_b) < TOLERANCE and abs(i_b - i_c) < TOLERANCE:
        return "spherical top"
    if atomicity > 0 and i_a < i_b and abs(i_b - i_c) < TOLERANCE:
        return "prolate symmetric top"
    if atomicity > 0 and i_a < i_c and abs(i_a - i_b) < TOLERANCE:
        return "oblate symmetric top"
    if atomicity > 0 and i_a < i_b < i_c and (i_b - i_a) > 0.4 and (i_c - i_a) > 0.4:
        return "asymmetric top"
    return "unclassified top"


def atomicity_name(atomicity):
    if atomicity == 0:
        return "Monatomic"
    if atomicity == 1:
        return "Diatomic"
    return "Polyatomic"


def read_species(reader, temp):
    atomicity = int(reader.next_data_line()[0])

    num_atom, masses, mass_tot, mass_tot_kg = calc_mass(atomicity, reader)
    com, i_a, i_b, i_c, atoms = calc_iner(atomicity, num_atom, mass_tot, masses, reader)

    vol_fac = float(reader.next_data_line()[0])
    sym_fac = int(reader.next_data_line()[0])

    trans_val = trans_parti(mass_tot_kg, vol_fac, temp)
    rot_val, i_a_unit, i_b_unit, i_c_unit = rot_parti(atomicity, sym_fac, temp, i_a, i_b, i_c)
    vib_val, freqs = vib_parti(atomicity, temp, reader)
    parti_val = parti(trans_val, rot_val, vib_val)

    return {
        "atomicity": atomicity,
        "num_atom": num_atom,
        "masses": masses,
        "mass_tot": mass_tot,
        "mass_tot_kg": mass_tot_kg,
        "atoms": atoms,
        "com": com,
        "I_A": i_a,
        "I_B": i_b,
        "I_C": i_c,
        "I_A_unit": i_a_unit,
        "I_B_unit": i_b_unit,
        "I_C_unit": i_c_unit,
        "vol_fac": vol_fac,
        "sym_fac": sym_fac,
        "trans_val": trans_val,
        "rot_val": rot_val,
        "vib_val": vib_val,
        "freqs": freqs,
        "parti_val": parti_val,
        "top_type": top_type(atomicity, i_a, i_b, i_c),
    }


def write_species(out, label, species):
    out.write(f"________________________________Information about {label}____________________________________________________\n")
    out.write(f"{label} is {atomicity_name(species['atomicity'])}.\n")
    out.write(f"Total mass of {label} is {species['mass_tot_kg']:.8e} kg.\n")

    if species["atomicity"] == 0:
        out.write(f"Again, {label} is monatomic.\n")
    else:
        out.write(
            f"Moments of inertia (kg m2) of {label} are: "
            f"I_A={species['I_A_unit']:.8e}, "
            f"I_B={species['I_B_unit']:.8e}, "
            f"I_C={species['I_C_unit']:.8e}\n"
        )

    out.write(f"Volume factor of {label} is {species['vol_fac']} m3.\n")
    out.write(f"Symmetry factor of {label} is {species['sym_fac']}.\n")
    out.write(f"Translational partition function of {label} is {species['trans_val']:.8e} m-3.\n")

    if species["top_type"] == "single atom":
        out.write(f"Your {label} has only a single atom !!!\n")
    else:
        out.write(f"{label} is {species['top_type']}.\n")
        out.write(f"Rotational partition function of {label} is {species['rot_val']:.8e}.\n")

    if species["atomicity"] == 0:
        out.write(f"Again, your {label} has only a single atom !!!\n")
    else:
        freqs = species["freqs"]
        word = "frequency" if len(freqs) == 1 else "frequencies"
        out.write(f"{label} has {len(freqs)} vibrational {word}.\n")
        for i, freq in enumerate(freqs, start=1):
            out.write(f"{i} frequency is {freq} cm-1.\n")
        out.write(f"Vibrational partition function of {label} is {species['vib_val']:.8e}.\n")

    out.write(f"Total partition function of {label} is {species['parti_val']:.8e} m-3.\n\n")


def run(input_file):
    input_file = Path(input_file)
    reader = InputReader(input_file)

    reaction = " ".join(reader.next_data_line())
    temp = float(reader.next_data_line()[0])

    r1 = read_species(reader, temp)
    r2 = read_species(reader, temp)
    ac = read_species(reader, temp)

    e_0 = float(reader.next_data_line()[0])
    parti_fac = ac["parti_val"] / (r1["parti_val"] * r2["parti_val"])
    rate_val = rate(parti_fac, e_0, temp)

    output_file = input_file.with_suffix(".out")

    with output_file.open("w") as out:
        out.write("=" * 112 + "\n")
        out.write("CTSTpy: Python code for Rate Theory, A CTST rate calculator\n")
        out.write("=" * 112 + "\n")
        out.write(f"The reaction is {reaction}\n\n")

        write_species(out, "Reactant-1", r1)
        write_species(out, "Reactant-2", r2)
        write_species(out, "Activated complex", ac)

        out.write("________________________Information about Temperature, Energy barrier and Rate constant_________________________\n")
        out.write(f"Temperature is {temp} K.\n")
        out.write(f"The energy barrier of the reaction is {e_0} kj mol-1.\n")
        out.write(f"The CTST rate constant of the reaction is {rate_val:.8e} m3 s-1.\n")
        out.write("________________________________________________________________________________________________________________\n")
        out.write("| Hurray! A successful calculation! |\n")
        out.write("------------------------------------------------End of the file-------------------------------------------------\n")

    return output_file


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="CTST rate calculator converted from Fortran to Python.")
    parser.add_argument("input_file", nargs="?", default="H_HBr.inp", help="Input .inp file")
    args = parser.parse_args()

    output = run(args.input_file)
    print(f"Wrote {output}")
