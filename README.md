
# Conventional Transition-State Theory (CTST) in Python
<img width="1829" height="949" alt="Untitled 001" src="https://github.com/user-attachments/assets/d24a0710-5160-439e-97ed-6ef32e5d0327" />

This project is a part of the final assignment of a Python based coursework offered by Weizmann Institute [Basic programming skills (Python)-Summer, 2026].
## 1. Motivation
Conventional Transition State Theory (CTST) provides a theoretical framework for calculating chemical reaction rate constants by assuming a quasi-equilibrium between reactants and the transition state. CTST relates the reaction rate to the statistical properties of molecules and the activated complex located at the saddle point of the potential energy surface. Implementing CTST in Python enables automated and reproducible calculations of reaction kinetics from quantum chemical data. Such computational tools facilitate the investigation of reaction mechanisms and the prediction of temperature-dependent rate constants for complex chemical systems.

## 2. Theory
The theory of reaction rates, published almost simultaneously by Henry Eyring and M. G. Evans, and M. Polanyi in 1935, is now referred to as conventional transition-state theory (CTST). This theory involves the same assumptions and approximations used in the calculations of equilibrium constants using statistical mechanics.\
So after merging the statistical mechanics and chemical equilibrium and to calculate the CTST rate constant, the partition functins of both the reactants and the activated complex are needed. The total partition function of a molecule can be written as:
<p align="center">


$$
q = \sum_i g_i e^{-\frac{ε_i}{k_BT}} \quad...... (1)
$$


</p>

Where the summation is taken over all energy levels. The energy ε<sub>i</sub> is the energy of the i<sup>th</sup> state relative to the zero-point energy, and g<sub>i</sub> is the degeneracy. Usually, it is assume that the various types of energy- electronic, vibrational, rotational, and translational- are independent of one another. So, the total energy corresponding to the i<sup>th</sup> energy state is thus expressed as the sum of the different types of energies:
<p align="center">


$$
ε_i = e_i + ν_i + r_i + t_i  \quad...... (2)
$$


</p>

Where, e<sub>i</sub> = electronic energy, ν<sub>i</sub> = vibrational energy, r<sub>i</sub> = rotational energy, and t<sub>i</sub> = translational energy. So, using equation (2), equation (1) can be written as:
<p align="center">


$$
q = \sum_i g_{e_i} e^{-\frac{e_i}{k_BT}}  g_{ν_i} e^{-\frac{ν_i}{k_BT}}  g_{r_i} e^{-\frac{r_i}{k_BT}}  g_{t_i} e^{-\frac{t_i}{k_BT}}  \quad ...... (3)
$$


</p>

Omitting the degeneracy, equation (3) may be written as:
<p align="center">


$$
q = q_{e} q_{ν} q_{r} q_{t}  \quad...... (4)
$$


</p>

Where, q<sub>e</sub> (electronic), q<sub>ν</sub> (vibrational), q<sub>r</sub> (rotational), and q<sub>t</sub> (translational) are separate partition functions, each referring to one type of energy. At ordinary temperatures, the excited electronic levels of an atom or molecule are
usually too high to make a significant contribution to the partition function, and equation (4) becomes,
<p align="center">


$$
q = q_{ν} q_{r} q_{t}  \quad...... (5)
$$


</p>

The formula of separate partition functions are:
<p align="center">


$$
q_{ν} (per\ normal\ mode) =  \frac{1}{(1 - e^{-\frac{hν}{k_BT}})}  \quad...... (6)
$$

$$
q_{r} (linear\ molecule) =  \frac{8π^2Ik_BT}{σh^2}  \quad...... (7)
$$

$$
q_{r} (nonlinear\ molecule) =  \frac{8π^2(8π^3I_AI_BI_C)^{1/2}(k_BT)^{3/2}}{σh^3}  \quad...... (8)
$$

$$
q_{t} =  \frac{(2πmk_BT)^{3/2}}{h^2}  \quad...... (9)
$$
</p>

Where, h = Planck constant, ν = normal-mode vibrational frequency, k<sub>B</sub> = Boltzmann constant, T = absolute temperature, I = moment of inertia for the linear molecule, (I<sub>A</sub>,I<sub>B</sub>, and I<sub>C</sub>) = moment of inertia for a nonlinear molecule about three axes at right angles, σ = symmetry number, m = mass of molecule. For solving the moment of inertia, the moment of inertia tensor [3x3] is generated first using the equations,
<p align="center">


$$
I_{xx} = \sum_im_i(y_i^2 + z_i^2)  \quad......(11)
$$


$$
I_{yy} = \sum_im_i(z_i^2 + x_i^2)   \quad......(12)
$$


$$
I_{zz} = \sum_im_i(x_i^2 + y_i^2)  \quad......(13)
$$


$$
I_{xy} = -\sum_im_ix_iy_i = I_{yx}  \quad......(14)
$$


$$
I_{xz} = -\sum_im_ix_iz_i = I_{zx}  \quad......(15)
$$


$$
I_{yz} = -\sum_im_iy_iz_i = I_{zy}  \quad......(16)
$$


</p>

Finally, using equations (11)-(16), the moment of inertia tensor is generated.
  | I<sub>xx</sub> | I<sub>xy</sub> | I<sub>xz</sub> |
  |----------------|----------------|----------------|
  | I<sub>yx</sub> | I<sub>yy</sub> | I<sub>yz</sub> |
  | I<sub>zx</sub> | I<sub>zy</sub> | I<sub>zz</sub> |

And to obtain the values of moment of inertia(s), the tensor matrix needed to be diagonalized.

Now, if a reaction is
<p align="center">


$$
A + B => a^{‡} => C + D 
$$


</p> 
then using all the formulas, finally the CSTS rate constant equation becomes,
<p align="center">


$$
k = \frac{k_BT}{h} \frac{q_{‡}}{q_A q_B} e^{-\frac{E_0}{RT}}  \quad...... (10)
$$


</p>

Where, q<sub>‡</sub> is the total partition function of activated complex, q<sub>A</sub> and q<sub>B</sub> are the total partition function of reactant A and B respectively, E<sub>0</sub> = Barrier hight, and R = Gas constant.

## 3. File(s) management
In this project, the CTST rate constant is calculated by solving equation (10) using Python programming. Three Python files are generated: **requirement_scratch.py**, **partition.py**, and **rate_calc.py**. 
- The file **requirement_scratch.py** contains the required molecular data and calculations needed to solve the partition functions, such as total mass and moments of inertia. To calculate the moment of inertia, a [3×3] moment of inertia tensor is first constructed using equations (11)–(16). The tensor matrix is then diagonalized using Python scientific libraries such as NumPy and SciPy. The resulting principal moments of inertia are used according to the molecular type, such as linear, spherical top, prolate symmetric top, oblate symmetric top, or asymmetric top molecules.

- The second file, **partition.py**, contains the calculation of translational, rotational, vibrational, and electronic partition functions using the values obtained from requirement_scratch.py. This file also includes the calculation of the CTST rate constant.

- Finally, **rate_calc.py** acts as the main Python file. It reads all required information from the input file, ```H_HBr.inp```, connects the functions from requirement_scratch.py and partition.py, and generates the final output file, ```H_HBr.out```.
 
## 4. Requirements
(i) NumPy\
(ii) Scipy

- a requirement.txt file will be added with the other files in directory.

## 5. Necessary installation
- ```pip install -r requirements.txt```

## 6. How to run
- ```python rate_calc.py H_HBr.py``` 

After running this command, the CTST rate constant and related partition function results are generated in the output file.


**N.B.**
- Here the reaction H + HBr = H<sub>2</sub> + Br is used as example for this code from the book Chemical Kinetics by K. J. Laidler.
- For more detailed information about CTST and example reaction, please check: Laidler, K. J., & Keith, J. **(1965)**. Chemical kinetics (Vol. 2). New York: McGraw-Hill.
- For more information on moment of inertia tensor, please visit: Bernath, P. F. **(2020)**. Spectra of atoms and molecules. Oxford university press.
