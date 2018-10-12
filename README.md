# LDOS_localized_QWZ

## Introduction

This program was a part of my summer internship (2018) at the Hungarian Academy of Sciences funded by the (Hungarian) National Research, Development and Innovation Office. The code uses the [Kernel Polynomial Method] (KPM) on the toy model introduced by Qi, Wu and Zhang (QWZ model), equation (6.6) in [this book] or the original paper can be found [here]. We investigated the Anderson localization on the lattice with a random onsite disturbance for each lattice site drawn from a box distribution. A statistical ensemble is required of the Local Density of States (LDOS) to deduce localization, e.g. a possible [numerical approach]. 

## General usage

You should run the program from the command line. The onsite potential for this model should be either u=1.0 or u=3.0, as these are the two qualitatively different cases, for the topologically non-trivial and trivial cases respectively. I calculated some ideal shrinking factors for both cases for the disorder weight in the range $[0.0,13.0]$ for the program to run. There are a few input parameters, which are:

+ -s, sets random seed,
+ -L, linear lattice size, the program creates an L by L square shaped two-dimensional lattice,
+ -u sets the onsite potential,
+ -d sets the disroder strength,
+ -m sets the number of moments to be calculated with the KPM.

So for example you could run the code with: python3 qwz_kpm_ldos.py -s 12345 -L 200 -u 3.0 -d 3.0 > data_rseed12345.

Note: due to the (over) shrinking at the edges we get LDOS values out of our actual energy range.

[Kernel Polynomial Method]: https://arxiv.org/abs/cond-mat/0504627
[this book]: https://arxiv.org/abs/1509.02295
[here]: https://journals.aps.org/prb/abstract/10.1103/PhysRevB.74.085308
[numerical approach]: https://journals.aps.org/prb/abstract/10.1103/PhysRevB.81.155106
