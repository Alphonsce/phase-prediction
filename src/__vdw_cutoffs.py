import itertools

vdw_radii = {
    "H": 1.20,
    "C": 1.70,
    "N": 1.55,
    "O": 1.52,
    "F": 1.47,
    "P": 1.80,
    "S": 1.80,
    "Cl": 1.75,
    "Br": 1.85
}

vdw_cutoffs = {
    (a, b): vdw_radii[a] + vdw_radii[b]
    for a, b in itertools.combinations_with_replacement(vdw_radii, 2)
}