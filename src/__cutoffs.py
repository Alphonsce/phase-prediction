# cutoffs = {
#     ("O", "H"): 1.27,
#     ("O", "C"): 1.72,
#     ("C", "H"): 1.37,
#     ("C", "N"): 1.60,
#     ("C", "C"): 1.70,
#     ("O", "O"): 1.50,
#     # ("C", "S"): 1.90,
#     # ("O", "O"): 1.50
# }

cutoffs_gpt = {
    ("C", "H"): 1.20,   # C–H bond length ≈ 1.09 Å + buffer
    ("C", "C"): 1.70,   # C–C single bond ≈ 1.54 Å + buffer
    ("C", "N"): 1.60,   # C–N single bond ≈ 1.47 Å + buffer
    ("C", "O"): 1.60,   # C–O single bond ≈ 1.43 Å + buffer
    ("C", "S"): 1.90,   # C–S single bond ≈ 1.82 Å + buffer
    ("C", "Cl"): 1.90,  # C–Cl bond length ≈ 1.76 Å + buffer
    ("C", "Br"): 2.00,  # C–Br bond length ≈ 1.94 Å + buffer
    ("C", "I"): 2.20,   # C–I bond length ≈ 2.14 Å + buffer
    ("N", "H"): 1.10,   # N–H bond length ≈ 1.01 Å + buffer
    ("N", "O"): 1.50,   # N–O single bond ≈ 1.45 Å + buffer
    ("N", "N"): 1.50,   # N–N single bond ≈ 1.45 Å + buffer
    ("O", "H"): 1.10,   # O–H bond length ≈ 0.97 Å + buffer
    ("O", "O"): 1.50,   # O–O single bond ≈ 1.45 Å + buffer
    ("S", "H"): 1.40,   # S–H bond length ≈ 1.34 Å + buffer
    ("S", "O"): 1.60,   # S–O single bond ≈ 1.48 Å + buffer
    ("Cl", "H"): 1.30,  # H–Cl bond length ≈ 1.27 Å + buffer
    ("Br", "H"): 1.50,  # H–Br bond length ≈ 1.41 Å + buffer
    ("I", "H"): 1.70,   # H–I bond length ≈ 1.61 Å + buffer
    ("H",  "H"): 0.84,   # H–H bond ≈0.74 Å +0.10 Å buffer :contentReference[oaicite:0]{index=0}
    ("S",  "S"): 2.15,   # S–S bond ≈2.03 Å +0.12 Å buffer :contentReference[oaicite:1]{index=1}
    ("Cl", "Cl"): 2.09,  # Cl–Cl bond ≈1.99 Å +0.10 Å buffer :contentReference[oaicite:2]{index=2}
    ("Br", "Br"): 2.38,  # Br–Br bond ≈2.28 Å +0.10 Å buffer :contentReference[oaicite:3]{index=3}
    ("I",  "I"): 2.76,    # I–I bond ≈2.66 Å +0.10 Å buffer :contentReference[oaicite:4]{index=4}

    ### Update: ###

    # Homoatomic bonds
    ("Ag", "Ag"): 2*1.45 + 0.10,   # 3.00 Å; r_cov(Ag)=1.45 Å :contentReference[oaicite:1]{index=1}
    ("Al", "Al"): 2*1.21 + 0.10,   # 2.52 Å; r_cov(Al)=1.21 Å :contentReference[oaicite:2]{index=2}
    ("As", "As"): 2*1.19 + 0.10,   # 2.48 Å; r_cov(As)=1.19 Å :contentReference[oaicite:3]{index=3}
    ("Au", "Au"): 2*1.36 + 0.10,   # 2.82 Å; r_cov(Au)=1.36 Å :contentReference[oaicite:4]{index=4}
    ("B",  "B"):  2*0.84 + 0.10,   # 1.78 Å; r_cov(B)=0.84 Å :contentReference[oaicite:5]{index=5}
    ("Ba", "Ba"): 2*2.15 + 0.10,   # 4.40 Å; r_cov(Ba)=2.15 Å :contentReference[oaicite:6]{index=6}
    ("Bi", "Bi"): 2*1.48 + 0.10,   # 3.06 Å; r_cov(Bi)=1.48 Å :contentReference[oaicite:7]{index=7}
    ("Br", "Br"): 2*1.20 + 0.10,   # 2.50 Å; r_cov(Br)=1.20 Å :contentReference[oaicite:8]{index=8}
    ("C",  "C"):  2*0.76 + 0.10,   # 1.62 Å; r_cov(C)=0.76 Å :contentReference[oaicite:9]{index=9}
    ("Ca", "Ca"): 2*1.76 + 0.10,   # 3.62 Å; r_cov(Ca)=1.76 Å :contentReference[oaicite:10]{index=10}
    ("Cl", "Cl"): 2*1.02 + 0.10,   # 2.14 Å; r_cov(Cl)=1.02 Å :contentReference[oaicite:11]{index=11}
    ("Co", "Co"): 2*1.26 + 0.10,   # 2.62 Å; r_cov(Co)=1.26 Å :contentReference[oaicite:12]{index=12}
    ("Cr", "Cr"): 2*1.39 + 0.10,   # 2.88 Å; r_cov(Cr)=1.39 Å :contentReference[oaicite:13]{index=13}
    ("Cs", "Cs"): 2*2.44 + 0.10,   # 4.98 Å; r_cov(Cs)=2.44 Å :contentReference[oaicite:14]{index=14}
    ("Cu", "Cu"): 2*1.39 + 0.10,   # 2.88 Å; r_cov(Cu)=1.39 Å :contentReference[oaicite:15]{index=15}
    ("F",  "F"):  2*0.57 + 0.10,   # 1.24 Å; r_cov(F)=0.57 Å :contentReference[oaicite:16]{index=16}
    ("Fe", "Fe"): 2*1.32 + 0.10,   # 2.74 Å; r_cov(Fe)=1.32 Å :contentReference[oaicite:17]{index=17}
    ("Ga", "Ga"): 2*1.22 + 0.10,   # 2.54 Å; r_cov(Ga)=1.22 Å :contentReference[oaicite:18]{index=18}
    ("Ge", "Ge"): 2*1.22 + 0.10,   # 2.54 Å; r_cov(Ge)=1.22 Å :contentReference[oaicite:19]{index=19}
    ("H",  "H"):  0.84,           # fixed H–H cutoff   
    ("Hg", "Hg"): 2*1.36 + 0.10,   # 2.82 Å; r_cov(Hg)=1.36 Å :contentReference[oaicite:20]{index=20}
    ("I",  "I"):  2*1.39 + 0.10,   # 2.88 Å; r_cov(I)=1.39 Å :contentReference[oaicite:21]{index=21}
    ("In", "In"): 2*1.42 + 0.10,   # 2.94 Å; r_cov(In)=1.42 Å :contentReference[oaicite:22]{index=22}
    ("K",  "K"):  2*2.03 + 0.10,   # 4.16 Å; r_cov(K)=2.03 Å :contentReference[oaicite:23]{index=23}
    ("Li", "Li"): 2*1.28 + 0.10,   # 2.66 Å; r_cov(Li)=1.28 Å :contentReference[oaicite:24]{index=24}
    ("Mg", "Mg"): 2*1.41 + 0.10,   # 2.92 Å; r_cov(Mg)=1.41 Å :contentReference[oaicite:25]{index=25}
    ("Mn", "Mn"): 2*1.39 + 0.10,   # 2.88 Å; r_cov(Mn)=1.39 Å :contentReference[oaicite:26]{index=26}
    ("Na", "Na"): 2*1.66 + 0.10,   # 3.42 Å; r_cov(Na)=1.66 Å :contentReference[oaicite:27]{index=27}
    ("Ni", "Ni"): 2*1.24 + 0.10,   # 2.58 Å; r_cov(Ni)=1.24 Å :contentReference[oaicite:28]{index=28}
    ("O",  "O"):  2*0.66 + 0.10,   # 1.42 Å; r_cov(O)=0.66 Å :contentReference[oaicite:29]{index=29}
    ("P",  "P"):  2*1.07 + 0.10,   # 2.24 Å; r_cov(P)=1.07 Å :contentReference[oaicite:30]{index=30}
    ("Pb", "Pb"): 2*1.36 + 0.10,   # 2.82 Å; r_cov(Pb)=1.36 Å :contentReference[oaicite:31]{index=31}
    ("Pt", "Pt"): 2*1.36 + 0.10,   # 2.82 Å; r_cov(Pt)=1.36 Å :contentReference[oaicite:32]{index=32}
    ("Rb", "Rb"): 2*2.20 + 0.10,   # 4.50 Å; r_cov(Rb)=2.20 Å :contentReference[oaicite:33]{index=33}
    ("S",  "S"):  2*1.05 + 0.10,   # 2.20 Å; r_cov(S)=1.05 Å :contentReference[oaicite:34]{index=34}
    ("Sb", "Sb"): 2*1.39 + 0.10,   # 2.88 Å; r_cov(Sb)=1.39 Å :contentReference[oaicite:35]{index=35}
    ("Sc", "Sc"): 2*1.70 + 0.10,   # 3.50 Å; r_cov(Sc)=1.70 Å :contentReference[oaicite:36]{index=36}
    ("Se", "Se"): 2*1.20 + 0.10,   # 2.50 Å; r_cov(Se)=1.20 Å :contentReference[oaicite:37]{index=37}
    ("Si", "Si"): 2*1.11 + 0.10,   # 2.32 Å; r_cov(Si)=1.11 Å :contentReference[oaicite:38]{index=38}
    ("Sn", "Sn"): 2*1.39 + 0.10,   # 2.88 Å; r_cov(Sn)=1.39 Å :contentReference[oaicite:39]{index=39}
    ("Sr", "Sr"): 2*1.95 + 0.10,   # 4.00 Å; r_cov(Sr)=1.95 Å :contentReference[oaicite:40]{index=40}
    ("Te", "Te"): 2*1.38 + 0.10,   # 2.86 Å; r_cov(Te)=1.38 Å :contentReference[oaicite:41]{index=41}
    ("Ti", "Ti"): 2*1.60 + 0.10,   # 3.30 Å; r_cov(Ti)=1.60 Å :contentReference[oaicite:42]{index=42}
    ("Tl", "Tl"): 2*1.45 + 0.10,   # 3.00 Å; r_cov(Tl)=1.45 Å :contentReference[oaicite:43]{index=43}
    ("Zn", "Zn"): 2*1.39 + 0.10,   # 2.88 Å; r_cov(Zn)=1.39 Å :contentReference[oaicite:44]{index=44}

    # Heteroatomic bonds (element–H)
    # for each X: (X, "H"): r_cov(X)+r_cov(H)+0.10**
    # e.g. ("Ag","H"): 1.45+0.31+0.10 = 1.86 Å; r_cov(H)=0.31 Å :contentReference[oaicite:45]{index=45}**
    ("Ag","H"): 1.86,  ("Al","H"): 1.62, ("As","H"): 1.60, ("Au","H"): 1.77,
    ("B","H"):  1.25,  ("Ba","H"): 2.56, ("Bi","H"): 1.89, ("Br","H"): 1.61,
    ("C","H"):  1.17,  ("Ca","H"): 2.17, ("Cl","H"): 1.43, ("Co","H"): 1.67,
    ("Cr","H"): 1.80,  ("Cs","H"): 2.85, ("Cu","H"): 1.80, ("F","H"):  0.98,
    ("Fe","H"): 1.73,  ("Ga","H"): 1.63, ("Ge","H"): 1.63, ("Hg","H"): 1.77,
    ("I","H"):  1.80,  ("In","H"): 1.83, ("K","H"):  2.44, ("Li","H"): 1.69,
    ("Mg","H"): 1.82,  ("Mn","H"): 1.80, ("Na","H"): 2.07, ("Ni","H"): 1.65,
    ("O","H"):  1.07,  ("P","H"):  1.48, ("Pb","H"): 1.77, ("Pt","H"): 1.77,
    ("Rb","H"): 2.61,  ("S","H"):  1.46, ("Sb","H"): 1.80, ("Sc","H"): 2.11,
    ("Se","H"): 1.61,  ("Si","H"): 1.52, ("Sn","H"): 1.80, ("Sr","H"): 2.36,
    ("Te","H"): 1.79,  ("Ti","H"): 2.01, ("Tl","H"): 1.86, ("Zn","H"): 1.80
}

cutoffs_grok = {
    ("C", "H"): 1.20,
    ("C", "C"): 1.70,
    ("C", "N"): 1.60,
    ("C", "O"): 1.60,
    ("C", "S"): 1.90,
    ("C", "Cl"): 1.90,
    ("C", "Br"): 2.00,
    ("C", "I"): 2.20,
    ("N", "H"): 1.10,
    ("N", "O"): 1.50,
    ("N", "N"): 1.50,
    ("O", "H"): 1.10,
    ("O", "O"): 1.50,
    ("S", "H"): 1.40,
    ("S", "O"): 1.60,
    ("Cl", "H"): 1.30,
    ("Br", "H"): 1.50,
    ("I", "H"): 1.70,
    ("H", "H"): 0.84,
    ("S", "S"): 2.15,
    ("Cl", "Cl"): 2.09,
    ("Br", "Br"): 2.38,
    ("I", "I"): 2.76,
    ("Ag", "Ag"): 2.68,
    ("Al", "Al"): 2.64,
    ("As", "As"): 2.54,
    ("Au", "Au"): 2.84,
    ("B", "B"): 1.82,
    ("Ba", "Ba"): 4.42,
    ("Bi", "Bi"): 3.08,
    ("Ca", "Ca"): 3.54,
    ("Co", "Co"): 2.34,
    ("Cr", "Cr"): 2.56,
    ("Cs", "Cs"): 5.00,
    ("Cu", "Cu"): 2.36,
    ("F", "F"): 1.40,
    ("Fe", "Fe"): 2.44,
    ("Ga", "Ga"): 2.60,
    ("Ge", "Ge"): 2.54,
    ("Hg", "Hg"): 2.76,
    ("In", "In"): 2.96,
    ("K", "K"): 4.04,
    ("Li", "Li"): 2.78,
    ("Mg", "Mg"): 2.90,
    ("Mn", "Mn"): 2.50,
    ("Na", "Na"): 3.22,
    ("Ni", "Ni"): 2.32,
    ("P", "P"): 2.34,
    ("Pb", "Pb"): 3.04,
    ("Pt", "Pt"): 2.84,
    ("Rb", "Rb"): 4.32,
    ("Sb", "Sb"): 2.92,
    ("Sc", "Sc"): 3.08,
    ("Se", "Se"): 2.44,
    ("Si", "Si"): 2.44,
    ("Sn", "Sn"): 2.92,
    ("Sr", "Sr"): 3.82,
    ("Te", "Te"): 2.84,
    ("Ti", "Ti"): 2.84,
    ("Tl", "Tl"): 3.02,
    ("Zn", "Zn"): 2.48
}