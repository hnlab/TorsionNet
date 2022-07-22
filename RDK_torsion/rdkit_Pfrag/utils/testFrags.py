# %%
"""
The OEGetFuncGroupFragments function uses the following heuristics to fragment a molecule:

1. rings are left intact i.e. considered as a unit
2. exo-ring double bonds are not cleaved from adjacent rings
3. hetero atoms next to each other are kept together
4. sp and sp2 atoms next to each other are kept together
5. in order to avoid generating one atom fragments isolated atoms are attached to the smallest neighbor fragment

    Functional group inclusion criteria:
    - <= 5 heavy atoms
    - must contain at least one hetero atom
    - non-ring
"""

RDK_NCOPS_group_SMARTS_NOS_simplified = [
    "c=O", # for enolate, added by qcxia
    "[CX3]=[OX1]",   # fr_C=O
    "[C!$(C=O)]-[OH]",  # fr_Al_OH
    "c[OH1]", # fr_Ar_OH
    "[OX2](-[#6])-[CH3]",  # fr_methoxy, keep because of capping CH3
    "[CX3]=[NX2]-[OX2]",  # fr_oxime C=N-O-
    "[CX3](=O)[OX2H0]",  # C(=O)O added by qcxia
    # "[#6][CX3](=O)[OX2H0][#6]",  # fr_ester CC=OOC
    "C(=O)[O;H1,-]", # COO, COO-, COOH, added by qcxia
    # "C-C(=O)[O;H1,-]",  # fr_Al_COO C-COOH, C-COO-
    # "c-C(=O)[O;H1,-]",  # fr_Ar_COO c-COOH, c-COO-
    # "[#6]C(=O)[O;H,-1]",  # fr_COO  C-COOH, C-COO-
    "[CX3](=O)[OX1H0-,OX2H1]",  # fr_COO2 -COO-， 
    "[#6][CX3](=O)[#6]",  # fr_ketone, -C(=O)-C
    # "[OD2]([#6])[#6]", # ether  # fr_ether -O-, include unexpected
    "[CX3H1](=O)[#6]",  # fr_aldehyde H-C(=O)-
    "c-[NX3]", # fr_Ar_NH aniline 
    "[Nv3](=C)",  # C=N, added by qcxia
    # "[Nv3](=C)-[#6]",  # fr_Ar_NH imine
    "[NX1]#[CX2]", # fr_nitrile nitrile
    "[NX3]-[NX3]", # fr_hdrzine -N-N-
    "C=N-[NX3]", # fr_hdrzone C=N-N
    "[N!$(N-O)]=O", # fr_nitroso N=O
    "[N!$(N=O)](-O)", # N-(O) added by qcxia
    # "[N!$(N=O)](-O)-C", # fr_N N-(O)-C
    "[$([NX3](=O)=O),$([NX3+](=O)[O-])]", # fr_nitro, added by qcxia
    # "[$([NX3](=O)=O),$([NX3+](=O)[O-])][!#8]", # fr_nitro 
    "N(=O)(O)", # nitro, added by qcxia
    "N=N", # N=N, added by qcxia
    # "N(=O)(O)[#6]", # #fr_nitro
    # "[#6]-N=N-[#6]", # fr_azo -N=N-
    "[N+]#N", # fr_diazo
    "[$(*-[NX2-]-[NX2+]#[NX1]),$(*-[NX2]=[NX2+]=[NX1-])]", # fr_azide
    "C(=O)-N", # fr_amide
    "C(=O)-[NH2]", # fr_priamide
    "C(=N)(-N)" # C(=N)(-N), added by qcxia
    # "C(=N)(-N)-[!#7]", # fr_amidine
    "C(=N)(N)N", # fr_guanido
    "N(-C(=O))-C=O", # fr_imide
    "N=C=O", # fr_isocyan
    "N=C=S", # fr_isothiocyan
    "S-C#N", # fr_thiocyan
    # "[SX2](-[#6])-C", # fr_sulfide
    "C=[SX1]", # fr_C=S
    "S(=,-[OX1;+0,-1])(=,-[OX1;+0,-1])", # O=S=O, added by qcxia
    # "S(=,-[OX1;+0,-1])(=,-[OX1;+0,-1])(-[#6])-[#6]", #fr_sulfone
    "N-S(=,-[OX1;+0,-1])(=,-[OX1;+0,-1])", #N-S(=O)=O, added by qcxia
    # "N-S(=,-[OX1;+0,-1])(=,-[OX1;+0,-1])-[#6]", #fr_sulfonamd
    "[NH2]-S(=,-[OX1;+0;-1])(=,-[OX1;+0;-1])", # NH2-S(=O)=O, added by qcxia 
    # "[NH2]-S(=,-[OX1;+0;-1])(=,-[OX1;+0;-1])-[#6]", #fr_prisulfonamd
    "[P](=[#8])(-[!#1])(-[!#1])-[!#1]", # phosphate, added by qcxia
    # "[$(P(=[OX1])([$([OX2H]),$([OX1-]),$([OX2]P)])([$([OX2H]),$([OX1-]),$([OX2]P)])[$([OX2H]),$([OX1-]),$([OX2]P)]),$([P+]([OX1-])([$([OX2H]),$([OX1-]),$([OX2]P)])([$([OX2H]),$([OX1-]),$([OX2]P)])[$([OX2H]),$([OX1-]),$([OX2]P)])]", #fr_phos_acid
    # "[$(P(=[OX1])([OX2][#6])([$([OX2H]),$([OX1-]),$([OX2][#6])])[$([OX2H]),$([OX1-]),$([OX2][#6]),$([OX2]P)]),$([P+]([OX1-])([OX2][#6])([$([OX2H]),$([OX1-]),$([OX2][#6])])[$([OX2H]),$([OX1-]),$([OX2][#6]),$([OX2]P)])]", #fr_phos_ester
]




RDK_NCOPS_group_SMARTS_most = [
    "[CX3]=[OX1]",
    "[C!$(C=O)]-[OH]",
    "[OX2](-[#6])-[CH3]",
    "[CX3]=[NX2]-[OX2]",
    "[#6][CX3](=O)[OX2H0][#6]",
    "C-C(=O)[O;H1,-]",
    "c-C(=O)[O;H1,-]",
    "[#6]C(=O)[O;H,-1]",
    "[CX3](=O)[OX1H0-,OX2H1]",
    "[#6][CX3](=O)[#6]",
    # "[OD2]([#6])[#6]", # ether
    # "[OX2H]-c1ccccc1",  # fr_phenol include non-ring
    "[CX3H1](=O)[#6]",
    "[$([NX4+]),$([NX4]=*)]"
    "[NH2,nH2]",
    "[NH1,nH1]",
    "[NH0,nH0]",
    "[Nv3](=C)-[#6]",
    "[NX1]#[CX2]",
    "[NX3]-[NX3]",
    "C=N-[NX3]",
    "[N!$(N-O)]=O",
    "[N!$(N=O)](-O)-C",
    "[$([NX3](=O)=O),$([NX3+](=O)[O-])][!#8]",
    "N(=O)(O)[#6]",
    "[#6]-N=N-[#6]",
    "[N+]#N",
    "[$(*-[NX2-]-[NX2+]#[NX1]),$(*-[NX2]=[NX2+]=[NX1-])]",
    "C(=O)-N",
    "C(=O)-[NH2]",
    "C(=N)(-N)-[!#7]",
    "C(=N)(N)N",
    "[nH]",
    "N(-C(=O))-C=O",
    "N=C=O",
    "N=C=S",
    "S-C#N",
    "[SX2](-[#6])-C",
    "[SH]",
    "C=[SX1]",
    "S(=,-[OX1;+0,-1])(=,-[OX1;+0,-1])(-[#6])-[#6]",
    "S(=,-[OX1;+0,-1])(=,-[OX1;+0,-1])(-[#6])-[#6]",
    "N-S(=,-[OX1;+0,-1])(=,-[OX1;+0,-1])-[#6]",
    "[NH2]-S(=,-[OX1;+0;-1])(=,-[OX1;+0;-1])-[#6]",
    "[$(P(=[OX1])([$([OX2H]),$([OX1-]),$([OX2]P)])([$([OX2H]),$([OX1-]),$([OX2]P)])[$([OX2H]),$([OX1-]),$([OX2]P)]),$([P+]([OX1-])([$([OX2H]),$([OX1-]),$([OX2]P)])([$([OX2H]),$([OX1-]),$([OX2]P)])[$([OX2H]),$([OX1-]),$([OX2]P)])]",
    "[$(P(=[OX1])([OX2][#6])([$([OX2H]),$([OX1-]),$([OX2][#6])])[$([OX2H]),$([OX1-]),$([OX2][#6]),$([OX2]P)]),$([P+]([OX1-])([OX2][#6])([$([OX2H]),$([OX1-]),$([OX2][#6])])[$([OX2H]),$([OX1-]),$([OX2][#6]),$([OX2]P)])]",
]


# Mainly COPS acid, including N-related groups
# Some cannot be read may due to format incompatibility between RDKit and Omega
raw_ChayaSt_fragmentation = [
    "[NX3][NX3]",   # hydrazine
    "[NX3][NX2]",   # hydrazone
    "[N]-[O]",   # nitric oxide
    "[#7][#6](=[#8])", # amide
    "[#7][#6](-[O-])",   # amide
    "[NX3][CX3](=[OX1])[NX3]",   # urea
    "[CX3H1](=O)[#6]",   # aldehyde
    "[#16X3]=[OX1]", # sulfoxide
    "[#16X3+][OX1-]", # sulfoxide
    "[#16X4](=[OX1])=([OX1])",   # sulfonyl
    "[#16X3](=[OX1])[OX2H,OX1H0-]",   # sulfinic acid
    "[#16X4](=[OX1])=([OX1])([NX3R0])",   # sulfinamide
    "[#16X4](=[OX1])(=[OX1])[OX2H,OX1H0-]",   # sulfonic acid
    "[PX4](=[OX1])([#6])([#6])([#6])",   # phosphine oxide
    "P(=[OX1])([OX2H,OX1-])([OX2H,OX1-])",   # phosphonate
    "[PX4](=[OX1])([#8])([#8])([#8])",   # phosphate
    "[CX3](=O)[OX1H0-,OX2H1]",   # carboxylic acid
    "([NX3+](=O)[O-])", # nitro
    "([NX3](=O)=O)",   # nitro
    "[CX3](=O)[OX2H0]",   # ester
    "[#6]((([F,Cl,I,Br])[F,Cl,I,Br])[F,Cl,I,Br])",  # tri-halides
]


ChayaSt_fragmentation = [
    "[NX3]-[NX3]",   # hydrazine
    "[NX3]-[NX2]",   # hydrazone
    "[N]-[O]",   # nitric oxide
    "[#7]-[#6](=[#8])", # amide
    "[#7]=[#6](-[O-])",   # amide
    "[NX3]-[CX3](=[OX1])-[NX3]",   # urea
    "[CX3H1](=O)-[#6]",   # aldehyde
    "[#16X3]=[OX1]", # sulfoxide
    "[#16X3+]-[OX1-]", # sulfoxide
    "[#16X4](=[OX1])=[OX1]",   # sulfonyl
    "[#16X3](=[OX1])-[OX2H]",   # sulfinic acid
    "[#16X3](=[OX1])-[OX1H0-]",   # sulfinic acid
    "[NX3R0]-[#16X4](=[OX1])=[OX1]",   # sulfinamide
    "[#16X4](=[OX1])(=[OX1])-[OX2H]",   # sulfonic acid
    "[#16X4](=[OX1])(=[OX1])-[OX1H0-]",   # sulfonic acid
    "[PX4](=[OX1])(-[#6])(-[#6])(-[#6])",   # phosphine oxide
    "P(=[OX1])(-[OX2H])(-[OX2H])",   # phosphonate
    "P(=[OX1])(-[OX2H])(-[OX1-])",   # phosphonate
    "P(=[OX1])(-[OX1-])(-[OX2H])",   # phosphonate
    "P(=[OX1])(-[OX1-])(-[OX1-])",   # phosphonate
    "[PX4](=[OX1])(-[#8])(-[#8])(-[#8])",   # phosphate
    "[CX3](=O)-[OX1H0-]",   # carboxylic acid
    "[CX3](=O)-[OX2H1]",   # carboxylic acid
    "[NX3+](=O)-[O-]", # nitro
    "[NX3](=O)=O",   # nitro
    "[CX3](=O)-[OX2H0]",   # ester
    "[CX3]=[OX1]",   # fr_C=O, added by qcxia at 2022/07/17 to solve aromaticity problem
    "c=O",   # Ar c=O, also added by qcxia at 2022/07/17 to solve aromaticity problem

] + [ f"[#6](-[{i}])(-[{j}])-[{k}]" for i in ["F", "Cl", "Br", "I"] for j in ["F", "Cl", "Br", "I"] for k in ["F", "Cl", "Br", "I"] ] #  tri-halides 
