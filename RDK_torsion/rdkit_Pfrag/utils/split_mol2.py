import sys
"""
Used to separate mol2 mols in DOCK's output .mol2.gz
"""
# Read mol2 molecule once per time
def next_mol2_lines(infile):
    """Method to return one mol2 block once."""
    lines = list()

    for line in open(infile):
        if "@<TRIPOS>MOLECULE" in line:
            if len(lines) == 0:
                lines.append(line)
            else: # in case there are multiple mol2blocks in infile
                yield lines
                lines = list()
                lines.append(line)
        else:
            lines.append(line)

    yield lines

if __name__ == "__main__":
    mol2file = sys.argv[1]
    # mol2file = "/pubhome/qcxia02/Downloads/dataset/ShuoGu/AmpC.mol2"
    # 1. for mol2gz
    multi_mol2lines = next_mol2_lines(mol2file)
    multi_mol2lines = list(multi_mol2lines)[1:]
    for mol2lines in multi_mol2lines:
        end = -1
        for i,mol2line in enumerate(mol2lines):
            if mol2line.startswith("##########"):
                end = i
                break
        mol2name = mol2lines[1].split()[0]
        mol2lines = mol2lines[:end]
        # with open(mol2name + ".mol2",'w') as f:
            # f.write(''.join(mol2lines))
        