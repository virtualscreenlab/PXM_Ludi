from pymol import cmd
import numpy as np

POCKET      = [-8.257, 85.181, 19.050]
PDB_CRYSTAL = r"C:/Users/PC1/Downloads/8C7Y_TXV_protein.pdb"
PDB_MODEL   = r"C:/Users/PC1/Downloads/8c7y_corrected.pdb"
OUT_ITER    = r"C:/Users/PC1/Downloads/8c7y_corrected_aligned.pdb"


def top_cas(sel, pocket, n=5):
    cas = cmd.get_model(sel + " and name CA").atom
    ref = np.array(pocket, dtype=float)
    out = []
    for a in cas:
        d = float(np.linalg.norm(np.array(a.coord) - ref))
        out.append((a.resn, a.resi, a.chain, d))
    out.sort(key=lambda t: t[3])
    return out[:n]


def show(title, rows):
    print("")
    print("=== " + title + " ===")
    for resn, resi, chain, d in rows:
        print("  %3s %-4s %s  %6.2f A" % (chain, resn, resi, d))


cmd.delete("all")
cmd.reinitialize()
cmd.load(PDB_CRYSTAL, "crystal")
cmd.load(PDB_MODEL,   "mdl")
cmd.pseudoatom("ref", pos=POCKET)

cmd.align("mdl", "crystal")
show("ITERATIVE", top_cas("mdl", POCKET, 5))
cmd.save(OUT_ITER, "mdl")
print("saved -> " + OUT_ITER)

cmd.align("mdl", "crystal", cycles=0, cutoff=100.0)
show("NON-ITERATIVE", top_cas("mdl", POCKET, 5))

print("DONE.")