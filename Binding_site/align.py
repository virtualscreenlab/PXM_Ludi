from pymol import cmd
import numpy as np

POCKET      = [-8.257, 85.181, 19.050]
PDB_CRYSTAL = r"C:/Users/PC1/Downloads/8C7Y_TXV_protein.pdb"
PDB_MODEL   = r"C:/Users/PC1/Downloads/8c7y_corrected.pdb"
OUT_ITER    = r"C:/Users/PC1/Downloads/8c7y_corrected_aligned.pdb"
CHAIN       = "B"
OFFSET      = 448


def ca_coord(sel):
    return np.array(cmd.get_atom_coords(sel))


def top_cas(sel, pocket, n=5):
    cas = cmd.get_model(sel + " and name CA").atom
    ref = np.array(pocket, dtype=float)
    rows = [(a.resn, a.resi, a.chain,
             float(np.linalg.norm(np.array(a.coord) - ref))) for a in cas]
    rows.sort(key=lambda t: t[3])
    return rows[:n]


def show_fingerprint(title, rows):
    print("")
    print("=== " + title + " : top-5 CAs to pocket point ===")
    print("  chain  resn  resi   dist(A)")
    for resn, resi, chain, d in rows:
        print("  %3s  %-4s  %s   %6.2f" % (chain, resn, resi, d))


def ca_shift_table(pairs):
    print("")
    print("=== CA shift: crystal vs mdl (same residue after alignment) ===")
    print("  resn   crystal#   mdl#    |dCA| (A)")
    dsq = 0.0
    n = 0
    for resi_m, resn in pairs:
        resi_c = str(int(resi_m) + OFFSET)
        resi_m = str(resi_m)
        try:
            p_c = ca_coord("crystal and chain %s and resi %s and name CA" % (CHAIN, resi_c))
            p_m = ca_coord("mdl and chain %s and resi %s and name CA" % (CHAIN, resi_m))
            d = float(np.linalg.norm(p_c - p_m))
            print("  %-4s     %4s      %4s    %6.2f" % (resn, resi_c, resi_m, d))
            dsq += d * d
            n += 1
        except Exception as e:
            print("  %-4s     %4s      %4s    ERROR: %s" % (resn, resi_c, resi_m, e))
    if n > 0:
        rms = np.sqrt(dsq / n)
        print("")
        print("  pocket CA RMSD (crystal vs mdl) = %.2f A  (n=%d)" % (rms, n))


# ---- fresh state ----
cmd.delete("all")
cmd.reinitialize()

cmd.load(PDB_CRYSTAL, "crystal")
cmd.load(PDB_MODEL,   "mdl")
cmd.pseudoatom("ref", pos=POCKET)

# ---- iterative ----
cmd.align("mdl", "crystal")
rows_iter = top_cas("mdl", POCKET, 5)
show_fingerprint("ITERATIVE", rows_iter)
cmd.save(OUT_ITER, "mdl")
print("")
print("saved -> " + OUT_ITER)

pairs = [(resi, resn) for resn, resi, ch, d in rows_iter]
ca_shift_table(pairs)

# ---- non-iterative ----
cmd.align("mdl", "crystal", cycles=0, cutoff=100.0)
rows_non = top_cas("mdl", POCKET, 5)
show_fingerprint("NON-ITERATIVE", rows_non)

pairs_n = [(resi, resn) for resn, resi, ch, d in rows_non]
ca_shift_table(pairs_n)

print("")
print("DONE.")
