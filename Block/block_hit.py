import numpy as np
import matplotlib.pyplot as plt

# Load protein RMSD
data = np.loadtxt('hit.dat', skiprows=1)
rmsd = data[:, 1]

# Analysis parameters
cutoff = 0       # Excluded frames
dt = 0.3           # Sampling interval (ns)

# Production trajectory
x = rmsd[cutoff:]
time = np.arange(len(rmsd)) * dt

# Block averaging
block_sizes = []
errors = []

print("\nPROTEIN RMSD BLOCK ANALYSIS")
print("=" * 55)
print(f"Total frames:      {len(rmsd)}")
print(f"Excluded frames:   {cutoff}")
print(f"Production frames: {len(x)}")
print(f"Mean RMSD:         {np.mean(x):.4f} Å")
print(f"Standard deviation: {np.std(x, ddof=1):.4f} Å")
print("-" * 55)
print("Blocks\tBlock (ns)\tBSE (Å)")

for n in [128, 64, 32, 16, 8, 4, 2]:

    size = len(x) // n

    if size < 1:
        continue

    blocks = x[:n * size].reshape(n, size).mean(axis=1)

    error = np.std(blocks, ddof=1) / np.sqrt(n)

    block_sizes.append(size * dt)
    errors.append(error)

    print(f"{n:6d}\t{size * dt:10.2f}\t{error:.4f}")

# Cumulative mean
cumulative = np.cumsum(x) / np.arange(1, len(x) + 1)

# Plot results
fig, ax = plt.subplots(1, 2, figsize=(10, 4))

# A. Block averaging
ax[0].plot(block_sizes, errors, 'bo-')
ax[0].set_xscale('log')
ax[0].set_xlabel('Block duration (ns)')
ax[0].set_ylabel('Block standard error (Å)')
ax[0].set_title('(A) Block averaging')

# B. Cumulative RMSD
ax[1].plot(time[cutoff:], cumulative, 'b-')
ax[1].axhline(np.mean(x), color='red', linestyle='--')
ax[1].set_xlabel('Simulation time (ns)')
ax[1].set_ylabel('Cumulative mean RMSD (Å)')
ax[1].set_title('(B) RMSD convergence')

plt.tight_layout()
plt.savefig('block_analysis.png', dpi=600)
plt.savefig('block_analysis.pdf')
plt.show()

umulative = np.cumsum(x) / np.arange(1, len(x) + 1)

np.savetxt(
    'cumulative_rmsd.dat',
    np.column_stack((time[cutoff:], cumulative)),
    fmt='%.4f',
    header='Time(ns) Cumulative_RMSD(Angstrom)',
    comments='# '
)
