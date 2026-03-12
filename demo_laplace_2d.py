import matplotlib.pyplot as plt
import numpy as np
import scipy.sparse as sp


def solve_laplace_equation(mask, peak_coords):
    """
    Solve Laplace equation to generate smooth reference field.

    Parameters:
        mask (2D array): Boolean matrix, True indicates interior region, False indicates boundary or exterior.
        peak_coords (tuple): Peak point coordinates in (y, x) format.

    Returns:
        u (2D array): Solved field with values in range [0, 1].
    """
    H, W = mask.shape
    N = H * W

    # 1. Prepare sparse matrix construction data (COO format: row, col, value)
    rows = []
    cols = []
    data = []
    rhs = np.zeros(N)  # Right-hand side vector b of equation

    peak_idx = peak_coords[0] * W + peak_coords[1]

    # 2. Iterate over grid to build equations
    # Note: For large grids, Python loops are slower, but at 100x100 scale usually only takes hundreds of milliseconds
    for r in range(H):
        for c in range(W):
            k = r * W + c

            # --- Case A: Peak point (enforce constraint u = 1) ---
            if k == peak_idx:
                rows.append(k)
                cols.append(k)
                data.append(1.0)
                rhs[k] = 1.0
                continue

            # --- Case B: Boundary or exterior point (enforce constraint u = 0) ---
            if not mask[r, c]:
                rows.append(k)
                cols.append(k)
                data.append(1.0)
                rhs[k] = 0.0
                continue

            # --- Case C: Interior point (discrete Laplace equation) ---
            # Equation: 4*u_center - u_up - u_down - u_left - u_right = 0

            # Main diagonal term (coefficient 4)
            rows.append(k)
            cols.append(k)
            data.append(4.0)

            # Neighbor terms (coefficient -1)
            # Check four directions to ensure no out-of-bounds
            # (Note: Even if a neighbor is a boundary point, it exists as a variable in the equation, just with value 0)
            neighbors = []
            if r > 0:
                neighbors.append((r - 1) * W + c)  # Up
            if r < H - 1:
                neighbors.append((r + 1) * W + c)  # Down
            if c > 0:
                neighbors.append(r * W + (c - 1))  # Left
            if c < W - 1:
                neighbors.append(r * W + (c + 1))  # Right

            for n_idx in neighbors:
                rows.append(k)
                cols.append(n_idx)
                data.append(-1.0)

    # 3. Construct sparse matrix A
    A = sp.csr_matrix((data, (rows, cols)), shape=(N, N))

    # 4. Solve linear system Au = b
    print(f"Solving linear system with {N} variables...")
    u_flat = sp.linalg.spsolve(A, rhs)

    return u_flat.reshape(H, W)


# --- Example: Solve within a circular region ---

# 1. Define grid and shape
size = 100
Y, X = np.ogrid[:size, :size]
center = (size // 2, size // 2)
radius = size / 2
# Create circular mask (interior is True)
mask = (X - center[1]) ** 2 + (Y - center[0]) ** 2 <= radius**2

# 2. Set peak point (slightly offset from center to observe asymmetric diffusion)
peak_point = (center[0] - 20, center[1] + 10)

# 3. Run solver
u_field = solve_laplace_equation(mask, peak_point)

# rescale the values to higher contrast
# u_field = (u_field - u_field.min()) / (u_field.max() - u_field.min())

# 4. Visualize results
plt.figure(figsize=(12, 5))

plt.subplot(1, 2, 1)
plt.imshow(u_field, origin="lower", cmap="winter")
plt.colorbar(label="Potential Value u(x)")
plt.scatter([peak_point[1]], [peak_point[0]], c="cyan", marker="x", label="Peak")
plt.title("Solution of Laplace Equation")
plt.legend()

plt.subplot(1, 2, 2)
# Plot horizontal cross-section passing through the peak point
plt.plot(u_field[peak_point[0], :])
plt.title(f"Cross Section at y={peak_point[0]}")
plt.xlabel("x")
plt.ylabel("u(x)")
plt.grid(True, linestyle="--", alpha=0.6)

plt.tight_layout()
plt.show()
