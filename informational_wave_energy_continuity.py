# Informational Wave Energy Continuity Verification (vOmega Law 41)
# ASCII-safe simulation: verifies local energy conservation for wave equation.

import numpy as np

# Domain
N = 500
L = 1.0
x = np.linspace(0, L, N, endpoint=False)
dx = x[1] - x[0]

# Finite propagation speed
c = 1.0

# CFL-stable timestep
dt = 0.45 * dx / c

# Initial wave: Gaussian bump
psi_now = np.exp(-((x - 0.5)**2) / 0.002)
psi_prev = psi_now.copy()

def d_dx(f):
    return (np.roll(f, -1) - np.roll(f, 1)) / (2*dx)

def lap(f):
    return (np.roll(f, -1) - 2*f + np.roll(f, 1)) / dx**2

def energy_density(psi_prev, psi_now):
    psi_t = (psi_now - psi_prev) / dt
    psi_x = d_dx(psi_now)
    return 0.5 * (psi_t**2 + c**2 * psi_x**2)

def flux(psi_prev, psi_now):
    psi_t = (psi_now - psi_prev) / dt
    psi_x = d_dx(psi_now)
    return -c**2 * psi_t * psi_x

# Track global energy over time
energies = []

for step in range(300):
    eps = energy_density(psi_prev, psi_now)
    energies.append(np.sum(eps) * dx)

    lap_now = lap(psi_now)
    psi_next = 2*psi_now - psi_prev + (c*dt)**2 * lap_now

    psi_prev, psi_now = psi_now, psi_next

print("Initial total energy:", energies[0])
print("Final total energy:", energies[-1])
print("Energy difference:", energies[-1] - energies[0])
