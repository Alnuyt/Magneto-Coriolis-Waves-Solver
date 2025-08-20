#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 18 15:52:07 2025

@author: alexandrenuyt
"""

# === USEFUL PACKAGES ===
import numpy as np
import mpmath as mp
import matplotlib.pyplot as plt
from scipy.linalg import eig
from scipy.special import yv, jv, jvp
from scipy.optimize import root_scalar
from scipy.interpolate import griddata

# === PHYSICAL SITUATION ===
solve_mhd = True # True : MHD, False : fluid (B0 = 0)
solve_Lundquist = False # True : m = 0 

if solve_mhd == True : 
    if solve_Lundquist == True :
        print("-" * 67)
        print("                    Cylindrical-Alfvèn Waves           ")
        print("-" * 67)
    else : 
        print("-" * 67)
        print("                    Magneto-Coriolis Waves           ")
        print("-" * 67)
else : 
    print("-" * 67)
    print("                          Inertial Waves           ")
    print("-" * 67)
    
# === PHYSICAL PARAMETERS ===
h = 1
a = 1
n_axial = 1
Omega_mass = 1
Omega_rot = 0 if solve_Lundquist else Omega_mass
mu0 = 1
rho = 1
m_az = 1
B0 = 1 if solve_mhd else 0
m = 0 if solve_Lundquist else m_az
k = (n_axial * np.pi) / h
C = 1 # pour le champ extérieur
phi_indices = [0,50]
times = [1]

# === Solution Greenspan ===
# Fonction à résoudre
def equation(xi):
    Jm = jv(abs(m), xi)
    dJm_dxi = jvp(abs(m), xi, n=1)
    factor = np.sqrt(1 + (xi**2 * h**2) / (n_axial**2 * np.pi**2 * a**2))
    return xi * dJm_dxi + m * factor * Jm

# Recherche des racines
xi_vals = np.linspace(0.0001, 200, 1000)
f_vals = [equation(xi) for xi in xi_vals]

# Roots of the function 
roots = []
for i in range(len(xi_vals) - 1):
    if f_vals[i] * f_vals[i+1] < 0:
        try:
            sol = root_scalar(equation, bracket=[xi_vals[i], xi_vals[i+1]], method='brentq')
            roots.append(sol.root)
        except:
            pass

# lambda_{knm}
def lambda_knm(xi):
    return 2 / np.sqrt(1 + (xi**2 * h**2) / (n_axial**2 * np.pi**2 * a**2))

roots = roots[:50]  # 50 premières valeurs
lambdas = [lambda_knm(xi) for xi in roots]

# === Solution de Lundquist ===
if solve_Lundquist:
    N_l = 10
    
    # Alfvén speed
    V_a = B0 / np.sqrt(mu0 * rho)

    # Zeros of J1
    j1_zeros = np.array([float(mp.besseljzero(1, n)) for n in range(1, N_l + 1)])

    # Axial wavenumber
    k_z = (np.pi / h) * np.arange(1, N_l + 1, dtype=float)

    # Eigenfrequencies
    omega_L_rad = V_a * j1_zeros / a      # radial Lundquist
    omega_L_ax  = V_a * k_z               # axial Alfvén
else : 
    pass

# === CHEBYSHEV MATRIX ===
def cheb_diff(N):
    if N == 0:
        return np.array([[0]]), np.array([1])
    x = -np.cos(np.pi * np.arange(N + 1) / N) # points sur [-1, 1] (ordre croissant, d'où le signe -)
    
    # Chebyshev coefficients 
    c = np.ones(N + 1)
    c[0], c[-1] = 2, 2
    c *= (-1) ** np.arange(N + 1)
    
    # Computation of matrix D
    X = np.tile(x, (N + 1, 1)).T
    dX = X - X.T + np.eye(N + 1)
    D = (np.outer(c, 1 / c)) / dX
    D -= np.diag(np.sum(D, axis=1))
    
    # s s on [0, 1]
    s = 0.5 * (x + 1)
    
    return 2 * D, s

# === DISCRÉTISATION ===
N = 100
D, s = cheb_diff(N)
D2 = D @ D
n_p = N + 1
size = 4 * n_p

T_idx = np.arange(0, n_p)
P_idx = np.arange(n_p, 2 * n_p)
G_idx = np.arange(2 * n_p, 3 * n_p)
F_idx = np.arange(3 * n_p, 4 * n_p)

A = np.zeros((size, size), dtype=complex)
B = np.zeros((size, size), dtype=complex)

# === SYSTEM ===
for j, sj in enumerate(s):
    
    # Equation 1 (curl of the motion eq., radial projection)
    row = j
    A[row, P_idx[j]] = -2j * Omega_rot * (k**3 * sj**2 + k * m**2)
    A[row, G_idx[j]] = -1j * B0 / (mu0 * rho) * (k**3 * sj**2 + k * m**2)
    A[row, F_idx[j]] = -2j * B0 / (mu0 * rho) * (k**2 * m)
    B[row, T_idx[j]] = 1j * Omega_mass * (k**2 * sj**2 + m**2)
    B[row, P_idx[j]] = 2j * Omega_mass * k * m
    
    # Equation 2 (double curl of the motion eq., radial projection)
    row = n_p + j
    A[row, T_idx[j]] = 2j * Omega_rot * (k**3 * sj**4 + k * m**2 * sj**2)
    A[row, P_idx[j]] = 4j * Omega_rot * (k**2 * m * sj**2)
    A[row, G_idx[j]] = 2j * B0 / (mu0 * rho) * (k**2 * m * sj**2)
    coef0_F = 1j * B0 / (mu0 * rho) * (k**5 * sj**4 + k**3 * (2 * m**2 + 1) * sj**2 + k * (m**4 - m**2))
    coef1_F = 1j * B0 / (mu0 * rho) * (-k**3 * sj**3 + k * m**2 * sj)
    coef2_F = -1j * B0 / (mu0 * rho) * (k**3 * sj**4 + k * m**2 * sj**2)
    A[row, F_idx] += coef0_F * np.eye(n_p)[j] + coef1_F * D[j] + coef2_F * D2[j]
    B[row, T_idx[j]] = -2j * Omega_mass * (k * m * sj**2)
    coef0_P = 1j * Omega_mass * ((m**2 - m**4) - (2 * m**2 + 1) * k**2 * sj**2 - k**4 * sj**4)
    coef1_P = 1j * Omega_mass * (k**2 * sj**3 - m**2 * sj)
    coef2_P = 1j * Omega_mass * (k**2 * sj**4 + m**2 * sj**2)
    B[row, P_idx] += coef0_P * np.eye(n_p)[j] + coef1_P * D[j] + coef2_P * D2[j]
    
    # Equation 3 (induction equation, radial projection)
    row = 2 * n_p + j
    A[row, P_idx[j]] = -1j * B0 / (mu0 * rho) * (k**3 * sj**2 + k * m**2)
    B[row, F_idx[j]] = 1j * (k**2 * sj**2 + m**2)

    # Equation 4 (curl of induction equation, radial projection)
    row = 3 * n_p + j
    A[row, T_idx[j]] = 1j * B0 / (mu0 * rho) * (k**3 * sj**2 + k * m**2)
    A[row, P_idx[j]] = 2j * B0 / (mu0 * rho) * (k**2 * m)
    B[row, G_idx[j]] = -1j * (k**2 * sj**2 + m**2)
    B[row, F_idx[j]] = -2j * k * m

# === BOUNDARY CONDITIONS ===
# Match with B_ext
def compute_F1_dF1_G1(m, k, C):
    F1 = ((-1j*C*k)/2*(k**2 + m**2)) * (yv(m + 1, -1j * k) - yv(m - 1, -1j * k))
    dF1 = F1 - C * yv(m, -1j * k)
    G1 = 0
    return F1, dF1, G1

F1_val, dF1_val, G1_val = compute_F1_dF1_G1(m, k, C)

# Condition in s=0 : P(0)=0 or dP(0)=0 regarding m
if m == 0 or m % 2 == 1:
    A[P_idx[0], :] = 0
    A[P_idx[0], P_idx[0]] = 1
    B[P_idx[0], :] = 0
else:
    A[P_idx[0], :] = 0
    A[P_idx[0], P_idx] = D[0, :]
    B[P_idx[0], :] = 0
    
A[P_idx[-1], :] = 0
A[P_idx[-1], P_idx[-1]] = 1
B[P_idx[-1], :] = 0

if solve_mhd :
    A[G_idx[-1], :] = 0
    A[G_idx[-1], G_idx[-1]] = 1
    B[G_idx[-1], :] = 0

    A[F_idx[-1], :] = 0
    A[F_idx[-1], F_idx[-1]] = 1
    B[F_idx[-1], :] = 0
    B[F_idx[-1], F_idx[-1]] = F1_val
else : 
    pass

# === SYSTEM SOLVING ===
omega_all, VEP_all = eig(A, B)
mask = np.isfinite(omega_all) & (np.abs(omega_all.imag) < 1e-4)
omega = omega_all[mask]
VEP = VEP_all[:, mask]

# === SOLUTION DOMINANTE ===
idc_dom = np.argsort(-omega.real)
omega_phys = omega[idc_dom]
top_omegas = omega_phys[:50]
u_dom = VEP[:, idc_dom[0]]
u_total = u_dom

T_vec = u_total[T_idx]
P_vec = u_total[P_idx]
G_vec = u_total[G_idx]
F_vec = u_total[F_idx]

# === AFFICHAGE DES RÉSULTATS ===
#### Greenspan if Fluid case ####
if solve_mhd == False : 
    # Affichage des résultats
    print("Greenspan's eigenvalues vs MCWS : ")
    print(" " * 67)
    for i, (xi_val, lambda_val, omega_val) in enumerate(zip(roots[:5], lambdas[:5], top_omegas[:5]) ):
        print(rf"ξ_{i+1} = {xi_val:.5f}   ->   λ_{i+1} = {lambda_val:.5f} | Re(ω)_{i+1} = {omega_val.real:.5f}")
    print(" " * 67)    
    # Comparaison with Greenspan's eigenvalues for Hydro
    plt.figure(figsize=(8, 3))
    plt.scatter(np.arange(len(roots)), lambdas, marker='_',label='Greenspan')
    plt.scatter(np.arange(len(top_omegas)), np.real(top_omegas), marker='.',label='Chebychev')
    plt.xlabel(r"$i$", fontsize=12)
    plt.ylabel(r"$\lambda_{knm}$ - Re$(\omega_{phys})$", fontsize=12)
    plt.title(rf" Greenspan's $\lambda$ vs MCWS's $\omega$ pour $n={n_axial}$, $m={m}$")
    plt.legend()
    plt.grid(True)
    plt.show()

    # Error 
    error = np.abs(lambdas - top_omegas)
    MSE = np.real(np.square(np.subtract(lambdas,top_omegas)).mean())
    print(f"MSE = {MSE:.4}")
    print(" " * 67)
    plt.figure(figsize=(10, 6))
    plt.subplot(2, 2, 1)
    plt.plot(s, np.abs(T_vec),color = 'royalblue'); plt.title("|T(s)|"); plt.grid()
    plt.subplot(2, 2, 2)
    plt.plot(s, np.abs(P_vec),color = 'royalblue'); plt.title("|P(s)|"); plt.grid()
else : 
    #### Lundquist ####
    if solve_Lundquist :
        # previous code results
        omega_MCWS_rad = np.array([3.83492, 7.19234, 10.39528, 13.55841, 16.70897,
                                   19.85460, 22.99808, 26.14051, 29.28247, 32.42410], dtype=float)
        omega_MCWS_ax  = np.array([3.14159, 6.28319, 9.42478, 12.56637, 15.70796,
                                   18.84957, 21.99115, 25.13277, 28.27434, 31.41593], dtype=float)

        print(f"Ideal MHD (Lundquist) eigenfrequencies :")
        print(" " * 67)
        print(f"Alfvén speed V_A = {float(V_a):.1f}\n")
        print(" " * 67)
        print(r"|   n |  omega_L_rad |omega_MCWS_rad|   omega_L_ax | omega_MCWS_ax|")
        print("-" * 67)
        for n, olr, omr, ola, oma in zip(range(1, N_l + 1), omega_L_rad, omega_MCWS_rad, omega_L_ax, omega_MCWS_ax):
            print(f"| {n:3d} | {float(olr):12.6f} | {float(omr):12.6f} | {float(ola):12.6f} | {float(oma):12.6f} |")

        # Check d’erreur
        abs_err_rad = omega_MCWS_rad - omega_L_rad
        abs_err_ax  = omega_MCWS_ax  - omega_L_ax
        print("\nMAE_rad = {:.3e}, MAE_ax = {:.3e}".format(np.mean(np.abs(abs_err_rad)), np.mean(np.abs(abs_err_ax))))
    
    #### MHD GENERAL ####
    else : 
        print(f"Top 5 eigenvalues for (m, n) = ({m}, {n_axial}) :")
        for i, val in enumerate(top_omegas[:5]):
            print(rf"{i+1}. Re(ω) = {val.real:.5f}, Im(ω) = {val.imag:.5e}")
        plt.figure(figsize=(8, 4))
        plt.plot(range(1,51), np.real(omega_phys[:50]), '.-', 
                 linewidth='1', markeredgewidth=2, color = 'orange' if solve_mhd else 'royalblue')
        plt.xlabel(r"$i$")
        plt.ylabel(r"Re$(\omega_{phys})$")
        plt.title(rf"Eigenvalues $\omega$ for (m, n) = ({m}, {n_axial})")
        plt.grid(True)
        plt.show()
        
        # Eigenvectors
        plt.figure(figsize=(10, 6))
        plt.subplot(2, 2, 1)
        plt.plot(s, np.abs(T_vec),color = 'royalblue'); plt.title("|T(s)|"); plt.grid()
        plt.subplot(2, 2, 2)
        plt.plot(s, np.abs(P_vec),color = 'royalblue'); plt.title("|P(s)|"); plt.grid()
        plt.subplot(2, 2, 3)
        plt.plot(s, np.abs(G_vec),color = 'orange'); plt.title("|G(s)|"); plt.grid()
        plt.subplot(2, 2, 4)
        plt.plot(s, np.abs(F_vec),color = 'orange'); plt.title("|F(s)|"); plt.grid()
        plt.tight_layout()
        plt.show()
