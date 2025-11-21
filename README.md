# Linear Programming Control (LPC)

A **computationally efficient control framework** for **constrained input-affine nonlinear systems**, developed by **Ali Safi** as part of the research paper:

> *A Computationally Efficient Linear Programming Control Framework for Constrained Input-Affine Nonlinear Systems*  
> by Ali Safi, Ali Taghavian, Esmaeel Khanmirza and Fateme Namdarpour (2025)

---

## 🧠 Overview

**Linear Programming Control (LPC)** is an optimization-based control method that computes control inputs for nonlinear systems using a **linear programming (LP)** formulation.  
The approach operates directly on the nonlinear dynamics and handles **input, state, and sampling-time constraints** efficiently.

The LPC framework:
- Uses a **linear optimization problem** at each iteration (not iterative MPC).  
- Achieves **real-time feasibility** even on low-cost processors.  
- Provides **Lyapunov-based stability guarantees**.  
- Includes an **optional variable sampling-time mechanism** for adaptive control frequency.  

LPC was conceptually inspired by the *Fine-Tuner* section of the authors’ previous **Computational Hybrid Controller (CHC)**, later developed into a **standalone, generalized controller** for continuous nonlinear systems.

---

## ⚙️ Features
- Handles **constrained input-affine nonlinear systems** directly.  
- **Lightweight** and **easy to implement** in MATLAB or Python.  
- Compatible with **SISO**, **MIMO**, and **underactuated systems**.  
- Includes example simulations for:
  - Single-Link Robot Arm  
  - Revolute–Prismatic Robot  
  - Rotary Inverted Pendulum  

---

## 🧩 Installation

Clone the repository:

```bash
git clone https://github.com/A-Safi/LPC.git
cd LPC

### Rotary Inverted Pendulum (RIP) — LPC Demo
Run `rip_lpc_demo.m` to reproduce the LPC results on the RIP benchmark with adaptive sampling.
- Constraints: $\dot{\alpha}, \dot{\theta} \in [-5,5]$ rad/s, $u \in [-1.2,1.2]$ N·m
- Sampling: $t \in [0.01, 0.02]$ s (adaptive)
- Plots: states, input torque, and sampling interval

