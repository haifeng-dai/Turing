# Impact of Noise on Turing Patterns of Multiplex Networks

This repository provides the official implementation and source code for the simulations and manuscript of the research paper: **"Impact of Noise on Turing Patterns of Multiplex Networks"**.

## Abstract

This study investigates the stochastic dynamics of Turing instability in activator-inhibitor systems (e.g., FitzHugh-Nagumo and Mimura-Murray models) residing on multiplex networks. By incorporating multiplicative noise into both intra-layer and inter-layer interactions, we explore how the network's super-Laplacian structure and stochastic perturbations jointly determine the threshold and stability of pattern formation. The simulations are conducted using a robust Euler-Maruyama scheme to handle the high-dimensional stochastic differential equations (SDEs) inherent in multiplex networked systems.

## Core Mathematical Model

The system is modeled as a set of coupled SDEs on a multiplex network with $K$ layers and $N$ nodes per layer:

$$ \dot{\mathbf{x}} = h(\mathbf{x}) - \alpha (\mathcal{L} \otimes \varsigma) \mathbf{x} - \eta \alpha (\mathcal{L} \otimes \varsigma) \mathbf{x} \xi $$

where:
- $\mathbf{x} = [u, v]^\top$ represents the state vector of activator and inhibitor densities.
- $\mathcal{L}$ is the multiplex super-Laplacian matrix.
- $\alpha$, $\beta$ are the intra-layer and inter-layer diffusion coefficients.
- $\eta$ denotes the intensity of multiplicative noise.
- $\xi$ represents independent Gaussian white noise.

## Repository Contents

### 📄 Manuscript
- `manuscript/main.tex`: LaTeX source code of the paper.
- [Scientific Manuscript Style]: Prepared according to the *Science China Information Sciences* template.

### 💻 Source Code (MATLAB)
The code is organized to ensure full reproducibility of the results presented in the paper.

#### 1. Simulation Engines (`code/simulations/`)
- `solve_multiplex.m`: The core numerical solver implementing the stochastic Euler-Maruyama algorithm.
- `find_thresholds.m`: Analyzes the shifts in Turing thresholds under varying noise intensities.
- `sweep_param.m`: Simulates phase transitions and hysteresis loops in pattern formation.
- `simulate_sde.m`: Simulates the time-evolution of patterns using Euler-Maruyama.

#### 2. Analysis and Visualization (`code/`)
The scripts in the root of the `code` directory correspond to the figures in the manuscript:
- `phase_diagram.m`: Generates threshold analysis plots (Corresponds to Fig. 3c).
- `compare_connectivity.m`: Generates network connectivity influence plots (for ER networks).
- `pattern_evolution.m`: Visualizes the dynamic emergence of patterns across layers (Fig. 1b-e).
- `hysteresis_loop.m`: Plots transition curves/hysteresis loop (Fig. 1a).
- `threshold_shifts.m`: Detailed envelopes of critical values across parameters.
- `noise_analysis.m`: Detailed dissection of noise impact on bifurcation boundaries.

#### 3. Network Generation (`code/networks/`)
Implementations for generating different multiplex topologies:
- Barabási–Albert (BA), Erdős–Rényi (ER), Watts–Strogatz (WS), and Scale-Free (SF) with controllable degree exponents.

## Simulation Pipeline

To reproduce the study's results:
1. **Initial Setup**: Ensure the system parameters (N, K, diffusion constants) in the simulation scripts match those specified in Section 4 of the paper.
2. **Execute Simulations**: Run the simulation scripts in `code/simulations/`. Results are saved to `.mat` files in `code/results/`.
3. **Data Analysis**: Run the corresponding `*.m` scripts to process the raw output and generate the manuscript figures.

## Requirements
- MATLAB R2021a or higher.
- (Optional) LaTeX environment for compiling the manuscript.

## Citation

```bibtex
@article{Dai2025Turing,
  title={Impact of noise on Turing patterns of multiplex networks},
  author={Dai, Haifeng and Sun, Yongzheng and Wen, Guanghui and Chen, Guanrong},
  journal={Science China Information Sciences},
  year={2025}
}
```

---
**Author Correspondence**: Haifeng Dai (Southeast University / City University of Hong Kong)
