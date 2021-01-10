# Physics-aware registration based auto-encoder for convection dominated PDEs

#### [[project website]](https://arxiv.org/abs/2006.15655)
<img src="data/schematic.png" width="250">

## Table of contents
* [Introduction](#introduction)
* [Method](#Method)
    * [Video-description](#Video-description)
* [FAQ](#FAQ)
* [Requirements](#Requirements)
* [Experiments](#Experiments)
    * [Rotating A](#Rotating-A)
    * [Two-dimensional fluid flows](#Two-dimensional-fluid-flows)
    * [Physics aware auto-encoder in an LSTM architecture](#Physics-aware-auto-encoder-in-an-LSTM-architecture)
* [Citation](#citation)

## Introduction
We design a physics-aware auto-encoder to specifically reduce the dimensionality of solutions arising from convection-dominated nonlinear physical systems. Although existing nonlinear manifold learning methods seem to be compelling tools to reduce the dimensionality of data characterized by a large Kolmogorov n-width, they typically lack a straightforward mapping from the latent space to the high-dimensional physical space. Moreover, the realized latent variables are often hard to interpret. Therefore, many of these methods are often dismissed in the reduced order modeling of dynamical systems governed by the partial differential equations (PDEs). Accordingly, we propose an auto-encoder type nonlinear dimensionality reduction algorithm. The unsupervised learning problem trains a diffeomorphic spatio-temporal grid, that registers the output sequence of the PDEs on a non-uniform parameter/time-varying grid, such that the Kolmogorov n-width of the mapped data on the learned grid is minimized. We demonstrate the efficacy and interpretability of our approach to separate convection/advection from diffusion/scaling on various manufactured and physical systems.

## Method

### Video-description 
- [Low-rank registeration based manifolds](https://youtu.be/fDYPAj9WAbk)

## FAQ
- How does the method handle noisy data?
	- SVD (singular value decomposition) and truncation at the heart of the algorithm acts as a filter removing the low energy-containing features of the data, i.e. the noise is filtered as a result of SVD-truncate.

- Why the method is considered *physics aware*?
	- The existance of a low-rank grid that minimize the Kolmogrov n-width is an assumption based on the physics of many of the convection dominated flows read sec 3.1 of [this](https://arxiv.org/abs/1701.04343).

- What interpolation scheme to use?
	- Any of the off-the-shelf interpolation schemes, e.g. (bi-)linear, cubic, spline. High order interpolation schemes only becomes advantageous in higher rank of reconstruction, i.e. higher \\(k_r\\). This is due to the local aliasing of high wave-number bases (features) on the coarsened grid.

## Requirements
- Matlab R2016+
- python 3.6
	- [scipy](https://pypi.org/project/scipy/)
	- [numpy](https://pypi.org/project/numpy/)
- [TensorFlow 2](https://www.tensorflow.org/install)
- [Keras 2.3.1](https://pypi.org/project/Keras/)

## Experiments
### Rotating A
Rotating A [Rotated A Location](./Experiments/rotatedA) 

open matlab
```
matlab -nodisplay -nosplash
```

Run the manifold learning
```
run main_rotatingA.m
```

Evaluate the snapshots on the learned manifold
```
run post_process.m
```

### Two-dimensional fluid flows
2D Fluid Flows (2D Riemann) [2D Fluid Flows](./Experiments/2DRiemann) 

open matlab
```
matlab -nodisplay -nosplash
```

Run the manifold learning
```
run main_opt_config03.m
run main_opt_config12.m
```
Evaluate the snapshots on the learned manifold
```
main_solve_config03.m
main_solve_config12.m
```


### Physics aware auto-encoder in an LSTM architecture
LSTM architecture [LSTM architecture](./Experiments/LSTM) 
* The notebooks are Google colab ready, make sure to have .py, *.pkl, *.h5 files at the same directory of the notebooks (*.ipynb)

#### Burgers' equation [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](./Experiments/LSTM/Burgers/main_manifold_burgers.ipynb)

Run the manifold learning and train the LSTM with Physics-aware registration based auto-encoder
```
main_manifold_burgers.ipynb
```

Run the manifold learning and train the LSTM with Neural Network based auto-encoder
```
main_eulerian_burgers.ipynb
```

#### Wave equation [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](./Experiments/LSTM/Wave/main_manifold_wave_small.ipynb)

Run the manifold learning and train the LSTM with Physics-aware registration based auto-encoder
```
 main_eulerian_wave_small.ipynb
```

Run the manifold learning and train the LSTM with Neural Network based auto-encoder
```
main_manifold_wave_small.ipynb
```


## Citation
To appear in AAAI-21; Meanwhile, read more on [[arXiv]](https://arxiv.org/abs/2006.15655)
```
@article{Mojgani_arxiv_2020,
	author = {Mojgani, Rambod and Balajewicz, Maciej},
	title={Physics-aware registration based auto-encoder for convection dominated {PDE}s},
	journal={arXiv preprint arXiv:2006.15655},
	archivePrefix = "arXiv",
	eprint = {2006.15655},
	year = 2020,
}
```
