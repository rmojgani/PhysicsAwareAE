# Physics-aware registration based auto-encoder for convection dominated PDEs
#### [[project website]](http://www.rmojgani.com)
<img src="data/schematic.png" width="250">

**Requirements**
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

### Two-dimensional fluid flows (2D Riemann Flows)

### Physics based auto-encoder in an LSTM architecture

#### Burgers' equation

Run the manifold learning
```
python 
```

Run post process
```
python 
```

#### Wave equation
