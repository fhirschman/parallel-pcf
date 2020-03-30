# parallel-pcf
Python script for parallel calculation of pair correlation functions of simulation trajectories


### Prerequisites
Several python packages are needed:
```
MDAnalysis
matplotlib
multiprocessing
numpy
ctypes
```

### Running

An exemplary Brownian Dynamics simulation trajectory of
343 coarse-grained, 6-sited molecules containing 100 frames
is provided in 'traj.gsd'. By running

```
$ python main.py
```
the script will automatically use all available CPU threads
to calculate the pair-correlation function of centers-of-mass
of all molecules, more precisely:

<img src="https://render.githubusercontent.com/render/math?math=
g(r) =  \frac{1}{\rho N} \left\langle \sum_{m=1}^{N}  \sum_{n\neq m}^{N} \delta \left[ \boldsymbol{r} - \left( \boldsymbol{r}_m - \boldsymbol{r}_n \right) \right] \right\rangle ">

<img src="https://render.githubusercontent.com/render/math?math=%3D%20%20%5Cfrac%7B1%7D%7B%5Cfrac%7B4%20%5Cpi%7D%7B3%7D%20%5Cleft((r%2B%5CDelta%20r)%5E3-r%5E3%5Cright)%7D%20%20%5Cfrac%7B1%7D%7B%5Crho_i%20N_i%7D%20%5Csum_%7Bm%3D1%7D%5E%7BN%7D%20%5Csum_%7Bn%5Cneq%20m%7D%5E%7BN%7D%20%5Cleft%5C%7B%20%5Cbegin%7Barray%7D%7Blr%7D%201%20%26%20%20%5Ctext%7Bif%7D%20%20%5C%2C%20%5C%2C%20%5C%2C%20%7C%20%5Cboldsymbol%7Br%7D_m%20-%20%5Cboldsymbol%7Br%7D_n%7C%20%5Cin%20%5Br%2Cr%2B%5CDelta%20r%5D%20%5C%5C%200%20%26%20%20%5Ctext%7Belse%7D%20%5C%2C%20%5C%2C%20%5C%2C%20%5C%2C%20%20%5C%2C%20%20%5C%2C%20%5C%2C%20%5C%2C%20%5C%2C%20%5C%2C%20%5C%2C%20%5C%2C%20%5C%2C%20%5C%2C%20%5C%2C%20%5C%2C%20%5C%2C%20%5C%2C%20%20%20%5C%2C%20%5C%2C%20%5C%2C%20%5C%2C%20%5C%2C%20%5C%2C%20%5C%2C%20%5C%2C%20%5C%2C%20%5C%2C%20%5C%2C%20%5C%2C%20%5Cend%7Barray%7D%20%5Cright.">

\
It will output the result in 'pcf.dat' and also a produce a graph 'pcf.pdf':

![pcf](./pcf.png)
