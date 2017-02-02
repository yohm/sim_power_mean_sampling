# A model of network sampling by communication channel

A simulation code for the model for sampling social network by communication channel.
The model is proposed in paper entitled "What Big Data tells: Sampling the social network by communication channels" published in Phys.Rev.E [link](https://journals.aps.org/pre/abstract/10.1103/PhysRevE.94.052319).
The paper is also on [arXiv](https://arxiv.org/abs/1511.08749).

# Model definition

First, a Erdos-Renyi random network of size `N` and average degree `k0` is constructed.

Then, a scalar value , called "affinity" of node i, is assigned randomly to each node.
The affinity values are taken from the Weibull distribution with exponent `alpha` and a characteristic scale `f0`.

After affinity values are assigned to each node, then sample the links with the probability dependent on the affinity of the nodes in both ends.
The sampling probability is given by `p_{ij} = p(fi,fj)`, where p(x,y) is the generalized mean with exponent `beta`.
See section V of the paper for more accurate description.

# How to Run

## Compiling

To copmile the program, you need boost C++ library.
You do not have to compile boost since header-only libraries are required.
Put the boost library in an appropriate path. Then, run `make` in "simulator/" and "network_analysis" directories.

```sh
cd simulator && make && cd -
cd analyzer && make && cd -
```

If you would like to specify the compiler and/or include path explicitly, set `CXX` and `INCLUDE` environment variables.

```
env CXX=g++ INCLUDE=-I/usr/local/include make
```

## Running

To run the simulation, run the following code and you'll find an output file named `sampled.edg`.

```
./sampling.out <input net file> <alpha> <beta> <f0> <q> <seed>
```

Instead of specifying `f0`, you can also specify the sampled degree. In this case, f0 is tuned so that the sampled degree becomes `<k>+-<dk>`.

```
./sampling_targeted_k.out <input net file> <alpha> <beta> <k> <dk> <q> <seed>
```

Note: When beta < -5, we use min(fi,fj) instead of the generalized mean for numerical stability. For the same reason, we use max(fi,fj) when beta > 5.

