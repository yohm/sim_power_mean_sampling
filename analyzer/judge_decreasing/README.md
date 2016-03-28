To prepare python environment, introduce miniconda.

Then, create an environment called network as follows:

```
conda create -n network numpy networkx
```

run.sh loads this environment and then run the python script.

```
./run.sh sampled.edg > _output.json
```

