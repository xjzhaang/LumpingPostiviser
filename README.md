# LumpingPositiviser

This repository contains the implementation of the algorithm designed in the paper *Interpretable exact linear reductions via positivity* by Gleb Pogudin and Xingjian Zhang.
The code takes as input an ODE model and a set of observables, performs exact linear reduction using [CLUE](https://github.com/pogudingleb/CLUE), and then performs re-parametrization of the reduction aiming at improving the interpretability of the model.

## How to use:

The input ODE model and observables should be stored in \*.ode and \*.obs files, respectively.

### The format of \*.ode files

The detailed description of the format can be found [here](https://doi.org/10.1007/978-3-030-31304-3_13).
Briefly, one can use the format to write either an ODE system directly or to specify a chemical reaction network.

An example of defining an ODE system directly:

```
begin model ExampleODE
  begin ODE
    d(x) = x + y
    d(y) = x^2 - 3 * y
  end ODE
end model
```

An example of defining an ODE system as CRN:

```
begin model ExampleCRN
  begin reactions
    A + B -> C, k1
    B + 2 * C -> A, k2
  end reactions
end model
```

### The format of \*.obs files

An \*.obs file contains the observables (linear combinations of the variables of the system) one by line.

### How to run

```sh 
python3 run.py model.ode model.obs [path_to_julia_binary]
``` 

If the path to the [julia](https://julialang.org/) binary is not given, a shell command `julia` will be used.
As the result of the computation, there will be two files model_final_macrovariables.txt and model_final_ode.txt containing the new variables and the reduced ODE system, respectively.
