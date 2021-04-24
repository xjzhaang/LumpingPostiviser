
# LumpingPositiviser

## How to use:

**Require**: *.ode file and *.obs file. 
*.obs file is a file containing lists of observables with one observable per line. For example a file:

   ["S0"]
   
   ["S0 + S1"]
   
   ["S0", "S1"]
   
**Input**:
```sh 
python3 run.py ../model.ode ../model.obs
``` 
will output three files in the directory: model_CLUE_result.txt, model_final_macrovariables.txt, model_final_ode.txt
