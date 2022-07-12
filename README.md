# Sbox_LinIneq_Reduction_DM
Divide and Merge approach for reduction of linear inequalities corresponding to the DDT of S-box\
**_This repository is a part of the research paper to reduce number of linear inequalties corresponding to the DDT of S-box.\
This repository contains four files._** 
* _Inequalities_Reduction_Divide_and_Merge.py_
* _CPLEX_Problem_WARP_sbox_-_10_2.lp_
* _CPLEX_Solution_WARP_sbox_-_10_2.sol_
* _ineq_WARP_sbox_-_10_2.txt_

## Reduction of Linear Inequalities 
_Inequalities_Reduction_Divide_and_Merge.py_ contains the source code to 
* Compute DDT of S-box
* To generate linear inequalities corresponding using H-representation of convex hull (pycddlib library used)
* To apply divide and merge appraoch to introduce new set of linear inequalities
* To reduce the number of linear ineuqalities by constructing MILP problem whcih is solved using GUROBI/CPLEX solver

## Parameters 
1 - name of cipher (PRESENT/GIFT/WARP/TWINE/ASCON/FIDES-5/SC2000-5)\
2 - sbox or prob   (This research uses sbox)\
3 - DDT points (2/4/6), '-' for all  (This research uses '-')\
4 - If do not want to introduce new inequalities then use '-' otherwise specify batch size(beta) in which inequalities are to be added e.g.'10' or '50', if batch size is '0' then every inequality will added with every other inequality (specifying '0' is same as len(inequalities), 0 must be used to get results equivalent to Boura and Coggia's approach(https://hal.inria.fr/hal-03046211/document) of inequalities reduction)\
5 - k(=2/3/4) inequalities to be added (it should not be 1), if parameter 4 is '-' then use the parameter 5 as '-'

## Examples
```python Inequalities_Reduction_Divide_and_Merge.py WARP sbox - 10 2```\
```python Inequalities_Reduction_Divide_and_Merge.py WARP sbox - 20 3```\
The .lp and .sol files are created corresponding to MILP problem and solution.\
The minimized set of inequalities will be printed on the screen and will be written to .txt file.

## Acknowledgement
1. https://github.com/stephane-caron/pypoman
2. https://pypi.org/project/pycddlib