# smtcutpartial
Compute tropical partial equilibrations and associated truncated systems.

To install the program you need to install the following packages: sympy, gmpy2, pysmt (you will need to install a solver via the command pysmt-install --smat (for the solver msat, else you can check the list solvers via pysmt-install --check). 
Under Windows, there is a linkage problem with msat, then you should make a pyproj project under visual sutdio and run the program with the debugger.

The program take take as argument the folder when you can find four files : conservation_laws_q.txt, parameters_q.txt, slow.txt, vector_field.txt.
You can use the option -p <INTEGER> to change the precision, -e <INTEGER> to change the epsilon (it will be of the form 1/<INTEGER>).
  
The code have been tested with python 3.7 and the solver msat.

To analyse results, we provide some maple files. 
