Preliminary implementation of the ONA model (Ghobadi and Jayaraman, Soft Matter, 2016, 12, 2276) on MoSDeF platform (https://mosdef.org).

Installation Guide:
1. Create a new conda env, and install mbuild and foyer. 
    * In linux shell:
      * `conda create -n mosdef37 -c mosdef -c conda-forge -c omnia mbuild foyer python=3.7 openbabel`
1. Install mbuild_ONA.
1. Compile mbuild_ONA.
    * In linux shell:
      * `pip install -e .`

The topology built in this tool can either be saved out as a bare-bone untyped file to memory, in desired format supported by mbuild (as shown in examples), or fed to other mosdef backend tools, such as foyer and gmso, for atomtyping and force field application. Please refer to https://github.com/zijiewu3/forcefield_ona or https://github.com/mosdef-hub/mosdef-workflows for more details. 

Thank you very much!
