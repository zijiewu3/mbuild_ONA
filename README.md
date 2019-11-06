A preliminary implementation of the ONA model (Ghobadi and Jayaraman, Soft Matter, 2016, 12, 2276) on MoSDeF platform (https://mosdef.org).

Guide for testers:
1.Create a new conda env, and install mbuild and foyer. In linux shell,
conda create -n mosdef35 -c mosdef -c conda-forge -c omnia mbuild foyer python=3.7 openbabel
2.Install mbuild_ONA. Download or clone this repo and run 
pip install -e .
3.Refer to example.py or example_from_mbuild.py

Thank you for helping us and trying this out!
