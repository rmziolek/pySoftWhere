Basic Use
=====

.. _importing:

Importing pySoftWhere
----------------------

.. code-block:: python3

    import MDAnalysis as mda
    import matplotlib.pyplot as plt
    import numpy as np
    import sys


Set the directory where pySoftWhere is located (if a full install has not been completed)

.. code-block:: python3

    sys.path.append('<path to>/pysoftwhere') 

Import the pySoftWhere functions

.. code-block:: python3

    import pysoftwhere.clustering.cluster
    import pysoftwhere.ICSI.icsi
    import pysoftwhere.IPI.ipi



.. _ipi_simple:

Intrinsic Planar Interface (IPI)
----------------------------------

Firstly, define an MDAnalysis Universe

.. code-block:: python3

   u = mda.Universe('system.tpr','system.xtc')


Now call ``pysoftwhere.IPI.ipi()`` to calculate the intrinsic density for a single frame

.. code-block:: python3
    
   ''' 
   u:           an MDAnalysis Universe that contains bond information
   core_sel:    an MDAnalysis atom selection of the atoms used to construct the intrinsic interface
   density_sel: an MDAnalysis atom selection of the atoms whose intrinsic positions will be calculated
   interface"  
   '''

    intrinsic_z,interface_vals  =  pysw_ipi(u,
                                            frame=-1,
                                            core_sel='resname UNK and prop mass >2',
                                            density_sel='name OW',
                                            interface='Lower',                        #Lower or Upper...                                 
                                            no_bins=51,
                                            recombine=True,                           ###set true if you're using a slab/bilayer (not seperated monolayers!)
                                            cluster=True,                             #if you want to check if anything has escaped slab/bilayer so it's not included
                                            interpolate_interface=False)()             ###can use this if the layer is patchy


The ``pysoftwhere.IPI.ipi()`` function returns two results when called. These are ``intrinsic_z``,  a list of the intrinsic positions of atoms selected in ``density_sel``, and ``interface```, an array containing the z-positions of the intrinsic interface. 
