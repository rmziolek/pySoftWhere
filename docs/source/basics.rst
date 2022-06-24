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

    import pysoftwhere.clustering.cluster as pysw_cluster
    import pysoftwhere.ICSI.icsi as pysw_isci
    import pysoftwhere.IPI.ipi as pysw_ipi



.. _ipi_simple:

Intrinsic Planar Interface (IPI)
----------------------------------

Firstly, define an MDAnalysis Universe

.. code-block:: python3

   u = mda.Universe('system.tpr','system.xtc')


Now call pysw_ipi to calculate the intrinsic density for a single frame

.. code-block:: python3

    intrinsic_z,interface = pySW_IPI.ipi(u,
                                         frame=-1,
                                         core_sel='resname UNK and prop mass >2',
                                         density_sel='name OW',
                                         interface='Lower',                        #Lower or Upper... 
                                         no_bins=51,
                                         recombine=True,                           ###set true if you're using a slab/bilayer (not seperated monolayers!)
                                         cluster=True,                             #if you want to check if anything has escaped slab/bilayer so it's not included
                                         interpolate_interface=False)()             ###can use this if the layer is patchy

