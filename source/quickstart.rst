.. _qstart:

==============
Quickstart
==============

-------------------
Directory structure
-------------------

---------------------
Case files
---------------------

-------------------
Scripts
-------------------
In this section, some useful batch scripts for compiling and running a NekIBM cases. The folloing steps are listed for running a standard Nek5000/CMTnek case. 

preprocessing
- ``genbox <case>`` prepare flow domain and spectral element discretization
- ``genmap <case>`` mapping

compling NekIBM
- ``makenek <case>`` compiles your case

running a parallel case

- ``nekmpi/nekbmpi <case> <number of ranks>`` runs a parallel job

postprocessing
- ``visnek <case>`` creates metadata file required by `VisIt <https://wci.llnl.gov/simulation/computer-codes/visit/>`_ and `ParaView <https://www.paraview.org/>`_. 


-------------------
Preprocessing
-------------------

Nek5000 is mainly a solver. However, simple box type meshes can be generated with the ``genbox`` tool. For more complex meshes please consider using ``PRENEK`` and the meshing tools ``nekmerge`` and ``n2to3``. We provide mesh converters like ``exo2nek`` and ``msh2nek`` which are quite handy if you want to use your favorite mesh generator. Also check our 
`Bazaar <https://github.com/Nek5000/NekBazaar>`_ for 3rd party meshing tools.

.. _qstart_vis:

-------------------
Postprocessing 
-------------------

- ``Flow visualization`` 
  
Nek5000 output (``.fld`` or ``0.f%05d``) files can be read by `VisIt <https://wci.llnl.gov/simulation/computer-codes/visit/>`_ or 
`ParaView <https://www.paraview.org/>`_. This requires using ``visnek`` to generate a metadata file.  
There is also an build-in X-Window based postprocessor called ``POSTNEK`` located in tools.

- ``Particle visualization`` 
  
