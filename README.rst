#################
Simple VOF Solver
#################

THIS IS AN ON-GOING PROJECT.

|License|_ |LastCommit|_

.. |License| image:: https://img.shields.io/github/license/NaokiHori/SimpleVOFSolver
.. _License: https://opensource.org/licenses/MIT

.. |LastCommit| image:: https://img.shields.io/github/last-commit/NaokiHori/SimpleVOFSolver/main
.. _LastCommit: https://github.com/NaokiHori/SimpleVOFSolver/commits/main

.. image:: https://github.com/NaokiHori/SimpleVOFSolver/blob/main/docs/source/snapshot2d.png
   :width: 100%

.. image:: https://github.com/NaokiHori/SimpleVOFSolver/blob/main/docs/source/snapshot3d.png
   :width: 50%

********
Overview
********

This library numerically solves the motion of sharp interfacial structures in two- and three-dimensional Cartesian domains using finite-difference and volume-of-fluid methods.
This is built on top of `SimpleNSSolver <https://github.com/NaokiHori/SimpleNavierStokesSolver>`_.

**********
Dependency
**********

* `C compiler <https://gcc.gnu.org>`_
* `MPI <https://www.open-mpi.org>`_
* `FFTW3 <https://www.fftw.org>`_

Docker image is available.

.. code-block:: console

   $ mkdir /path/to/your/working/directory
   $ cd    /path/to/your/working/directory
   $ git clone https://github.com/NaokiHori/SimpleVOFSolver
   $ cd SimpleVOFSolver
   $ docker build -t simplenavierstokessolver:latest .

***********
Quick start
***********

Please check the README of `SimpleNavierStokesSolver <https://github.com/NaokiHori/SimpleNavierStokesSolver>`_.

****************
Acknowledgements
****************

The volume-of-fluid method is based on `THINC/QQ scheme <https://www.sciencedirect.com/science/article/pii/S0021999117305995?via%3Dihub>`_ with some modifications.

