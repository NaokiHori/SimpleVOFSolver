#################
Simple VOF Solver
#################

THIS IS AN ON-GOING PROJECT.

|License|_ |CI|_ |DOC|_ |LastCommit|_

.. |License| image:: https://img.shields.io/github/license/NaokiHori/SimpleVOFSolver
.. _License: https://opensource.org/licenses/MIT

.. |CI| image:: https://github.com/NaokiHori/SimpleVOFSolver/actions/workflows/ci.yml/badge.svg
.. _CI: https://github.com/NaokiHori/SimpleVOFSolver/actions/workflows/ci.yml

.. |DOC| image:: https://github.com/NaokiHori/SimpleVOFSolver/actions/workflows/documentation.yml/badge.svg
.. _DOC: https://naokihori.github.io/SimpleVOFSolver/index.html

.. |LastCommit| image:: https://img.shields.io/github/last-commit/NaokiHori/SimpleVOFSolver/main
.. _LastCommit: https://github.com/NaokiHori/SimpleVOFSolver/commits/main

.. image:: https://github.com/NaokiHori/SimpleVOFSolver/blob/main/docs/source/snapshot2d.png
   :width: 100%

.. image:: https://github.com/NaokiHori/SimpleVOFSolver/blob/main/docs/source/snapshot3d.png
   :width: 75%

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

***********
Quick start
***********

I guess people who are interested in this project already have installed the above dependencies.
(If not, please check the README of `SimpleNSSolver <https://github.com/NaokiHori/SimpleNavierStokesSolver>`_.)

Just in case, Docker image is available:

.. code-block::

   $ docker pull naokihori/snss:latest
   $ docker run -it --rm --cpuset-cpus="0-1" -u runner -v ${PWD}:/home/runner naokihori/snss:latest

By default, two-dimensional version is enabled.
If you are interested in 3D cases, please replace ``NDIMS=2`` with ``NDIMS=3`` in ``Makefile``.

.. code-block::

   $ make clean  # do not forget to do this after NDIMS is changed
   $ make output # just in case, to create directories for artifacts
   $ make all    # compile
   $ sh exec.sh  # change parameters if needed

For three-dimensional cases, you need to give ``glksize`` and ``lz`` as additional environmental variables (see ``exec.sh``).

*************
Documentation
*************

Documentation can be found `here <https://naokihori.github.io/SimpleVOFSolver>`_ (still on half way).
Please refer to the `documentation of SimpleNSSolver <https://naokihori.github.io/SimpleNavierStokesSolver>`_ for other details.

****************
Acknowledgements
****************

The volume-of-fluid method is based on `THINC/QQ scheme <https://www.sciencedirect.com/science/article/pii/S0021999117305995?via%3Dihub>`_ with some modifications.

