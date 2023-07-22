#################
Simple VOF Solver
#################

|License|_ |CI|_ |DOC|_ |LastCommit|_

.. |License| image:: https://img.shields.io/github/license/NaokiHori/SimpleVOFSolver
.. _License: https://opensource.org/licenses/MIT

.. |CI| image:: https://github.com/NaokiHori/SimpleVOFSolver/actions/workflows/ci.yml/badge.svg
.. _CI: https://github.com/NaokiHori/SimpleVOFSolver/actions/workflows/ci.yml

.. |DOC| image:: https://github.com/NaokiHori/SimpleVOFSolver/actions/workflows/documentation.yml/badge.svg
.. _DOC: https://naokihori.github.io/SimpleVOFSolver/index.html

.. |LastCommit| image:: https://img.shields.io/github/last-commit/NaokiHori/SimpleVOFSolver/main
.. _LastCommit: https://github.com/NaokiHori/SimpleVOFSolver/commits/main

.. image:: https://github.com/NaokiHori/SimpleVOFSolver/blob/main/docs/source/thumbnail.gif
   :target: https://youtu.be/CpqQJxSkm7Q
   :width: 100%

.. image:: https://github.com/NaokiHori/SimpleVOFSolver/blob/main/docs/source/snapshot3d.png
   :width: 50%

********
Overview
********

This library numerically solves the motion of free-surface and deformable droplets in two- and three-dimensional Cartesian domains using the finite-difference and volume-of-fluid methods.
This is built on top of `SimpleNSSolver <https://github.com/NaokiHori/SimpleNSSolver>`_.

**********
Dependency
**********

Same dependency as the `SimpleNSSolver <https://github.com/NaokiHori/SimpleNSSolver>`_.

***********
Quick start
***********

#. Prepare workplace

   .. code-block:: console

      mkdir -p /path/to/your/directory
      cd       /path/to/your/directory

#. Get source

   .. code-block:: console

      git clone --recurse-submodules https://github.com/NaokiHori/SimpleVOFSolver
      cd SimpleVOFSolver

#. Set initial condition

   Here ``Python3`` is used to initialise the flow fields conveniently.
   One can give ``NPY`` files in different way under ``initial_condition/output/``.

   .. code-block:: console

      cd initial_condition
      make output
      bash exec.sh
      cd ..

#. Build solver

   .. code-block:: console

      make output
      make all

#. Run

   .. code-block:: console

      bash exec.sh

*************
Documentation
*************

Free-surface treatment is briefly documented `here <https://naokihori.github.io/SimpleVOFSolver>`_.
Please refer to the `documentation of SimpleNSSolver <https://naokihori.github.io/SimpleNSSolver>`_ for other details.

****************
Acknowledgements
****************

The volume-of-fluid method is based on `THINC/QQ scheme <https://www.sciencedirect.com/science/article/pii/S0021999117305995?via%3Dihub>`_ with some modifications.

