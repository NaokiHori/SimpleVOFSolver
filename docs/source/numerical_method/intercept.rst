
.. _intercept:

.. include:: /references.txt

##########################
Intercept and THINC scheme
##########################

********
Overview
********

In the :ref:`previous section <normal>`, the surface normal is computed at each cell center.
This information is, however, not enough to uniquely define the surface position uniquely:

.. figure:: images/intercept0.png
   :width: 100%

   Three different interfacial positions having the same surface normal.

In this image, for example, although all three surface patterns have the same surface normal, they are obviously different.

Since the three cases have the different occupations, I use the information to fix the surface position:

.. math::

   \frac{1}{V_{ij}} \int_{V_{ij}} H dV_{ij} = \phi_{ij},

which is the definition of :math:`\phi`.

.. note::

   In :ref:`the previous section <normal>`, I defined local coordinate systems such that the size of the computational cells :math:`V_{ij}` is :math:`1`.
   The origin of the local coordinates locates at the cell center, and as a result the cell edges exist at :math:`x = \pm \frac{1}{2}` and :math:`y = \pm \frac{1}{2}`.

Mathematically or physically speaking, since the free surface is discontinuous, :math:`H` is a step function.
Conventional VOF methods solve the equation based on this assumption.

In this project, I treat this problem in a different way here: instead of :math:`H` (step function), I consider a slightly smoothed function :math:`\hat{H}`:

.. math::

   H \left( n_X X + n_Y Y + d \right)
   & \approx
   \hat{H} \left( n_X X + n_Y Y + d \right) \\
   & \equiv
   \frac{1}{2} \left[ 1 + \tanh \left\{ \beta \left( n_X X + n_Y Y + d \right) \right\} \right].

.. figure:: images/intercept1.png
   :width: 75%

   Sharp (left, :math:`H \left( x, y \right)`) and smoothed (right, :math:`\hat{H} \left( x, y \right)`) interfacial representations.

:math:`d` is a parameter to be defined here, giving the surface position.
:math:`\beta` is a sharpness parameter which specifies the thickness, which corresponds to the Cahn number :math:`Cn` in the phase-field methods:

.. figure:: images/intercept2.png
   :width: 75%

   Interface sharpness as a function of the control parameter :math:`\beta`.
   :math:`\beta \rightarrow \infty` leads to :math:`H`.
   :math:`\Delta x = 1` corresponds to the size of the grid.

.. note::

   The larger the :math:`\beta` is, the sharper the surface becomes and the more realistic it is.
   However, there are two drawbacks for large :math:`\beta`:

   * Inaccurate surface representation

      As discussed below, the numerical integral of :math:`\hat{H}` is justified by the fact that :math:`\hat{H}` is smooth.
      Larger :math:`\beta` violates this assumption, and higher-order Gaussian quadratures are needed to accurately compute the surface position, which is computationally demanding.

   * Large spurious currents

      THINC schemes tend to suppress the spurious currents well compared to the conventional VOF method, which comes from the diffused nature of the surface.
      Larger :math:`\beta` diminishes this advantage.

Recall that here I want to decide :math:`d` to uniquely decide the position of the local surface.
Although the above equation has an analytical solution for two-dimensional cases, it contains Eulerian logarithmic integrals, which is computationally expensive.
Thanks to the diffused assumption of the surface, however, numerical integral gives a good approximation, namely:

.. math::

   \int_{V} \hat{H} \left( n_X X + n_Y Y + d \right) dV
   \approx
   \sum_{i}^N \sum_{j}^N w_i w_j \hat{H} \left( n_X X_i + n_Y Y_j + d \right)
   =
   \phi,

where :math:`w_i` and :math:`X_i` are the weights and the points of :math:`N`-th-order Gaussian quadrature defined in the interval :math:`\left[ -\frac{1}{2}, +\frac{1}{2} \right]`.
Notice that the volume integral is replaced with the summations in all dimensions.

To obtain :math:`d`, I need to minimise the following function :math:`f`:

.. math::

   f \left( d \right)
   \equiv
   -
   \phi
   +
   \sum_{i}^N \sum_{j}^N w_i w_j \hat{H} \left( n_X X_i + n_Y Y_j + d \right).

For simplicity, I introduce

.. math::

   \psi
   \equiv
   \psi \left( X, Y \right)
   \equiv
   n_X X + n_Y Y.

here.
Also I use an identity:

.. math::

   \hat{H} \left( x \right)
   =
   \frac{1}{2} \left\{ 1 + \tanh \left( x \right) \right\}
   \equiv
   \frac{1}{1 + \exp{\left( - 2 \beta x \right)}},

to yield

.. math::

   f \left( d \right)
   \equiv
   -
   \phi
   +
   \sum_{i}^N \sum_{j}^N \frac{w_i w_j}{
      1 + \exp{
        \left\{
          - 2 \beta \left( \psi + d \right)
        \right\}
      }
   }.

Since :math:`\exp{\left( x \right)}` is a monotonic function, finding :math:`d` is equivalent to finding the corresponding :math:`D \equiv \exp{\left(-2 \beta d \right)}`:

.. math::

   f \left( D \right)
   \equiv
   -
   \phi
   +
   \sum_{i}^N \sum_{j}^N \frac{w_i w_j}{
      1 + P_{ij} D
   },

and its derivative (with respect to :math:`D` instead of :math:`d`) is

.. math::

   f^{\prime} \left( D \right)
   \equiv
   \sum_{i}^N \sum_{j}^N \frac{- w_i w_j P_{ij}}{
      \left( 1 + P_{ij} D \right)^2
   },

where

.. math::

   P_{ij} \equiv \exp{\left( - 2 \beta \psi_{ij} \right)}.

*********************
Newton-Raphson method
*********************

Roughly speaking, to obtain :math:`d` (or :math:`D`), I move the hyper-surface (see below) and find :math:`d`, whose (numerical) integral coincides with the given local volume fraction :math:`\phi`.

.. figure:: images/intercept3.png
   :width: 75%

   Hyper-surfaces :math:`\hat{H} \left( x, y, d \right)` as a function of :math:`d`, while the normal vector is fixed.

which is achieved by `the Newton-Raphson method <https://en.wikipedia.org/wiki/Newton's_method>`_ following |XIE2017|.

Notice that :math:`P_{ij} = \exp{\left( -2\beta \psi_{ij} \right)}` is independent of :math:`D` and thus I compute them before starting the iteration:

.. myliteralinclude:: /../../src/interface/curvature_tensor.c
   :language: c
   :tag: compute constants a priori

Also, to start the iteration, I need the initial :math:`D_0`.
As a guess, I take the solution of the first-order Gaussian quadrature:

.. math::

   P_{00} = \exp{\left\{ \psi \left( X = 0, Y = 0 \right) \right\}} = 1

and

.. math::

   w_0 = 1,

giving

.. math::

   f \left( D_0 \right)
   =
   - \phi
   +
   \frac{1}{1 + D_0}
   =
   0

and thus

.. math::

   D_0 = \frac{1}{\phi} - 1:

.. myliteralinclude:: /../../src/interface/curvature_tensor.c
   :language: c
   :tag: initial guess

Also the function

.. math::

   f \left( D \right)
   \equiv
   -
   \phi
   +
   \sum_{i}^N \sum_{j}^N \frac{w_i w_j}{
      1 + P_{ij} D
   }

and its derivative

.. math::

   f^{\prime} \left( D \right)
   \equiv
   \sum_{i}^N \sum_{j}^N \frac{- w_i w_j P_{ij}}{
      \left( 1 + P_{ij} D \right)^2
   }

are computed here:

.. myliteralinclude:: /../../src/interface/curvature_tensor.c
   :language: c
   :tag: sum up

:math:`D` is updated using the computed :math:`f` and :math:`f^{\prime}`:

.. math::

   D_{n+1} \leftarrow D_{n} - \frac{f}{f^{\prime}},

which is implemented as

.. myliteralinclude:: /../../src/interface/curvature_tensor.c
   :language: c
   :tag: update D

This iteration continues until the residual :math:`| f |` becomes smaller than the given threshold or until the number of iterations exceeds the specified maximum:

.. myliteralinclude:: /../../src/interface/curvature_tensor.c
   :language: c
   :tag: Newton-Raphson method, loop terminating conditions

Finally I convert :math:`D` to :math:`d` following :math:`D = \exp{\left( -2\beta d \right)} \Leftrightarrow d = -\frac{1}{2\beta} \log{\left( D \right)}`:

.. myliteralinclude:: /../../src/interface/curvature_tensor.c
   :language: c
   :tag: convert D to d

In the :ref:`next step <flux>`, I use this uniquely-determined surface to evaluate the flux to integrate the advection equation.

