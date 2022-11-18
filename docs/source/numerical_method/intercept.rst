
.. _intercept:

.. include:: /references.txt

##########################
Intercept and THINC scheme
##########################

Now I have obtained local surface normals.
This information is, however, not enough to decide the piecewise interfacial position uniquely, since a translational freedom still exists:

.. figure:: images/intercept0.png
   :width: 100%

   Three different interfacial positions having the same surface normal.

In this image, for example, although all three surface patterns have the same surface normal, they are obviously different from case to case.

Notice that these cases have different phase volumes: In order to *fix* the surface position, the definition of :math:`\phi`:

.. math::

   \frac{1}{V_{ij}} \int_{V_{ij}} H dV_{ij} = \phi_{ij},

i.e., the volume integral of the indicator function :math:`H` gives the local volume fraction :math:`\phi_{ij}`, can be used.
Since :math:`H` is a function of the intercept (say :math:`d`), I can uniquely determine the surface position by solving this integral equation with respect to :math:`d`.

.. note::

   In :ref:`the previous step <normal>`, I defined a local coordinate system which attaches to each computational cell.
   Thus the size of the computational cell size :math:`V_{ij}` is :math:`1`.
   Also I consider that, for each cell, the origin of the local coordinate system locates at the cell center, and the cell edges locate at :math:`x = \pm \frac{1}{2}` and :math:`y = \pm \frac{1}{2}`.

Mathematically speaking, :math:`H` is a step function since I regard the interface as a discontinuity.
Also, it is possible to solve the above equation for this sharp surface, which is what the conventional VOF methods do.

However, I treat this problem differently here.
Namely, instead of handling a discontinuous step function :math:`H`, I decide to consider a bit smoothened version :math:`\hat{H}`:

.. math::

   H \left( n_X X + n_Y Y + d \right)
   & \approx
   \hat{H} \left( n_X X + n_Y Y + d \right) \\
   & \equiv
   \frac{1}{2} \left[ 1 + \tanh \left\{ \beta \left( n_X X + n_Y Y + d \right) \right\} \right].

.. figure:: images/intercept1.png
   :width: 100%

   Sharp (left, :math:`H \left( x, y \right)`) and smoothed (right, :math:`\hat{H} \left( x, y \right)`) interfacial representations.

Here :math:`\beta` is a sharpness parameter specifying the thickness of the diffusion, whose corresponding concept is Cahn number :math:`Cn` in phase-field methods, and it is unfortunately dependent on the user's choice.

.. figure:: images/intercept2.png
   :width: 100%

   Interface sharpness as a function of a control parameter :math:`\beta`, which are distinguished by different coloured lines.
   :math:`\beta \rightarrow \infty` corresponds to the analytical discontinuous step function.

Note that :math:`\Delta x = 1` corresponds to the size of the grid.

.. note::

   The larger :math:`\beta` is, the sharper the surface becomes and the closer to the reality.
   However, there are two drawbacks for large :math:`\beta`:

   * Inaccurate surface representation

      As being discussed above, the numerical integral of :math:`\hat{H}` is justified by the fact that :math:`\hat{H}` is smooth.
      Larger :math:`\beta` violates this assumption, and higher-order Gaussian quadratures are needed to accurately compute the surface position, which is computationally demanding.

   * Large spurious currents

      THINC schemes tend to suppress the spurious current more than conventional VOF method, which comes from the diffused nature of the surface.
      Larger :math:`\beta` would diminish this advantage and thus the spurious currents tend to get larger.

Recall that my task is to decide :math:`d` to uniquely decide the position of the local surface.
Although the above equation has an analytical solution for two-dimensional case, it needs Eulerian logarithmic integral which is computationally too expensive and thus is impractical.
Thanks to the diffused assumption of the surface, however, numerical integral gives a good approximation, namely:

.. math::

   \int_{V} \hat{H} \left( n_X X + n_Y Y + d \right) dV
   \approx
   \sum_{i}^N \sum_{j}^N w_i w_j \hat{H} \left( n_X X_i + n_Y Y_j + d \right)
   =
   \phi,

where :math:`w_i` and :math:`X_i` are weights and points of :math:`N`-th-order Gaussian quadrature defined in the interval of :math:`\left[ -\frac{1}{2}, +\frac{1}{2} \right]`.
Notice that the volume integral is replaced with summations in each dimension.

In order to obtain :math:`d`, I need to minimise the residual of

.. math::

   f \left( d \right)
   \equiv
   -
   \phi
   +
   \sum_{i}^N \sum_{j}^N w_i w_j \hat{H} \left( n_X X_i + n_Y Y_j + d \right),

which is achieved by `the Newton-Raphson iteration <https://en.wikipedia.org/wiki/Newton's_method>`_ following |XIE2017|.

So roughly speaking, what I will do is to move the hyper-surface (see below) and to find :math:`d`, whose (numerical) integral coincides with the given local volume fraction :math:`\phi`.

.. figure:: images/intercept3.png
   :width: 100%

   :math:`\hat{H} \left( x, y \right)` for various :math:`d`

Before explaining the details, I introduce

.. math::

   \psi
   \equiv
   \psi \left( X, Y \right)
   \equiv
   n_X X + n_Y Y

for notational simplicity.
Also I use an identity:

.. math::

   \hat{H} \left( x \right)
   =
   \frac{1}{2} \left\{ 1 + \tanh \left( x \right) \right\}
   \equiv
   \frac{1}{1 + \exp{\left( - 2 \beta x \right)}}.

Thus, the function whose residual should be minimised leads

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

Since :math:`\exp{\left( x \right)}` is a monotonic function, finding :math:`d` which minimises :math:`f` is equivalent to finding the corresponding :math:`\exp{\left(-2\beta d \right)}`.

Thus, I am able to make this a bit simpler

.. math::

   f \left( D \right)
   \equiv
   -
   \phi
   +
   \sum_{i}^N \sum_{j}^N \frac{w_i w_j}{
      1 + P_{ij} D
   },

and its derivative (with respect to :math:`D` instead of :math:`d`) leads

.. math::

   f^{\prime} \left( D \right)
   \equiv
   \sum_{i}^N \sum_{j}^N \frac{- w_i w_j P_{ij}}{
      \left( 1 + P_{ij} D \right)^2
   },

where

.. math::

   & P_{ij} \equiv \exp{\left( - 2 \beta \psi_{ij} \right)}, \\
   & D      \equiv \exp{\left( - 2 \beta d         \right)}.

Now it's time to iterate to find :math:`D` which satisfies :math:`f \left( D \right) = 0`.

Notice that :math:`P_{ij} = \exp{\left( -2\beta \psi_{ij} \right)}` is only functions of Gaussian quadrature and surface normal, and independent of :math:`D`.
Thus I can compute them a priori:

.. myliteralinclude:: /../../src/interface/tensor.c
   :language: c
   :tag: compute constants a priori

In order to iterate, I need initial :math:`D_0`.
As a initial guess, I adopt a solution of the first-order Gaussian quadrature; namely, from

.. math::

   P_{00} = \exp{\left\{ \psi \left( X = 0, Y = 0 \right) \right\}} = 1

and

.. math::

   w_0 = 1,

I have

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

.. myliteralinclude:: /../../src/interface/tensor.c
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

.. myliteralinclude:: /../../src/interface/tensor.c
   :language: c
   :tag: sum up

:math:`D` is updated using the computed :math:`f` and :math:`f^{\prime}`:

.. math::

   D_{n+1} \leftarrow D_{n} - \frac{f}{f^{\prime}},

which is implemented as

.. myliteralinclude:: /../../src/interface/tensor.c
   :language: c
   :tag: update D

This iteration continues until the residual :math:`| f |` becomes smaller than the threshold or the number of iterations exceeds the maximum value:

.. myliteralinclude:: /../../src/interface/tensor.c
   :language: c
   :tag: Newton-Raphson method, loop terminating conditions

Finally I convert :math:`D` to :math:`d` following :math:`D = \exp{\left( -2\beta d \right)} \Leftrightarrow d = -\frac{1}{2\beta} \log{\left( D \right)}`:

.. myliteralinclude:: /../../src/interface/tensor.c
   :language: c
   :tag: convert D to d

The resulting normal vector and the intercept are stored to ``NORMAL``:

.. myliteralinclude:: /../../src/interface/tensor.c
   :language: c
   :tag: store normal and intercept

Now I have uniquely determined the location of the piecewise (diffused) surface, whose information will be used in :ref:`the next step <flux>` to compute numerical fluxes in each direction.

