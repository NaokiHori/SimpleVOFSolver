
.. _intercept:

##########################
Intercept and THINC scheme
##########################

Now I have obtained local surface normals.
This information is, however, insufficient to decide surface uniquely since translational freedom still exists.
In order to *fix* the surface position, I use the relation

.. math::

   \frac{1}{V_{ij}} \int_{V_{ij}} H dV_{ij} = \phi_{ij},

namely, the volume integral of the indicator function :math:`H` gives the local volume fraction :math:`\phi_{ij}`, which is nothing else but the definition of :math:`\phi`.

.. note::

   In :ref:`the previous step <normal>`, I defined a local coordinate system which attaches to each computational cell.
   Thus the size of the computational cell size :math:`V_{ij}` is :math:`1`.
   Also I consider that, for each cell, the cell center locates at the origin and the cell edges locate at :math:`\pm \frac{1}{2}`.

Mathematically speaking, :math:`H` is a step function since the interface is discontinuous from the continuum point of view.
However, I introduce a big approximation here: instead of considering a discontinuous step function, I decide to use a continuous approximation:

.. math::

   H \left( n_X X + n_Y Y + d \right)
   & \approx
   \hat{H} \left( n_X X + n_Y Y + d \right) \\
   & \equiv
   \frac{1}{2} \left[ 1 + \tanh \left\{ \beta \left( n_X X + n_Y Y + d \right) \right\} \right].

Here :math:`\beta` is a sharpness parameter specifying the thickness of the diffusion, whose corresponding concept is Cahn number :math:`Cn` in phase-field method.

So my task here is to decide :math:`d`, which uniquely decides the position of surface.
Although the above equation does have an analytical solution, it uses Eulerian logarithmic integral, which is computationally expensive.
Thanks to the diffused assumption of the surface, however, numerical integral gives a good approximation, namely:

.. math::

   \int_{V} \hat{H} \left( n_X X + n_Y Y + d \right) dV
   \approx
   \sum_{i}^N \sum_{j}^N w_i w_j \hat{H} \left( n_X X_i + n_Y Y_j + d \right)
   =
   \phi,

where :math:`w_i` and :math:`X_i` are weights and points of :math:`N`-th-order Gaussian quadrature defined in :math:`\left[ -\frac{1}{2}, +\frac{1}{2} \right]`.
Notice that the volume integral is replaced with summations in each dimension.

In order to obtain :math:`d`, I need to minimise the residual of

.. math::

   f \left( d \right)
   \equiv
   -
   \phi
   +
   \sum_{i}^N \sum_{j}^N w_i w_j \hat{H} \left( n_X X_i + n_Y Y_j + d \right),

which is achieved by the Newton-Raphson iteration.

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

Now it's time to iterate to find :math:`D`.

Notice that :math:`P_{ij} = \exp{\left( -2\beta \psi_{ij} \right)}` is only functions of Gaussian quadrature and surface normal, and independent of :math:`D`.
Thus I can compute them a priori:

.. myliteralinclude:: /../../src/interface/tensor.c
   :language: c
   :tag: compute constants a priori

In order to iterate, I need initial :math:`D_0`.
As a initial guess, I adopt a solution of :math:`0`-th-order Gaussian quadrature; namely, from

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

Finally I convert :math:`D` to :math:`d` following :math:`D = \exp{\left( -2\beta d \right)}` or :math:`d = -\frac{1}{2\beta} \log{\left( D \right)}`:

.. myliteralinclude:: /../../src/interface/tensor.c
   :language: c
   :tag: convert D to d

The resulting normal vector and the intercept are stored to ``NORMAL``:

.. myliteralinclude:: /../../src/interface/tensor.c
   :language: c
   :tag: store normal and intercept

Now I have uniquely determined the location of (diffused) surface, whose information will be used in :ref:`the next step <flux>` to compute numerical fluxes in each direction.

.. note::

   The larger :math:`\beta` is, the sharper the surface becomes and closer to the reality.
   However, there are two drawbacks for large :math:`\beta`:

   * Inaccurate surface representation

      As being discussed above, the numerical integral of :math:`\hat{H}` is justified by the fact that :math:`\hat{H}` is smooth.
      Larger :math:`\beta` violates this assumption, and higher-order Gaussian quadratures are needed to accurately compute the surface position, which is computationally demanding.

   * Large spurious currents

      THINC schemes tend to suppress the spurious current more than conventional VOF method, which comes from the diffused nature of the surface.
      Larger :math:`\beta` would diminish this advantage and thus the spurious currents tend to get larger.

