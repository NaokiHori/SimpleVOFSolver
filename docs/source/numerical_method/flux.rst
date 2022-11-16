
.. _flux:

################
Numerical fluxes
################

In :ref:`the previous step <intercept>`, I have computed the surface normal and its intercept to uniquely determine the surface function:

.. math::

   \psi \left( X, Y \right) + d
   \equiv
   n_X X
   +
   n_Y Y
   +
   d

and its diffused representation:

.. math::

   H \left( X, Y \right)
   \equiv
   \frac{1}{
      1 + \exp{
        \left\{
          -2 \beta
          \left(
            \psi + d
          \right)
        \right\}
      }
   }

for each cell center :math:`\left( \pic, \pjc \right)`.

Now I am ready for computing the numerical fluxes

.. math::

   \frac{1}{V_{ij}} \int_{\partial V} u_j H n_j dS

to update :math:`\phi`.
(Recall that my final objective is to integrate the advection equation of :math:`\phi` in time.)

Explicitly, on a two-dimensional Cartesian domain, this integral can be written as

.. math::

   \int_{x, \pim} \ux H dy
   -
   \int_{x, \pip} \ux H dy
   +
   \int_{y, \pjm} \uy H dx
   -
   \int_{y, \pjp} \uy H dx

divided by the cell size

.. math::

   V_{ij}
   =
   \Delta x_{\pic} \Delta y.

Now I only focus on :math:`\frac{1}{V_{ij}} \int_{x_{\pim}} \ux H dy`, since other terms can be evaluated in the same manner.

First of all, in the current framework, :math:`\ux` is constant on the local cell face, and thus I can take it out of the integral:

.. math::

   \frac{1}{\Delta x_i \Delta y} \ux \int_{x, \pim} H dy.

The integral is evaluated at the cell face, in which information of :math:`H` is needed.
Because of a stability limitation, I need to use information from upward, i.e.,

.. math::

   \vat{H}{\pimm, \pjc} \,\, & \text{if} \,\, \vat{u}{\pim, \pjc} \ge 0, \\
   \vat{H}{\pic,  \pjc} \,\, & \text{if} \,\, \vat{u}{\pim, \pjc} <   0.

Same holds for the other fluxes, which are implemented as

.. myliteralinclude:: /../../src/interface/update.c
   :language: c
   :tag: use upwind information

Recall that :math:`H` is defined on local coordinate system attached to each cell center.
Thus, for the left cell :math:`\left( \pimm, \pjc \right)`, the integrand should be evaluated at :math:`x = \frac{1}{2}`, while :math:`x = -\frac{1}{2}` for the right cell :math:`\left( \pic, \pjc \right)`.
Also notice that

.. math::

   dy
   =
   \der{y}{Y} dY
   \approx
   \Delta y dY,

and thus

.. math::

   \frac{1}{\Delta x_i} \ux \int_{x, \pim} H dY,

whose integral is again approximated by :math:`N`-th-order Gaussian quadrature:

.. math::

   \frac{1}{\Delta x_i} \ux \sum_{j}^N w_j H \left( x, y_j \right),

which is implemented here:

.. myliteralinclude:: /../../src/interface/update.c
   :language: c
   :tag: evaluate flux

All fluxes evaluated at cell faces are stored in ``VOFFLUXX``, ``VOFFLUXY``, and ``VOFFLUXZ``, which are used to update :math:`\phi`:

.. myliteralinclude:: /../../src/interface/update.c
   :language: c
   :tag: compute right-hand-side of advection equation

Notice that fluxes are normalised by the pre-factors here for simplicity.

A conventional three-step Runge-Kutta scheme is adopted to integrate the equation in time.

