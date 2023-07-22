
.. _flux:

##############
Numerical flux
##############

In the :ref:`previous step <intercept>`, I have obtained the surface function

.. math::

   \psi \left( X, Y \right) + d
   \equiv
   n_X X
   +
   n_Y Y
   +
   d,

and its diffused representation

.. math::

   \hat{H} \left( X, Y \right)
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

In this section, I compute the numerical flux

.. math::

   \frac{1}{V_{ij}} \int_{\partial V} u_j \hat{H} n_j dS

to integrate the avection equation and to update :math:`\phi`.

Explicitly, on a two-dimensional Cartesian domain, this integral is consisted of four terms:

.. math::

   \int_{x, \pim} \ux \hat{H} dy
   -
   \int_{x, \pip} \ux \hat{H} dy
   +
   \int_{y, \pjm} \uy \hat{H} dx
   -
   \int_{y, \pjp} \uy \hat{H} dx

divided by the cell size

.. math::

   V_{ij}
   =
   \Delta x_{\pic} \Delta y.

Now I focus the first term, and the other terms can be evaluated in the same manner.

First of all, since :math:`\ux` is constant on the local cell face, I can take it out of the integral:

.. math::

   \frac{1}{\Delta x_i \Delta y} \ux \int_{x, \pim} \hat{H} dy.

To stabilise the cell-face integral, I need to use the upward scheme: namely information of the upwind cell should be used, i.e.

.. math::

   \vat{\hat{H}}{\pimm, \pjc} \,\, & \text{if} \,\, \vat{\ux}{\pim, \pjc} \ge 0, \\
   \vat{\hat{H}}{\pic,  \pjc} \,\, & \text{if} \,\, \vat{\ux}{\pim, \pjc} <   0.

.. figure:: images/flux0.png
   :width: 75%

   :math:`\vat{\ux}{\pim, \pjc} \ge 0`

.. figure:: images/flux1.png
   :width: 75%

   :math:`\vat{\ux}{\pim, \pjc}  <  0`

Recall that :math:`\hat{H}` is defined on the local coordinate system attached to each cell center.
Thus, for the left cell :math:`\left( \pimm, \pjc \right)`, the integrand should be evaluated at :math:`x = \frac{1}{2}`, while :math:`x = -\frac{1}{2}` for the right cell :math:`\left( \pic, \pjc \right)`.

Same applied to the all fluxes.

.. myliteralinclude:: /../../src/interface/update/flxx.c
   :language: c
   :tag: use upwind information

.. myliteralinclude:: /../../src/interface/update/flxy.c
   :language: c
   :tag: use upwind information

.. myliteralinclude:: /../../src/interface/update/flxz.c
   :language: c
   :tag: use upwind information

Also notice that

.. math::

   dy
   =
   \der{y}{Y} dY
   \approx
   \Delta y dY,

and thus

.. math::

   \frac{1}{\Delta x_i \Delta y} \ux \int_{x, \pim} \hat{H} dy
   =
   \begin{cases}
      \ux \ge 0 & \frac{1}{\Delta x_i} \ux \int_{x = +\frac{1}{2}} \hat{H} dY_{i-1, j}, \\
      \ux  <  0 & \frac{1}{\Delta x_i} \ux \int_{x = -\frac{1}{2}} \hat{H} dY_{i  , j},
   \end{cases}

whose integral is again approximated by the :math:`N`-th-order Gaussian quadrature:

.. math::

   \frac{1}{\Delta x_i} \ux \sum_{j}^N w_j \hat{H} \left( x, y_j \right),

which is implemented here:

.. myliteralinclude:: /../../src/interface/update/flxx.c
   :language: c
   :tag: evaluate flux

.. myliteralinclude:: /../../src/interface/update/flxy.c
   :language: c
   :tag: evaluate flux

.. myliteralinclude:: /../../src/interface/update/flxz.c
   :language: c
   :tag: evaluate flux

.. note::

   For very small (:math:`\approx 0`) or very large (:math:`\approx 1`) volume fractions, i.e. single-phase regions,

   .. math::

      \hat{H} \approx \phi = const.

   gives a good approximation and thus is directly used.

All the fluxes evaluated at the cell faces are stored in ``FLXX`` and ``FLXY``, which are used to update :math:`\phi`:

.. myliteralinclude:: /../../src/interface/update/main.c
   :language: c
   :tag: compute right-hand-side of advection equation, x flux

.. myliteralinclude:: /../../src/interface/update/main.c
   :language: c
   :tag: compute right-hand-side of advection equation, y flux

.. myliteralinclude:: /../../src/interface/update/main.c
   :language: c
   :tag: compute right-hand-side of advection equation, z flux

A conventional three-step Runge-Kutta scheme is adopted to integrate the equation in time.

.. myliteralinclude:: /../../src/interface/update/main.c
   :language: c
   :tag: update vof, alpha contribution

.. myliteralinclude:: /../../src/interface/update/main.c
   :language: c
   :tag: update vof, beta contribution

