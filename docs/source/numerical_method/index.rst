################
Numerical method
################

It would be very easy if I could use naive second-order central difference scheme to discretise the evolution of :math:`H`:

.. math::

   u_j \frac{\partial H}{\partial x_j}
   \approx
   \overline{u_j \frac{\delta H}{\delta x_j}}^j,

which is, as is well-known, impractical because of numerical instability.
Thus it is inevitable to use up-wind schemes in general.
This treatment, however, suffers from numerical diffusions, which tend to eliminate the sharp nature of the surface.

One solution is to control the diffusivity of :math:`H` locally to keep the thickness of the surface, which is categorised as phase-field methods.
The other solution is to reconstruct the surface every (or one in a few) time steps, which is categorised as volume-of-fluid methods.

Here I adopt an intermediate solution (very roughly speaking) which is named as `THINC <https://onlinelibrary.wiley.com/doi/abs/10.1002/fld.975>`_ scheme, in which I allow the surface to have a finite thickness, and also conduct surface reconstruction.

First of all, instead of discretising the above equation directly, I integrate it in a control volume:

.. math::

   \int_V \der{H}{t} dV
   +
   \int_V u_j \der{H}{x_j} dV
   =
   0.

By assuming the temporal and spatial treatments commutable, the first term leads

.. math::

   \der{}{t} \int_V H dV.

By adopting the incompressibility constraint, the second terms leads

.. math::

   \int_V u_j \der{H}{x_j} dV
   =
   \int_V u_j \der{H}{x_j} dV
   +
   \int_V H \der{u_j}{x_j} dV
   =
   \int_V \der{u_j H}{x_j} dV.

Thanks to the Gauss theorem, this yields

.. math::

   \int_{\partial V} u_j H n_j dS.

Note that :math:`n_j` is *not* the surface normal of the free surface, but surface normal to the surface of the control volume.

Dividing the equation by the volume

.. math::

   \int_{V} dV,

gives me

.. math::

   \der{\phi}{t}
   +
   \frac{1}{V} \int_{\partial V} u_j H n_j dS
   =
   0,

which is the equation playing a key role in this project.

Here, :math:`\phi` is called ``volume-of-fluid``, whose definition is explained by the name.
The amount of increase or decrease of this quantity :math:`\phi` is given by the second term, which is called flux.
In order to numerically integrate this equation in time and simulate the evolution of free surface, I need to describe :math:`H` and the integrand, which will be extensively discussed below, as well as the implementations.

.. toctree::
   :maxdepth: 2

   normal
   intercept
   flux
   force

