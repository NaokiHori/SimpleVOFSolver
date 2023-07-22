################
Numerical method
################

It would be very easy if I could use the following naive second-order central-difference scheme to discretise the evolution of :math:`H`:

.. math::

   u_j \frac{\partial H}{\partial x_j}
   \approx
   \overline{u_j \frac{\delta H}{\delta x_j}}^j,

which is unconditionally unstable if treated explicitly.
Since the implicit treatment of the advective terms are too computationally-demanding, it is inevitable to adopt an upwind scheme for turbulent simulations.
This treatment, however, suffers from numerical diffusions which eliminates the sharp nature of the free surface.

One solution is to control the diffusivity of :math:`H` *locally* to maintain the thickness of the surface, which is categorised as a phase-field method.
The other solution is to reconstruct the surface every (or one in a few) time steps, which is categorised as a volume-of-fluid method.
In this project I adopt an intermediate solution (very roughly speaking) which is known as the `THINC scheme <https://onlinelibrary.wiley.com/doi/abs/10.1002/fld.975>`_, where I assume the surface is slightly diffused and also conduct the surface reconstruction.

To begin with, instead of discretising the advection equation directly, I integrate it in a control volume:

.. math::

   \int_V \der{H}{t} dV
   +
   \int_V u_j \der{H}{x_j} dV
   =
   0.

By assuming the temporal and the spatial treatments commutative, the first term leads to

.. math::

   \der{}{t} \int_V H dV.

The second term yields

.. math::

   \int_V u_j \der{H}{x_j} dV
   =
   \int_V u_j \der{H}{x_j} dV
   +
   \int_V H \der{u_j}{x_j} dV
   =
   \int_V \der{u_j H}{x_j} dV,

where the incompressibility constraint is used.

Thanks to the Gauss theorem, this is equal to

.. math::

   \int_{\partial V} u_j H n_j dS.

Note that :math:`n_j` is *not* the normal vector with respect to the free surface, but is the normal vector to the surface of the control volume.

Dividing the equation by the volume

.. math::

   \int_{V} dV,

yields

.. math::

   \der{\phi}{t}
   +
   \frac{1}{V} \int_{\partial V} u_j H n_j dS
   =
   0,

which is the central equation in this project.

Here, :math:`\phi` is called ``volume-of-fluid``, whose definition is explained by the name.
The rate of increase or decrease of this quantity :math:`\phi` is given by the second term, called the flux.
In order to numerically integrate this equation in time and to simulate the motion of the free surface, I need to describe the integrand :math:`u_j H n_j`, which will be extensively discussed in the following sections.

.. toctree::
   :maxdepth: 1

   normal
   intercept
   flux
   force

