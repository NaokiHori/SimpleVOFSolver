
.. _normal:

.. include:: /references.txt

##############
Surface normal
##############

I assume that I know the profile of :math:`\phi` at :math:`n` step:

.. image:: images/normal0.png
   :width: 100%

Note that the left-bottom and right-top cells are fully occupied by one of the two phases, while other cells are *shared* by the two phases.

Interface reconstruction starts by computing the local surface normal.
So a tentative goal here is know the local normal vector :math:`n_{\pic,\pjc}` at each cell center.
Note that I need to evaluate it from the surrounding discrete :math:`\phi` values.

There are lots of possibilities, and each way has different accuracy and computational expense.
Here, I adopt one of the most simple approach using the fact that the normal can be described by

.. math::

   n_i
   =
   \frac{m_i}{\left| m_i \right|},

where :math:`m_i` is the local gradient approximated as

.. math::

   m_i
   \approx
   \der{\phi}{x_i}.

To compute this quantity at the cell center :math:`\left( \pic, \pjc \right)`, one of the simplest ways would be to use second-order accurate central-difference scheme:

.. math::

   \der{\phi}{x} \approx \frac{\phi_{\pipp, \pjc} - \phi_{\pimm, \pjc}}{2 \Delta x},

which tends to be unstable since this does not use :math:`\phi_{\pic, \pjc}` (i.e., not central-dominant).

Instead I adopt Youngs approach (|YOUNGS1982|):

#. Compute gradients at cell corners

   I use 4 (or 9 in 3D) surrounding cells to compute gradients at cell corners, i.e.,

   Gradients in :math:`x`:

   .. math::

      \vat{m_x}{\pim, \pjm} & \approx \vat{\der{\phi}{x}}{\pim, \pjm}, \\
      \vat{m_x}{\pip, \pjm} & \approx \vat{\der{\phi}{x}}{\pip, \pjm}, \\
      \vat{m_x}{\pim, \pjp} & \approx \vat{\der{\phi}{x}}{\pim, \pjp}, \\
      \vat{m_x}{\pip, \pjp} & \approx \vat{\der{\phi}{x}}{\pip, \pjp}.

   .. myliteralinclude:: /../../src/interface/tensor.c
      :language: c
      :tag: x gradient

   Gradients in :math:`y`:

   .. math::

      \vat{m_y}{\pim, \pjm} & \approx \vat{\der{\phi}{y}}{\pim, \pjm}, \\
      \vat{m_y}{\pip, \pjm} & \approx \vat{\der{\phi}{y}}{\pip, \pjm}, \\
      \vat{m_y}{\pim, \pjp} & \approx \vat{\der{\phi}{y}}{\pim, \pjp}, \\
      \vat{m_y}{\pip, \pjp} & \approx \vat{\der{\phi}{y}}{\pip, \pjp}.

   .. myliteralinclude:: /../../src/interface/tensor.c
      :language: c
      :tag: y gradient

   Gradients in :math:`z`:

   .. myliteralinclude:: /../../src/interface/tensor.c
      :language: c
      :tag: z gradient

   They are gradients and not normals.
   I here normalise to obtain normal, which are assigned to the resulting array ``DVOF``.

   .. myliteralinclude:: /../../src/interface/tensor.c
      :language: c
      :tag: normalise and obtain corner normals

   .. image:: images/normal1.png
      :width: 100%

   .. note::

      In this image, lengths of arrows are adjusted for the sake of appearances and thus may be inaccurate.

#. Compute normals at cell centers

   Normals at cell centers are computed by averaging the normal of the surrouding cell corners:

   .. myliteralinclude:: /../../src/interface/tensor.c
      :language: c
      :tag: average nx

   .. myliteralinclude:: /../../src/interface/tensor.c
      :language: c
      :tag: average ny

   .. myliteralinclude:: /../../src/interface/tensor.c
      :language: c
      :tag: average nz

   Notice that, because of the averaging, the norm of this central vector might not be unity.

   .. image:: images/normal2.png
      :width: 100%

   .. note::

      In this image, lengths of arrows are adjusted for the sake of appearances and thus may be inaccurate.

#. Coordinate transformation

   These normal vectors at cell centers are defined on the *global* coordinate system (e.g., coordinate system used to describe domain), whose basis is :math:`\left( \underline{e}_x, \underline{e}_y \right)`.
   Later on, it is convenient to consider thing on *local* coordinate system, which is a coordinate system attached to each cell :math:`\left( i, j \right)`, whose basis is :math:`\left( \underline{e}_{X}, \underline{e}_{Y} \right)`.
   These two coordinate systems are tied with the relation:

   .. math::

      \begin{pmatrix}
         \underline{e}_X \\
         \underline{e}_Y
      \end{pmatrix}
      =
      \begin{pmatrix}
         \Delta x_i & 0        \\
                  0 & \Delta y
      \end{pmatrix}
      \begin{pmatrix}
         \underline{e}_x \\
         \underline{e}_y
      \end{pmatrix},

   or apparently

   .. math::

      \begin{pmatrix}
         \underline{e}_x \\
         \underline{e}_y
      \end{pmatrix}
      =
      \begin{pmatrix}
         \frac{1}{\Delta x_i} & 0                 \\
                            0 & \frac{1}{\Delta y}
      \end{pmatrix}
      \begin{pmatrix}
         \underline{e}_X \\
         \underline{e}_Y
      \end{pmatrix}.

   .. image:: images/normal3.png
      :width: 100%

   .. note::

      In this image, lengths of arrows are adjusted for the sake of appearances and thus may be inaccurate.

   Based on these formula, the above normal vector attached to the *global* coordinate system

   .. math::

      n_i^{global}
      =
      \underline{e}_x n_x
      +
      \underline{e}_y n_y

   is written as

   .. math::

      \underline{e}_X \frac{1}{\Delta x_i} n_x
      +
      \underline{e}_Y \frac{1}{\Delta y  } n_y

   on the *local* coordinate system.

   Obviously the length of this vector is not unity, which should be rescaled so that it is normalised with respect to the *local* coordinate system:

   .. math::

      n_i^{local}
      \equiv
      \underline{e}_X n_X
      +
      \underline{e}_Y n_Y,

   which is

   .. math::

      n_X &= \frac{1}{k} \frac{n_x}{\Delta x_{i}}, \\
      n_Y &= \frac{1}{k} \frac{n_y}{\Delta y    },

   where

   .. math::

      k
      \equiv
      \sqrt{
         \left( \frac{n_x}{\Delta x_{i}} \right)^2
         +
         \left( \frac{n_y}{\Delta y    } \right)^2
      }.

   whose results are assigned to the resulting array ``NORMAL``.

   .. myliteralinclude:: /../../src/interface/tensor.c
      :language: c
      :tag: normalise and obtain center normals

Now, on the new coordinate system, each computational grid is a unit square (or a unit cube in 3D), which greatly simplifies :ref:`the following discussion <intercept>`.

