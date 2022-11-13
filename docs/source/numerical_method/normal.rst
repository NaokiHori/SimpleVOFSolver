
.. _normal:

##############
Surface normal
##############

I assume that I know the profile of :math:`\phi` at :math:`n` step:

+---+---+---+---+
|   |i-1| i |i+1|
+---+---+---+---+
|j+1|0.1|0.2|0.3|
+---+---+---+---+
| j |0.2|0.4|0.6|
+---+---+---+---+
|j-1|0.4|0.7|1.0|
+---+---+---+---+

Here this ugly table denotes the profile of :math:`\phi`, where e.g., the left-top, center, left-bottom value represents :math:`\phi_{i-1, j+1}`, :math:`\phi_{i, j}`, and :math:`\phi_{i-1, j-1}`, respectively.

Note that the right-bottom cell is fully occupied by one phase, while other cells are *shared* by the two phases.

Interface reconstruction starts by computing the local surface normal.
So a tentative goal here is know the local normal vector :math:`n_{\pic,\pjc}`.
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

To compute this quantity at the cell center :math:`\left( \pic, \pjc \right)`, one of the simplest way is to use second-order accurate central difference:

.. math::

   \der{\phi}{x} \approx \frac{\phi_{\pipp, \pjc} - \phi_{\pimm, \pjc}}{2 \Delta x},

which tends to be unstable since this does not use :math:`\phi_{\pic, \pjc}` (i.e., not central-diminant).

Instead I adopt Youngs approach (TODO: ref):

#. Compute gradients at cell corners

   I use 4 (or 9 in 3D) surrounding cells to compute gradients at cell corners.

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

   Since averaging operations might break the unitary nature, I normalise them again here.

   These normal vectors are defined on the *global* coordinate system (physical coordinate system), whose basis is :math:`\left( \underline{e}_x, \underline{e}_y \right)`.
   Later on, it is convenient to consider thing on *local* coordinate system, which is a coordinate system attached to the cell :math:`\left( i, j \right)`, whose basis is :math:`\left( \underline{e}_{X}, \underline{e}_{Y} \right)`.
   These two coordinate systems are tied with

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

   or equivalently

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

   Based on these relation, the above normal vector attached to the *global* coordinate system

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
      \underline{e}_Y \frac{1}{\Delta y  } n_y.

   on the *local* coordinate system.
   Thus I can conclude that the surface normal on the local coordinate system

   .. math::

      n_i^{local}
      \equiv
      \underline{e}_X n_X
      +
      \underline{e}_Y n_Y

   leads

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

Now, on the new coordinate system, the cell sizes in all directions are :math:`1`, which greatly simplifies :ref:`the following discussion <intercept>`.

