**NUMA-TF (Numerical Analysis of Trusses and Frames)**
 
 Welcome to the NUMA-TF program, an open-source MATLAB program to perform
 numerical analyses of reticulated structural models.
 
 GitLab repository:
 https://gitlab.com/rafaelrangel/numa-tf
 
 Matlab File Exchange with release versions:
 https://mathworks.com/matlabcentral/fileexchange/79406-numa-tf
 
 Types of analysis:
 * Linear-elastic.
 * Geometrically nonlinear.
 * Stability (buckling modes and critical loads).
 * Free vibration.
 * Natural vibration (vibration modes and natural frequencies).
 
 Structural model options:
 * 2D and 3D trusses and frames.
 * Euler-Bernoulli and Timoshenko theories for bending behavior.
 * Concentrated nodal loads, uniform/linear distributed loads, and thermal loads.

 Run the main.m file to select the input file with model and analysis information
 and perform the analysis.

 If the input file is valid, the program can print the analysis feedback in
 the MATLAB command window, write an analysis report file with a .pos
 extension in the same directory of the input file, and plot the requested
 results in MATLAB figure windows.

 If the input file is not valid, the progam will send a warning or an error
 message.
 
 In case of bad behavior of the program or interest in collaboration,
 please contact the author.

--------------------------------------------------------------------------------
**Authorship**

 The program (originally called NLframe2D) was created for the MSc thesis
 "Educational Tool for Structural Analysis of Plane Frame Models with Geometric
 Nonlinearity", 2019, avaiable in: 
 webserver2.tecgraf.puc-rio.br/~lfm/teses/RafaelRangel-Mestrado-2019.pdf

 Developed at Tecgraf Institute of Technical-Scientific Software Development of
 PUC-Rio (Tecgraf/PUC-Rio)

 Author:
 * Rafael Lopez Rangel - CIMNE/UPC BarcelonaTech (rafaelrangel@tecgraf.puc-rio.br)

 Collaborators:
 * Luiz Fernando Martha - PUC-Rio (lfm@tecgraf.puc-rio.br)
 * Marcos Antonio Rodrigues - UFES (rodriguesma.civil@gmail.com)
 * Cláudio Horta Barbosa de Resende - PUC-Rio (claudio.horta@aluno.puc-rio.br)

--------------------------------------------------------------------------------
**Analysis Model Types**

 Reticulated models can be frames or trusses, which are made up of uniaxial
 members with one dimension much larger than the others, such as bars
 (members with axial behavior) or beams (members with flexural behavior).
 In the present context, members of any type are generically called elements.

 Truss models

 * Truss elements are bars connected by friction-less pins.
 * A truss element does not present any secondary bending moment or torsion
   moment induced by rotation continuity at joints.
 * Therefore, there is only one type of internal force in a truss element:
   axial internal force, which may be tension or compression.
 * A 2D truss model is considered to be laid in the XY-plane, with only
   in-plane behavior, that is, there is no displacement transversal to the
   truss plane.
 * Each node of a 2D truss model has two d.o.f.'s (degrees of freedom):
   a horizontal displacement in X direction and a vertical displacement in
   Y direction.
 * Each node of a 3D truss model has three d.o.f.'s (degrees of freedom):
   displacements in X, Y and Z directions.

 Frame models

 * Frame elements are usually rigidly connected at joints.
   However, a frame element might have a hinge (rotation liberation) at
   an end or hinges at both ends.
 * A 2D frame model is considered to be laid in the XY-plane, with only
   in-plane behavior, that is, there is no displacement transversal to
   the frame plane.
 * Each node of a 2D frame model has three d.o.f.'s: a horizontal
   displacement in X direction, a vertical displacement in Y direction,
   and a rotation about the Z direction.
 * Each node of a 3D frame model has six d.o.f.'s: displacements in
   X, Y and Z directions, and rotations about X, Y and Z directions.
 * It is considered that a hinge realeses all the rotations of an element
   end. In a 3D frame, it means that the rotations about the three
   direction are free.
 * Internal forces at any cross-section of a 2D frame element are:
   axial force (local x direction), shear force (local y direction),
   and bending moment (about local z direction).
 * Internal forces at any cross-section of a 3D frame element are:
   axial force, shear forces (local y direction and local z direction),
   bending moments (about local y direction and local z direction),
   and torsion moment (about x direction).

 Attention: When setting the structural model properties in the input file,
 all six d.o.f's of a node are considered. The program will manage to work only
 with the d.o.f's related to the selected analysis model.
 
--------------------------------------------------------------------------------
**Element Types**

 For frame models, whose elements have bending effects, two types of
 beam elements are considered:

 * Euler-Bernoulli beam element.
 * Timoshenko beam element

 In Euler-Bernoulli flexural behavior, it is assumed that there is no
 shear deformation. As a consequence, bending of an element is such that
 its cross-section remains plane and normal to the element longitudinal axis.

 In Timoshenko flexural behavior, shear deformation is considered in an
 approximated manner. Bending of an element is such that its cross-section
 remains plane but it is not normal to the element longitudinal axis.

 In truss models, the two types of elements may be used indistinguishably,
 since there is no bending behavior of a truss element, and
 Euler-Bernoulli elements and Timoshenko elements are equivalent for the
 axial behavior.

--------------------------------------------------------------------------------
**Element Local Axes**

 In 2D models, the local axes of an element are defined uniquely in the
 following manner:

 * The local x-axis of an element is its longitudinal axis, from its
   initial node to its final node.
 * The local z-axis of an element is always in the direction of the
   global Z-axis, which is perpendicular to the model plane and its
   positive direction points out of the screen.
 * The local y-axis of an element lays on the global XY-plane and is
   perpendicular to the element x-axis in such a way that the
   cross-product (x-axis * y-axis) results in a vector in the
   global Z direction.

 In 3D models, the local y-axis and z-axis are defined by an auxiliary
 vector vz = (vzx, vzy, vzz), which is an element property and should
 be specified as an input data of each element:

 * The local x-axis of an element is its longitudinal axis, from its
   initial node to its final node.
 * The auxiliary vector lays in the local xz-plane of an element, and the
   cross-product (vz * x-axis) defines the local y-axis vector.
 * The direction of the local z-axis is then calculated with the
   cross-product (x-axis * y-axis).

 In 2D models, the auxiliary vector vz must be (0,0,1).
 
 In 3D models, it is important that the auxiliary vector is not parallel
 to the local x-axis; otherwise, the cross-product (vz * x-axis) is zero
 and local y-axis and z-axis will not be defined.

--------------------------------------------------------------------------------
**Load Types**

 There are four types of loads that can be considered:

 * Concentrated nodal force and moment in global axes directions.
 * Uniform distributed force on an element, spanning its entire length,
   in local or in global axes directions.
 * Linearly distributed force on an element, spanning its entire length,
   in local or in global axes directions.
 * Temperature variation along elements local axes (temperature gradient).

 Distributed loads applied in local element direction act as follower
 loads, moving with the element local axes.

 Prescribed displacements (known support settlement values) can also be
 specified to fixed degrees-of-freedom.

 It is assumed that there is a single load case, and all loads are
 proportionally applied during the incremental process of a nonlinear
 anlysis.

 Element loads (distributed or thermal) are ignored on truss models.
 No type of load is considered in vibration analyses.

--------------------------------------------------------------------------------
**Materials**

 All materials are considered to have a linear-elastic bahavior.
 In adition, homogeneous and isotropic properties are also considered,
 that is, all materials have the same properties at every point and in
 all directions.

--------------------------------------------------------------------------------
**Nonlinear Analysis Options**

 Tangent Stiffness Matrix Formulations:
 * Corotational
 * Updated Lagrangian

 Solution Methods:
 * Euler's pure incremental method
 * Load control
 * Work control
 * Linear Arc-Length (fixed normal plane)
 * Linear Arc-Length (updated normal plane)
 * Constant Arc-Length (cylindrical)
 * Constant Arc-Length (spherical)
 * Minimum norm
 * Orthogonal residue
 * Generalized displacement
 
 Increment Adjustment:
 * Constant increment in all steps
 * Adjusted increment in each step
 
 Iteration Strategies:
 * Standard Newton-Raphson
 * Modified Newton-Raphson

 Nonlinear Analysis Parameters:
 * Initial increment: Increment of load ratio in the first iteration of the first step.
 * Maximum load ratio: Maximum value of the load ratio to stop analysis.
 * Maximum steps: Maximum number of steps to stop analysis.
 * Maximum iterations: Maximum number of iterations in each step.
 * Target iterations: Desired number of iterations in each step to adjust the increment.
 * Tolerance: Tolerance to assume that equilibrium has been reached.

--------------------------------------------------------------------------------
**Vibration Analysis Options**

 Mass Matrix Formulations:
 * Consistent
 * Lumped
 
 Damping:
 * Undamped
 * Proportional Damping (Rayleigh damping)
 
 Solution Methods:
 * Runge-Kutta (4th order)
 * Newmark
 * Wilson-Theta
 * Adams-Moulton (3rd order)
 
 Result domain:
 * Time
 * Frequency

 Vibration Analysis Parameters:
 * Time increment: Increment of time in each step.
 * Number of steps: Total number of incremental time steps.
