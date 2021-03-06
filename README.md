# FEA Matlab Code for Slope Stability
## Finite Element Code for Slope Stability using Elasto-plastic Analysis 
### Description of Packages, Modules and Functions
The main package have modules such as Elastic_Master_Code, Elastoplastic_Master_Code and setpath. The module Elastic_Master_Code sets input parameters for slope stability analysis using elastic finite element method. The input parameters such as dimensions of the geometry, mesh type, material properties, initial stress conditions, boundary conditions, and number of steps for the solver are defined in this module. Similarly, Elastoplastic_Master_Code sets input parameters for slope stability analysis using elasto-plastic finite element method.  

The setpath module sets path for other modules to be called from intergratio, meshing, plotting and subfunctions packages.      

The package integrations contains modules with functions named shape function, B matrix, and gauss point. The shape function takes type, coordinate and dimension of element and returns a shape function and its derivative with respect to coordinates. The function gauss_rule provides the weight vector and gauss point coordinate matrix of the element based on the integration order selected. The function gauss_wt_pt provides the weights and coordinates of the gauss integration points. Finally, based on the shape functions and derivatives of shape functions the Bmatrix returns strain displacement matrix for each element.     

The package meshing creates a finite element mesh based on the type and dimensions of the elements provided by the user. it contains make_elem, mesh_region, mesh_t6_elem, q4totq9, q4toto8, square_node_array, structured_q8_mesh, structured_q9_mesh. The module make_elem creates a connectivity lest of primary nodes in Q4 and T3 element. The function mesh_region generates an array of nodal connectivity (coordinates of each node and element connectivity (nodes of each element with counterclockwise ordering.) given the four corners points of the domain,number of elements in each direction (numx,numy)and the element type (Q4,Q8,Q9 and T3). The module mesh_t6_elem forms the element and node matrices for T6 element, with the element and node matrices arranged in a nodal counterclockwise order (primary nodes then secondary nodes). The module q4toto8 forms the element and node matrices for eight node rectangular element, from the elemet and node matrices of a four node rectangular element with the element and node matrices arranged in a counterclockwish order Inputs:-element,node,numx and numy. The module q4totq9 forms the element and node matrices for nine node rectangular element given the elemet and node matrices of a four node rectangular element and the number of elements in x and y direction. The module square_node_array generates a quadratleral array of nodes between the counterclockwise ordering of nodes pt1 - pt4,given number of elements in x and y direction. The module structured_q8_mesh %Forms the element and node matrices for Q8 element, with the element and node matrices arranged in a nodal counterclockwise order (primary nodes then secondary nodes). The module structured_q9_mesh forms the element and node matrices for Q9 element, with the element and node matrices arranged in a nodal counterclockwise order (primary nodes-> secondary nodes then internal nodes).   

The package plotting is used to plot geometry of the model, mesh, contours of displacement, and contours of stress. This package contains plot_defo, pot_field, plot_m, plot_mesh, plot_sig, plot_strain modules. The module plot_defo % plots the color coded displacement intensity in the finite element region along with color bar scale. The module plot_field forms a color dependent finite element mesh for plotting outputs. The module plot_mesh plots the finite element mesh. The module plot_field plots the color coded stress distribution in the finite element region along with color bar scale. The module plots the color coded strain distribution in the finite element region along with color bar scale.        

The package subfunctions evluates different variabes, functions, used during solving the the finite element method. This module contains displacements, elementdof, force_matrix, formdg, formm, internalrxn, invariants, makedfds, makedgds, plastic_mat, principal_stress, selfwt_matrix, stiffness_matrix, stress_calulation, supportcond. The module discplacements evaluates the unknown degree of freedom (displacements) at the nodes. The module elementdof forms the global degree of freedom nodes in each element from the node identification number. The module force_matrix generates the force matrix due to externally applied loads for Q4 and T3 elements and the location of these loads. the module formdg calculates the partial derivatives of the plastic potential with respect to p, J2 and J3. The module formm generates the partial derivatives of the p, J2 and J3 with respect to stress for use in the plastic potential derivative. The module internalrxn generates the vector of nodal force reactions due to internal stress. it can be used out side the iteration loop if the stresses at each gauss point is given. The module invariants calculates the stress invariants for a single gauss point. The module makedfds forms the partial derivatives of the yield function with respect to stress at each gauss point. The module plastic_mat generates the plastic constitutive matix CP at each gauss point. The module principal_stress calcualtes the principal stresses from the stress invariants. The module selfwt_matirx generates the force matrix due to self weight. The module stiffness_matrix generates the element and global stiffness matrix. The module stress_calculation calculates the element strains and stresses at the nodes in x, y and xy directions. The module supportcond forms a vector of nodes which are restrained in x direction or in both x and y direction which are the supports of the domain. It also forms a matrix of nodes for load application at the top of the domain.       

References 
The code could have not been written without the help of Peeyush Garg and MASc thesis from NTNU

### Output 

Finite Element Mesh    

![image](https://user-images.githubusercontent.com/61520478/142275103-b1595245-bdbe-434a-a029-b4bdf20a7c5c.png)     

Horizontal and Vertical Displacement (m)   

![image](https://user-images.githubusercontent.com/61520478/142275329-0e100c40-de54-4e25-a568-707dd1ab527f.png)    

Horizontal and Inplane stress stress (KPa)  

![image](https://user-images.githubusercontent.com/61520478/142275944-44482f91-408a-4be2-b25a-38bd28b8978e.png)

Vertical and Shear stress (KPa)

![image](https://user-images.githubusercontent.com/61520478/142275796-4150dc8d-f557-46d1-b3f8-41c87a70606b.png)
