%  This script file is a simple FE code for a Rectangual element with 
%  4,8 or 9 nodes or for a triangular element with 3 or 6 nodes 
%  under uniform loading.the input parameters are the dimensions,
%  material properties and selection of element types.
%  The out put will be the deformation and stress of the region.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
clc
tic;   
feval('setpath')

%   GLOBAL VARIABLES
global node element elemType E nu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp([num2str(toc),'    INPUT PARAMETERS '])

E  = 60e6 ;      % modulus of elasticity, KPa
nu = 0.25 ;       % poissons ratio  
L = 50;           % width of domain
D = 30;           % depth of domain
numx=100;          % number of elements in x-direction
numy=60;          % number of elements in y-direction
%gamma=14;         % unit weight of soil in KN/m3 
%%%% Modification in expression of gamma
gamma = [10 20  28;
           20 30 26;
          30 40 24];
    % each row has values:: y1 y2  gamma
    % where y = y1 represents the lower boundary line of segment 
    %       y = y2 represents the upper boundary line of segment
%%%%%
sigmatoy=0;     % vertical load
sigmatox =0;     % horizontal load
load_edge1=0;    % starting x-coordinate point of load, SHOULD BE AT A NODE
load_edge2=10;    % ending x-coordinate point of load, SHOULD BE AT A NODE
elemType ='Q4';
normal_order = 2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% elastic compliance matrix 
 %C= E/(1+nu)/(1-2*nu)*[ 1-nu nu 0;nu 1-nu 0;0 0 0.5-nu];
 C= E/(1+nu)/(1-2*nu)*[ 1-nu nu 0 nu;nu 1-nu 0 nu;0 0 0.5-nu 0;nu nu 0 1-nu];
% Four corner points
pt1 = [0 10];
pt2 = [50  10];
pt3 = [ 50   18];
%%%% Additional
%%%% Conrner points at top layer matrix pte1 is introduced
pte1 = [10 40; 25 25; 35  25; 45 20; 50 18 ];% corner points at topside
% % % n1=15; % Number of point need for curve 
% % %     x1= linspace(pte1(1),pt3(1),n1);
% % %     for i1= 1:n1
% % %         y1(i1) = fun1(x1(i1));
% % %     end
% % %     x1= x1(2:n1);
% % %     y1=y1(2:n1);
% % % pte1= [pte1; [x1' y1']];
%%%%%%%%%%%%%%%%%%%%%
pt4 = [ 0  40];

%%%% Additional :: Plotting of surface after addition of corner pts
xkk=[pt1(:,1)  pt4(:,1)  [pte1(:,1)]'  pt3(:,1)   pt2(:,1)  pt1(:,1) ];
ykk=[pt1(:,2)  pt4(:,2)  [pte1(:,2)]'  pt3(:,2)   pt2(:,2)  pt1(:,2)];
plot([xkk],[ykk],'LineWidth',3);
hold on;
%%%%%%%

disp([num2str(toc),'    MESH GENERATION....'])

 switch elemType
     case {'Q9','Q8','Q4','T3'}
         %%% Adding an argument in existing function : pte1
         %%% Definition of function also adjusted/changed accordingly
      [node,element] = mesh_region(pt1, pt2, pt3,pte1,pt4,numx,numy,elemType);
     otherwise 
      [node,element] =mesh_t6_elem(L,D,numx,numy);
 end
 [topEdge,topEdge1,dispNodes,dispNodes1,leftNodes1]=supportcond(elemType,numx,numy);  

 
disp([num2str(toc),'    STIFFNESS MATRIX COMPUTATION....'])

K=stiffness_matrix(node,element,elemType,normal_order,C);
numnode = size(node,1);
numelem = size(element,1);
total_unknown = numnode*2;
selfwt=selfwt_matrix(elemType,normal_order,gamma,node,element);
  switch elemType
     case {'Q9','Q8','T6'}          
            [f,sctry]=force_matrix689(node,topEdge,topEdge1,sigmatoy,sigmatox,load_edge1,load_edge2);
     case{'Q4','T3'}
            [f,sctry]=force_matrix(node,topEdge,sigmatoy,sigmatox,load_edge1,load_edge2);
  end


disp([num2str(toc),'    CALCULATING DISPLACEMENTS'])

[U,u_x,u_y] =displacements(dispNodes,dispNodes1,numnode,K,f,selfwt);


disp([num2str(toc),'    STRESS COMPUTATION....'])
[stress,strain] =stress_calculation(node,element,elemType,U,normal_order,C);
[r] = internalrxn( node,element,elemType,normal_order,stress);
 s=permute(stress,[2 1 3]);  sigma=permute(stress,[3 2 1]);  
 s1=permute(strain,[2 1 3]);  sigma1=permute(strain,[3 2 1]); %additinal line 
% Plot the FEM mesh 
plot_m(elemType,dispNodes,dispNodes1)
title('Undeformed FE mesh')
%Plot numerical deformed configuration
dispnorm=L/max(sqrt(u_x.^2+u_y.^2));
fac =dispnorm*0.05;                    %magnification factor
plot_def(fac,u_x,u_y,elemType,dispNodes,dispNodes1);
%plot stress and deformation intensity with a colormap
plot_defo(fac,u_x,u_y,elemType)
plot_sig(fac,u_x,u_y,elemType,sigma)
plot_strain(fac,u_x,u_y,elemType,sigma1)% Additional Line
% end of the script file


