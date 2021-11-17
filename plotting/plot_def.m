function plot_def(fac,u_x,u_y,elemType,dispNodes,dispNodes1)

% Plots the deformed finite element mesh with the support condition

global node element 

figure
clf
hold on
plot_mesh(node+fac*[u_x u_y],element,elemType,'b-');
title(' Numerical deformed mesh ')
plot(node(dispNodes,1),node(dispNodes,2),'ks');
plot(node(dispNodes1,1),node(dispNodes1,2),'ko');

end % end of function