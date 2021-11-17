function plot_defo(fac,u_x,u_y,elemType)

% plots the color coded displacement intensity in the finite element 
% region along with color bar scale.

global node element

figure
clf
subplot(2,1,1);
plot_field(node+fac*[u_x u_y],element,elemType,u_x);
colorbar
title('Deformation plot, U_X')

subplot(2,1,2);
plot_field(node+fac*[u_x u_y],element,elemType,u_y);
colorbar
title('Deformation plot, U_Y')
end  % end of function


