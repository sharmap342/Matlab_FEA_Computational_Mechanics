function plot_strain(fac,u_x,u_y,elemType,strain)

% plots the color coded strain distribution in the finite element
% region along with color bar scale.
global node element
figure
clf
subplot(2,1,1);
plot_field(node+fac*[u_x u_y],element,elemType,strain(:,:,1));
colorbar
title('Strain plot, sigma__xx')

subplot(2,1,2);
plot_field(node+fac*[u_x u_y],element,elemType,strain(:,:,2));
colorbar
title('strain plot, sigma_yy')

figure
clf
subplot(2,1,1);
plot_field(node+fac*[u_x u_y],element,elemType,strain(:,:,3));
colorbar
title('strain plot, sigma__xy')

subplot(2,1,2);
plot_field(node+fac*[u_x u_y],element,elemType,strain(:,:,4));
colorbar
title('strain plot, sigma__zz')
end % end of function



