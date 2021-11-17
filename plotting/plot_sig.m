function plot_sig(fac,u_x,u_y,elemType,stress)

% plots the color coded stress distribution in the finite element
% region along with color bar scale.
global node element
figure
clf
subplot(2,1,1);
plot_field(node+fac*[u_x u_y],element,elemType,stress(:,:,1));
colorbar
title('Stress plot, sigma__xx')

subplot(2,1,2);
plot_field(node+fac*[u_x u_y],element,elemType,stress(:,:,2));
colorbar
title('Stress plot, sigma_yy')

figure
clf
subplot(2,1,1);
plot_field(node+fac*[u_x u_y],element,elemType,stress(:,:,3));
colorbar
title('Stress plot, sigma__xy')

subplot(2,1,2);
plot_field(node+fac*[u_x u_y],element,elemType,stress(:,:,4));
colorbar
title('Stress plot, sigma__zz')

figure
strength=30;
plot_field(node+fac*[u_x u_y],element,elemType,strength\stress(:,:,3));
colorbar
title('Factor of Safety plot, strength \Shear Stress(sigma__xy) ')

end % end of function




