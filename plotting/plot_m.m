function plot_m(elemType,dispNodes,dispNodes1)

% Plots the finite element mesh with the support condition
global node element

v=get(0,'ScreenSize');
figure('Color',[1 1 1])
hold on
plot_mesh(node,element,elemType,'b-');
plot(node(dispNodes,1),node(dispNodes,2),'ks');
plot(node(dispNodes1,1),node(dispNodes1,2),'ko');
axis off

end   % end of function