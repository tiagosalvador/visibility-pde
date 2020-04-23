function plot_visibility_set(xstar,ystar,X,Y,g,U,v,nametitle)
% Plots the visibility set as determined by the v level set of U together
% with the obstacles g.

cmap = [0 0 0;1 0 0]; % color map
contourf(X,Y,-U,[-v -v],'r');
hold on;
contourf(X,Y,g,[v v],'k');
colormap(cmap);
plot(xstar,ystar,'*k')
axis([-2 2 -2 2]);
hold off
xticks(-2:0.5:2)
yticks(-2:0.5:2)
title(nametitle)
daspect([1 1 1])

end