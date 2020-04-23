%% 3D Example
% Computing the visibility set from a single viewpoint to two buildings

%% Setup Grid

n = 256;
x = linspace(-5,5,n);
y = linspace(-5,5,n);
z = linspace(-5,5,n);
[X,Y,Z] = meshgrid(x,y,z);

%% Setup Obstacles

a1 = -2;
b1 = 0;
c1 = 0;
a2 = 3;
b2 = 4;
c2 = 0;
g1 = max(abs(X-a1),max(abs(Y-b1),Z-c1))-1;
g2 = max(abs(X-a2),max(abs(Y-b2),Z-c2))-2;
g = -min(g1,g2); 

%% Setup viewpoint
xstar = 0;
ystar = 0;
zstar = 0;
gstar = -min(max(abs(xstar-a1),max(abs(ystar-b1),zstar-c1))-1,max(abs(xstar-a2),max(abs(ystar-b2),zstar-c2))-2);
v = 0;

%% Solve
W = visibility_solver3D_single(g,x,y,z,xstar,ystar,zstar,gstar);

%% Plotting

figure(1);
p1 = patch(isosurface(X,Y,Z,g,v));
isonormals(X,Y,Z,g,p1)
p1.FaceColor = [0 0.1 0.9];
p1.EdgeColor = 'none';
daspect([1 1 1])
view(27,11);
lightangle(27,21)
axis off
axis([-5 5 -5 5 -5 5])
alpha(0.7)
hold on
scatter3(xstar,ystar,zstar,100,'k','Filled')
hold off
title('Obstacles')

figure(2)
p1 = patch(isosurface(X,Y,Z,g,v));
isonormals(X,Y,Z,g,p1)
p1.FaceColor = [0 0.1 0.9];
p1.EdgeColor = 'none';
hold on
p2 = patch(isosurface(X,Y,Z,W,v));
isonormals(X,Y,Z,W,p2)
p2.FaceColor = [1 0 0];
p2.EdgeColor = 'none';
daspect([1 1 1])
view(27,11);
lightangle(27,21)
axis off
axis([-5 5 -5 5 -5 5])
alpha(.8)
scatter3(xstar,ystar,zstar,100,'k','Filled')
hold off
title('Visibility')