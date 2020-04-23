%% Example Bunny - computing the visibility set to a bunny shape

%% Setup Bunny
clear
addpath(genpath(pwd))
n = 128; %grid size
filename_levelset = ['./bunny reconstruction/bunny_level_set',num2str(n),'.mat'];
if exist(filename_levelset, 'file') == 0
    filename = './bunny reconstruction/bunny_points.mat';
    if exist(filename,'file') == 0
        compute_bunny_points
    end
    compute_bunny_level_set(n)
end
load(filename_levelset)     

%% Setup grid
x = linspace(-0.2,0.4,n);
y = linspace(-0.2,0.4,n);
z = linspace(-0.2,0.4,n);
[X,Y,Z] = meshgrid(x,y,z);

%% Define viewpoint
istar = find(0.15>x,1,'last');
xstar = x(istar);
jstar = find(0.1>y,1,'last');
ystar = y(jstar);
kstar = find(-0.15>z,1,'last');
zstar = z(kstar);
gstar = g(istar,jstar,kstar);

v = 0;

%%  Solve
W = visibility_solver3D_single(g,x,y,z,xstar,ystar,zstar,gstar);


%% Plotting

figure(1);
p1 = patch(isosurface(X,Y,Z,g,v));
isonormals(X,Y,Z,g,p1)
p1.FaceColor = [0 0.1 0.9];
p1.EdgeColor = 'none';
daspect([1 1 1])
view(83,-79);
lightangle(83,-79)
axis off
axis([-0.2 0.2 -0.2 0.2 -0.2 0.2])
alpha(0.7)
hold on
scatter3(xstar,ystar,zstar,100,'k','Filled')
hold off
title('Bunny obstacle')

%%
figure(2);
p2 = patch(isosurface(X,Y,Z,W,v-1.0e-6));
isonormals(X,Y,Z,W,p2)
p2.FaceColor = [1 0 0];
p2.EdgeColor = 'none';
hold on
p1 = patch(isosurface(X,Y,Z,g,v));
isonormals(X,Y,Z,g,p1)
p1.FaceColor = [0 0.1 0.9];
p1.EdgeColor = 'none';
daspect([1 1 1])
view(83,-79);
lightangle(83,-79)
axis off
axis([-0.2 0.2 -0.2 0.2 -0.2 0.2])
alpha(0.8)
scatter3(xstar,ystar,zstar,100,'k','Filled')
hold off
title('Visibility')
