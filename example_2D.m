%% 2D Example
% Computing the visibility set from a single viewpoint
% and multiple viewpoints

%% Setup (general)

n = 256;
x = linspace(-2,2,n);
y = linspace(-2,2,n);
[X,Y] = meshgrid(x,y);

a1 = -1.5;
b1 = -0.2;
g1 = (1/0.25)/2*max(abs(X-a1),abs(Y-b1));

a2 = 0;
b2 = 0.3;
g2 = max(abs(X-a2),abs(Y-b2));

a3 = -0.3;
b3 = 1.5;
g3 = sqrt((X-a3).^2+(Y-b3).^2);

a4 = -0.3;
b4 = -1.4;
g4 = sqrt((X-a4).^2+(Y-b4).^2);

g = min(g1,g2); 
g = min(g,g3);
g = min(g,g4);
g = -g;

xstar = [-1.5 1.5];
ystar = [-1.4 -0.3];

gstar1 = (1/0.25)/2*max(abs(xstar-a1),abs(ystar-b1));
gstar2 = max(abs(xstar-a2),abs(ystar-b2));
gstar3 = sqrt((xstar-a3).^2+(ystar-b3).^2);
gstar4 = sqrt((xstar-a4).^2+(ystar-b4).^2);
gstar = min(gstar1,gstar2);
gstar = min(gstar,gstar3);
gstar = min(gstar,gstar4);
gstar = -gstar;


%%  Solve (single viewpoint)
k = 1;
U1 = visibility_solver2D_single(g,x,y,xstar(k),ystar(k),gstar(k));
k = 2;
U2 = visibility_solver2D_single(g,x,y,xstar(k),ystar(k),gstar(k));

%% Plot visibility set (single viewpoint)
v = -0.5;
subplot(2,2,1)
k = 1;
plot_visibility_set(xstar(k),ystar(k),X,Y,g,U1,v,'Visibility to a single point')

subplot(2,2,2)
k = 2;
plot_visibility_set(xstar(k),ystar(k),X,Y,g,U2,v,'Visibility to a single point')

%%  Solve
W1 = visibility_solver2D_alo(g,x,y,xstar,ystar,gstar);
W2 = visibility_solver2D_all(g,x,y,xstar,ystar,gstar);

%% Plot visibility set (multiple viewpoints)
subplot(2,2,3)
plot_visibility_set(xstar,ystar,X,Y,g,W2,v,'Visibility by all viewpoints')
subplot(2,2,4)
plot_visibility_set(xstar,ystar,X,Y,g,W1,v,'Visibility by at least one viewpoint')


