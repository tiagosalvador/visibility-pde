function W = visibility_solver2D_alo(g,x,y,xstar,ystar,gstar)
% Solver to compute the visibility from multiple viewpoints (xstar,ystar)
% where visibility is defined as being seen by all the viewpoints.
% This solves the PDE min(u-g,max_i min_t(du/dni(x)) = 0
% where du/dni(x) = \grad u \cdot (x-xi*)
% In order to do this, we solve min(u-g,du/dni(x)) = 0
% for each viewpoint and take the minimum of the solutions
% Args:
%    g: obstacle function evaluated at the grid
%    x,y: x, y are the coordinates of the grid points, respectively
%    gstar: g at (xstar,ystar)
% Returns:
%    Solution U.


W = Inf*ones(size(g));

for i = 1:length(xstar)
    W = min(W,visibility_solver2D_single(g,x,y,xstar(i),ystar(i),gstar(i)));
end