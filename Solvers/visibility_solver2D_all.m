function W = visibility_solver2D_all(g,x,y,xstar,ystar,gstar)
% Solver to compute the visibility from multiple viewpoints (xstar,ystar)
% where visibility is defined as being seen by all the viewpoints.
% Args:
%    g: obstacle function evaluated at the grid
%    x,y: x, y are the coordinates of the grid points, respectively
%    gstar: g at (xstar,ystar)
% Returns:
%    Solution W.

W = -Inf*ones(size(g));

for i = 1:length(xstar)
    W = max(W,visibility_solver2D_single(g,x,y,xstar(i),ystar(i),gstar(i)));
end