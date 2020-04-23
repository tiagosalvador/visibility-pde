function U = visibility_solver1D_single(g,x,xstar,gstar)
% Solver to compute the 1D upper SS envelope (visibility) of g.
% This solves the PDE min(u-g,du/dn(x)) = 0 where
% du/dn(x) = \grad u \cdot (x-x*) with x* = xstar.
% We march out from the point x* and we take the max
% Args:
%    g: obstacle function evaluated at the grid
%    x: coordinates of the grid points
%    gstar: g at xstar
% Returns:
%    Solution U.

% We are given the function g, which is a rectangular matrix
n = length(g);
U = g;

ix = max(find(x>xstar,1,'first')-1,1);

% The vector field goes out from x*.
% We are going to evaluate the vector field  (i-ix,j-jx)*(Wx,Wy)
% and compare the upwind solve with the current value of so that
% Wnew(i,j) = max(W(i,j), W(upwind(i,j)))


U(ix) = max(g(ix),gstar);
U(ix+1) = max(g(ix+1),gstar);

for i = (ix+2):n
    U(i) = max(g(i),U(i-1));
end

for i = (ix-1):-1:1
    U(i) = max(g(i),U(i+1));
end