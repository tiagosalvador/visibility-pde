function W = upper_star_shapped_solver2D(g,ix,jx)
% Solver to compute the 2D upper SS envelope (visibility) of g
% This solves the PDE max(u-g,-du/dn(x)) = 0
% where du/dn(x) = \grad u \cdot (x-x*) with x* = x_(ix,jx).
% We march toward the point x* and we take the min
% Args:
%    g: obstacle function evaluated at the grid
%    x,y: x, y are the coordinates of the grid points, respectively
%    gstar: g at (xstar,ystar)
% Returns:
%    Solution W.

% We are given the function g, which is a rectangular matrix
[n,m] = size(g);
W = g;
% The vector field goes in to x*.
% We are going to evaluate the vector field  (i-ix,j-jx)*(Ux,Uy)
% and compare the upwind solve with the current value of so that
% Unew(i,j) = min(U(i,j), U(upwind(i,j)))

%% First over the bottom part.
for i = n-1:-1:ix
    % Loop over the bottom left square
    for j = 2:jx
        val = ValueUpdate(ix,jx,W,i,j);
        % Update
        W(i,j) =  min(val, W(i,j));
    end
    % Loop over the bottom right square
    for j = m-1:-1:jx+1
        val = ValueUpdate(ix,jx,W,i,j);
        % Update
        W(i,j) =  min(val, W(i,j));
    end
    
end

%% Next loop up on the top half,
for i = 2:ix
    % Loop over the top left square
    for j = 2:jx
        val = ValueUpdate(ix,jx,W,i,j);
        % Update
        W(i,j) =  min(val, W(i,j));
    end

    % Loop over the top right square
    for j = m-1:-1:jx+1
        val = ValueUpdate(ix,jx,W,i,j);
        % Update
        W(i,j) =  min(val, W(i,j));
    end    
end
end


function val = ValueUpdate(ix,jx,U,i,j)
% this function interpolates the values of the vector field
% m1,m2 intersected with the nine point grid values of U.
% the line itersects a segment with two of the neighbors.
m1 = i-ix;
m2 = j-jx;
s1 = sign(m1);
s2 = sign(m2);
if abs(m2) < abs(m1)
    % |m2| < |m1| then uS is neighbor.
    % m2 < 0  other nbr is USW
    t = abs(m2)/abs(m1);
    val = (1-t)*U(i+s1,j) + t*U(i+s1,j+s2);
else
    % in this case, m2 < 0, the neighbors are uW, uSW
    t = abs(m1)/abs(m2);
    val = (1-t)*U(i,j+s2) + t*U(i+s1,j+s2);
end
end

