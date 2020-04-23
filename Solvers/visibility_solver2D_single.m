function U = visibility_solver2D_single(g,x,y,xstar,ystar,gstar)
% Solver to compute the 2D upper SS envelope (visibility) of g.
% This solves the PDE min(u-g,du/dn(x)) = 0 where
% du/dn(x) = \grad u \cdot (x-x*) with x* = (xstar,ystar).
% We march out from the point x* and we take the max
% Args:
%    g: obstacle function evaluated at the grid
%    x,y: x, y are the coordinates of the grid points, respectively
%    gstar: g at (xstar,ystar)
% Returns:
%    Solution U.

% We are given the function g, which is a rectangular matrix
[n,m] = size(g);
U = g;

ix = max(find(y>ystar,1,'first')-1,1);
jx = max(find(x>xstar,1,'first')-1,1);

% The vector field goes out from x*.
% We are going to evaluate the vector field  (i-ix,j-jx)*(Wx,Wy)
% and compare the upwind solve with the current value of so that
% Wnew(i,j) = max(W(i,j), W(upwind(i,j)))


U(ix,jx) = max(g(ix,jx),gstar);
U(ix+1,jx) = max(g(ix+1,jx),gstar);
U(ix,jx+1) = max(g(ix,jx+1),gstar);
U(ix+1,jx+1) = max(g(ix+1,jx+1),gstar);


%% Loops over the strips

for j = jx-1:-1:1
    for i=ix:1:ix+1
        val = ValueUpdateFlipped(xstar,ystar,x,y,U,i,j);
        % Update
        U(i,j) =  max(val, g(i,j));
    end
end

for j = jx+2:1:m
    for i=ix:1:ix+1
        val = ValueUpdateFlipped(xstar,ystar,x,y,U,i,j);
        % Update
        U(i,j) =  max(val, g(i,j));
    end
end

for i = ix-1:-1:1
    for j = jx:1:jx+1
        val = ValueUpdateFlipped(xstar,ystar,x,y,U,i,j);
        % Update
        U(i,j) =  max(val, g(i,j));
    end
end
for i = ix+2:1:n
    for j = jx:1:jx+1
        val = ValueUpdateFlipped(xstar,ystar,x,y,U,i,j);
        % Update
        U(i,j) =  max(val, g(i,j));
    end
end


%% Loops over bottom part
for i = ix+2:1:n
    % Loop over the bottom left square
    for j = jx-1:-1:1
        val = ValueUpdateFlipped(xstar,ystar,x,y,U,i,j);
        % Update
        U(i,j) =  max(val, g(i,j));
    end
    % Loop over the bottom right square
    for j = jx+2:1:m
        val = ValueUpdateFlipped(xstar,ystar,x,y,U,i,j);
        % Update
        U(i,j) =  max(val, g(i,j));
    end
end

%% Loops over top half
for i = ix-1:-1:1
    % Loop over the top left square
    for j = jx-1:-1:1
        val = ValueUpdateFlipped(xstar,ystar,x,y,U,i,j);
        % Update
        U(i,j) =  max(val, g(i,j));
    end
    
    % Loop over the top right square
    for j = jx+2:1:m
        val = ValueUpdateFlipped(xstar,ystar,x,y,U,i,j);
        % Update
        U(i,j) =  max(val, g(i,j));
    end
end

end

function val = ValueUpdateFlipped(xstar,ystar,x,y,U,i,j)
% this function interpolates the values of the vector field 
% m1,m2 intersected with the nine point grid values of U.
% the line itersects a segment with two of the neighbors.
%
% We are flipping the sign of the vector field from the other case.
m1 = -(y(i)-ystar);
m2 = -(x(j)-xstar);
s1 = sign(m1);
s2 = sign(m2);
if abs(m2) < abs(m1)
    %   |m2| < |m1| then uS is neighbor.
    %   m2 < 0  other neighbor is USW
    t = abs(m2)/abs(m1);
    val = (1-t)*U(i+s1,j) + t*U(i+s1,j+s2);
else
    % in this case, when m2 < 0, the neighbors are uW, uSW
    t = abs(m1)/abs(m2);
    val = (1-t)*U(i,j+s2) + t*U(i+s1,j+s2);
end

end