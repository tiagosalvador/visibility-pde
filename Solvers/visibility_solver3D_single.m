function U = visibility_solver3D_single(g,x,y,z,xstar,ystar,zstar,gstar)
% Solver to compute the 3D upper SS envelope (visibility) of g.
% This solves the PDE min(u-g,du/dn(x)) = 0 where
% du/dn(x) = \grad u \cdot (x-x*) with x* = (xstar,ystar,zstar).
% We march out from the point x* and we take the max
% Args:
%    g: obstacle function evaluated at the grid
%    x,y,z: x, y, z are the coordinates of the grid points, respectively
%    gstar: g at (xstar,ystar,zstar)
% Returns:
%    Solution W.

% We are given the function g, which is a rectangular matrix
[ni,nj,nk] = size(g);
U = g;
ix = max(find(y>ystar,1,'first')-1,1);
jx = max(find(x>xstar,1,'first')-1,1);
kx = max(find(z>zstar,1,'first')-1,1);

% The vector field goes out from x*.
% We are going to evaluate the vector field  (i-ix,j-jx,k-kx)*(Wx,Wy,Wz)
% and compare the upwind solve with the current value of so that
% Wnew(i,j,k) = max(W(i,j,k), W(upwind(i,j,k)))

% 
U(ix,jx,kx) = max(g(ix,jx,kx),gstar);
U(ix+1,jx,kx) = max(g(ix+1,jx,kx),gstar);
U(ix,jx+1,kx) = max(g(ix,jx+1,kx),gstar);
U(ix+1,jx+1,kx) = max(g(ix+1,jx+1,kx),gstar);

U(ix,jx,kx+1) = max(g(ix,jx,kx+1),gstar);
U(ix+1,jx,kx+1) = max(g(ix+1,jx,kx+1),gstar);
U(ix,jx+1,kx+1) = max(g(ix,jx+1,kx+1),gstar);
U(ix+1,jx+1,kx+1) = max(g(ix+1,jx+1,kx+1),gstar);
 
%% Loop over tubes

% along y
for i = ix+2:1:ni
    for j = jx:1:jx+1
        for k = kx:1:kx+1
            val = ValueUpdateFlipped(xstar,ystar,zstar,x,y,z,U,i,j,k);
            U(i,j,k) =  max(val, g(i,j,k));
        end
    end
end
for i = ix-1:-1:1
    for j = jx:1:jx+1
        for k = kx:1:kx+1
            val = ValueUpdateFlipped(xstar,ystar,zstar,x,y,z,U,i,j,k);
            U(i,j,k) =  max(val, g(i,j,k));
        end
    end
end

% along x
for j = jx+2:1:nj
    for i = ix:1:ix+1
        for k = kx:1:kx+1
            val = ValueUpdateFlipped(xstar,ystar,zstar,x,y,z,U,i,j,k);
            U(i,j,k) =  max(val, g(i,j,k));
        end
    end
end
for j = jx-1:-1:1
    for i = ix:1:ix+1
        for k = kx:1:kx+1
            val = ValueUpdateFlipped(xstar,ystar,zstar,x,y,z,U,i,j,k);
            U(i,j,k) =  max(val, g(i,j,k));
        end
    end
end

% along z
for k = kx+2:1:nk
    for i = ix:1:ix+1
        for j = jx:1:jx+1
            val = ValueUpdateFlipped(xstar,ystar,zstar,x,y,z,U,i,j,k);
            U(i,j,k) =  max(val, g(i,j,k));
        end
    end
end
for k = kx-1:-1:1
    for i = ix:1:ix+1
        for j = jx:1:jx+1
            val = ValueUpdateFlipped(xstar,ystar,zstar,x,y,z,U,i,j,k);
            U(i,j,k) =  max(val, g(i,j,k));
        end
    end
end
%% Loops over the plane strips

% y cross section
for j = jx+2:1:nj
    for k = kx+2:1:nk
        for i = ix:1:ix+1
            val = ValueUpdateFlipped(xstar,ystar,zstar,x,y,z,U,i,j,k);
            U(i,j,k) =  max(val, g(i,j,k));
        end
    end
    for k = kx-1:-1:1
        for i = ix:1:ix+1
            val = ValueUpdateFlipped(xstar,ystar,zstar,x,y,z,U,i,j,k);
            U(i,j,k) =  max(val, g(i,j,k));
        end
    end
end
for j = jx-1:-1:1
    for k = kx+2:1:nk
        for i = ix:1:ix+1
            val = ValueUpdateFlipped(xstar,ystar,zstar,x,y,z,U,i,j,k);
            U(i,j,k) =  max(val, g(i,j,k));
        end
    end
    for k = kx-1:-1:1
        for i = ix:1:ix+1
            val = ValueUpdateFlipped(xstar,ystar,zstar,x,y,z,U,i,j,k);
            U(i,j,k) =  max(val, g(i,j,k));
        end
    end
end

% x cross section
for i = ix+2:1:ni
    for k = kx+2:1:nk
        for j = jx:1:jx+1
            val = ValueUpdateFlipped(xstar,ystar,zstar,x,y,z,U,i,j,k);
            U(i,j,k) =  max(val, g(i,j,k));
        end
    end
    for k = kx-1:-1:1
        for j = jx:1:jx+1
            val = ValueUpdateFlipped(xstar,ystar,zstar,x,y,z,U,i,j,k);
            U(i,j,k) =  max(val, g(i,j,k));
        end
    end
end
for i = ix-1:-1:1
    for k = kx+2:1:nk
        for j = jx:1:jx+1
            val = ValueUpdateFlipped(xstar,ystar,zstar,x,y,z,U,i,j,k);
            U(i,j,k) =  max(val, g(i,j,k));
        end
    end
    for k = kx-1:-1:1
        for j = jx:1:jx+1
            val = ValueUpdateFlipped(xstar,ystar,zstar,x,y,z,U,i,j,k);
            U(i,j,k) =  max(val, g(i,j,k));
        end
    end
end

% z cross section
for i = ix+2:1:ni
    for j = jx+2:1:nj
        for k = kx:1:kx+1
            val = ValueUpdateFlipped(xstar,ystar,zstar,x,y,z,U,i,j,k);
            U(i,j,k) =  max(val, g(i,j,k));
        end
    end
    for j = jx-1:-1:1
        for k = kx:1:kx+1
            val = ValueUpdateFlipped(xstar,ystar,zstar,x,y,z,U,i,j,k);
            U(i,j,k) =  max(val, g(i,j,k));
        end
    end
end
for i = ix-1:-1:1
    for j = jx+2:1:nj
        for k = kx:1:kx+1
            val = ValueUpdateFlipped(xstar,ystar,zstar,x,y,z,U,i,j,k);
            U(i,j,k) =  max(val, g(i,j,k));
        end
    end
    for j = jx-1:-1:1
        for k = kx:1:kx+1
            val = ValueUpdateFlipped(xstar,ystar,zstar,x,y,z,U,i,j,k);
            U(i,j,k) =  max(val, g(i,j,k));
        end
    end
end

%% loop over the 8 octants
 
for i = ix+2:1:ni
    for j = jx+2:1:nj
        for k = kx+2:1:nk
            val = ValueUpdateFlipped(xstar,ystar,zstar,x,y,z,U,i,j,k);
            U(i,j,k) =  max(val, g(i,j,k));
        end
        for k = kx-1:-1:1
            val = ValueUpdateFlipped(xstar,ystar,zstar,x,y,z,U,i,j,k);
            U(i,j,k) =  max(val, g(i,j,k));
        end
    end
    for j = jx-1:-1:1
        for k = kx+2:1:nk
            val = ValueUpdateFlipped(xstar,ystar,zstar,x,y,z,U,i,j,k);
            U(i,j,k) =  max(val, g(i,j,k));
        end
        for k = kx-1:-1:1
            val = ValueUpdateFlipped(xstar,ystar,zstar,x,y,z,U,i,j,k);
            U(i,j,k) =  max(val, g(i,j,k));
        end
    end
end
for i = ix-1:-1:1
    for j = jx+2:1:nj
        for k = kx+2:1:nk
            val = ValueUpdateFlipped(xstar,ystar,zstar,x,y,z,U,i,j,k);
            U(i,j,k) =  max(val, g(i,j,k));
        end
        for k = kx-1:-1:1
            val = ValueUpdateFlipped(xstar,ystar,zstar,x,y,z,U,i,j,k);
            U(i,j,k) =  max(val, g(i,j,k));
        end
    end
    for j = jx-1:-1:1
        for k = kx+2:1:nk
            val = ValueUpdateFlipped(xstar,ystar,zstar,x,y,z,U,i,j,k);
            U(i,j,k) =  max(val, g(i,j,k));
        end
        for k = kx-1:-1:1
            val = ValueUpdateFlipped(xstar,ystar,zstar,x,y,z,U,i,j,k);
            U(i,j,k) =  max(val, g(i,j,k));
        end
    end
end

end

function val = ValueUpdateFlipped(xstar,ystar,zstar,x,y,z,U,i,j,k)
% this function interpolates the values of the vector field 
% m1,m2,m3 intersected with the 27 point grid values of U.
% the line itersects a plane with four of the neighbors.
%
% We are flipping the sign of the vector field from the other case.
m1 = -(y(i)-ystar);
m2 = -(x(j)-xstar);
m3 = -(z(k)-zstar);
s1 = sign(m1);
s2 = sign(m2);
s3 = sign(m3);
if abs(m1) > max(abs(m2),abs(m3))
    t = abs(m2)/abs(m1);
    s = abs(m3)/abs(m1);
    U1 = U(i+s1,j,   k);
    U2 = U(i+s1,j,   k+s3);
    U3 = U(i+s1,j+s2,k);
    U4 = U(i+s1,j+s2,k+s3);
elseif abs(m2) > abs(m3)
    t = abs(m1)/abs(m2);
    s = abs(m3)/abs(m2);
    U1 = U(i,   j+s2,k);
    U2 = U(i,   j+s2,k+s3);
    U3 = U(i+s1,j+s2,k);
    U4 = U(i+s1,j+s2,k+s3);
else
    t = abs(m1)/abs(m3);
    s = abs(m2)/abs(m3);
    U1 = U(i,   j,   k+s3);
    U2 = U(i,   j+s2,k+s3);
    U3 = U(i+s1,j,   k+s3);
    U4 = U(i+s1,j+s2,k+s3);
end
val = s*(t*U4+(1-t)*U2)+(1-s)*(t*U3+(1-t)*U1);
end