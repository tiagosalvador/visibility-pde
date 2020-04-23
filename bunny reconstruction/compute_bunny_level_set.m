function compute_bunny_level_set(n)
% This scrip computes a level-set function whose
% zero level-set function is a bunny

load('bunny_points.mat','points')
z = double(points(:,1));
y = double(points(:,2));
x = double(points(:,3));

npoints = length(x);

xgrid = linspace(-0.2,0.4,n);
ygrid = xgrid;
zgrid = xgrid;
[Xgrid,Ygrid,Zgrid] = meshgrid(xgrid,ygrid,zgrid);

start_progress(' - Determining enclosing grid points')
b = zeros(size(Xgrid));
for i=1:npoints
    ix = find(x(i)>xgrid,1,'last');
    jx = find(y(i)>ygrid,1,'last');
    kx = find(z(i)>zgrid,1,'last');
    b(ix,jx,kx) = b(ix,jx,kx)+1;
    display_progress(i,npoints,1);
end

i = 1:n-1;
j = 1:n-1;
k = 1:n-1;

midIn = b(i,j,k)+b(i+1,j,k)+b(i+1,j+1,k)+b(i+1,j,k+1)+b(i+1,j+1,k)+b(i+1,j,k+1)+b(i,j+1,k+1)+b(i+1,j+1,k+1);

dims = [n n n];
midIn(n,n,n) = 0;

start_progress(' - Determining grid points outside')
ind_outside = 1;
adding = 1;
while true
    display_progress(floor(100*length(ind_outside)/prod(dims)),100,1);
    [I,J,K] = ind2sub(dims,adding);
    [indNup,indNEup,indEup,indSEup,indSup,indSWup,indWup,indNWup,indup,...
        indNdown,indNEdown,indEdown,indSEdown,indSdown,indSWdown,indWdown,indNWdown,inddown,...
        indN,indNE,indE,indSE,indS,indSW,indW,indNW] = find_ind_neighbors3d(I,J,K,dims);
    ind_neighbors = [indNup;indNEup;indEup;indSEup;indSup;indSWup;indWup;indNWup;indup;...
        indNdown;indNEdown;indEdown;indSEdown;indSdown;indSWdown;indWdown;indNWdown;inddown;...
        indN;indNE;indE;indSE;indS;indSW;indW;indNW];
    
    maybe_adding = setdiff(ind_neighbors,ind_outside);
    
    adding = maybe_adding(midIn(maybe_adding)==0);
    
    ind_outside = [adding;ind_outside];
    
    if isempty(adding)
        break
    end
end
display_progress(100,100,1);

inside = ones(dims);
inside(ind_outside) = 0;

start_progress(' - Computing distance function')
% Computing distance function
g = inf(size(Xgrid));
for i=1:npoints
    d = (x(i)-Ygrid).^2+(y(i)-Xgrid).^2+(z(i)-Zgrid).^2;
    g = min(g,d);
	display_progress(i,npoints,1);
end
%%

g(inside<=0) = -g(inside<=0);

save(['./bunny reconstruction/bunny_level_set',num2str(n)],'g')

%% Plotting results

% p1 = patch(isosurface(Xgrid,Ygrid,Zgrid,g,0));
% isonormals(Xgrid,Ygrid,Zgrid,g,p1)
% p1.FaceColor = [0 0.1 0.9];
% p1.EdgeColor = 'black';
% xlabel('x')
% ylabel('y')
% zlabel('z')
end