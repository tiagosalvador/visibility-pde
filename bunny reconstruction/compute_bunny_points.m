% data from http://graphics.stanford.edu/data/3Dscanrep/
% 
ptCloud = pcread('bun_zipper.ply');
pcshow(ptCloud)
points = ptCloud.Location;
save('./bunny reconstruction/bunny_points','points')
