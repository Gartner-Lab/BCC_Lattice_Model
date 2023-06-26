
function [fighandle] = plot_tissue(tissue)
c = tissue.const;
colorL = [255, 185, 15] / 255; colorM = [100, 30, 100] / 255;

% make a Delaunay triangulation from the points
DT = delaunayTriangulation(tissue.p);
% they will be made in the same order so nicely labeled by tissue.is
idxC = find(tissue.is(:,c.C));
% compute Voronoi diagram
[vtx, regions] = voronoiDiagram(DT);

% Plot the whole tissue
colors = colorM .* tissue.is(:,c.M) + colorL .* tissue.is(:,c.L);
figure(1); plot_idx3D(idxC, colors, vtx, regions);
fighandle = gcf;

end


