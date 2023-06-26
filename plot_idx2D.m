%% [] = plot_idx2D(tissue)
% 23 March 2023
% Vasudha Srivastava
%
%   Plots 2D Voronoi diagram. Required inputs:
%   - tissue: tissue object (2D hexagonal lattice typically)
%   - colors: [nx3] matrix of color assignments to all points
%  
%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

function [] = plot_idx2D(tissue)
    c = tissue.const;
    colorL = [255, 185, 15] / 255; colorM = [100, 30, 100] / 255;

    % make a Delaunay triangulation from the points
    DT = delaunayTriangulation(tissue.p);
    % they will be made in the same order so nicely labeled by tissue.is
    idxC = find(tissue.is(:,c.C));
    % compute Voronoi diagram
    [vtx, regions] = voronoiDiagram(DT);
    colors = colorM .* tissue.is(:,c.M) + colorL .* tissue.is(:,c.L);
    figure(1);
    clf; pbaspect([1 1 1]); daspect([1 1 1]);
    set(gcf,'position',[20,20,100,100])
    ax = gca; ax.XColor = 'none'; ax.YColor = 'none';
    hold on;
    for i = 1:length(idxC)
        id = idxC(i);
        pts = vtx(regions{id},:);
        region = convhull(pts);
        color = colors(id,:);
        p = fill(pts(:,1),pts(:,2),color) ;
        p.LineWidth = 0.25;
        p.EdgeColor = [1 1 1];
        %fill(v(c{i},1),v(c{i},2),char(randsample(color,1))) ;
    end
    hold off;
    fighandle = gcf;
end