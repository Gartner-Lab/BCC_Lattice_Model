%% [] = plot_idx3D(idx, colors, vtx, regions)
% 8 February 2021
% Jennifer Hu
%
%   Plots various subsets of a Voronoi diagram. Required inputs:
%   - idx: numeric indices of points to plot
%   - colors: [nx3] matrix of color assignments to all points
%   - vtx, regions: voronoiDiagram()
%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

function [] = plot_idx3D(idx, colors, vtx, regions)
    clf; pbaspect([1 1 1]); daspect([1 1 1]); grid on;
    ax = gca; ax.XColor = 'none'; ax.YColor = 'none'; ax.ZColor = 'none';
    ax.Clipping = 'off'; hold on;
    for i = 1:length(idx)
        id = idx(i);
        pts = vtx(regions{id},:);
        region = convhull(pts);
        color = colors(id,:);
        trisurf(region, pts(:,1), pts(:,2), pts(:,3), ...
            'FaceColor', color, 'FaceAlpha', 0.6, ...
            'LineStyle', 'none');
    end
    hold off; view([-65,20]); camlight left; camlight right; lighting gouraud
end