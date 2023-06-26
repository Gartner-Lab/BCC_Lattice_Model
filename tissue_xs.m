%% [points, edges] = tissue_xs(tissue, dim, n)
% 8 February 2021
% Jennifer Hu
%
%   Given a tissue, the dimension to slice, and the number of layers,
%   returns the sets of points and edges that correspond to
%   the cross-section of the tissue taken at and around tissue.radius.
%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

function [points, edges] = tissue_xs(tissue, dim, n)
    % make 2D cross section by excluding all edges with points outside XS
    points = tissue.p(:,dim) == tissue.radius; 
    if (n > 1)
        for i = 2:n
            height = floor(i/2); if mod(i,2) == 1, height = -height; end
            points = points | (tissue.p(:,dim) == tissue.radius + height);
        end
    end
    % valid edges are only those between points both within the XS
    edges = tissue.edges.all & (points & points');
end