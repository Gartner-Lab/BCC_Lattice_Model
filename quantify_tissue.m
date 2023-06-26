%% [q] = quantify_tissue(tissue)
% 4 February 2021
% Jennifer Hu
%
%   Given a tissue, returns 13x16 quantities or 1x10 depending on andXS
%   Interface types will be sorted according to their order in const. 
%   Returns metrics for cross-sections taken at and around radius height.
%       [XSdim, XSlayers, XSn, XSnC, XSnL, XSnB, XSnBL, areas(x10)]
% 16 March 2022
% Vasudha Srivastava
% Also return the number of edge LEP
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

function [q] = quantify_tissue(tissue, andXS, is2D)
    if nargin == 1 | nargin == 3
        andXS = false;
    end
    if nargin<3
        is2d = false;
    end
    c = tissue.const;
    interfaces = [c.LL, c.LM, c.LX, c.MM, c.MX];
    if ~is2D
        q = zeros(1+3*4, 7+10);
    else
        q = zeros(1+3*4, 7+5);
    end
    isEL = tissue.is(:,c.EL);
    q(1,7) = sum(isEL);
    if ~is2D
        q(1, 8:end) = cell2mat(arrayfun(@(i) ...
            [sum(triu(tissue.edges.sq,1) .* tissue.edges.types(:,:,i),'all'), ...
             sum(triu(tissue.edges.hex,1) .* tissue.edges.types(:,:,i),'all')], ...
             interfaces, 'UniformOutput', false));
    else
        q(1, 8:end) = cell2mat(arrayfun(@(i) ...
            [sum(triu(tissue.edges.areas,1) .* tissue.edges.types(:,:,i),'all')],... 
             interfaces, 'UniformOutput', false));
    end
    if andXS
        % make 2D cross section by excluding all edges with points outside XS
        row = 2;
        for dim = 1:3
            for layers = 1:4
                [points, edges] = tissue_xs(tissue, dim, layers);
                nC = sum(points & tissue.is(:,c.C));
                L = points & tissue.is(:,c.L); 
                B = points & tissue.is(:,c.B); BL = B & L;
                qXS = [dim, layers, nC, sum(L), sum(B), sum(BL),c.EL, ... % metrics about the cross sections
                    cell2mat(arrayfun(@(i) ...
                    [sum(triu(tissue.edges.sq,1) .* edges .* tissue.edges.types(:,:,i),'all'), ...
                     sum(triu(tissue.edges.hex,1) .* edges .* tissue.edges.types(:,:,i),'all')], ...
                     interfaces, 'UniformOutput', false))];
                q(row,:) = qXS; row = row+1;
            end
        end
    else % dump those extra columns/rows
        q = q(1,7:end);
    end
end