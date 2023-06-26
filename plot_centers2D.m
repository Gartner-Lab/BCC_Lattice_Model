%% [] = plot_centers2D(tissue)
% 23 March 2023
% Vasudha Srivastava
%
%   Scatter plot for cell centers. Required inputs:
%   - tissue: tissue object (2D hexagonal lattice typically)
%   - colors: [nx3] matrix of color assignments to all points
%  
%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

function [] = plot_centers2D(tissue)
    c = tissue.const; 
    colorL = [255, 185, 15] / 255; colorM = [100, 30, 100] / 255; colorX = [255, 255, 255] / 255;
    colors = colorM .* tissue.is(:,c.M) + colorL .* tissue.is(:,c.L) + colorX .* tissue.is(:,c.X);
    colors = colors(~tissue.is(:,c.X),:);
    
    figure(1);
    clf; 
    ax = gca; ax.XColor = 'none'; ax.YColor = 'none'; ax.ZColor = 'none';
        ax.Clipping = 'off';
    set(gcf,'position',[20,20,100,100])
    cells = tissue.p(~tissue.is(:,c.X),:);
    scatter(cells(:,1), cells(:,2), 60,colors,"filled");
    pbaspect([1 1 1]); daspect([1 1 1]); 
    axis off;
    fighandle = gcf;
end