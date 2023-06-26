%% scripts

%datadir = find_datadir();
datadir = fullfile('/Users','Vasudha','Library','CloudStorage','Box-Box', ...
                    'Gartnerlab Data','Individual Folders','Jennifer Hu', ...
                    'Analysis','DAH Modeling','DAH Data');
%datadir = '/Users/Vasudha/Library/CloudStorage/Box-Box/Gartnerlab Data/Individual Folders/Jennifer Hu/Analysis/DAH Modeling/DAH Data/BCC/tissues/'
% small test tissue
tissuename = 'r6_BL30_k60';
figname = fullfile(datadir, 'BCC', tissuename);
load(fullfile(datadir, 'BCC', 'tissues', [tissuename, '.mat']));
const = struct('L',1,'M',2,'C',3,'B',4,'X',5, ...
    'LL',1,'LM',2,'LX',3,'MM',4,'MX',5);
colorL = [255, 185, 15] / 255; colorM = [100, 30, 100] / 255;
isM = tissue.is(:,const.M); isL = tissue.is(:,const.L);

%% Show the centers of each point and the edges
G = graph(tissue.edges.all); 
npoints = size(tissue.p,1);
Marker = repmat("o",npoints,1);
Marker(tissue.is(:,const.B)) = "^"; 
Marker(tissue.is(:,const.X)) = "+";
NodeColor = zeros(npoints,3);
NodeColor(isM,:) = repmat(colorM, sum(isM),1);
NodeColor(isL,:) = repmat(colorL, sum(isL),1);

plot(G, 'XData',tissue.p(:,1),'YData',tissue.p(:,2),'ZData',tissue.p(:,3), ...
    'LineWidth',2,'MarkerSize',10, ...
    'NodeColor', NodeColor, ...
    'Marker', Marker)

%% Make a 3D Voronoi diagram
% make a Delaunay triangulation from the points
DT = delaunayTriangulation(tissue.p);
% they will be made in the same order so nicely labeled by tissue.is
idxC = find(tissue.is(:,const.C)); idxL = find(isL); idxM = find(isM);
% compute Voronoi diagram
[vtx, regions] = voronoiDiagram(DT);

%% Plot the whole tissue
colors = colorM .* isM + colorL .* isL;
figure(1); plot_idx3D(idxC, colors, vtx, regions);
exportgraphics(gcf, [figname,'.pdf'], 'ContentType', 'vector');
% Just LEP and just MEP
figure(2); plot_idx3D(idxM, colors, vtx, regions);
exportgraphics(gcf, [figname,'_M.pdf'], 'ContentType', 'vector');
figure(3); plot_idx3D(idxL, colors, vtx, regions);
exportgraphics(gcf, [figname,'_L.pdf'], 'ContentType', 'vector');
% plot cross sections
for i = 1:3
    [points, edges] = tissue_xs(tissue,3,i);
    idxXS = find(points & tissue.is(:,const.C));
    figure(4); plot_idx3D(idxXS, colors, vtx, regions);
    exportgraphics(gcf, sprintf('%s_XS%d.pdf', figname, i), 'ContentType', 'vector');
end


