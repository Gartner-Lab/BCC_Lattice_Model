%% 3 February 2021
% Jennifer Hu
%
%   Make a truncated octahedral foam (Voronoi tesselation of BCC lattice)
%   in 3D. Return a tissue struct containing
%   - p, [nx3] matrix of coordinates px, py, pz
%   - edges, struct containing
%       - all, boolean of edges that exist
%       - types, [nxnx5] ordered by interface type (set in const)
%       - areas, in units of R^2 (typical: 100 um^2)
%   - is, [nx5] matrix of point types (columns set in const)
%   - n, [1x5] matrix of point counts
%   - gamma, to fill in later with interfacial tensions
%   - E, to fill in later with current tissue energy
%   - const
%   - radius
%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

function [tissue] = make_tissue(radius, LF, const)

%% Set up
diameter = round(radius * 2); center = repmat((diameter+1)/2, 1,3);

fprintf('Radius %g... ', radius);

% In the BCC, points have either all odd or all even coordinates
odds = 1:2:diameter; evens = 2:2:diameter;
% coordinates of every point (very fast b/c vectorized). linearize after.
[xo, yo, zo] = meshgrid(odds); 
[xe, ye, ze] = meshgrid(evens);
xo = xo(:); yo = yo(:); zo = zo(:);
xe = xe(:); ye = ye(:); ze = ze(:);
px = [xo; xe]; py = [yo; ye]; pz = [zo; ze];
% remove points that are more than radius away from center
on = sqrt(sum(([px, py, pz] - center).^2, 2)) <= radius;
px = px(on,:); py = py(on,:); pz = pz(on,:);
fprintf('%d points: ', length(px));

%% Use pairwise distances to find edges
% distances between each pt and every other in x, y, z *separately*
dx = pdist(px); dy = pdist(py); dz = pdist(pz);
% edge exists within odd/even groups if 2 dists are 0 and 1 is 2
edgewithin = squareform(...
    (dx == 0 & dy == 0 & dz == 2) | ...
    (dx == 2 & dy == 0 & dz == 0) | ...
    (dx == 0 & dy == 2 & dz == 0));
% edge exists between groups if dz = 1 and dx+dy = 2
edgebetween = squareform((dz == 1 & dx+dy == 2));
% ensure matches BCC connectivity (or less for edge points)
assert(all(sum(edgebetween, 2) <= 8), 'At most 8 between edges for any point.');
assert(all(sum(edgewithin, 2) <= 6), 'At most 6 within edges for any point.');

%% Assign cell and ECM identities
isC = ones(length(px),1);
% all edges 
edges = edgewithin | edgebetween; edgecounts = sum(edges, 2);
assert(~any(edgecounts > 14), '14 is the maximum connectivity in BCC.');
% ECM points have fewer than max edges. all others are cell points
isC(edgecounts < 14) = 0; isX = ~isC;
% remove edges between ECM points
edgeXX = isX' & isX; edges(edgeXX) = 0;
edgebetween(~edges) = 0; edgewithin(~edges) = 0;
% every point attached to an ECM (by design none are ECM)
isB = isneighbor(edges, isX);
% assign cell ID randomly based on LEP fraction LF
nC = sum(isC); nL = round(LF*nC);
% pick nL random cells
Cidx = find(isC);
Lidx = Cidx(randperm(nC, nL));
isL = zeros(size(isC)); isL(Lidx) = 1; isM = isC & ~isL;
assert(sum(isL) == nL && sum(isM) == nC-nL);

%% All the important information in the struct
% compile the point info into the 'is' matrix
is = zeros(length(px),5);
is(:,const.C) = isC;
is(:,const.L) = isL;
is(:,const.M) = isM;
is(:,const.X) = isX;
is(:,const.B) = isB;
is = logical(is);

n = zeros(1,5);
n(const.C) = nC;
n(const.L) = nL;
n(const.M) = sum(isM);
n(const.X) = sum(isX);
n(const.B) = sum(isB);
fprintf('%d ECM, %d cells (%d boundary, %d core)\n', ...
    n(const.X), n(const.C), n(const.B), n(const.C)-n(const.B));

gammas = zeros(1,5);
energy = 0;

% Weight the edges according to their area:
% In units of cell-radius R, faces have sides of length a,
% where a^3 = pi/(6*sqrt(2)) * R^3. a^2 = (pi/(6*sqrt(2)))^(2/3)
% - squares (within groups) have area a^2 = (pi/(6*sqrt(2)))^(2/3) * R^2
% - hexagons (between) have area 3*sqrt(3)/2 * (pi/(6*sqrt(2)))^(2/3) * R^2
edgeareas = (edgewithin * (pi/(6*sqrt(2)))^(2/3) + ...
    edgebetween * 3*sqrt(3)/2 * (pi/(6*sqrt(2)))^(2/3)) * 100;
% multiplied by 100 um^2, assuming cell radius of 10 um

tissue = struct('p', [px py pz], ...
    'edges', struct('all', edges, 'areas', edgeareas, ...
                    'sq', edgewithin, 'hex', edgebetween, ...
                    'types', edgetypes(is, edges, const)), ...
    'is', is, 'n', n, 'radius', radius, ...
    'gammas', gammas, 'E', energy, ...
    'const', const);