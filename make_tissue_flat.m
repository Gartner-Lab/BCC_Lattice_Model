%% 25 May 2021
% Jennifer Hu
%
%   Make a truncated octahedral foam (Voronoi tesselation of BCC lattice)
%   in 3D, but shaped like a bilayer sheet. Return a tissue struct containing
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

function [tissue] = make_tissue_flat(radius, LF, const)

%% Set up
diameter = round(radius * 2); center = repmat((diameter+1)/2, 1,2);

fprintf('Radius %g... ', radius);

% In the BCC, points have either all odd or all even coordinates. we'll use
% all odds for this
xys = 1:2:(diameter + 1); zs = 1:2:5; 
% coordinates of every point (very fast b/c vectorized). linearize after.
[xs, ys, zs] = meshgrid(xys, xys, zs);
px = xs(:); py = ys(:); pz = zs(:);
% remove those with xy distance > radius from center
on = sqrt(sum(([px, py] - center).^2, 2)) <= radius;
px = px(on); py = py(on); pz = pz(on);
fprintf('%d points: ', length(px));

% Use pairwise distances to find edges
% distances between each pt and every other in x, y, z *separately*
dx = pdist(px); dy = pdist(py); dz = pdist(pz);
% edge exists within odd/even groups if 2 dists are 0 and 1 is 2
edgenear = squareform(...
    (dx == 0 & dy == 0 & dz == 2) | ...
    (dx == 2 & dy == 0 & dz == 0) | ...
    (dx == 0 & dy == 2 & dz == 0));
% adding the 12 next-nearest neighbors
edgenext = squareform(...
    (dx == 0 & dy == 2 & dz == 2) | ...
    (dx == 2 & dy == 0 & dz == 2) | ...
    (dx == 2 & dy == 2 & dz == 0));
% ensure matches BCC connectivity (or less for edge points)
assert(all(sum(edgenext, 2) <= 12), 'At most 12 next-nearest edges for any point.');
assert(all(sum(edgenear, 2) <= 6), 'At most 6 within edges for any point.');

%% Assign cell and ECM identities
isC = (pz ~= 5); isX = (pz == 5);
% all edges 
edges = edgenear | edgenext; edgecounts = sum(edges, 2);
assert(~any(edgecounts > 18), '18 is the maximum connectivity for cubes.');
% remove edges between ECM points
edgeXX = isX' & isX; edges(edgeXX) = 0; edgenext(~edges) = 0; edgenear(~edges) = 0;
% every point attached to an ECM (by design none are ECM)
isB = isneighbor(edges, isX);
% assign cell ID randomly based on LEP fraction LF
nC = sum(isC); nL = round(LF*nC);
% pick nL random cells
Cidx = find(isC); Lidx = Cidx(randperm(nC, nL));
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

% Weight the edges according to their fake area:
% In units of cell-radius R, faces have sides of length a,
% where a^3 = pi/(6*sqrt(2)) * R^3. a^2 = (pi/(6*sqrt(2)))^(2/3)
% - squares (next-nearest) have area a^2 = (pi/(6*sqrt(2)))^(2/3) * R^2
% - hexagons (nearest) have area 3*sqrt(3)/2 * (pi/(6*sqrt(2)))^(2/3) * R^2
edgeareas = (edgenext * (pi/(6*sqrt(2)))^(2/3) + ...
    edgenear * 3*sqrt(3)/2 * (pi/(6*sqrt(2)))^(2/3)) * 100;
% multiplied by 100 um^2, assuming cell radius of 10 um

tissue = struct('p', [px py pz], ...
    'edges', struct('all', edges, 'areas', edgeareas, ...
                    'sq', edgenext, 'hex', edgenear, ...
                    'types', edgetypes(is, edges, const)), ...
    'is', is, 'n', n, 'radius', radius, ...
    'gammas', gammas, 'E', energy, ...
    'const', const);