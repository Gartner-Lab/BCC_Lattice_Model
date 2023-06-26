%% March 16, 2022
% Vasudha Srivastava
% - Edit Jennifer's code to have a single cell at the tissue center. 
%   Hence, choose odd value of tissue radius. 
% - As the edge sites are assigned as ECM, make Diameter = 2*Radius+2 and
% Center = Radius+1
% - A truncated octahedron has 6 square faces (nearest neighbor) and 8
% hexagonal faces (next nearest neighbor) (JH used 12??)
% - Export the number of boundary sites
% - Add option to set a target EL
%% 3 February 2021
% Jennifer Hu
%
%   Make a truncated octahedral foam (Voronoi tesselation of BCC lattice)
%   in 3D, but only with cubic grid. Return a tissue struct containing
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

function [tissue] = make_tissue_cube(radius, LF, const, targetEL)

%% Set up
% In the BCC with cell radius = 1 unit, all cell centers will have odd coordinates. 
% The cells on the edge are ECM. Hence, we initialize a lattice with the
% diameter = 2*(radius+2). Hence the tissue center is at radius+2.

diameter = round(radius*2+4); center = repmat(radius+2, 1,3); %%Why is center not diameter/2???
fprintf('Radius %g... ', radius);
odds = 1:2:(diameter+1);
% coordinates of every point (very fast b/c vectorized). linearize after.
[xs, ys, zs] = meshgrid(odds);
px = xs(:); py = ys(:); pz = zs(:);

% remove points that are more than radius away from center
on = sqrt(sum(([px, py, pz] - center).^2, 2)) <= radius+2;
px = px(on); py = py(on); pz = pz(on);
fprintf('%d points: ', length(px));

%% Use pairwise distances to find edges
% distances between each pt and every other in x, y, z *separately*
dx = pdist(px); dy = pdist(py); dz = pdist(pz);
% edge exists if 2 dists are 0 and 1 is 2
edgenear = squareform(...
    (dx == 0 & dy == 0 & dz == 2) | ...
    (dx == 2 & dy == 0 & dz == 0) | ...
    (dx == 0 & dy == 2 & dz == 0));
% adding the 12 next-nearest neighbors
edgenext = squareform(...
    (dx == 2 & dy == 2 & dz == 2) );
% ensure matches BCC connectivity (or less for edge points)
assert(all(sum(edgenext, 2) <= 8), 'At most 8 next-nearest edges for any point.');
assert(all(sum(edgenear, 2) <= 6), 'At most 6 within edges for any point.');

%% Assign cell and ECM identities
isC = ones(length(px),1);
% all edges 
edges = edgenear | edgenext; edgecounts = sum(edges, 2);
assert(~any(edgecounts > 14), '14 is the maximum connectivity in cube lattice.');
% ECM points have fewer than max edges. all others are cell points
isC(edgecounts < 14) = 0; isX = ~isC;
% remove edges between ECM points
edgeXX = isX' & isX; edges(edgeXX) = 0;
edgenext(~edges) = 0; edgenear(~edges) = 0;
% every point attached to an ECM (by design none are ECM)
isB = isneighbor(edgenear, isX); % B are cells with >1 X as nearest neighbor

nC = sum(isC); nL = round(LF*nC); nB = sum(isB);
% Either assign LEP identity randomly or based on a target boundary LEP
% occupancy
if(nargin==4) 
    % assign based on targetEL
    nBL = round(targetEL*nB);
    % check that 0< #core LEP < #core cell
    assert(nL-nBL <= nC-nB && nL-nBL >= 0, 'Invalid target EdgeLEP fraction.');
    Bidx = find(isB);
    BLidx = Bidx(randperm(nB, nBL));
    core_idx = find(isC & ~isB);
    coreL_idx = core_idx(randperm(nC-nB, nL-nBL));
    isL = zeros(size(isC)); isL(BLidx) = 1; isL(coreL_idx) = 1; isM = isC & ~isL;
    assert(sum(isL) == nL && sum(isM) == nC-nL);
else
    % assign cell ID randomly based on LEP fraction LF
    % pick nL random cells
    Cidx = find(isC);
    Lidx = Cidx(randperm(nC, nL));
    isL = zeros(size(isC)); isL(Lidx) = 1; isM = isC & ~isL;
    assert(sum(isL) == nL && sum(isM) == nC-nL);
end

%% All the important information in the struct
% compile the point info into the 'is' matrix
is = zeros(length(px),6);
is(:,const.C) = isC;
is(:,const.L) = isL;
is(:,const.M) = isM;
is(:,const.X) = isX;
is(:,const.B) = isB;
is(:,const.EL) = isB & isL;
is = logical(is);
isEL = isB & isL;
nEL = sum(isEL);
n = zeros(1,6);
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
edgeareas = (edgenext * (pi/(6*sqrt(2)))^(2/3) + ...
    edgenear * 3*sqrt(3)/2 * (pi/(6*sqrt(2)))^(2/3)) * 100;
% multiplied by 100 um^2, assuming cell radius of 10 um

tissue = struct('p', [px py pz], ...
    'edges', struct('all', edges, 'areas', edgeareas, ...
                    'sq', edgenext, 'hex', edgenear, ...
                    'types', edgetypes(is, edges, const)), ...
    'is', is, 'n', n, 'radius', radius, ...
    'gammas', gammas, 'E', energy, ...
    'const', const, 'nEL',nEL);