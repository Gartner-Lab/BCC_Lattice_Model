%% 23 March 2023
% Vasudha Srivastava
%
%   Make a 2D hexagonal lattice. Return a tissue struct containing
%   - p, [nx2] matrix of coordinates px, py
%   - edges, struct containing
%       - all, boolean of edges that exist
%       - types, [nxnx5] ordered by interface type (set in const)
%       - interface lengths, in units of edge lengths of hexagon.
%   - is, [nx5] matrix of point types (columns set in const)
%   - n, [1x5] matrix of point counts
%   - gamma, to fill in later with interfacial tensions
%   - E, to fill in later with current tissue energy
%   - const
%   - radius
%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

function [tissue] = make_tissue_2Dhex(radius, LF, const)

%% Set up
% Radius is in the units of edge lengths. 
% Make grid of radius R+1 (extra boundary cell layer.
% Use odd values of R.
diameter = round(radius * 2 + 2); 
C = [round(radius + 1)*1.5,round(radius + 3)* sqrt(3)/2];
center = repmat(C, diameter*2*diameter, 1);

fprintf('Radius %g... ', radius);
% For hexagons with flat-top orientations, width is 2*a and height is
% sqrt(3)*a. 
% Horizontal spacing between cell centers is 3/2*a.
% vertical spacing between cell centers is sqrt(3)*a.
% Odd rows are offset by sqrt(3)/2*a

xys = zeros(diameter*2*diameter, 2);
i = 1;
for row = 1:diameter*2
    for col = 1:diameter
        xys(i,1) = row * 1.5;
        xys(i,2) = (col + mod(row, 2)*0.5)*sqrt(3);
        i = i + 1;
    end
end
px = xys(:,1); py = xys(:,2); 
% remove points that are more than radius away from center
on = sqrt(sum(([px, py] - center).^2, 2)) <= radius+1;
px = px(on); py = py(on);
%scatter(px, py);
        
fprintf('%d points: ', length(px));

%% Use pairwise distances to find edges
% distances between each pt and every other in x, y, z *separately*
dx = pdist(px); dy = pdist(py);
% edge exists if distance is sqrt(3)
edgenear = squareform(...
    (dx == 0 & round(dy,3) == round(sqrt(3),3) ) | ...
    (dx == 1.5 & round(dy,3) == round(sqrt(3)*0.5,3) ) );

assert(all(sum(edgenear, 2) <= 6), 'At most 6 within edges for any point.');

%% Assign cell and ECM identities
isC = ones(length(px),1);
% all edges 
edges = edgenear;
edgecounts = sum(edges, 2);
assert(~any(edgecounts > 6), '6 is the maximum connectivity in hexagonal lattice.');
% ECM points have fewer than max edges. all others are cell points
isC(edgecounts < 6) = 0; isX = ~isC;
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

edgeareas = (edgenear * 100);
% multiplied by 100 um^2, assuming edge length and height of 10 um

tissue = struct('p', [px py], ...
    'edges', struct('all', edges, 'areas', edgeareas, ...
                    'types', edgetypes(is, edges, const)), ...
    'is', is, 'n', n, 'radius', radius, ...
    'gammas', gammas, 'E', energy, ...
    'const', const, 'nEL',nEL);