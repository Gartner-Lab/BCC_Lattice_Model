%% [tissue] = reset_tissue(tissue)
% Instead of making a fresh new one, just reassign all the cell identities
% in an existing tissue. This saves time in the connectivity calculation.
% Use the same number of LEP and MEP as in tissue.n by default, or pass a
% second argument specifying LF.
function [tissue] = reset_tissue(tissue, LF, targetEL)
    c = tissue.const;
    Cidx = find(tissue.is(:,c.C)); 
    if nargin == 1 % keep the same number of LEPs
        nL = tissue.n(c.L);
    else % change the number of LEPs
        assert(LF <= 1 && LF >= 0, 'LEP fraction must be in [0,1].');
        nL = round(LF * tissue.n(c.C));
        if tissue.n(:,c.L) ~= nL
            tissue.n(:,c.L) = nL;
            tissue.n(:,c.M) = tissue.n(c.C)-nL;
        end
    end
    if nargin == 3 % pick nEL and nL-nEL LEP in boundayr and core
        % assign based on targetEL
        nC = tissue.n(c.C);
        nB = tissue.n(c.B);
        nBL = round(targetEL*nB);
        % check that 0 < #core LEP < #core cell
        assert(nL-nBL <= nC-nB && nL-nBL >= 0, 'Invalid target EdgeLEP fraction.');
        Bidx = find(tissue.is(:,c.B));
        BLidx = Bidx(randperm(nB, nBL));
        core_idx = find(tissue.is(:,c.C) & ~tissue.is(:,c.B));
        coreL_idx = core_idx(randperm(nC-nB, nL-nBL));
        tissue.is(:,c.L) = 0; tissue.is(BLidx,c.L) = 1; tissue.is(coreL_idx,c.L) = 1;
        tissue.is(:,c.M) = tissue.is(:,c.C) & ~tissue.is(:,c.L);
        assert(sum(tissue.is(:,c.L)) == nL && sum(tissue.is(:,c.M)) == tissue.n(c.M));
        tissue.is(:,c.EL) = tissue.is(:,c.L) & tissue.is(:,c.B);
    else % pick n.L random cells to assign as LEP, rest are MEP
        Lidx = Cidx(randperm(tissue.n(c.C), nL));
        tissue.is(:,c.L) = 0; tissue.is(Lidx,c.L) = 1;
        tissue.is(:,c.M) = tissue.is(:,c.C) & ~tissue.is(:,c.L);
        assert(sum(tissue.is(:,c.L)) == nL && sum(tissue.is(:,c.M)) == tissue.n(c.M));
        tissue.is(:,c.EL) = tissue.is(:,c.L) & tissue.is(:,c.B);
    end
    % recalculate edges
    tissue.edges.types = edgetypes(tissue.is, tissue.edges.all, c);
end