%% evolve_tissue(tissue, temperature)
% Uses the acceptance function 1/(1+e^(d/T)).
function [tissue] = evolve_tissue(tissue, temperature)
    c = tissue.const;
    % randomly pick one of the cells to swap with a neighbor
    rn = randi(tissue.n(c.C)); i1 = find(tissue.is(:,c.C)); i1 = i1(rn);
    % all its neighbors that are cells
    is_neighbor = tissue.edges.all(:, i1) & tissue.is(:,c.C);
    n_neighbors = sum(is_neighbor);
    rn = randi(n_neighbors); i2 = find(is_neighbor); i2 = i2(rn);
    assert(all(tissue.is([i1,i2],c.C)), 'Should be cells.');
    % if the cell types were different this vector would be [1 1]
    if ~all(tissue.is(i1,[c.M, c.L]) + tissue.is(i2,[c.M, c.L]) == 1)
        % if both same cell type, skip
        return
    end
    % create the new swapped tissue by changing tissue.is
    new_tissue = tissue;
    new_tissue.is(i1,[c.M, c.L]) = tissue.is(i2,[c.M, c.L]);
    new_tissue.is(i2,[c.M, c.L]) = tissue.is(i1,[c.M, c.L]);
    % update edge types:
    newtypes = tissue.edges.types;
    if tissue.is(i1, c.L)
        iL = i1; iM = i2;
    else
        iL = i2; iM = i1;
    end
    % LL becomes LM, LM becomes MM, LX becomes MX
    LL = tissue.edges.types(iL,:,c.LL);
    newtypes(iL,LL,c.LL) = 0; newtypes(LL,iL,c.LL) = 0;
    newtypes(iL,LL,c.LM) = 1; newtypes(LL,iL,c.LM) = 1;
    LM = tissue.edges.types(iL,:,c.LM);
    LM(iM) = 0; % except the cell being swapped with!
    newtypes(iL,LM,c.LM) = 0; newtypes(LM,iL,c.LM) = 0;
    newtypes(iL,LM,c.MM) = 1; newtypes(LM,iL,c.MM) = 1;
    LX = tissue.edges.types(iL,:,c.LX);
    newtypes(iL,LX,c.LX) = 0; newtypes(LX,iL,c.LX) = 0;
    newtypes(iL,LX,c.MX) = 1; newtypes(LX,iL,c.MX) = 1;
    % MM becomes LM, LM becomes LL, MX becomes LX
    MM = tissue.edges.types(iM,:,c.MM);
    newtypes(iM,MM,c.MM) = 0; newtypes(MM,iM,c.MM) = 0;
    newtypes(iM,MM,c.LM) = 1; newtypes(MM,iM,c.LM) = 1;
    ML = tissue.edges.types(iM,:,c.LM); ML(iL) = 0;
    newtypes(iM,ML,c.LM) = 0; newtypes(ML,iM,c.LM) = 0;
    newtypes(iM,ML,c.LL) = 1; newtypes(ML,iM,c.LL) = 1;
    MX = tissue.edges.types(iM,:,c.MX);
    newtypes(iM,MX,c.MX) = 0; newtypes(MX,iM,c.MX) = 0;
    newtypes(iM,MX,c.LX) = 1; newtypes(MX,iM,c.LX) = 1;
    % quickly calculate energy change instead of redoing whole tissue E
    deltaE = 100 * ( ... don't forget the 100 um^2
        sum(tissue.edges.areas(iL,LL))*(tissue.gammas(c.LM)-tissue.gammas(c.LL)) + ...
        sum(tissue.edges.areas(iL,LM))*(tissue.gammas(c.MM)-tissue.gammas(c.LM)) + ...
        sum(tissue.edges.areas(iL,LX))*(tissue.gammas(c.MX)-tissue.gammas(c.LX)) + ...
        sum(tissue.edges.areas(iM,MM))*(tissue.gammas(c.LM)-tissue.gammas(c.MM)) + ...
        sum(tissue.edges.areas(iM,ML))*(tissue.gammas(c.LL)-tissue.gammas(c.LM)) + ...
        sum(tissue.edges.areas(iM,MX))*(tissue.gammas(c.LX)-tissue.gammas(c.MX)));
    % if new energy is higher
    if (deltaE > 0)
        % acceptance function is between 0 and 0.5
        p_accept = 1/(1+exp(deltaE/temperature));
        % do not swap; return
        if (rand() > p_accept)
            return
        end
    end
    % calculate new energy
    new_tissue.edges.types = newtypes;
    new_tissue.E = tissue.E + deltaE; % tissue_energy(new_tissue);
    % assert(abs(new_tissue.E - tissue.E - deltaE) < 1e-7, 'Delta E calculation check (to 6 decimal places).');
    % swap the tissues
    tissue = new_tissue;
end