function [E] = tissue_energy(tissue)
    c = tissue.const;
    assert(~all(tissue.gammas == 0), 'Energy parameters not loaded.');
    % tissue.gammas is a row vector.
    gammamatrix = zeros(size(tissue.edges.all));
    for interface = [c.LL, c.LM, c.LX, c.MM, c.MX]
        % populate the interface sites with gammas
        gammamatrix(tissue.edges.types(:,:,interface)) = tissue.gammas(interface);
    end
    % triu to avoid double counting
    E = sum(triu(tissue.edges.areas,1) .* gammamatrix, 'all');
end