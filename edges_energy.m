%% edges_energy(edges, const)
% 3 February 2021
% Jennifer Hu
%
%   Calculates tissue energy given a tissue's edges and a set of energies
%   whose elements are ordered by the values in const.
%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

function [energy] = edges_energy(edges, const)
    % use only the upper triangle of the areas matrix
    energy = triu(edges.areas) .* ...
    	(const.IE(const.LL) * edges.types(:,:,const.LL) + ...
		const.IE(const.LM) * edges.types(:,:,const.LM) + ...
		const.IE(const.LX) * edges.types(:,:,const.LX) + ...
		const.IE(const.MM) * edges.types(:,:,const.MM) + ...
		const.IE(const.MX) * edges.types(:,:,const.MX));
end