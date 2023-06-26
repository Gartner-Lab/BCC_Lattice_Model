%% 3 February 2021
% Jennifer Hu
% edgetypes.m
%
%   Given is (point assignments), edges (adjacency matrix), and const 
%   with type indices, returns adjacency matrix ordered by interface type.
%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

function [types] = edgetypes(is, edges, const)
	types = zeros([size(edges),5]);
	% LL, LM, LX, MM, MX (and the reverse) by truth tables. reverse by switching the row/column
	types(:,:,const.LL) = edges & (is(:,const.L)' & is(:,const.L));
	types(:,:,const.LM) = edges & ((is(:,const.L)' & is(:,const.M)) | (is(:,const.L) & is(:,const.M)'));
	types(:,:,const.LX) = edges & ((is(:,const.L)' & is(:,const.X)) | (is(:,const.L) & is(:,const.X)'));
	types(:,:,const.MM) = edges & (is(:,const.M)' & is(:,const.M));
	types(:,:,const.MX) = edges & ((is(:,const.M)' & is(:,const.X)) | (is(:,const.M) & is(:,const.X)'));
	assert(~any(sum(types,3) > 1, 'all'), 'Edge types should not overlap.');
    types = logical(types);
end