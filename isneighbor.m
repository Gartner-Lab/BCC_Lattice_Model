%% 3 February 2021
% Jennifer Hu
% isneighbor.m
%
%   Given an adjacency matrix and a px id(s), returns the px ids of all
%   neighboring points as a column vector of booleans. The input idx can be
%   either the numeric indices or a vector of booleans.
%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

function [is] = isneighbor(edges, idx)
    is = any(edges(:, idx), 2);
end