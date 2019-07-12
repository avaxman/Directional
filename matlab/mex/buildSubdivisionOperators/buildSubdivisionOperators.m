% buildSubdivisionOperators.m Help file for buildSubdivisionOperators MEX-file.
%  buildSubdivisionOperators.cpp - Builds subdivision operators
% 
%  Takes the input face matrix and produces edge data for the coarse and fine level, along with
%  subdivision operators for all spaces.
% 
%  Function signature:
%  [ED0_Struct, EDK_Struct, S_Struct] = buildSubdivisionOperators(F, level)
%       Here the ED structs have keys EF, sFE, EI and F, corresponding to the topology data.
%       The S struct contains the subdivision operators under the keys F, C, E, V and Gamma for the spaces