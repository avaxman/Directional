% buildDirectionalSubdivisionOperators.m Help file for buildDirectionalSubdivisionOperators MEX-file.
%  buildDirectionalSubdivisionOperators.cpp - Builds subdivision operators for directional fields
% 
%  Takes the input face matrix and produces edge data for the coarse and fine level, along with
%  subdivision operators for all spaces.
% 
%  Function signature:
%  [ED0, EDK, SubdivisionStruct, Matching] = buildDirectionalSubdivisionOperators(F, level, matching, N, singularityMat)
%       Here the ED structs have keys EF, sFE, EI and F, corresponding to the topology data.
%       The S struct contains the subdivision operators under the keys F, C, E, V and Gamma for the spaces