% curlFreeParametrization.m Help file for curlFreeParametrization MEX-file.
%  curlFreeParametrization.cpp - Solves for a polyvector field to become more curl free. Only for framefield elements.
% 
%  Takes the input face matrix and produces edge data for the coarse and fine level, along with
%  subdivision operators for all spaces.
% 
%  Function signature:
%  [CutV, CutF, UV, CombedField, EF, SingVertices, SingIndices] = curlFreeParametrization(F, V, RawField)
%       Here the ED structs have keys EF, sFE, EI and F, corresponding to the topology data.
%       The S struct contains the subdivision operators under the keys F, C, E, V and Gamma for the spaces