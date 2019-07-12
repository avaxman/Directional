% findMatchingAndSings.m Help file for findMatchingAndSings MEX-file.
%  findMatchingAndSings.cpp - Finds a matching and singularities with indices by using principal matching
% 
%  Function signature:
%  [Matching, SingularityVertices, SingularityIndices] = curlFreeParametrization(F, V, RawField)
% Here matching is the left to right face matching per edge. SingularityVertices is a list of vertex indices that are similar, 
% where SingularityIndices gives the corresponding index.