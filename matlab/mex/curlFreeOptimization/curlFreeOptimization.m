% curlFreeOptimization.m Help file for curlFreeOptimization MEX-file.
%  curlFreeOptimization.cpp - Solves for a polyvector field to become more curl free. Only for framefield elements.
% 
% 
%  Function signature:
%  [CF_Field, CF_Combed_Field, SingVertices, SingIndices, MaxCurl, MaxCurlCombed] = curlFreeOptimization(F, V, RawField, N, Batches)
%       