% edgeTopology.m Help file for edgeTopology MEX-file.
%  edgeTopology.cpp - Constructs edge topology data for a mesh
% 
%  Takes the input face matrix and produces edge data for the coarse and fine level, along with
%  subdivision operators for all spaces.
% 
%  Function signature:
%  [ED_struct] = edgeTopology(F, V)
%       Here the ED struct has keys EF, sFE, EI, E and F, corresponding to the topology data.