% exportSubdivisionOperators.m Help file for exportSubdivisionOperators MEX-file.
%  exportSubdivisionOperators.cpp - Exports incremental subdivision operators up to a specified level
% 
%  Takes the input face matrix and produces edge data for the coarse and fine level, along with
%  subdivision operators for all spaces.
% 
%  Function signature:
%  buildSubdivisionOperators(F, level, pathPrefix)
%       Exports three structs to per level, ED0, EDK and S. ED0 and EDK contain topology (EF, EV, sFE, EI) matrices. The S struct has all the subdivision operators (F,E,C,V,Gamma).
%       The gamma subdivision operator acts on Gamma2 elements.