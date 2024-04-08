//This file is part of Directional, a library for directional field processing.
// Copyright (C) 2017 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.

#ifndef DIRECTIONAL_EFFORT_TO_INDICES_H
#define DIRECTIONAL_EFFORT_TO_INDICES_H

#include <vector>
#include <cmath>
#include <Eigen/Core>
#include <directional/definitions.h>
#include <directional/dual_cycles.h>
#include <directional/CartesianField.h>


namespace directional
{
    // Computes cycle-based indices from adjaced-space efforts of a directional field.
    // Note: input is effort (sum of rotation angles), and not individual rotation angles
    // Input:
    //  basisCycles:    #c by #iE (inner edges of the mesh) the oriented basis cycles around which the indices are measured
    //  effort:         #iE the effort (sum of rotation angles) of matched vectors across the dual edge. Equal to N*rotation angles for N-RoSy fields.
    //  cycleCurvature: #c the cycle curvature (for instance, from directional::dual_cycles)
    //  N:              The degree of the field
    // Output:
    //  indices:     #c the index of the cycle x N (always an integer).
    inline void effort_to_indices(const Eigen::SparseMatrix<double>& basisCycles,
                                      const Eigen::VectorXd& effort,
                                      const Eigen::VectorXd& cycleCurvature,
                                      const int N,
                                      Eigen::VectorXi& indices)
    {
        using namespace std;
        Eigen::VectorXd dIndices = ((basisCycles * effort + (double)N*cycleCurvature).array() / (2.0*directional::PI));  //this should already be an integer up to numerical precision

        indices.conservativeResize(dIndices.size());
        for (int i=0;i<indices.size();i++) {
            assert("Indices are not naturally integer!" && fabs(std::round(dIndices(i))-dIndices(i))<1e-6);
            indices(i) = std::round(dIndices(i));
        }

    }


    // version that accepts a cartesian field object and operates on it as input and output.
    inline void effort_to_indices(directional::CartesianField& field)
    {
        //field.effort = Eigen::VectorXd::Zero(field.adjSpaces.rows());
        Eigen::VectorXd effortInner(field.tb->innerAdjacencies.size());
        for (int i=0;i<field.tb->innerAdjacencies.size();i++)
            effortInner(i)=field.effort(field.tb->innerAdjacencies(i));
        Eigen::VectorXi fullIndices;
        directional::effort_to_indices(field.tb->cycles, effortInner, field.tb->cycleCurvatures, field.N, fullIndices);

        Eigen::VectorXi indices(field.tb->local2Cycle.size());
        for (int i=0;i<field.tb->local2Cycle.size();i++)
            indices(i)=fullIndices(field.tb->local2Cycle(i));

        std::vector<int> singCyclesList;
        std::vector<int> singIndicesList;
        for (int i=0;i<field.tb->local2Cycle.size();i++)
            if (indices(i)!=0){
                singCyclesList.push_back(i);
                singIndicesList.push_back(indices(i));
            }

        Eigen::VectorXi singCycles(singCyclesList.size());
        Eigen::VectorXi singIndices(singIndicesList.size());
        for (int i=0;i<singCyclesList.size();i++){
            singCycles(i)=singCyclesList[i];
            singIndices(i)=singIndicesList[i];
        }
        field.set_singularities(singCycles, singIndices);
    }
}

#endif


