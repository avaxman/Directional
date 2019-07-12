#ifndef DIRECTIONAL_VERTEXFIELDSUBDIVIDER_H
#define DIRECTIONAL_VERTEXFIELDSUBDIVIDER_H

#include <directional/Subdivision/SubdivisionBuilder.h>

template<typename CoefficientProvider>
class VertexFieldSubdivider
{
	LeveledSparseConstructor builder;
	CoefficientProvider cf;
	EdgeData* ED;
	Eigen::VectorXi seenVs;
public:
	Eigen::SparseMatrix<double> getMatrix() const
	{
		return builder.matrix;
	}
	VertexFieldSubdivider() : cf({}) {}

	void setup(EdgeData& ED)
	{
		builder = LeveledSparseConstructor(ED.vertexCount(), ED.vertexCount());
		builder.makeId();
		this->ED = &ED;
	}

	template<typename...T>
	void prepareNext(SubdivisionBuilder<T...>& subdivider)
	{
		builder.cols = builder.rows;
		builder.rows += subdivider.ED.edgeCount();
		seenVs.setConstant(builder.rows, 0);
	}

	void handleBoundaryRing(const std::vector<int>& edges, const std::vector<int>& edgeOrients)
	{
		VertFinder vfind(edges, edgeOrients, ED);
		//Central vertex in loop
		const int vCentral = vfind.vertex(0);

		const int valence = vfind.boundaryValence();
		const int vCount = builder.cols;

		std::vector<int> inds;
		std::vector<double> coeffs;
		cf.getEvenBoundaryStencil(valence, inds, coeffs);
		for (int i = 0; i < inds.size(); i++) builder.addCoeff(vCentral, vfind.vertex(inds[i]), coeffs[i]);
		seenVs(vfind.vertex(0)) += 1;
		// Odd vertices around boundary vertex
		for(int v = 1 ; v <= valence; v++)
		{
			// Avoid duplicating stencils here
			if (vfind.vertexEdgeOrient(v) == 1) continue;
			
			inds.clear();
			coeffs.clear();
			cf.getOddBoundaryStencil(valence, v, inds, coeffs);
			const int target = vCount + vfind.vertexEdge(v);
			seenVs(target) += 1;
			for(int i = 0; i < inds.size(); i++)
			{
				builder.addCoeff(target, vfind.vertex(inds[i]), coeffs[i]);
			}
		}
	}

	void finalize()
	{
		builder.finalize();

		for(int i = 0; i < seenVs.rows(); i++)
		{
			DIR_ASSERT(seenVs(i) == 1);
		}
	}

	void handleRegularRing(const std::vector<int>& edges, const std::vector<int>& edgeOrients)
	{
		VertFinder vfind(edges, edgeOrients, ED);
		//Central vertex in loop
		const int vCentral = vfind.vertex(0);

		const int valence = vfind.boundaryValence();

		const int vCount = builder.cols;

		std::vector<int> inds;
		std::vector<double> coeffs;
		cf.getEvenRegularStencil(valence, inds, coeffs);
		for (int i = 0; i < inds.size(); i++) builder.addCoeff(vCentral, vfind.vertex(inds[i]), coeffs[i]);
		seenVs(vfind.vertex(0)) += 1;
		// Odd vertices around boundary vertex
		for (int v = 1; v <= valence; v++)
		{
			// Only handle odd element once
			if (vfind.vertexEdgeOrient(v) == 1) continue;

			inds.clear();
			coeffs.clear();
			cf.getOddRegularStencil(valence, v, inds, coeffs);
			const int target = vCount + vfind.vertexEdge(v);
			seenVs(target) += 1;
			for (int i = 0; i < inds.size(); i++)
			{
				builder.addCoeff(target, vfind.vertex(inds[i]), coeffs[i]);
			}
		}
	}
};
#endif