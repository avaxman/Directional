#ifndef DIRECTIONAL_SUBDIVIDER_H
#define DIRECTIONAL_SUBDIVIDER_H
#include <vector>
#include <Eigen/Eigen>
#include <directional/iterate_rings.h>
#include <directional/Subdivision/SubdivisionBuilder.h>
/**
 * Check a local ring for a singularity: follow the matching, if the final matched field is not the same, a singularity is present
 */ 
void getSingularity(int vertex, const std::vector<int>& edges, const std::vector<int>& edgeSides, int N, const Eigen::MatrixXi& Matching, std::vector<std::tuple<int,int>>& output){
    // Active ''level'': field that is currently matched.
    int currentLevel = 0;
    for(int e = 2; e < edgeSides.size(); e+=2){
        // Always positive?
        const int match = Matching(edges[e], 1-edgeSides[e]);
        currentLevel += match;
    }
    currentLevel += Matching(edges[0],1-edgeSides[0]);
    currentLevel = currentLevel % N;
    if(currentLevel != 0) output.emplace_back(vertex, currentLevel);
}

void circShift_L(const Eigen::VectorXi& input, unsigned int shift, Eigen::VectorXi& output)
{
	if (shift == 0) {
		output = input;
		return;
	}
	output.resize(input.size());
	output.head(input.size() - shift) = input.tail(input.size() - shift);
	output.tail(shift) = input.head(shift);
}
void circShift_R(const Eigen::VectorXi& input, unsigned int shift, Eigen::VectorXi& output)
{
	if (shift == 0) {
		output = input;
		return;
	}
	output.resize(input.size());
	output.tail(input.size() - shift) = input.head(input.size() - shift);
	output.head(shift) = input.tail(shift);
}

/**
 * Circularly shifts elements in the vector. A positive shift cyclically moves elements to higher
 * indices, negative shifts in the opposite direction.
 */
void circShift(const Eigen::VectorXi& input, int shift, Eigen::VectorXi& output)
{
	if (shift < 0) circShift_L(input, -shift, output);
	else circShift_R(input, shift, output);
}

Eigen::VectorXi circShift(const Eigen::VectorXi& input, int shift)
{
	Eigen::VectorXi out;
	circShift(input, shift, out);
	return out;
}

int lcm(int left, int right, int maxIt = 100)
{
	int currL = left;
	int currR = right;
	for(int i = 0; i < maxIt; i++)
	{
		if (currL == currR) return currL;
		if (currL > currR) currR += right;
		else currL += currL;
	}
	return -1;
}

void unwrapField(const std::vector<int>& edges, 
	const std::vector<int>& edgeSides, 
	int numWraps, 
	int N, 
	const Eigen::MatrixXi& Matching,
	Eigen::MatrixXi& faceLevels, 
	Eigen::MatrixXi& edgeLevels)
{
	const int originalValence = edges.size() / 2;
	// ''New'' valence: valence of the locally unfolded geometry (branched covering space?)
    const int singularValence = numWraps * originalValence;

    faceLevels.setConstant(N, singularValence, 0);
	faceLevels.col(0) = Eigen::VectorXi::LinSpaced(N, 0, N - 1);
	int meshF = 1;
    for(int i = 1; i < singularValence; i++){
		const int e = 2 * meshF;
		// Update new levels for target face by applying matching
		faceLevels.col(i) = circShift(faceLevels.col(i - 1), -Matching(edges[e], 1 - edgeSides[e]));
		// Next face. Same face will be seen multiple times due to unwrapping.
		meshF = (meshF + 1) % originalValence;
    }
	// SANITY CHECK
	{
		for(int i = 0; i < singularValence; i++)
		{
			int seen[] = { 0,0,0,0 };
			int start = faceLevels(0, i);
			bool fail = false;
			for(int j = 1; j < N ;j++)
			{
				if(faceLevels(j,i) != ((start + 1)%N))
				{
					fail = true; break;
				}
				start = faceLevels(j, i);
			}
			if(fail) std::cout << "Invalid circshift occureed" << Eigen::VectorXi(faceLevels.col(i)) << std::endl;
		}
		/*Eigen::VectorXi last = circShift(faceLevels.col(singularValence - 1), -Matching(edges[0], 1 - edges[0]));
		if (last(0) != faceLevels(0, 0)) std::cout << "!!Face did not fully wrap to same after " << numWraps << " wraps: " <<
			last(0) << ", vs level " << faceLevels(0,0)
		<< std::endl;*/
	}

	// Setup edge levels
	edgeLevels.setConstant(N, edges.size() * numWraps, 0);

	for(int w = 0; w < numWraps; w++)
	{
		const int offset = w * edges.size();
		const int faceOffset = w * originalValence;
		// Fix levels due to orientation of edges. 
		for(int i = 0; i < originalValence; i++)
		{
			const int e0 = edges[2 * i], e1 = edges[2 * i + 1];
			const int s0 = edgeSides[2 * i], s1 = edgeSides[2 * i + 1];
			// Only apply extra compensation if the (half)edge is not clockwise wrt the face normal. So if s0/s1 = 0,
			// no compensation is applied.
			edgeLevels.col(offset + 2 * i) = circShift(faceLevels.col(faceOffset + i), -s0 * Matching(e0, s0) );
			edgeLevels.col(offset + 2 * i + 1) = circShift(faceLevels.col(faceOffset + i), -s1 * Matching(e1, s1) );
		}
	}
}

template<typename CoeffProvider>
struct DirectionalSubdivider
{
	CoeffProvider cp;
	int N;
	// Matching over edges. Both from left to right face (column 0) and right to left face (column 1).
	Eigen::MatrixXi* Matching;
	// Map from coarse level edge to fine level edges. Each row corresponds to an edge e and the columns 
	// correspond to : 0,1 -> even elements, 0 starts at start vertex of e and ends at a new odd vertex, 
	// 1 starts at the new odd vertex and ends at the end vertex of 1.
	// Column 2 corresponds to the odd edge created in the left face, which is kept in the same orientation (locally) as the original
	// edge. 
	Eigen::MatrixXi* E0ToEk = nullptr;

	const Eigen::MatrixXi* Singularities;

	EdgeData* ED = nullptr;

	LeveledSparseConstructor builder;

	bool isSigned;

	Eigen::SparseMatrix<double> getMatrix() const
	{
		return builder.matrix;
	}

	void setup(EdgeData& ED)
	{
		builder = LeveledSparseConstructor(ED.edgeCount() * N, ED.edgeCount() * N);
		builder.makeId();
		this->ED = &ED;
	}

	void updateMatching(const Eigen::MatrixXi& matching, Eigen::MatrixXi& output)
	{
		Eigen::MatrixXi newMat;
		newMat.setConstant(ED->edgeCount() * 4, 2, 0);
		for(int i = 0; i < E0ToEk->rows(); i++)
		{
			newMat(E0ToEk->coeff(i, 0), 0) = matching(i, 0);
			newMat(E0ToEk->coeff(i, 0), 1) = matching(i, 1);
			newMat(E0ToEk->coeff(i, 1), 0) = matching(i, 0);
			newMat(E0ToEk->coeff(i, 1), 1) = matching(i, 1);
		}
		output = newMat;
	}

	DirectionalSubdivider(int N, const Eigen::MatrixXi* singularities, Eigen::MatrixXi* Matching, bool isSigned):
		Singularities(singularities), 
		N(N), 
		isSigned(isSigned),
		Matching(Matching)
	{}

	template<typename...T>
	void prepareNext(SubdivisionBuilder<T...>& subdivider)
	{
		E0ToEk = &subdivider.E0ToEk;
		builder.cols = builder.rows;
		builder.rows += N * (subdivider.ED.edgeCount() + 3 * subdivider.ED.faceCount());

		std::cout << "New size " << builder.rows << "," << builder.cols << std::endl;
	}

	void handleBoundaryRing(
		const std::vector<int>& edges,
		const std::vector<int>& edgeSides
	)
	{
		//TODO
		//Ignored for now
	}

	void handleRegularRing(
		const std::vector<int>& edges,
		const std::vector<int>& edgeSides
	)
	{
		const int meshValence = edges.size() / 2;
		int centralV = ED->E(edges[0], edgeSides[0]);
		int wraps = 1;

		// We only check the first level singularities, since they don't move and 
		// initial level vertices retain their index.
		if(centralV < Singularities->rows() && Singularities->coeff(centralV,0) == 1)
		{
			wraps = lcm(Singularities->coeff(centralV,1), N) / Singularities->coeff(centralV, 1);
			std::cout << "Sing wraps:" << wraps << ", ind" << Singularities->coeff(centralV,1) << std::endl;
		}
		// Branched covering space valence?
		const int valence = meshValence * wraps;
		//std::cout << "Regular ring with wraps: " << wraps << std::endl;

		Eigen::MatrixXi faceLevels, edgeLevels;
		unwrapField(edges, edgeSides, wraps, N, *Matching, faceLevels, edgeLevels);
		//std::cout << "Field unwrapped" << std::endl;

		//SANITY CHECK
		{
			int wrapNum = 0;
			int startLevel = faceLevels(0, 0);
			int invalidCount = 0;
			for(int w = 0; w < wraps-1; w++)
			{
				for(int i = 0; i < meshValence; i++)
				{
					if (faceLevels(0, w * meshValence + i) == faceLevels(0, (w+1) * meshValence + i))
					{
						invalidCount++;
					}
				}	
			}
			if (invalidCount > 0)
			{
				std::cout << "--- Checking facelevels failed " << invalidCount << " times at wraps " << wraps << std::endl;
			}
		}
		
		// Old and new edge count
		const int oldECount = builder.cols / N;
		const int newECount = builder.rows / N;

		// Assume regular stencil is used everywhere along the ring.
		std::vector<int> inds;
		std::vector<double> coeffs;

		std::vector<double> signs(edges.size());
		for(int i = 0; i < edges.size(); i++)
		{
			if (isSigned) signs[i] = 1.0 - 2.0 * edgeSides[i];
			else signs[i] = 1.0;
		}

		// Even elements
		for(int e = 0; e < edges.size(); e += 2)
		{
			inds.clear();
			coeffs.clear();
			cp.getEvenRegularStencil(valence, e, inds, coeffs);
			for(int n = 0; n < N; n++)
			{
				const int target = (*E0ToEk)(edges[e], edgeSides[e]) + edgeLevels(n, e) * newECount;
				const double sign = signs[e];
				for(int i = 0; i < coeffs.size(); i++)
				{
					const int iW = inds[i] % edges.size();
					builder.addCoeff(target, edges[iW] + edgeLevels(n, inds[i]) * oldECount, coeffs[i] * sign * signs[iW]);
				}
			}
			
		}
		// Odd elements
		for (int e = 1; e < edges.size(); e += 2)
		{
			inds.clear();
			coeffs.clear();
			cp.getOddRegularStencil(valence, e, inds, coeffs);
			for (int n = 0; n < N; n++)
			{
				const int target = (*E0ToEk)(edges[e], 2 + edgeSides[e]) + faceLevels(n, e/2) * newECount;
				const double sign = signs[e];
				for (int i = 0; i < coeffs.size(); i++)
				{
					const int iW = inds[i] % edges.size();
					builder.addCoeff(target, edges[iW] + edgeLevels(n, inds[i]) * oldECount, coeffs[i] * sign * signs[iW]);
				}
			}
		}
	}

	void finalize()
	{
		std::cout << "DirBuilder finalizing" << std::endl;
		builder.finalize();
		//Matching has not been updated yet.
		if (Matching->rows() == builder.cols / N)
		{
			std::cout << "Updating matching" << std::endl;
			updateMatching(*Matching, *Matching);
		}
		std::cout << "DirBuilder finalized" << std::endl;
	}
};



#endif