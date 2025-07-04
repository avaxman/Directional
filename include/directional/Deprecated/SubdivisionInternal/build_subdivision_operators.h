#ifndef DIRECTIONAL_BUILD_SUBDIVISION_OPERATORS_H
#define DIRECTIONAL_BUILD_SUBDIVISION_OPERATORS_H
#include <Eigen/Eigen>
#include "quadrisect.h"
#include "iterate_rings.h"

namespace directional{

	// Type of a triplet provider
	template<typename CoeffProvider>
	using triplet_provider_t = void(*)(const int&, //Vertex count
			const Eigen::MatrixXi&, // F0
			const Eigen::MatrixXi&, // SFE0
			const Eigen::MatrixXi&, // E0 
			const Eigen::MatrixXi&, // EI0
			const Eigen::MatrixXi&, // EF0
			const Eigen::MatrixXi&, // E0ToEk
			const std::vector<int>&, // edges
			const std::vector<int>&, // edge orientations
			const CoeffProvider&, //Coefficient provider
			std::vector<Eigen::Triplet<double>>&, // Output triplets
			int&); //Maximum row number assigned

	//(bool isBoundary, bool isEven, int valence, int location, std::vector<int>& inds, std::vector<double>& coeffs)
	using coefficient_provider_t = void(*)(bool, bool, int, int, std::vector<int>&, std::vector<double>&);

	template<typename CoeffProvider>
	auto triplet_provider_wrapper(
			const CoeffProvider& coeffProvider,
			triplet_provider_t<CoeffProvider> provider
	)
	{
		return [coeffProvider, provider](const int& vertexCount,
			const Eigen::MatrixXi& F0,
			const Eigen::MatrixXi& SFE0,
			const Eigen::MatrixXi& E0,
			const Eigen::MatrixXi& EI0,
			const Eigen::MatrixXi& EF0,
			const Eigen::MatrixXi& E0ToEk,
			const std::vector<int>& edges,
			const std::vector<int>& edgeOrients,
			std::vector<Eigen::Triplet<double>>& output,
			int& maxRowSize){
				provider(vertexCount, F0, SFE0, E0, EI0, EF0, E0ToEk, edges, edgeOrients, coeffProvider, output, maxRowSize);
		};
	}

	template<typename...TripletProviders,std::size_t...Is>
	void handleRing(
		const int& vertexCount,
		const Eigen::MatrixXi& F0,
		const Eigen::MatrixXi& SFE0,
		const Eigen::MatrixXi& E0,
		const Eigen::MatrixXi& EI0,
		const Eigen::MatrixXi& EF0,
		const Eigen::MatrixXi& E0ToEk,
		const std::vector<int>& edges, 
		const std::vector<int>& edgeSides, 
		std::vector<std::vector<Eigen::Triplet<double>>>& output,
		std::vector<int>& rowSizes,
		std::tuple<TripletProviders...> tripletProviders, 
		std::index_sequence<Is...>
	)
	{
		// Hack to apply the triplet providers as functor on the given edge data.
		using exp = int[];
		(void) exp{
			(std::get<Is>(tripletProviders)(vertexCount, F0, SFE0, E0, EI0, EF0, E0ToEk, edges,edgeSides, output[Is], rowSizes[Is]), 0)...
		};
	}

	template<typename...TripletProviders>
    void build_subdivision_operators(
        const Eigen::MatrixXd& V0,
		const Eigen::MatrixXi& F0,
		const Eigen::MatrixXi& E0,
		const Eigen::MatrixXi& EF0,
		const Eigen::MatrixXi& EI0,
		const Eigen::MatrixXi& SFE0,
		const std::vector<int>& initialSizes,
		int level,
		Eigen::MatrixXi& FK,
		Eigen::MatrixXi& EK,
		Eigen::MatrixXi& EFK,
		Eigen::MatrixXi& EIK,
		Eigen::MatrixXi& SFEK,
		std::vector<Eigen::SparseMatrix<double>>& output,
		TripletProviders...tripletProviders
	)
	{
		constexpr int N = sizeof...(TripletProviders);
        

		// Matrices representing geometry connectivity.
		Eigen::MatrixXi Fs[2] = { F0,{} };
		Eigen::MatrixXi Es[2] = { E0, {} };
		Eigen::MatrixXi SFEs[2] = { SFE0,{} };
		Eigen::MatrixXi EFs[2] = { EF0,{} };
		Eigen::MatrixXi EIs[2] = { EI0,{} };
		Eigen::MatrixXi E0ToEK; // Map edges from level to next level
		
		// The current vertex count
		int currentVCount = V0.rows();

		// The index of the matrices to fill next
		int toFill = 1;

		// Tuple of subdivision constructors
		std::tuple<TripletProviders...> constructors = std::make_tuple(tripletProviders...);

		// The triplets
		std::vector<std::vector<Eigen::Triplet<double>>> triplets(N, std::vector<Eigen::Triplet<double>>{});

		// The row sizes for the jump level subdivision operator
		std::vector<int> rowSizes(N, 0);

		// Function to handle a new ring
		auto ringHandler = [&Fs,&SFEs,&Es,&EFs, &EIs, &E0ToEK, &triplets,&constructors, &currentVCount, &toFill, &rowSizes](const std::vector<int>& edges, const std::vector<int>& edgeSides)
		{
			const int filled = 1- toFill;
			/*const int& vertexCount,
				const Eigen::MatrixXi& F0,
				const Eigen::MatrixXi& SFE0,
				const Eigen::MatrixXi& E0,
				const Eigen::MatrixXi& EI0,
				const Eigen::MatrixXi& EF0,
				const Eigen::MatrixXi& E0ToEk,
				const std::vector<int>& edges,
				const std::vector<int>& edgeSides,
				std::vector<std::vector<Eigen::Triplet<double>>>& output,
				std::vector<int>& rowSizes,
				std::tuple<TripletProviders...> tripletProviders,
				std::index_sequence<Is...>*/
			handleRing(currentVCount, 
				Fs[filled],
				SFEs[filled],
				Es[filled],
				EIs[filled],
				EFs[filled],
				E0ToEK, 
				edges, 
				edgeSides, 
				triplets, 
				rowSizes, 
				constructors, std::index_sequence_for<TripletProviders...>{});
		};


		// Initialize the output to identity matrices initially. We are 
		// going to progressively build the subdivision operator by 
		// multiplying the operator for the different levels.
		for(int i = 0; i < N; i++)
		{
			output.emplace_back(initialSizes[i],initialSizes[i]);
			output.back().setIdentity();
		}

		// Construct subdivision per level
		for(int i = 0; i < level; i++)
		{
			// Quadrisect connectivity data first
			const int filled = 1 - toFill;
			quadrisect(Fs[filled], currentVCount, Es[filled], SFEs[filled], EFs[filled], EIs[filled], E0ToEK,
				Fs[toFill], Es[toFill], SFEs[toFill], EFs[toFill], EIs[toFill]);

			// Iterate over all rings in the mesh, apply the subdivision constructors to acquire
			// the triplets for every matrix.
			iterate_rings(currentVCount, Es[filled], EFs[filled], EIs[filled], SFEs[filled], ringHandler);

			// Construct subdivision operators
			for(int j = 0; j < N; j++)
			{
				// Construct the operator to move one subdivision level up
				Eigen::SparseMatrix<double> levelJumpMat(rowSizes[j], output[j].rows() );
				levelJumpMat.setFromTriplets(triplets[j].begin(), triplets[j].end());

				// Apply the operator to the previous subdivision operator
				output[j] = levelJumpMat * output[j];
				triplets[j].clear();
			}
			
			// Update target
			currentVCount += Es[filled].rows();
			toFill = filled;
		}

		// Setup output connectivity matrices
		const int outputInd = 1 - toFill;
		EFK = EFs[outputInd];
		SFEK = SFEs[outputInd];
		FK = Fs[outputInd];
		EK = Es[outputInd];
		EIK = EIs[outputInd];
    }
}

#endif