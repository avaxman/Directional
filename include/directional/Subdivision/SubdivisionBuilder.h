#ifndef DIRECTIONAL_SUBDIVISIONBUILDER_H
#define DIRECTIONAL_SUBDIVISIONBUILDER_H
#include <vector>
#include <Eigen/Eigen>
#include <Eigen/Sparse>
#include <directional/Subdivision/EdgeData.h>
#include <directional/Subdivision/SimpleTimer.h>
#include <directional/Subdivision/DCEL.h>
namespace Subdivision
{
	using CoefficientList = std::vector<double>;
	using IndexList = std::vector<int>;
}

template<typename T,template<typename> class Container>
struct CircularLookup
{
	const Container<T>* target;
	int len;
	CircularLookup(const Container<T>* target) : target(target), len(target->size()) {}

	void setTarget(const Container<T>* target)
	{
		this->target = target;
		len = target->size();
	}
	void range(int start, int endExcl, Container<T>& target)
	{
		target = Container<T>(endExcl - start);
		for(int i = 0; i < endExcl-start; i++)
		{
			target[i] = at(start + i);
		}
	}
	int at(int ind) const
	{
		if (ind < 0) return (*target)[len + ind];
		if (ind >= len) return (*target)[ind - len];
		return (*target)[ind];
	}
	int operator[](int ind) const
	{
		return at(ind);
	}
};

struct LeveledSparseConstructor
{
	std::vector<Eigen::Triplet<double>> m_coefficients;
	int rows;
	int cols;
	Eigen::SparseMatrix<double> matrix;

	LeveledSparseConstructor():rows(0),cols(0){}
	LeveledSparseConstructor(int rows, int cols) :rows(rows), cols(cols) {}
	LeveledSparseConstructor(int rowsAndCols) : rows(rowsAndCols),cols(rowsAndCols)
	{
		
	}
	void makeId()
	{
		std::vector<Eigen::Triplet<double>> trips;
		trips.reserve(rows);
		for (int i = 0; i < rows; i++) trips.emplace_back(i, i, 1.);
		matrix = Eigen::SparseMatrix<double>(rows, rows);
		matrix.setFromTriplets(trips.begin(), trips.end());
	}

	void finalize()
	{
		SimpleTimer s;
		s.start();
		Eigen::SparseMatrix<double,Eigen::RowMajor> newMat(rows, cols);
		newMat.setFromTriplets(m_coefficients.begin(), m_coefficients.end());
		s.stop();
		std::cout << "Construct from triplet subtime: " << s.elapsed() << " ms" << std::endl;

		// Premultiply original matrix to acquire new matrix.
		s.reset().start();
		matrix = newMat * matrix;
		s.stop();
		std::cout << "Mul subtime: " << s.elapsed() << " ms" << std::endl;
	}

	/**
	 * \brief Adds a coefficient for the subdivision matrix. Note that
	 * repeated values for the same (r,c) pair will be ADDED.
	 * \param r Target row
	 * \param c Target column
	 * \param v Value.
	 */
	void addCoeff(int r, int c, double v)
	{
		assert(r >= 0 && r < rows && c >= 0 && c < cols);
		m_coefficients.emplace_back(r, c, v);
	}
};


template<typename...SubBuilders>
struct SubdivisionBuilder
{
	using expander = int[];
	// The subdivision builders to use
	std::tuple< SubBuilders&...> m_builders;
	// The current level
	int level = 0;

	// The active Edge data.
	EdgeData ED;

	// Mapping from the previous level to the level to construct.
	Eigen::MatrixXi E0ToEk;

	// The vertex count of the mesh.
	int vertexCount = 0;

	bool wasSetup = false;

	SubdivisionBuilder(Eigen::MatrixXi F, SubBuilders&... builders): m_builders(std::forward_as_tuple(builders...))
	{
		ED.construct(F);
	}

	template<std::size_t...Is>
	void prepareNext(std::index_sequence<Is...>)
	{
		(void) expander{0, (std::get<Is>(m_builders).prepareNext(*this),0)... };
	}

	/**
	 * \brief Allow the builders to do some post iteration handling 
	 * \tparam Is Indices of the builder types to handle
	 * Expander hack:https://stackoverflow.com/questions/28547456/call-void-function-for-each-template-type-in-a-variadic-templated-function
	 */
	template<std::size_t...Is>
	void finalize(std::index_sequence<Is...>)
	{
		(void)expander{0, (std::get<Is>(m_builders).finalize(),0)... };
	}

	template<std::size_t...Is>
	void handleRegularRing_impl(const std::vector<int>& edges, const std::vector<int>& edgeSides, std::index_sequence<Is...>)
	{
		(void)expander{ 0, (std::get<Is>(m_builders).handleRegularRing(edges,edgeSides),0)... };
	}

	/**
	 * \brief Handles a regular one ring, specified as a list of edges with accompanying orientations.
	 * The edges are given as alternating spoke-ring edges, where spokes are the edges emanating from the 
	 * central vertex of the one ring and ring edges are the edges opposite the central vertex. The elements
	 * are given in counterclockwise order with respect to the ''vertex'' normal. The sides denote whether the
	 * global orientation of each edge is the same as the canonical orientation of the edges. The canonical orientation
	 * here is the orientation where the spokes are pointing away from the central vertex and the ring edges are pointing
	 * in counterclockwise direction with respect to the vertex.
	 * \param edges Edge indices
	 * \param edgeSides Orientations of the edges. 0 = aligned with canonical orientation, 1 = opposite direction wrt canonical orientation.
	 */
	void handleRegularRing(const std::vector<int>& edges, const std::vector<int>& edgeSides)
	{
		handleRegularRing_impl(edges, edgeSides, std::index_sequence_for<SubBuilders...>{});
	}

	template<std::size_t...Is>
	void setup_impl(std::index_sequence<Is...>)
	{
		(void)expander {
		0, (std::get<Is>(m_builders).setup(ED), 0)...
	};
	}
	void setup()
	{
		setup_impl(std::index_sequence_for<SubBuilders...>{});
	}

	template<std::size_t...Is>
	void handleBoundaryRing_impl(const std::vector<int>& edges, const std::vector<int>& edgeSides, std::index_sequence<Is...>)
	{
		(void)expander {
			0, (std::get<Is>(m_builders).handleBoundaryRing(edges,edgeSides), 0)...
		};
	}

	/**
	 * \brief Handles a boundary ring. The edges and orientations are given in the same manner as the regular ring edges. The only
	 * difference is that there is one more spoke edge than ring edges. The same orientation as for the regular ring is assumed.
	 * \param edges Edge indices
	 * \param edgeSides Edge orientations
	 */
	void handleBoundaryRing(const std::vector<int>& edges, const std::vector<int>& edgeSides)
	{
		handleBoundaryRing_impl(edges, edgeSides, std::index_sequence_for<SubBuilders...>{});
	}

	void constructUpToLevel(int targetLevel)
	{
		if(!wasSetup)
		{
			wasSetup = true;
			setup();
		}
		for(int i = level; i < targetLevel; i++)
		{
			constructNextLevel();
		}
	}

	/**
	 * \brief Constructs the next level subdivision operators by invoking the 1-ring handling functions
	 * on the provided subdivision builders for each one ring.
	 */
	void constructNextLevel()
	{
		if (!wasSetup)
		{
			wasSetup = true;
			setup();
		}
		level++;
		std::cout << "Constructing level " << level << std::endl;;

		SimpleTimer st;
		st.start();
		int vCount = ED.vertexCount(); // TODO
		ED.quadrisect(vCount, E0ToEk);
		st.stop();
		std::cout << "Quadrisect: " << st.elapsed() << " ms" << std::endl;;
		assert(ED.isConsistent());

		// Prepare for next step
		prepareNext(std::index_sequence_for<SubBuilders...>{});

		DCEL dcel(ED);

		SimpleTimer st2;
		st2.start();
		// Iterate over all 1 rings. Calls handleRegularRing() and handleBoundaryRing() on this object
		dcel.iterateRings(*this);
		st2.stop();
		std::cout << "Ring iteration : " << st2.elapsed() << " ms" << std::endl;;

		// Finalize construction
		st2.reset().start();
		finalize(std::index_sequence_for<SubBuilders...>{});
		st2.stop();
		std::cout << "Finalizing in : " << st2.elapsed() << " ms" << std::endl;
	}
};

template<typename...SubBuilders>
SubdivisionBuilder<SubBuilders...> create_sub_builder(const Eigen::MatrixXi& F, SubBuilders&...builders)
{
	return SubdivisionBuilder<SubBuilders...>(F, builders...);
}
#endif