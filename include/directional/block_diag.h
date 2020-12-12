#ifndef DIRECTIONAL_BLOCKDIAG_H
#define DIRECTIONAL_BLOCKDIAG_H
#include <Eigen/Eigen>
#include <numeric>
namespace directional{
template <typename Scalar>
void block_diag(
    const std::vector<Eigen::SparseMatrix<Scalar>*>& matrices,
    Eigen::SparseMatrix<Scalar> & C)
{
  using namespace Eigen;
  using SMat = Eigen::SparseMatrix<Scalar>;

  // Compute total new rows
  const size_t newRows = std::accumulate(matrices.begin(), matrices.end(), 0, [](const size_t& s1, const SMat* s2)-> size_t
  {
	  return s1 + s2->rows();
  });
  const size_t newCols = std::accumulate(matrices.begin(), matrices.end(), 0, [](const size_t& s1, const SMat* s2)-> size_t
  {
	  return s1 + s2->cols();
  });
  std::vector<Triplet<Scalar>> trips;
  
   	int rOff = 0, cOff = 0;
	for (int i = 0; i < matrices.size(); i++)
	{
		for (int k = 0; k < matrices[i]->outerSize(); ++k)
		{
			for (typename SMat::InnerIterator it(*matrices[i], k); it; ++it)
			{
                trips.emplace_back(rOff + it.row(), cOff + k, it.value());
			}
		}
		rOff += matrices[i]->rows();
		cOff += matrices[i]->cols();
	}
  C = SMat(newRows, newCols);
  C.setFromTriplets(trips.begin(), trips.end());
  C.makeCompressed();
}
}

#endif