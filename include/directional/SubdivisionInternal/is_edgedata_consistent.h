#ifndef DIRECTIONAL_IS_EDGEDATA_CONSISTENT_H
#define DIRECTIONAL_IS_EDGEDATA_CONSISTENT_H
#include <iostream>
#include <Eigen/Eigen>
namespace directional{

    inline bool is_edgedata_consistent(
        const Eigen::MatrixXi& F,
        const Eigen::MatrixXi& E,
        const Eigen::MatrixXi& EF,
        const Eigen::MatrixXi& EI,
        const Eigen::MatrixXi& SFE,
        bool verbose = false)
	{
		bool fail = false;

		for (int e = 0; e < E.rows(); e++)
		{
			if (EF(e, 0) != -1)
			{
				const int f = EF(e, 0);
				const int c = EI(e, 0);
				if (SFE(f, c) != e) {
					if (verbose) {
						std::cout << "Edge and sFE not the same at e = " << e << ", f,c = (" << f << "," << c << ")" << std::endl;
						fail = true;
					}
					else  return false;
				}
				if (SFE(f, c + 3) != 0) {
					if (verbose) {
						std::cout << "Corrner of sFE not 0 at e = " << e << ", f,c = (" << f << "," << c << ")" << std::endl;
						fail = true;
					}
					else return false;
				}
				if (F(f, (c + 1) % 3) != E(e, 0) && F(f, (c + 2) % 3) != E(e, 1))
				{
					if (verbose) {
						std::cout << "Edge vertices do not concur at e = " << e << ", f,c = (" << f << "," << c << "), E:" << E.row(e) << ", sFE" << SFE.row(f) << std::endl;
						fail = true;
					}
					else return false;
				}
			}
			if (EF(e, 1) != -1)
			{
				const int f = EF(e, 1);
				const int c = EI(e, 1);
				if (SFE(f, c) != e) {
					if (verbose) {
						std::cout << "Edge and sFE not the same at e = " << e << ", f,c = (" << f << "," << c << ")" << std::endl;
						fail = true;
					}
					else return false;
				}
				if (SFE(f, c + 3) != 1) {
					if (verbose) {
						std::cout << "Corrner of sFE not 1 at e = " << e << ", f,c = (" << f << "," << c << ")" << std::endl;
						fail = true;
					}
					else return false;
				}
				if (F(f, (c + 1) % 3) != E(e, 1) && F(f, (c + 2) % 3) != E(e, 0))
				{
					if (verbose) {
						std::cout << "Edge vertices do not concur at e = " << e << ", f,c = (" << f << "," << c << "), E:" << E.row(e) << ", sFE" << SFE.row(f) << std::endl;
						fail = true;
					}
					else return false;
				}
			}

		}
		for(int i = 0; i < SFE.rows(); i++)
		{
			for(int j = 0; j < 3; j++)
			{
				const int e = SFE(i, j);
				const int s = SFE(i, j + 3);
				// Check all edges are set.
				if (e < 0 || e >= EF.rows()) {
					if (verbose) {
						std::cout << "Edge index out of bounds: e = " << e << ", at f  = " << i;
						fail = true;
					}
					else return false;
				}
				if (EF(e, s) != i) {
					if (verbose) {
						std::cout << "Edge flap with e,s = (" << e << ',' << s << ") incorrect face: expected " << i << ", got " << EF(e, s) << std::endl;
						fail = true;
					}
					else return false;
				}
				if (EI(e, s) != j) {
					if (verbose) {
						std::cout << "Edge flap with e,s = (" << e << ',' << s << ") incorrect corner: expected " << j << ", got " << EI(e, s) << std::endl;
						fail = true;
					}
					else return false;
				}
			}
		}
		return !fail;
	}
}
#endif