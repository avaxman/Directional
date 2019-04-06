
struct OneFormCoefficientProvider
{
	static void assign(Eigen::VectorXd& v, const std::vector<double>& data, double factor = 1.0)
	{
		v = Eigen::VectorXd(data.size());
		for (int i = 0; i < data.size(); i++) v(i) = data[i] * factor;
	}

	static void assign(Eigen::VectorXi& v, const std::vector<int>& data)
	{
		v = Eigen::VectorXi(data.size());
		for (int i = 0; i < data.size(); i++) v(i) = data[i];
	}
	static void assignRange(int start, int endExcl, Eigen::VectorXi& target)
	{
		for (int i = 0; i < endExcl - start; i++)target(i) = start + i;
	}
	static void assignDecreasingRange(int start, int endExcl, Eigen::VectorXi& target)
	{
		for (int v = endExcl - 1,i=0; v>= start; v--,i++)target(i) = v;
	}


	void getEvenBoundaryStencil(int valence, int location, Eigen::VectorXi& inds, Eigen::VectorXd& coeffs)
	{
		const int maxEdgeIndex = 2 * valence - 2;
		if(location == 0 || location == maxEdgeIndex)
		{
			coeffs = { 3. / 8., -1. / 8. };
			inds = { location == 0 ? 0 : maxEdgeIndex, location == 0 ? maxEdgeIndex : 0 };
		}
		else if(location == 2 || location == maxEdgeIndex - 2)
		{
			if(valence == 3)
			{
				assign(coeffs, {
					5./3.,
					5./3.,
					26./3.,
					-5./3.,
					5./3.
					},1./32.);
				if (location == 2) assign(inds, { 0,1,2,3,4 });
				else assign(inds, { maxEdgeIndex, maxEdgeIndex-1, maxEdgeIndex-2, maxEdgeIndex-3, maxEdgeIndex-4 });
			}
			else if(valence == 4)
			{
				assign(coeffs, {
					5./3.,
					5./3.,
					10,
					-1./3.,
					3.,
					-4./3.,
					-8./3.
					}, 1./32.);
				if (location == 2) int_range(0, 7, inds);
				else assign(inds, { maxEdgeIndex, maxEdgeIndex - 1, maxEdgeIndex - 2, maxEdgeIndex - 3, maxEdgeIndex - 4, maxEdgeIndex - 5, maxEdgeIndex - 6 });
			}
			else
			{
				assign(coeffs, {
					5./3.,
					5./3.,
					10.,
					-1./3.,
					3.,
					-4./3.,
					4./3.,
					-4.
					},1./32.);
				inds = Eigen::VectorXi(8);
				if (location == 2) {
					assignRange(location - 2, location + 5, inds);
					inds(7) = maxEdgeIndex;
				}
				else 
				{
					assignDecreasingRange(maxEdgeIndex - 6, maxEdgeIndex+1, inds);
					inds(7) = 0;
				}
			}
		}
		else
		{
			assign(coeffs, {
				1./32,
				1./32,
				4./32,
				1./32,
				10./32.,
				-1./32.,
				4./32.,
				-1./32.,
				1./32.,
				-4/32., // Boundary coefficient
				-4./32.// Boundary coefficient
			});
			int_range(location - 4, location + 5, inds);
			inds.conservativeResize(11);
			inds(9) = 0;
			inds(10) = maxEdgeIndex;
		}
	}
	void getOddBoundaryStencil(int valence, int location, Eigen::VectorXi& inds, Eigen::VectorXd& coeffs)
	{
		const int maxEdgeIndex = valence * 2 - 2;
		if (location == 1 || location == maxEdgeIndex - 1)
		{
			if(valence == 2)
			{
				assign(coeffs,
					{ -1. / 4., 1. / 4, 1. / 4. });
				if (location == 1) inds = { 0,1,2 };
				else inds = { maxEdgeIndex, maxEdgeIndex - 1, maxEdgeIndex - 2 };
			}
			else
			{
				assign(coeffs,
					{ -5. / 32., 7. / 32, 6./32, 1./32, 3./32.});
				if (location == 1) assign(inds, { 0, 1, 2, 3, 4 });
				else assign(inds, { maxEdgeIndex, maxEdgeIndex - 1, maxEdgeIndex - 2, maxEdgeIndex - 3, maxEdgeIndex - 4 });
			}
		}
		else
		{
			assign(coeffs, {
				-3./32.,
				1./32,
				-3./32.,
				6./32.,
				3./32.,
				1./32.,
				3./32.
				});
			assign(inds, {location-3, location -2, location -1, location, location + 1, location + 2, location + 3});
		}
	}
	void getEvenRegularStencil(int valence, int location, Eigen::VectorXi& inds, Eigen::VectorXd& coeffs)
	{
		// Loop subdivision factor
		const double alfa = LoopCoefficientProvider::factor(valence);
		// Halfbox spline subdivision factor
		const double beta = HalfboxSplineProvider::factor(valence);
		// Indices
		Eigen::VectorXi localInds(valence*2);
		for (int i = 0; i < valence * 2; i++) localInds(i) = i;
		// 
		CircularAccess ca(&localInds);
		switch(valence)
		{
		case 0:
		case 1:
		case 2:
			throw std::runtime_error("Invalid regular stencil valence");
		case 3:
			coeffs = Eigen::VectorXd(6);
			coeffs << 1. / 8. - alfa + beta / 8., 
				beta / 8., 
				3. / 8. - alfa - .25 * beta, 
				-beta / 8.,
				1. / 8. - alfa + beta / 8.;
			ca.range(location - 2, location + 3, inds);
			break;
		case 4:
			coeffs = Eigen::VectorXd(8);
			coeffs << 1. / 8. - alfa, 
				beta / 8.,
				3. / 8. - alfa - .25 * beta, // Target edge
				-beta / 8.,
				1. / 8. - alfa, 
				-beta / 8.,
				beta / 4. - alfa;
			ca.range(location - 2, location + 4, inds);
			break;
		default:
			// Higher valence
			break;
		}

	}
	void getOddRegularStencil(int valence, int location, Eigen::VectorXi& inds, Eigen::VectorXd& coeffs)
	{
		Eigen::VectorXi localInds;
		int_range(0, valence * 2, localInds);
		CircularAccess ca(&localInds);

		if(valence == 3)
		{
			assign(coeffs, 
				{
				1./32., 
				-3./32.,
				6./32.,
				3./32.,
				1./32
				});
			ca.range(location - 2, location + 3, inds);
		}
		else
		{
			assign(coeffs,
				{
				-3./32.,
				1. / 32.,
				-3. / 32.,
				6. / 32.,
				3. / 32.,
				1. / 32,
				3. /32.
				});
			ca.range(location - 3, location + 4, inds);
		}
	}
};
