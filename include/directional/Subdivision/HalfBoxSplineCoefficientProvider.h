struct HalfboxSplineCoefficientProvider
{
	static int factor(int valence)
	{
		if (valence == 3) return 1.0 / 12.0;
		if (valence == 4) return 1.0 / 8.0;
		if (valence == 5) return 1.0 / 4.0 - 1. / (16. * std::sin(igl::PI * 2 * 0.2) * std::sin(igl::PI * 2 * 0.2));
		return 0.25;
	}
	/**
	 * Stencils for faces that are newly created at the center of old faces.
	 */
	void getEvenRegularStencil(int vertexValence, int location, Eigen::VectorXi& inds, Eigen::VectorXd& coeffs )
	{
		// We will handle this separately
		inds = { 0,1,2,3 };
		coeffs = Eigen::VectorXd::Constant(4, 1. / 16.);
	}
	/**
	 * Stencils for faces that move towards original vertex on subdividing.
	 */
	void getOddRegularStencil(int vertexValence, int location, Eigen::VectorXi& inds, Eigen::VectorXd& coeffs)
	{
		const double beta = factor(3);
		switch(vertexValence)
		{
		case 3:
			coeffs = Eigen::VectorXd::Constant(3, 1./32. + beta * 1./8.);
			coeffs(location) = 3. / 16. - beta / 4.;
			inds = { 0,1,2 };
			break;
		case 4:
			coeffs = Eigen::VectorXd::Constant(3, 1. / 32.);
			coeffs(location) = 3. / 16. - beta / 4.;
			coeffs(location < 2 ? location + 2 : location - 2) = beta / 4.;
			inds = { 0,1,2,3 };
			break;
		case 5:
			coeffs = Eigen::VectorXd::Constant(beta / 8., 5);
			coeffs(location) = 3. / 16. - beta / 4.;
			coeffs(location == 0 ? vertexValence - 1 : location - 1) = 1. / 32.;
			coeffs(location == vertexValence - 1 ? 0 : location + 1) = 1. / 32.;
			int_range(0, vertexValence, inds);
			break;
		default:
			coeffs = Eigen::VectorXd(5);
			coeffs << beta * 0.5, 1. / 8., 3. / 4. - beta, 1. / 8., beta * 0.5;
			int_range_wrapped(location - 2, location + 3, vertexValence, inds);
		}
	}
	void getEvenBoundaryStencil(int neighboursCount,  Eigen::VectorXi& inds, Eigen::VectorXd& coeffs)
	{
		switch(neighboursCount)
		{
		case 0:
			coeffs = Eigen::VectorXd(1);
			coeffs(0) = 1. / 4. ;
			break;
		case 1:
			coeffs = { 5. / 16., -1. / 16. };
			break;
		case 2:
			coeffs = { 6. / 32., 1. / 32.,1. / 32. };
			break;
		case 3:
			coeffs = {}; // Not boundary
			break;
		default:
			break;
		}
	}
	void getOddBoundaryStencil(int vertexValence, int location, Eigen::VectorXi& inds, Eigen::VectorXd& coeffs)
	{
		switch(vertexValence)
		{
		case 2:
			coeffs = Eigen::VectorXd(1);
			coeffs(0) = 1. / 4.;
			inds = Eigen::VectorXi::Constant(1, 0);
			break;
		case 3:
			coeffs = { 2. / 12., 1. / 12. };
			inds = { location, 1 - location };
			break;
		case 4:
			inds = { 0,1,2 };
			coeffs = Eigen::VectorXd::Constant(3, 1. / 24.);
			coeffs(location) = 4. / 24.;
			break;
		default:
			if(location == 0 || location == vertexValence-2)
			{
				coeffs = { 4. / 24., 1. / 24., 1. / 24. };
				if (location == 0) inds = { 0,1,2 };
				else inds = { location, location - 1, location - 2 };
			}
			else if(location == 1 || location == vertexValence - 3)
			{
				coeffs = { 5. / 96., 14. / 96., 2. / 96., 3. / 96. };
				if(location == 1)
				{
					inds = { 0,1,2,3 };
				}
				else inds = { location + 1, location, location - 1, location - 2 };
			}
			else
			{
				int_range_wrapped(location - 2, location + 3, vertexValence, inds);
				coeffs = Eigen::VectorXd(5);
				coeffs << 1. / 32., 1. / 32., 4. / 32., 1. / .32, 1. / 32. ;
			}
			break;
		}
	}
};