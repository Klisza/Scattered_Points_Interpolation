#include <sparse_interp/energy.h>
#include <sparse_interp/surface.h> 
#include <igl/Timer.h>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>
#include <sparse_interp/Types.hpp>
#include <cmath>
#include <iostream>
namespace SIBSplines{
	igl::Timer timer;
	double time0 = 0, time1 = 0, time2 = 0, time3 = 0;
	
	// These functions handle the polynomials.
	template <typename Tp, typename knotT, typename valueT>
	std::vector<knotT> ply_operations<Tp, knotT, valueT>::polynomial_add(const std::vector<knotT>& poly1, const std::vector<knotT>& poly2) {
		int size = std::max(poly1.size(), poly2.size());
		std::vector<knotT> result(size);
		for (int i = 0; i < size; i++) {
			bool flag1 = i < poly1.size();
			bool flag2 = i < poly2.size();
			if (flag1 && flag2) {
				result[i] = poly1[i] + poly2[i];
			}
			else if (flag1) {
				result[i] = poly1[i];
			}
			else {
				result[i] = poly2[i];
			}
		}
		return result;//ply_operations<Tp, knotT, valueT>::polynomial_simplify(result);
	}
	// Cauchy product/Polynomial product
	template <typename Tp, typename knotT, typename valueT>
	std::vector<knotT> ply_operations<Tp, knotT, valueT>::polynomial_times(const std::vector<knotT>& poly1, const std::vector<knotT>& poly2) {
		int size = poly1.size() + poly2.size() - 1;
		std::vector<knotT> result(size);
		for (int i = 0; i < size; i++) {// initialize the result
			result[i] = 0;
		}

		for (int i = 0; i < poly1.size(); i++) {
			for (int j = 0; j < poly2.size(); j++) {
				result[i + j] += poly1[i] * poly2[j];
			}
		}
		return result; //ply_operations<Tp, knotT, valueT>::polynomial_simplify(result);
	}
	// Scalar multiplication
	template <typename Tp, typename knotT, typename valueT>
	std::vector<Tp> ply_operations<Tp, knotT, valueT>::polynomial_times(const std::vector<Tp>& poly1, const Tp& nbr) {
		std::vector<Tp> result;
		result = poly1;
		for (int i = 0; i < result.size(); i++) {
			result[i] *= nbr;
		}

		return result;
	}
	// This might need to be updated to work on TinyAD since it should use polynomial_times because of the AD type.
	template <typename Tp, typename knotT, typename valueT>
	Tp ply_operations<Tp, knotT, valueT>::power(const valueT &value, const int order)
	{
		if (order == 0)
		{
			return 1;
		}
		Tp result = value;
		for (int i = 1; i < order; i++)
		{
			result = result * value;
		}
		return result;
	}
	template <typename Tp, typename knotT, typename valueT>
	Tp ply_operations<Tp, knotT, valueT>::polynomial_value(const std::vector<knotT>& poly, const valueT para) {
		double result = 0;
		for (int i = 0; i < poly.size(); i++) {
			result += poly[i] * power(para, i);
		}
		return result;
	}
	template <typename Tp, typename knotT, typename valueT>
	// I dont know if divison works
	std::vector<Tp> ply_operations<Tp, knotT, valueT>::polynomial_integration(const std::vector<Tp>& poly) {
		std::vector<Tp> result(poly.size() + 1);
		result[0] = 0;
		for (int i = 1; i < result.size(); i++) {
			result[i] = poly[i - 1] / i;
		}
		return result;
	}
	template <typename Tp, typename knotT, typename valueT>
	Tp ply_operations<Tp, knotT, valueT>::polynomial_integration(const std::vector<Tp>& poly, const Tp lower, const Tp upper) {
		Tp up = ply_operations::polynomial_value(ply_operations::polynomial_integration(poly), upper);
		Tp lw = ply_operations::polynomial_value(ply_operations::polynomial_integration(poly), lower);
		return up - lw;
	}
	// order 1 differential. I guess obsolete with TinyAD? ////
	const std::vector<double> polynomial_differential(const std::vector<double>& func) {
		std::vector<double> result;
		if (func.size() == 1) {
			result.resize(1);
			result[0] = 0;
			return result;
		}
		result.resize(func.size() - 1);
		for (int i = 0; i < result.size(); i++) {
			result[i] = func[i + 1] * (i + 1);
		}
		return result;

	}

	std::vector<double> Ni0_func(const int i, const double u, const std::vector<double> &U) {
		std::vector<double> result(1);
		if (u >= U[i] && u < U[i + 1]) {
			result[0] = 1;
			return result;
		}
		result[0] = 0;
		return result;
	}
	
	// TODO!! ///////////// Rewrite it for templated functions.
	template <typename Tp, typename knotT, typename valueT>
	std::vector<Tp> Nip_func(const int i, const int p, const double u, const std::vector<double> &U) {
		ply_operations<Tp, knotT, valueT> PO;
		std::vector<Tp> result;
		if (p == 0) {
			return Ni0_func(i, u, U);
		}
		if (u == U.back()) {
			result.resize(1);
			if (i == U.size() - 2 - p) {
				result[0] = 1;
				return result;
			}
			else {
				result[0] = 0;
				return result;
			}
		}
		std::vector<double> v;
		v = { {-U[i],1} };// u - U[i]
		std::vector<double> result1 = ply_operations<Tp, knotT, valueT>::polynomial_times(handle_division_func(v, U[i + p] - U[i]), Nip_func(i, p - 1, u, U));

		v = { {U[i + p + 1],-1} };// U[i+p+1] - u 
		std::vector<double> result2 = ply_operations<Tp, knotT, valueT>::polynomial_times(handle_division_func(v, U[i + p + 1] - U[i + 1]), Nip_func(i + 1, p - 1, u, U));
		return ply_operations<Tp, knotT, valueT>::polynomial_add(result1, result2);
	}
	
	// Same thing here TinyAD should autodifferenciate this whole thing. ////
	const std::vector<double> polynomial_differential(const std::vector<double>& func, const int order) {
		std::vector<double> result = func;
		if (order == 0) return result;
		if (func.size() == 1 && order > 0) {
			result.resize(1);
			result[0] = 0;
			return result;
		}
		std::vector<double> tmp;
		for (int i = order; i > 0; i--) {
			tmp = polynomial_differential(result);
			result = tmp;
		}
		return result;
	}
	template <typename Tp, typename knotT, typename valueT>
	void PolynomialBasis<Tp, knotT, valueT>::init(Bsurface<Tp, knotT, valueT>& surface) {
		Uknot = surface.U;
		Vknot = surface.V;
		degree1 = surface.degree1;
		degree2 = surface.degree2;
		nu = surface.nu();
		nv = surface.nv();
		Ubasis = calculate(0);
		Vbasis = calculate(1);
		inited = true;
		return;
	}
	template <typename Tp, typename knotT, typename valueT>
	void PolynomialBasis<Tp, knotT, valueT>::clear() {
		Uknot.clear();
		Vknot.clear();
		Ubasis.clear();
		Vbasis.clear();
		inited = false;
	}
	template <typename Tp, typename knotT, typename valueT>
	PolynomialBasis<Tp, knotT, valueT>::PolynomialBasis(Bsurface<Tp, knotT, valueT>& surface) {
		init(surface);
	}
	template <typename Tp, typename knotT, typename valueT>
	PolynomialBasis<Tp, knotT, valueT>::PolynomialBasis() {
	}
	// Poly return the coefficents 
	template <typename Tp, typename knotT, typename valueT>
	std::vector<double> PolynomialBasis<Tp, knotT, valueT>::poly(const int id, const double value, const bool UVknot) {
		if (!inited) {
			std::cout << "WRONG USAGE OF CLASS PolynomialBasis, YOU SHOULD INITIALIZE IT BY CALLING init()" << std::endl;
		}

		std::vector<double> kv;
		int degree;

		if (UVknot) {
			kv = Vknot;
			degree = degree2;
		}
		else {
			kv = Uknot;
			degree = degree1;

		}

		int which = -1;
		for (int i = 0; i < kv.size() - 1; i++) {
			if (value >= kv[i] && value < kv[i + 1]) {
				which = i;
				break;
			}
		}
		if (which == -1) {
			std::cout << "ERROR: DON'T USE POLYNOMIAL WHEN VALUE = " << value << std::endl;
			exit(0);
		}
		// the value is in [U[i], U[i+1]), the Nip are from N(i-p) to N(i)
		std::vector<double> result;
		int index = id - (which - degree);
		if (UVknot) {// check v
			assert(which < Vbasis.size());
			assert(index < Vbasis[which].size());
			result = Vbasis[which][index];
		}
		else {
			assert(which < Ubasis.size());
			assert(index < Ubasis[which].size());
			result = Ubasis[which][index];
		}
		return result;
	}
	// 
	template<typename Tp, typename knotT, typename valueT>
	std::vector<std::vector<std::vector<Tp>>> PolynomialBasis<Tp, knotT, valueT>::calculate_single(const int degree, const std::vector<knotT> &knotVector)
	{
		std::vector<std::vector<std::vector<Tp>>> pl;
		int n = knotVector.size() - 2 - degree;
		pl.resize(n + 1);
		for (int i = degree; i < n + 1; i++) {// in interval [U[i], U[i+1])
			pl[i].resize(degree + 1);
			for (int j = 0; j < degree + 1; j++) {
				pl[i][j] = Nip_func(i - degree + j, degree, knotVector[i], knotVector);
			}
		}
		return pl;
	}
	// Should work
	template<typename Tp, typename knotT, typename valueT>
	std::vector<Tp> basisValues(const int whichItv, const int degree, const std::vector<std::vector<std::vector<double>>>&basis, const double param)
	{
		std::vector<Tp> result(degree + 1);
		for(int i=0;i<degree+1;i++)
		{
			if(whichItv>=basis.size())
			{
				std::cout<<"which itv "<<whichItv<<" basis size "<<basis.size()<<'\n';
			}
			if(basis[whichItv].size()==0)
			{
				std::cout<<"The basis whihcitv is empty, which itv "<<whichItv<<", the size of basis functions "<<basis.size()<<"\n";
			}
			result[i] = ply_operations<Tp, knotT, valueT>::polynomial_value(basis[whichItv][i], param);
		}
		return result;
	}
	template<typename Tp, typename knotT, typename valueT>
	std::vector<std::vector<std::vector<double>>> PolynomialBasis<Tp, knotT, valueT>::calculate(const bool uorv) {
		std::vector<double> kv;
		int degree;
		int n;
		if (uorv) {
			kv = Vknot;
			degree = degree2;
			n = nv;
		}
		else {
			kv = Uknot;
			degree = degree1;
			n = nu;
		}
		std::vector<std::vector<std::vector<double>>> pl;
		pl.resize(n + 1);// n+1 intervals;

		for (int i = degree; i < n + 1; i++) {// in interval [U[i], U[i+1])
			pl[i].resize(degree + 1);
			for (int j = 0; j < degree + 1; j++) {
				pl[i][j] = Nip_func(i - degree + j, degree, kv[i], kv);
			}
		}

		return pl;

	}
	template<typename Tp, typename knotT, typename valueT>
	std::vector<std::vector<std::vector<double>>> PartialBasis<Tp, knotT, valueT>::do_partial(const
		std::vector<std::vector<std::vector<double>>>&basis) {
		std::vector<std::vector<std::vector<double>>> result(basis.size());
		for (int i = 0; i < basis.size(); i++) {
			result[i].resize(basis[i].size());
			for (int j = 0; j < basis[i].size(); j++) {
				result[i][j] = polynomial_differential(basis[i][j]);
			}
		}
		return result;
	}

	template<typename Tp, typename knotT, typename valueT>
	PartialBasis<Tp, knotT, valueT>::PartialBasis(Bsurface<Tp, knotT, valueT>& surface) {
		init(surface);
	}
	template<typename Tp, typename knotT, typename valueT>
	PartialBasis<Tp, knotT, valueT>::PartialBasis() {
	}
	template<typename Tp, typename knotT, typename valueT>
	void PartialBasis<Tp, knotT, valueT>::init(Bsurface<Tp, knotT, valueT>& surface) {
		PolynomialBasis pb(surface);
		Ubasis = pb.Ubasis;
		Vbasis = pb.Vbasis;
		Ubasis_1 = do_partial(Ubasis); Vbasis_1 = do_partial(Vbasis);
		Ubasis_2 = do_partial(Ubasis_1); Vbasis_2 = do_partial(Vbasis_1);
		Uknot = surface.U;
		Vknot = surface.V;
		degree1 = surface.degree1;
		degree2 = surface.degree2;
	}
	template<typename Tp, typename knotT, typename valueT>
	void PartialBasis<Tp, knotT, valueT>::init(PolynomialBasis<Tp, knotT, valueT> &pb){
		Ubasis = pb.Ubasis;
		Vbasis = pb.Vbasis;
		Ubasis_1 = do_partial(Ubasis); Vbasis_1 = do_partial(Vbasis);
		Ubasis_2 = do_partial(Ubasis_1); Vbasis_2 = do_partial(Vbasis_1);
		Uknot = pb.Uknot;
		Vknot = pb.Vknot;
		degree1 = pb.degree1;
		degree2 = pb.degree2;
	}
	template<typename Tp, typename knotT, typename valueT>
	void PartialBasis<Tp, knotT, valueT>::clear() {
		Uknot.clear();
		Vknot.clear();
		Ubasis.clear();
		Vbasis.clear();
		Ubasis_1.clear();
		Vbasis_1.clear();
		Ubasis_2.clear();
		Vbasis_2.clear();
	}
	// What exact type should this function be
	template<typename Tp, typename knotT, typename valueT>
	std::vector<double> PartialBasis<Tp, knotT, valueT>::poly(const int id, const double value, const bool UVknot, int partial) {
		std::vector<double> kv;
		int degree;

		if (UVknot) {
			kv = Vknot;
			degree = degree2;
		}
		else {
			kv = Uknot;
			degree = degree1;

		}

		int which = -1;
		for (int i = 0; i < kv.size() - 1; i++) {
			if (value >= kv[i] && value < kv[i + 1]) {
				which = i;
				break;
			}
		}
		if (which == -1) {
			std::cout << "ERROR: DON'T USE POLYNOMIAL WHEN VALUE = " << value << std::endl;
			exit(0);
		}
		// the value is in [U[i], U[i+1]), the Nip are from N(i-p) to N(i)
		std::vector<double> result;
		int index = id - (which - degree);
		if (UVknot) {// check v
			if (partial == 0) {
				result = Vbasis[which][index];
			}
			if (partial == 1) {
				result = Vbasis_1[which][index];
			}
			if (partial == 2) {
				result = Vbasis_2[which][index];
			}

		}
		else {
			if (partial == 0) {
				result = Ubasis[which][index];
			}
			if (partial == 1) {
				result = Ubasis_1[which][index];
			}
			if (partial == 2) {
				result = Ubasis_2[which][index];
			}
		}
		return result;
	}

	
	// construct an integration of multiplication of two B-spline basis (intergration of partial(Ni1)*partial(Ni2))
	// the integration domain is [u1, u2]
	template<typename Tp, typename knotT, typename valueT>
	double construct_an_integration(const int degree, const std::vector<double>& U,
		const int partial1, const int partial2, const int i1, const int i2, const double u1, const double u2,
		PolynomialBasis<Tp, knotT, valueT> &basis, const bool uv) {
		timer.start();
		std::vector<double> func1 = basis.poly(i1, u1, uv);
		//Nip_func(i1, degree, u1, U);
		std::vector<double> func2 = basis.poly(i2, u1, uv);
		//Nip_func(i2, degree, u1, U);
		/*if (basis.poly(i1, u1, uv) != func1) {
			std::cout << "NOT EQUAL" << std::endl;
			exit(0);
		}
		if (basis.poly(i2, u1, uv) != func2) {
			std::cout << "NOT EQUAL" << std::endl;
			exit(0);
		}
		assert(basis.poly(i1, u1, uv) == func1);
		assert(basis.poly(i2, u1, uv) == func2);*/
		timer.stop();
		time0 += timer.getElapsedTimeInMilliSec();
		//std::cout << "degree, "<<degree << std::endl;
		//std::cout << "func1 and func2" << std::endl;
		//print_vector(func1);
		//print_vector(func2);
		timer.start();
		func1 = polynomial_differential(func1, partial1);
		func2 = polynomial_differential(func2, partial2);
		timer.stop();
		time1 += timer.getElapsedTimeInMilliSec();
		//std::cout << "differencial" << std::endl;
		//print_vector(func1);
		//print_vector(func2);
		timer.start();
		std::vector<double> func = ply_operations<Tp, knotT, valueT>::polynomial_times(func1, func2);
		//std::cout << "times" << std::endl;
		//print_vector(func);
		double upper = u2;
		if (u2 == U.back()) {
			upper = U.back() - SCALAR_ZERO;
		}
		double result = ply_operations<Tp, knotT, valueT>::polynomial_integration(func, u1, upper);
		timer.stop();
		time2 += timer.getElapsedTimeInMilliSec();

		return result;
	}
	template<typename Tp, typename knotT, typename valueT>
	double construct_an_integration(const int degree, const std::vector<double>& U,
		const int partial1, const int partial2, const int i1, const int i2, const double u1, const double u2,
		PartialBasis<Tp, knotT, valueT> &basis, const bool uv) {
		timer.start();
		//std::vector<double> func1 = basis.poly(i1, u1, uv);
		//std::vector<double> func2 = basis.poly(i2, u1, uv);	
		timer.stop();
		time0 += timer.getElapsedTimeInMilliSec();

		timer.start();
		std::vector<double> func1 = basis.poly(i1, u1, uv, partial1);
		std::vector<double> func2 = basis.poly(i2, u1, uv, partial2);
		timer.stop();
		time1 += timer.getElapsedTimeInMilliSec();

		timer.start();
		std::vector<double> func = ply_operations<Tp, knotT, valueT>::polynomial_times(func1, func2);

		double upper = u2;
		if (u2 == U.back()) {
			upper = U.back() - SCALAR_ZERO;
		}
		double result = ply_operations<Tp, knotT, valueT>::polynomial_integration(func, u1, upper);
		timer.stop();
		time2 += timer.getElapsedTimeInMilliSec();

		return result;
	}
	template<typename Tp, typename knotT, typename valueT>
	// construct an integration of multiplication of two B-spline basis (intergration of partial(Ni1)*partial(Ni2))
	// the integration domain is [u1, u2]
	double construct_an_integration(const int degree, const std::vector<double>& U,
		const int partial1, const int partial2, const int i1, const int i2, const double u1, const double u2) {
		timer.start();
		std::vector<Tp> func1 = Nip_func(i1, degree, u1, U);
		std::vector<Tp> func2 = Nip_func(i2, degree, u1, U);
		timer.stop();
		time0 += timer.getElapsedTimeInMilliSec();
		//std::cout << "degree, "<<degree << std::endl;
		//std::cout << "func1 and func2" << std::endl;
		//print_vector(func1);
		//print_vector(func2);
		timer.start();
		func1 = polynomial_differential(func1, partial1);
		func2 = polynomial_differential(func2, partial2);
		timer.stop();
		time1 += timer.getElapsedTimeInMilliSec();
		//std::cout << "differencial" << std::endl;
		//print_vector(func1);
		//print_vector(func2);
		timer.start();
		std::vector<double> func = ply_operations<Tp, knotT, valueT>::polynomial_times(func1, func2);
		//std::cout << "times" << std::endl;
		//print_vector(func);
		double upper = u2;
		if (u2 == U.back()) {
			upper = U.back() - SCALAR_ZERO;
		}
		double result = ply_operations<Tp, knotT, valueT>::polynomial_integration(func, u1, upper);
		timer.stop();
		time2 += timer.getElapsedTimeInMilliSec();

		return result;
	}
	
	// do partial difference to Pi, the cofficient of jth element Pj.
	template<typename Tp, typename knotT, typename valueT>
	double surface_energy_least_square(Bsurface<Tp, knotT, valueT>& surface, const int i, const int j, PartialBasis<Tp, knotT, valueT>& basis) {
		// figure out which Pij corresponding to the ith control point
		int partial_i = i / (surface.nv() + 1);
		int partial_j = i - partial_i * (surface.nv() + 1);

		// figure out which Pij corresponding to the jth control point
		int coff_i = j / (surface.nv() + 1);
		int coff_j = j - coff_i * (surface.nv() + 1);

		// if do partial Pij, the related other Pi1j1 will be: i1 = i-p~i+p, j1= j-q~j+q
		int degree1 = surface.degree1, degree2 = surface.degree2;
		if (coff_i<partial_i - degree1 || coff_i>partial_i + degree1) return 0;
		if (coff_j<partial_j - degree2 || coff_j>partial_j + degree2) return 0;


		// do partial Pij
		// for each block, (U_k1, U_(k1+1)), (V_k2, V_(k2+1)). the related control points are k1-p,...,k1 and k2-q,..., k2
		double result = 0;
		for (int k1 = partial_i; k1 < partial_i + degree1 + 1; k1++) {
			for (int k2 = partial_j; k2 < partial_j + degree2 + 1; k2++) {
				//std::cout << " k1, k2 " << k1 << " " << k2 << std::endl;
				if (coff_i<k1 - degree1 || coff_i>k1 || coff_j<k2 - degree2 || coff_j>k2) {
					continue;
				}
				assert(k1 + 1 < surface.U.size());
				assert(k2 + 1 < surface.V.size());
				if (surface.U[k1] == surface.U[k1 + 1] || surface.V[k2] == surface.V[k2 + 1]) {
					continue;// if the area is 0, then no need to compute
				}
				// Suu part
				double value1 = construct_an_integration(degree1, surface.U, 2, 2,
					partial_i, coff_i, surface.U[k1], surface.U[k1 + 1], basis, false);
				double value2 = construct_an_integration(degree2, surface.V, 0, 0,
					partial_j, coff_j, surface.V[k2], surface.V[k2 + 1], basis, true);
				// Suv part
				double value3 = construct_an_integration(degree1, surface.U, 1, 1,
					partial_i, coff_i, surface.U[k1], surface.U[k1 + 1], basis, false);
				double value4 = construct_an_integration(degree2, surface.V, 1, 1,
					partial_j, coff_j, surface.V[k2], surface.V[k2 + 1], basis, true);
				// Svv part
				double value5 = construct_an_integration(degree1, surface.U, 0, 0,
					partial_i, coff_i, surface.U[k1], surface.U[k1 + 1], basis, false);
				double value6 = construct_an_integration(degree2, surface.V, 2, 2,
					partial_j, coff_j, surface.V[k2], surface.V[k2 + 1], basis, true);
				result += 2 * value1*value2 + 4 * value3*value4 + 2 * value5*value6;
			}
		}
		return result;

	}

	// do partial difference to Pi, the cofficient of jth element Pj.
	template<typename Tp, typename knotT, typename valueT>
	double surface_energy_least_square_tripletes(Bsurface<Tp, knotT, valueT>& surface, const int partial_i, const int partial_j,
	const int coff_i, const int coff_j, PartialBasis<Tp, knotT, valueT>& basis) {
		// figure out which Pij corresponding to the ith control point
		

		// if do partial Pij, the related other Pi1j1 will be: i1 = i-p~i+p, j1= j-q~j+q
		int degree1 = surface.degree1, degree2 = surface.degree2;
		if (coff_i<partial_i - degree1 || coff_i>partial_i + degree1) return 0;
		if (coff_j<partial_j - degree2 || coff_j>partial_j + degree2) return 0;


		// do partial Pij
		// for each block, (U_k1, U_(k1+1)), (V_k2, V_(k2+1)). the related control points are k1-p,...,k1 and k2-q,..., k2
		double result = 0;
		for (int k1 = partial_i; k1 < partial_i + degree1 + 1; k1++) {
			for (int k2 = partial_j; k2 < partial_j + degree2 + 1; k2++) {
				//std::cout << " k1, k2 " << k1 << " " << k2 << std::endl;
				if (coff_i<k1 - degree1 || coff_i>k1 || coff_j<k2 - degree2 || coff_j>k2) {
					continue;
				}
				assert(k1 + 1 < surface.U.size());
				assert(k2 + 1 < surface.V.size());
				if (surface.U[k1] == surface.U[k1 + 1] || surface.V[k2] == surface.V[k2 + 1]) {
					continue;// if the area is 0, then no need to compute
				}
				// Suu part
				double value1 = construct_an_integration(degree1, surface.U, 2, 2,
					partial_i, coff_i, surface.U[k1], surface.U[k1 + 1], basis, false);
				double value2 = construct_an_integration(degree2, surface.V, 0, 0,
					partial_j, coff_j, surface.V[k2], surface.V[k2 + 1], basis, true);
				// Suv part
				double value3 = construct_an_integration(degree1, surface.U, 1, 1,
					partial_i, coff_i, surface.U[k1], surface.U[k1 + 1], basis, false);
				double value4 = construct_an_integration(degree2, surface.V, 1, 1,
					partial_j, coff_j, surface.V[k2], surface.V[k2 + 1], basis, true);
				// Svv part
				double value5 = construct_an_integration(degree1, surface.U, 0, 0,
					partial_i, coff_i, surface.U[k1], surface.U[k1 + 1], basis, false);
				double value6 = construct_an_integration(degree2, surface.V, 2, 2,
					partial_j, coff_j, surface.V[k2], surface.V[k2 + 1], basis, true);
				result += 2 * value1*value2 + 4 * value3*value4 + 2 * value5*value6;
			}
		} 
		return result;

	}
	// in interval [U[i], U[i+1])x[V[j], V[j+1])
	template<typename Tp, typename knotT, typename valueT>
	double discrete_surface_partial_value_squared(const int partial1, const int partial2,
		const int i, const int j, Bsurface<Tp, knotT, valueT>& surface,
		PartialBasis<Tp, knotT, valueT>& basis, const double u, const double v) {
		int p = surface.degree1;
		int q = surface.degree2;
		Eigen::VectorXd Nl(p + 1);
		Eigen::VectorXd Nr(q + 1);
		for (int k = 0; k < p + 1; k++) {
			Nl[k] = ply_operations<Tp, knotT, valueT>::polynomial_value(basis.poly(i - p + k, u, 0, partial1), u);
		}
		for (int k = 0; k < q + 1; k++) {
			Nr[k] = ply_operations<Tp, knotT, valueT>::polynomial_value(basis.poly(j - q + k, v, 1, partial2), v);
		}
		Eigen::MatrixXd px(p + 1, q + 1), py(p + 1, q + 1), pz(p + 1, q + 1);
		for (int k1 = 0; k1 < p + 1; k1++) {
			for (int k2 = 0; k2 < q + 1; k2++) {
				px(k1, k2) = surface.control_points[i - p + k1][j - q + k2][0];
				py(k1, k2) = surface.control_points[i - p + k1][j - q + k2][1];
				pz(k1, k2) = surface.control_points[i - p + k1][j - q + k2][2];
			}
		}
		double x = (Nl.transpose()*px*Nr);
		double y = (Nl.transpose()*py*Nr);
		double z = (Nl.transpose()*pz*Nr);
		return x * x + y * y + z * z;
	}

	// calculate thin-plate-energy in region [Ui, U(i+1)]x[Vj, V(j+1)]
	template<typename Tp, typename knotT, typename valueT>
	Eigen::MatrixXd Bsurface<Tp, knotT, valueT>::surface_energy_calculation(Bsurface<Tp, knotT, valueT>& surface, PartialBasis<Tp, knotT, valueT>& basis,
		const int discrete, Eigen::MatrixXd &energy_uu, Eigen::MatrixXd &energy_vv, Eigen::MatrixXd& energy_uv) {
		int p = surface.degree1;
		int q = surface.degree2;
		std::vector<double> U = surface.U;
		std::vector<double> V = surface.V;
		int nu = surface.nu(); // U[p]=0, U[nu+1]=1
		int nv = surface.nv();
		int uint = nu + 1 - p; //nbr of u intervals
		int vint = nv + 1 - q;
		int n_sample = discrete + 2;// in each interval, there are discrete+2 sample points
		energy_uu.resize(uint, vint);
		energy_vv.resize(uint, vint);
		energy_uv.resize(uint, vint);
		for (int i = 0; i < uint; i++) {
			for (int j = 0; j < vint; j++) {
				double u0 = U[i + p];
				double u1 = U[i + p + 1];
				double v0 = V[j + q];
				double v1 = V[j + q + 1];
				double delta_u = (u1 - u0) / (n_sample - 1);
				double delta_v = (v1 - v0) / (n_sample - 1);
				int Ni = i + p;// interval is [U[Ni], U[Ni+1])
				int Nj = j + q;
				Eigen::MatrixXd values_uu(n_sample, n_sample);
				Eigen::MatrixXd values_uv(n_sample, n_sample);
				Eigen::MatrixXd values_vv(n_sample, n_sample);
				//std::cout << "u0, u1, v0, v1 " << u0 << " " << u1 << " " << v0 << " " << v1 << std::endl;
				for (int k1 = 0; k1 < n_sample; k1++) {
					//std::cout << "XXXXXXXXXXXXXXXXXXXXXXXXXX" << std::endl;
					for (int k2 = 0; k2 < n_sample; k2++) {
						//int Ni = k1 == n_sample - 1 ? i + p + 1 : i + p;// if select the last point, 
						double uratio = k1 < n_sample - 1 ? 1 : 1 - SCALAR_ZERO;
						double vratio = k2 < n_sample - 1 ? 1 : 1 - SCALAR_ZERO;
						double uvalue = u0 + delta_u * k1*uratio;
						double vvalue = v0 + delta_v * k2*vratio;
						//std::cout << "uv values " << uvalue << " " << vvalue << std::endl;
						// Suu
						values_uu(k1, k2) = discrete_surface_partial_value_squared(2, 0, Ni, Nj, surface, basis, uvalue, vvalue);
						//Svv
						values_vv(k1, k2) = discrete_surface_partial_value_squared(0, 2, Ni, Nj, surface, basis, uvalue, vvalue);
						//Suv
						values_uv(k1, k2) = discrete_surface_partial_value_squared(1, 1, Ni, Nj, surface, basis, uvalue, vvalue);

					}
				}
				double uusum = 0;
				double vvsum = 0;
				double uvsum = 0;
				double single_area = delta_u * delta_v;
				for (int k1 = 0; k1 < n_sample - 1; k1++) {
					for (int k2 = 0; k2 < n_sample - 1; k2++) {
						uusum += (values_uu(k1, k2) + values_uu(k1, k2 + 1) + values_uu(k1 + 1, k2) + values_uu(k1 + 1, k2 + 1)) / 4;
						vvsum += (values_vv(k1, k2) + values_vv(k1, k2 + 1) + values_vv(k1 + 1, k2) + values_vv(k1 + 1, k2 + 1)) / 4;
						uvsum += (values_uv(k1, k2) + values_uv(k1, k2 + 1) + values_uv(k1 + 1, k2) + values_uv(k1 + 1, k2 + 1)) / 4;

					}
				}
				energy_uu(i, j) = uusum * single_area;
				energy_vv(i, j) = vvsum * single_area;
				energy_uv(i, j) = uvsum * single_area;
			}
		}
		Eigen::MatrixXd energy(uint, vint);
		energy = energy_uu + 2 * energy_uv + energy_vv;
		return energy;
		std::cout << "energy\n" << energy << std::endl;
	}

	// [U[which],U[which+1]) is the problematic one
	template<typename Tp, typename knotT, typename valueT>
	void Bsurface<Tp, knotT, valueT>::detect_max_energy_interval(Bsurface& surface, const Eigen::MatrixXd& energy, const Eigen::MatrixXd &energy_uu,
		const Eigen::MatrixXd & energy_vv, bool& uorv, int &which, double& em) {
		int k1, k2;
		em = 0;
		for (int i = 0; i < energy.rows(); i++) {
			for (int j = 0; j < energy.cols(); j++) {
				double evalue = energy(i, j);
				if (evalue > em) {
					em = evalue;
					k1 = i;
					k2 = j;
				}
			}
		}
		std::cout << "max energy " << em << std::endl;
		int p = surface.degree1;
		int q = surface.degree2;
		int u_interval = k1 + p;
		int v_interval = k2 + q;
		if (energy_uu(k1, k2) > energy_vv(k1, k2)) {
			uorv = 1;
			which = v_interval;
		}
		else {
			which = u_interval;
			uorv = 0;
		}
		return;
	}


	// which_part = 0: Suu; which_part = 1, Suv; which_part = 2, Svv.
	template<typename Tp, typename knotT, typename valueT>
	Eigen::MatrixXd energy_part_of_surface_least_square(Bsurface<Tp, knotT, valueT>& surface, PartialBasis<Tp, knotT, valueT>& basis) {
		int psize = (surface.nu() + 1)*(surface.nv() + 1);// total number of control points.
		Eigen::MatrixXd result(psize, psize);
		for (int i = 0; i < psize; i++) {
			//std::cout << "the ith row of matrix" << std::endl;
			for (int j = 0; j < psize; j++) {
				result(i, j) = surface_energy_least_square(surface, i, j, basis);
			}
		}
		std::cout << "energy matrix finish calculation" << std::endl;
		return result;
	}
	template<typename Tp, typename knotT, typename valueT>
	void energy_part_of_surface_least_square(Bsurface<Tp, knotT, valueT>& surface, PartialBasis<Tp, knotT, valueT>& basis, std::vector<Trip>& tripletes)
	{
		int psize = (surface.nu() + 1)*(surface.nv() + 1);// total number of control points.
		int degree1 = surface.degree1, degree2 = surface.degree2;
		for(int partial_i = 0;partial_i<surface.nu()+1; partial_i++)
		{
			for(int partial_j=0;partial_j<surface.nv()+1;partial_j++)
			{
				int i = partial_i * (surface.nv()+ 1) + partial_j;
				if(i>=psize)
				{
					std::cout<<"i out of range\n";
				}
				for(int coff_i = partial_i-degree1;coff_i<=partial_i+degree1;coff_i++)
				{
					if(coff_i<0||coff_i>=surface.nu()+1)
					{
						continue;
					}
					for(int coff_j = partial_j-degree2;coff_j<=partial_j+degree2;coff_j++)
					{
						if(coff_j<0||coff_j>=surface.nv()+1)
						{
							continue;
						}
						int j = coff_i*(surface.nv()+1) + coff_j;
						if(j>=psize)
				{
					std::cout<<"j out of range\n";
				}
						double value = surface_energy_least_square_tripletes(surface, partial_i, partial_j, coff_i, coff_j, basis);
						tripletes.push_back(Trip(i,j,value));
					}
				}
			}
		}

	}
	template<typename Tp, typename knotT, typename valueT>
	double surface_error_least_square(Bsurface<Tp, knotT, valueT>& surface, const int i, const int j,
		const Eigen::MatrixXd& paras) {
		// figure out which Pij corresponding to the ith control point
		int partial_i = i / (surface.nv() + 1);
		int partial_j = i - partial_i * (surface.nv() + 1);

		// figure out which Pij corresponding to the jth control point
		int coff_i = j / (surface.nv() + 1);
		int coff_j = j - coff_i * (surface.nv() + 1);

		// if do partial Pij, the related other Pi1j1 will be: i1 = i-p~i+p, j1= j-q~j+q
		int degree1 = surface.degree1, degree2 = surface.degree2;
		/*if (coff_i<partial_i - degree1 || coff_i>partial_i + degree1) return 0;
		if (coff_j<partial_j - degree2 || coff_j>partial_j + degree2) return 0;*/


		// do partial Pij
		// for each block, (U_k1, U_(k1+1)), (V_k2, V_(k2+1)). the related control points are k1-p,...,k1 and k2-q,..., k2
		double result = 0;
		for (int k1 = 0; k1 < paras.rows(); k1++) {
			double u = paras(k1, 0);
			double v = paras(k1, 1);
			double N1 = Nip(partial_i, degree1, u, surface.U);
			double N2 = Nip(partial_j, degree2, v, surface.V);
			double N3 = Nip(coff_i, degree1, u, surface.U);
			double N4 = Nip(coff_j, degree2, v, surface.V);
			if (N1 == 0 || N2 == 0 || N3 == 0 || N4 == 0) {
				continue;
			}
			result += N1 * N2*N3*N4;
		}
		return result;

	}
	template<typename Tp, typename knotT, typename valueT>
	Eigen::MatrixXd error_part_of_surface_least_square(Bsurface<Tp, knotT, valueT>& surface, const Eigen::MatrixXd& paras) {
		// figure out which Pij corr) {
		int psize = (surface.nu() + 1)*(surface.nv() + 1);// total number of control points.
		Eigen::MatrixXd result(psize, psize);
		for (int i = 0; i < psize; i++) {
			//std::cout << "the ith row of matrix" << std::endl;
			for (int j = 0; j < psize; j++) {
				result(i, j) = surface_error_least_square(surface, i, j, paras);
			}
		}
		std::cout << "error matrix finish calculation" << std::endl;
		return result;
	}
	template<typename Tp, typename knotT, typename valueT>
	Eigen::VectorXd right_part_of_least_square_approximation(Bsurface<Tp, knotT, valueT>& surface, const Eigen::MatrixXd& paras,
		const Eigen::MatrixXd& ver, const int dimension) {

		int psize = (surface.nu() + 1)*(surface.nv() + 1);// total number of control points.
		Eigen::VectorXd result(psize);
		for (int i = 0; i < psize; i++) {
			double res = 0;
			for (int j = 0; j < paras.rows(); j++) {
				int partial_i = i / (surface.nv() + 1);
				int partial_j = i - partial_i * (surface.nv() + 1);
				int degree1 = surface.degree1, degree2 = surface.degree2;

				double u = paras(j, 0);
				double v = paras(j, 1);
				double N1 = Nip(partial_i, degree1, u, surface.U);
				double N2 = Nip(partial_j, degree2, v, surface.V);
				if (N1 == 0 || N2 == 0) {
					continue;
				}
				res += N1 * N2*ver(j, dimension);
			}
			result(i) = res;
		}
		return result;
	}
	template<typename Tp, typename knotT, typename valueT>
	Eigen::MatrixXd eqality_part_of_surface_least_square(Bsurface<Tp, knotT, valueT>& surface, const Eigen::MatrixXd& paras) {
		int psize = (surface.nu() + 1)*(surface.nv() + 1);// total number of control points.
		Eigen::MatrixXd result(paras.rows(), psize);
		int degree1 = surface.degree1;
		int degree2 = surface.degree2;
		std::vector<double> U = surface.U;
		std::vector<double> V = surface.V;
		for (int i = 0; i < result.rows(); i++) {
			for (int j = 0; j < result.cols(); j++) {
				// figure out the jth control point corresponding to which Pij
				int coff_i = j / (surface.nv() + 1);
				int coff_j = j - coff_i * (surface.nv() + 1);
				double u = paras(i, 0);
				double v = paras(i, 1);
				// the corresponding cofficient should be N_coffi(u) and N_coffj(v)
				double N1 = Nip(coff_i, degree1, u, U);
				double N2 = Nip(coff_j, degree2, v, V);
				result(i, j) = N1 * N2;
			}
		}
		return result;
	}
	// this function generate ld and ru.
	template<typename Tp, typename knotT, typename valueT>
	void eqality_part_of_surface_least_square(Bsurface<Tp, knotT, valueT>& surface, const Eigen::MatrixXd& paras, int shifti, int shiftj, std::vector<Trip>& tripletes) {
		int psize = (surface.nu() + 1)*(surface.nv() + 1);// total number of control points.
		// Eigen::MatrixXd result(paras.rows(), psize);
		int degree1 = surface.degree1;
		int degree2 = surface.degree2;
		std::vector<double> U = surface.U;
		std::vector<double> V = surface.V;
		for (int i = 0; i < paras.rows(); i++) {
			for (int j = 0; j < psize + 1; j++) {
				// figure out the jth control point corresponding to which Pij
				int coff_i = j / (surface.nv() + 1);
				int coff_j = j - coff_i * (surface.nv() + 1);
				double u = paras(i, 0);
				double v = paras(i, 1);
				// the corresponding cofficient should be N_coffi(u) and N_coffj(v)
				double N1 = Nip(coff_i, degree1, u, U);
				double N2 = Nip(coff_j, degree2, v, V);
				double value = N1 * N2;
				tripletes.push_back(Trip(i+shifti, j+shiftj, value));
				tripletes.push_back(Trip(j+shiftj, i+shifti, value));
			}
		}
		return;
	}

	template<typename Tp, typename knotT, typename valueT>
	Eigen::MatrixXd lambda_part_of_surface_least_square(Bsurface<Tp, knotT, valueT>& surface, const Eigen::MatrixXd& paras) {
		Eigen::MatrixXd A = eqality_part_of_surface_least_square(surface, paras);
		Eigen::MatrixXd result = -A.transpose();
		return result;
	}
	template<typename Tp, typename knotT, typename valueT>
	Eigen::MatrixXd surface_least_square_lambda_multiplier_left_part(Bsurface<Tp, knotT, valueT>& surface,
		const Eigen::MatrixXd& paras, PartialBasis<Tp, knotT, valueT>& basis) {
		std::cout << "inside left part" << std::endl;
		int psize = (surface.nu() + 1)*(surface.nv() + 1);// total number of control points.
		int target_size = paras.rows();// nbr of target data points
		int size = psize + target_size;
		Eigen::MatrixXd result(size, size);
		Eigen::MatrixXd rd = Eigen::MatrixXd::Zero(target_size, target_size);// right down corner part
		std::cout << "finish rd" << std::endl;
		Eigen::MatrixXd lu = energy_part_of_surface_least_square(surface, basis);
		std::cout << "finish lu" << std::endl;
		Eigen::MatrixXd ld = eqality_part_of_surface_least_square(surface, paras);
		std::cout << "finish ld" << std::endl;
		Eigen::MatrixXd ru = -ld.transpose();
		std::cout << "finish ru" << std::endl;
		//lambda_part_of_surface_least_square(surface, paras);

		std::cout << "sizes" << std::endl;
		std::cout << "lu, " << lu.rows() << " " << lu.cols() << std::endl;
		std::cout << "ru, " << ru.rows() << " " << ru.cols() << std::endl;
		std::cout << "ld, " << ld.rows() << " " << ld.cols() << std::endl;
		std::cout << "rd, " << rd.rows() << " " << rd.cols() << std::endl;
		result << lu, ru,
			ld, rd;
		return result;
	}
	template<typename Tp, typename knotT, typename valueT>
	void surface_least_square_lambda_multiplier_left_part(Bsurface<Tp, knotT, valueT>& surface,
		const Eigen::MatrixXd& paras, PartialBasis<Tp, knotT, valueT>& basis, std::vector<Trip>& tripletes) {
		std::cout << "inside left part" << std::endl;
		tripletes.clear();
		
		int psize = (surface.nu() + 1)*(surface.nv() + 1);// total number of control points.
		tripletes.reserve(psize);
		int target_size = paras.rows();// nbr of target data points
		int size = psize + target_size;
		energy_part_of_surface_least_square(surface, basis, tripletes);
		std::cout<<"finished energy part\n";
		int shifti = psize, shiftj = 0;
		eqality_part_of_surface_least_square(surface, paras, shifti, shiftj, tripletes);
		std::cout<<"finished equality part\n";

		
		return;
	}
	template<typename Tp, typename knotT, typename valueT>
	Eigen::MatrixXd surface_least_square_lambda_multiplier_right_part(Bsurface<Tp, knotT, valueT>& surface, const Eigen::MatrixXd& paras,
		const Eigen::MatrixXd & points, const int dimension) {
		int psize = (surface.nu() + 1)*(surface.nv() + 1);// total number of control points.
		int target_size = paras.rows();// nbr of target data points
		int size = psize + target_size;
		Eigen::MatrixXd result(size, 1);
		for (int i = 0; i < psize; i++) {
			result(i, 0) = 0;
		}
		int counter = 0;
		for (int i = psize; i < size; i++) {
			result(i, 0) = points(counter, dimension);
			counter++;
		}
		assert(counter == target_size);
		return result;
	}
	template<typename Tp, typename knotT, typename valueT>
	void push_control_point_list_into_surface(Bsurface<Tp, knotT, valueT>& surface, const std::vector<Vector3d>& cps) {
		int id = 0;
		std::vector<std::vector<Vector3d>> control;
		control.resize(surface.nu() + 1);
		int vs = surface.nv() + 1;
		for (int i = 0; i < control.size(); i++) {
			control[i].resize(vs);
			for (int j = 0; j < vs; j++) {
				control[i][j] = cps[id];
				id++;
			}
		}
		surface.control_points = control;
		return;
	}
	template<typename Tp, typename knotT, typename valueT>
	void Bsurface<Tp, knotT, valueT>::solve_control_points_for_fairing_surface(Bsurface<Tp, knotT, valueT>& surface, const Eigen::MatrixXd& paras,
		const Eigen::MatrixXd & points, PartialBasis<Tp, knotT, valueT>& basis) {
		//using namespace Eigen;
		typedef Eigen::SparseMatrix<double> SparseMatrixXd;
		assert(paras.rows() == points.rows());
		int psize = (surface.nu() + 1)*(surface.nv() + 1);// total number of control points.
		std::vector<Vector3d> cps(psize);// control points
		// Eigen::MatrixXd A;
		// A = surface_least_square_lambda_multiplier_left_part(surface, paras, basis);
		// std::cout<<"A\n"<<A<<"\n";

		Eigen::FullPivLU<Eigen::DenseBase<Eigen::MatrixXd>::PlainMatrix> decomp;
		std::vector<Trip> tripletes;
		surface_least_square_lambda_multiplier_left_part(surface, paras, basis, tripletes);
		SparseMatrixXd matB;
		Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
		int size = psize + paras.rows();
		matB.resize(size, size);
		// std::cout<<"matrix is solved, nbr of tripletes "<<tripletes.size()<<"\n";
		matB.setFromTriplets(tripletes.begin(), tripletes.end());
		// std::cout<<"B\n"<<matB<<"\n";
		
		solver.compute(matB);
		if (solver.info() != Eigen::Success) {
			// decomposition failed
			std::cout << "solving failed" << std::endl;
			return;
		}

		for (int i = 0; i < 3; i++) {
			Eigen::MatrixXd b = surface_least_square_lambda_multiplier_right_part(surface, paras, points, i);
			double err = 0.0;

			// solve the matrix contains the p and lambda
			std::cout << "before solving" << std::endl;
			Eigen::MatrixXd p_lambda;

			p_lambda = solver.solve(b);
			if (solver.info() != Eigen::Success) {
				std::cout << "solving failed" << std::endl;
				return;
			}



			double relative_error = (matB*p_lambda - b).norm() / b.norm(); // norm() is L2 norm
			std::cout << "after solving, error is " << relative_error << std::endl;
			push_p_lambda_vector_to_control_points(p_lambda, i, cps);
		}
		push_control_point_list_into_surface(surface, cps);
		return;
	}
	void output_timing() {
		std::cout << "get basis time " << time0 << std::endl;
		std::cout << "get differential time " << time1 << std::endl;
		std::cout << "get integration time " << time2 << std::endl;
	}
}

