#pragma once
#include <sparse_interp/Types.hpp>
#include <vector>

namespace SIBSplines
{
double Nip(const int i, const int p, const double u, const std::vector<double> &U);
double kth_derivate_Nip(const int k, const int i, const int p, const double u,
                        const std::vector<double> &U);

// parameterization for curve fitting
std::vector<double> Centripetal_parameterization(std::vector<Vector3d> &pts);
std::vector<double> Chord_parameterization(std::vector<Vector3d> &pts);
bool read_csv_data_lbl(const std::string fname, std::vector<std::vector<double>> &data);
} // namespace SIBSplines