#include <Eigen/Core>
#include <igl/file_dialog_open.h>
#include <igl/readOBJ.h>
#include <polyscope/point_cloud.h>
#include <polyscope/surface_mesh.h>
#include <sparse_interp/Types.hpp>
#include <sparse_interp/energy.h>
#include <sparse_interp/functions.hpp>

#include "methods.hpp"
#include <cmath>
#include <igl/Timer.h>
#include <igl/harmonic.h>
#include <igl/write_triangle_mesh.h>

#include <TinyAD/ScalarFunction.hh>
#include <TinyAD/Utils/LineSearch.hh>
#include <TinyAD/Utils/NewtonDecrement.hh>
#include <TinyAD/Utils/NewtonDirection.hh>

#define PARAMETER true
#define CP false
#define UDIR true
#define VDIR false

namespace SIBSplines
{

void write_control_pts(std::vector<std::vector<Vector3d>> &cps, std::string file)
{
    std::ofstream fout;
    fout.open(file);
    for (int i = 0; i < cps.size(); i++)
    {
        for (int j = 0; j < cps[i].size(); j++)
        {
            fout << "v " << cps[i][j][0] << " " << cps[i][j][1] << " " << cps[i][j][2] << std::endl;
        }
    }
    fout.close();
}

/*void write_csv(const std::string &file, const std::vector<std::string> titles,
               const std::vector<double> data)
{
    std::ofstream fout;
    fout.open(file);
    for (int i = 0; i < titles.size() - 1; i++)
    {
        fout << titles[i] << ",";
    }
    fout << titles.back() << std::endl;
    for (int i = 0; i < data.size() - 1; i++)
    {
        fout << data[i] << ",";
    }
    fout << data.back() << std::endl;
    fout.close();
}*/
void write_csv(const std::string &file, const std::vector<std::string> titles,
               const std::vector<double> data)
{
    std::ofstream fout;
    fout.open(file);

    if (!fout.is_open())
    {
        std::cerr << "Failed to open file: " << file << std::endl;
        return;
    }

    // Write headers
    for (size_t i = 0; i < titles.size(); ++i)
    {
        if (i > 0)
            fout << ",";
        fout << titles[i];
    }
    fout << std::endl;

    // Write data - handle the case where we have multiple rows
    if (data.empty())
    {
        std::cerr << "Warning: No data to write to CSV" << std::endl;
        fout.close();
        return;
    }

    // If data size is a multiple of titles size, write multiple rows
    size_t cols = titles.size();
    size_t rows = data.size() / cols;

    for (size_t row = 0; row < rows; ++row)
    {
        for (size_t col = 0; col < cols; ++col)
        {
            if (col > 0)
                fout << ",";
            fout << data[row * cols + col];
        }
        fout << std::endl;
    }

    fout.close();
}

void write_svg_pts(const std::string &file, const Eigen::MatrixXd &param)
{
    std::ofstream fout;
    fout.open(file);
    double scale = 1000;
    double dot = 10;
    double displace = 30;
    std::string color = "#ED6E46";
    std::string black = "#000000";
    fout << "<svg>" << std::endl;
    for (int i = 0; i < param.rows(); i++)
    {
        double x = displace + param(i, 0) * scale;
        double y = displace + param(i, 1) * scale;
        fout << " <circle cx=\"" << x << "\" cy = \"" << y << "\" r = \"" << dot << "\" fill = \""
             << color << "\" />" << std::endl;
    }
    fout << "<polyline points=\"";

    std::vector<double> xlist = {{0, 0, 1, 1}};
    std::vector<double> ylist = {{0, 1, 1, 0}};

    for (int i = 0; i < 4; i++)
    {
        fout << displace + xlist[i] * scale << "," << displace + ylist[i] * scale << " ";
    }
    fout << displace + xlist[0] * scale << "," << displace + ylist[0] * scale;
    fout << "\" fill = \"white\" stroke = \"" + black + "\" stroke-width = \"" + std::to_string(5) +
                "\" /> "
         << std::endl;
    /*<svg>
        <rect width = "200" height = "100" fill = "#BBC42A" / >
        < / svg>*/
    fout << "</svg>" << std::endl;
    fout.close();
}

Eigen::VectorXd list_to_vec(std::vector<double> &vars)
{
    const Eigen::Index n = static_cast<Eigen::Index>(vars.size());
    Eigen::VectorXd result(n);
    for (Eigen::Index i = 0; i < n; ++i)
    {
        result[i] = vars[i];
    }
    return result;
}

std::vector<double> vec_to_list(Eigen::VectorXd globVar)
{
    std::vector<double> list(static_cast<std::size_t>(globVar.size()));
    for (Eigen::Index i = 0; i < globVar.size(); ++i)
    {
        list[static_cast<std::size_t>(i)] = globVar(i); // or eig[i]
    }
    return list;
}

void write_svg_knot_vectors(const std::string &file, const std::vector<double> &U,
                            const std::vector<double> &V)
{
    std::ofstream fout;
    fout.open(file);
    double scale = 1000;
    double width1 = 5;
    double width2 = 3;
    double displace = 30;
    std::string color = "#000000";
    std::string color1 = "#7AA20D";
    fout << "<svg>" << std::endl;
    std::vector<double> xlist = {{0, 0, 1, 1}};
    std::vector<double> ylist = {{0, 1, 1, 0}};
    fout << "<polyline points=\"";
    for (int i = 0; i < 4; i++)
    {
        fout << displace + xlist[i] * scale << "," << displace + ylist[i] * scale << " ";
    }
    fout << displace + xlist[0] * scale << "," << displace + ylist[0] * scale;
    fout << "\" fill = \"white\" stroke = \"" + color + "\" stroke-width = \"" + std::to_string(5) +
                "\" /> "
         << std::endl;

    for (int i = 0; i < U.size(); i++)
    {
        double x1 = displace + U[i] * scale;
        double y1 = displace + 0 * scale;
        double x2 = displace + U[i] * scale;
        double y2 = displace + 1 * scale;
        fout << "<line x1=\"" + std::to_string(x1) + "\" y1 = \"" + std::to_string(y1) +
                    "\" x2=\"" + std::to_string(x2) + "\" y2=\"" + std::to_string(y2) +
                    "\" stroke = \"" + color1 + "\" stroke-width = \"" + std::to_string(width2) +
                    "\" /> "
             << std::endl;
    }
    for (int i = 0; i < V.size(); i++)
    {
        double x1 = displace + 0 * scale;
        double y1 = displace + V[i] * scale;
        double x2 = displace + 1 * scale;
        double y2 = displace + V[i] * scale;
        fout << "<line x1=\"" + std::to_string(x1) + "\" y1 = \"" + std::to_string(y1) +
                    "\" x2=\"" + std::to_string(x2) + "\" y2=\"" + std::to_string(y2) +
                    "\" stroke = \"" + color1 + "\" stroke-width = \"" + std::to_string(width2) +
                    "\" /> "
             << std::endl;
    }

    fout << "</svg>" << std::endl;
    fout.close();
}

template <typename Tp, typename knotT, typename valueT>
std::vector<knotT>
ply_operationsT<Tp, knotT, valueT>::polynomial_add(const std::vector<knotT> &poly1,
                                                   const std::vector<knotT> &poly2)
{
    int size = std::max(poly1.size(), poly2.size());
    std::vector<knotT> result(size);
    for (int i = 0; i < size; i++)
    {
        bool flag1 = i < poly1.size();
        bool flag2 = i < poly2.size();
        if (flag1 && flag2)
        {
            result[i] = poly1[i] + poly2[i];
        }
        else if (flag1)
        {
            result[i] = poly1[i];
        }
        else
        {
            result[i] = poly2[i];
        }
    }
    return result;
}

template <typename Tp, typename knotT, typename valueT>
std::vector<knotT>
ply_operationsT<Tp, knotT, valueT>::polynomial_times(const std::vector<knotT> &poly1,
                                                     const std::vector<knotT> &poly2)
{
    int size = poly1.size() + poly2.size() - 1;
    std::vector<knotT> result(size);
    for (int i = 0; i < size; i++)
    { // initialize the result
        result[i] = 0;
    }

    for (int i = 0; i < poly1.size(); i++)
    {
        for (int j = 0; j < poly2.size(); j++)
        {
            result[i + j] += poly1[i] * poly2[j];
        }
    }
    return result;
}

template <typename Tp, typename knotT, typename valueT>
std::vector<Tp> ply_operationsT<Tp, knotT, valueT>::polynomial_times(const std::vector<Tp> &poly1,
                                                                     const Tp &nbr)
{
    std::vector<Tp> result;
    result = poly1;
    for (int i = 0; i < result.size(); i++)
    {
        result[i] *= nbr;
    }

    return result;
}

template <typename Tp, typename knotT, typename valueT>
std::vector<knotT>
ply_operationsT<Tp, knotT, valueT>::polynomial_times_double(const std::vector<knotT> &poly1,
                                                            const double &nbr)
{
    std::vector<knotT> result(poly1.size());
    for (int i = 0; i < result.size(); i++)
    {
        result[i] = poly1[i] * nbr;
    }

    return result;
}

template <typename Tp, typename knotT, typename valueT>
std::vector<knotT>
ply_operationsT<Tp, knotT, valueT>::polynomial_derivative(const std::vector<knotT> &poly1)
{
    std::vector<knotT> result;
    if (poly1.size() - 1 == 0)
    {
        result = {0};
        return result;
    }
    result.resize(poly1.size() - 1);
    for (int i = 0; i < result.size(); i++)
    {
        result[i] = poly1[i + 1] * (i + 1);
    }
    return result;
}

template <typename Tp, typename knotT, typename valueT>
Tp ply_operationsT<Tp, knotT, valueT>::power(const valueT &value, const int order)
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
Tp ply_operationsT<Tp, knotT, valueT>::polynomial_value(const std::vector<knotT> &poly,
                                                        const valueT &para)
{
    Tp result = 0;
    for (int i = 0; i < poly.size(); i++)
    {
        if (i == 0)
        {
            result += poly[i];
        }
        else
        {
            result += poly[i] * power(para, i);
        }
    }
    return result;
}

template <typename Tp, typename knotT, typename valueT>
std::vector<Tp>
ply_operationsT<Tp, knotT, valueT>::polynomial_integration(const std::vector<Tp> &poly)
{
    std::vector<Tp> result(poly.size() + 1);
    result[0] = 0;
    for (int i = 1; i < result.size(); i++)
    {
        result[i] = poly[i - 1] / i;
    }
    return result;
}

template <typename Tp, typename knotT, typename valueT>
Tp ply_operationsT<Tp, knotT, valueT>::polynomial_integration(const std::vector<Tp> &poly,
                                                              const Tp &lower, const Tp &upper)
{
    Tp up = ply_operations::polynomial_value(ply_operations::polynomial_integration(poly), upper);
    Tp lw = ply_operations::polynomial_value(ply_operations::polynomial_integration(poly), lower);
    return up - lw;
}

template <typename Tp, typename knotT, typename valueT>
double splineBasis<Tp, knotT, valueT>::Ni0_func(const int i, const int uId)
{
    if (i == uId)
    {
        return 1;
    }
    return 0;
}

// this version can be really used for treating U as variables, since the equality of two knots is
// not predicted using the values, but using the degree: original_p to predict the equality of the
// ends of the clampped B-spline knot vector.
template <typename Tp, typename knotT, typename valueT>
std::vector<knotT> splineBasis<Tp, knotT, valueT>::Nip_func(const int i, const int p, const int uId,
                                                            const std::vector<knotT> &U,
                                                            const int original_p)
{
    if (p == 0)
    {
        std::cout << "Error: This function cannot handle degree = 0\n";
        exit(0);
    }
    int pnbr = U.size() - p - 1; // the nbr of basis functions under p
    if (uId < p || uId >= pnbr)
    {
        std::cout << "error in Nip_func: interval id out of range, uid, " << uId << "\n";
        exit(0);
    }

    std::vector<knotT> v;
    // the first repeatative knot region is from U_0 to U_degree, the second is from
    // U.size()-1-degree to U.size()-1
    bool firstTerm0 =
        (i + p <= original_p && i <= original_p) ||
        (i + p >= U.size() - 1 - original_p && i >= U.size() - 1 - original_p); // U[i + p] <= U[i];
    bool lastTerm0 = (i + p + 1 <= original_p && i + 1 <= original_p) ||
                     (i + p + 1 >= U.size() - 1 - original_p &&
                      i + 1 >= U.size() - 1 - original_p); // U[i + p + 1] <= U[i + 1]; // if the
                                                           // first or the second term vanishes
    bool degreeIs0 = p == 1;                               // if the lower level degree is 0
    if (firstTerm0 && !lastTerm0)                          // only keep the last term
    {
        v = {{U[i + p + 1] / (U[i + p + 1] - U[i + 1]),
              -1 / (U[i + p + 1] - U[i + 1])}}; // U[i+p+1] - u
        if (!degreeIs0)
        {
            return PO.polynomial_times(v, Nip_func(i + 1, p - 1, uId, U, original_p));
        }
        else
        {
            return PO.polynomial_times_double(v, Ni0_func(i + 1, uId));
        }
    }
    if (lastTerm0 && !firstTerm0) // only keep the first term
    {
        v = {{-U[i] / (U[i + p] - U[i]), 1 / (U[i + p] - U[i])}}; // u - U[i]
        if (!degreeIs0)
        {
            return PO.polynomial_times(v, Nip_func(i, p - 1, uId, U, original_p));
        }
        else
        {
            return PO.polynomial_times_double(v, Ni0_func(i, uId));
        }
    }
    if (firstTerm0 && lastTerm0)
    {
        std::cout << "impossible case in Nip_func\n";
        exit(0);
    }
    // division can be properly handled, thus both terms are kept.
    if (!degreeIs0)
    {
        v = {{-U[i] / (U[i + p] - U[i]), 1 / (U[i + p] - U[i])}}; // u - U[i]
        std::vector<knotT> result1 = PO.polynomial_times(v, Nip_func(i, p - 1, uId, U, original_p));

        v = {{U[i + p + 1] / (U[i + p + 1] - U[i + 1]),
              -1 / (U[i + p + 1] - U[i + 1])}}; // U[i+p+1] - u
        std::vector<knotT> result2 =
            PO.polynomial_times(v, Nip_func(i + 1, p - 1, uId, U, original_p));
        return PO.polynomial_add(result1, result2);
    }
    v = {{-U[i] / (U[i + p] - U[i]), 1 / (U[i + p] - U[i])}}; // u - U[i]
    std::vector<knotT> result1 = PO.polynomial_times_double(v, Ni0_func(i, uId));

    v = {
        {U[i + p + 1] / (U[i + p + 1] - U[i + 1]), -1 / (U[i + p + 1] - U[i + 1])}}; // U[i+p+1] - u
    std::vector<knotT> result2 = PO.polynomial_times_double(v, Ni0_func(i + 1, uId));
    return PO.polynomial_add(result1, result2);
}

template <typename Tp, typename knotT, typename valueT>
std::vector<Tp>
splineBasis<Tp, knotT, valueT>::computeBasisFunctionValues(const valueT &value, const int uId,
                                                           const int p, const std::vector<knotT> &U)
{
    std::vector<Tp> result(p + 1, 0);
    // first compute the associated basis functions. Their type is the same as the knot vector U
    std::vector<std::vector<knotT>> basisFunctions(p + 1);
    for (int i = 0; i < p + 1; i++)
    {
        basisFunctions[i] = Nip_func(uId - p + i, p, uId, U, p);
        result[i] = PO.polynomial_value(basisFunctions[i], value);
    }
    return result;
}
template <typename Tp, typename knotT, typename valueT>
void splineBasis<Tp, knotT, valueT>::computeBasisFunctionDerivativeValues(
    const valueT &value, const int uId, const int p, const int order, const std::vector<knotT> &U,
    std::vector<Tp> &bvalues, std::vector<Tp> &dvalues)
{

    bvalues = std::vector<Tp>(p + 1, 0);
    dvalues = std::vector<Tp>(p + 1, 0);
    std::vector<std::vector<knotT>> basisFunctions(p + 1), bfds(p + 1);
    for (int i = 0; i < p + 1; i++)
    {
        basisFunctions[i] = Nip_func(uId - p + i, p, uId, U, p);
        bfds[i] = basisFunctions[i];
        for (int r = 0; r < order; r++)
            bfds[i] = PO.polynomial_derivative(bfds[i]);

        bvalues[i] = PO.polynomial_value(basisFunctions[i], value);
        dvalues[i] = PO.polynomial_value(bfds[i], value);
    }
    return;
}

template <typename Tp, typename knotT, typename valueT>
void splineBasis<Tp, knotT, valueT>::computeBasisFunctionDerivativeAndValues(
    const valueT &value, const int uId, const int p, const std::vector<knotT> &U,
    std::vector<std::vector<knotT>> &basisfuncs, std::vector<std::vector<knotT>> &d1funcs,
    std::vector<std::vector<knotT>> &d2funcs, std::vector<Tp> &bvalues, std::vector<Tp> &d1values,
    std::vector<Tp> &d2values)
{
    basisfuncs.resize(p + 1);
    d1funcs.resize(p + 1);
    d2funcs.resize(p + 1);
    bvalues.resize(p + 1);
    d1values.resize(p + 1);
    d2values.resize(p + 1);
    for (int i = 0; i < p + 1; i++)
    {
        basisfuncs[i] = Nip_func(uId - p + i, p, uId, U, p);
        d1funcs[i] = PO.polynomial_derivative(basisfuncs[i]);
        d2funcs[i] = PO.polynomial_derivative(d1funcs[i]);
        bvalues[i] = PO.polynomial_value(basisfuncs[i], value);
        d1values[i] = PO.polynomial_value(d1funcs[i], value);
        d2values[i] = PO.polynomial_value(d2funcs[i], value);
    }
    return;
}

// Rescales parameters from [a,b] -> [a+epsilon, b-epsilon]
void rescale_param(Eigen::MatrixXd &param)
{
    double epsilon = 1e-6;
    for (int i = 0; i < param.rows(); i++)
    {
        for (int j = 0; j < 2; j++)
        {
            param(i, j) = epsilon + (1 - 2 * epsilon) * param(i, j);
        }
    }
}

// Perform a binary search to find the lower bound of the knot vector that fits the parameter.
// Return the knot vector index
int return_closest_knot_index_to_param(const std::vector<double> &UV, const double param)
{
    int index;
    auto it = std::lower_bound(UV.begin(), UV.end(), param);
    if (it == UV.begin())
    {
        std::cerr << "Error: No interval found for " << param << std::endl;
    }
    index = std::distance(std::begin(UV), it);

    return index - 1;
}

int parameterLocation(const bool Udirection, const int pos, const Bsurface &surface)
{
    if (Udirection)
    {
        return surface.cpSize * 3 + pos;
    }
    else
    {
        return surface.cpSize * 3 + pos + surface.paramSize;
    }
}

// Gives you the varaible depending on U or V direction, variable type (control point or parameter)
int variableMap(const bool Udirection, const bool ParOrCp, const int pos, const int pos2,
                const int xyz, const Bsurface &surface)
{
    if (ParOrCp)
    {
        return parameterLocation(Udirection, pos, surface);
    }
    // Control point location
    else
    {
        return pos * surface.control_points[0].size() + pos2 + xyz * surface.cpSize;
    }
}

Eigen::SparseMatrix<double> make_block_diagonal(const Eigen::SparseMatrix<double> &A,
                                                const Bsurface &surface)
{
    int n = A.rows();
    int m = A.cols();
    int parSize = surface.paramSize;
    Eigen::SparseMatrix<double> B(3 * n + 2 * parSize, 3 * m + 2 * parSize);
    std::vector<Eigen::Triplet<double>> triplets;
    for (int k = 0; k < A.outerSize(); ++k)
    {
        for (Eigen::SparseMatrix<double>::InnerIterator it(A, k); it; ++it)
        {
            for (int i = 0; i < 3; ++i)
            {
                int row = i * n + it.row();
                int col = i * m + it.col();
                triplets.emplace_back(row, col, it.value());
            }
        }
    }

    B.setFromTriplets(triplets.begin(), triplets.end());
    return B;
}
// Generates the gradient and hessian for the fairing energy
std::tuple<double, Eigen::VectorXd, SparseMatrixXd>
calculate_fairing_energy(Bsurface &surface, PartialBasis &basis, Eigen::VectorXd x)
{
    int psize = (surface.nu() + 1) * (surface.nv() + 1); // total number of control points.
    std::vector<Trip> tripletes;
    energy_part_of_surface_least_square(surface, basis, tripletes);
    SparseMatrixXd matE, matEx30, hessE;
    Eigen::VectorXd gradE;
    matE.resize(psize, psize);
    matE.setFromTriplets(tripletes.begin(), tripletes.end());
    // thin plate energy hessian
    matEx30 = make_block_diagonal(matE, surface);
    hessE = matEx30;
    gradE = hessE * x;
    double f = 0.5 * x.dot(hessE * x);
    return std::tuple<double, Eigen::VectorXd, SparseMatrixXd>(f, gradE, hessE);
}
// Calculates the backtracking alpha rescaling
double calculate_alpha(const Eigen::VectorXd &d, const Eigen::VectorXd &x, const int i,
                       const double upper, const double lower)
{
    double alpha;
    double epsilon = 1e-3;
    double updateVar = d(i) + x(i);
    if (updateVar >= upper)
    {
        if (abs(d(i)) < 1e-6)
            return 0;
        alpha = (upper - x(i)) / d(i);
        for (int i1 = 0; i1 < 64; i1++)
        {
            if (alpha * d(i) + x(i) < upper)
            {
                break;
            }
            alpha *= 0.9;
        }
    }
    else if (updateVar < lower)
    {
        if (abs(d(i)) < 1e-6)
            return 0;
        alpha = (lower - x(i)) / d(i);
        for (int i1 = 0; i1 < 64; i1++)
        {
            if (alpha * d(i) + x(i) >= lower)
            {
                break;
            }
            alpha *= 0.9;
        }
    }
    else
    {
        return 1 - epsilon;
    }
    double updateNew = d(i) * alpha + x(i);
    if (updateNew < upper && updateNew >= lower)
    {
        return alpha * (1 - epsilon);
    }
    else
        return 0;
}
// Scale the direction vector d so that every parameter + direction stays in between two knots
// determined prior.
// Should work correctly
Eigen::VectorXd stepBacktracker(Eigen::VectorXd &d,
                                const std::vector<std::array<int, 2>> &paraInInterval,
                                Bsurface &surface)
{
    Eigen::VectorXd x = list_to_vec(surface.globVars);
    // U parameters
    for (int i = 0; i < surface.paramSize; ++i)
    {
        // Get the lower and upper bounds for the current parameter
        double upper = surface.U[paraInInterval[i][0] + 1];
        double lower = surface.U[paraInInterval[i][0]];
        double alpha = calculate_alpha(d, x, i + surface.cpSize * 3, upper, lower);
        d(i + surface.cpSize * 3) = alpha * d(i + surface.cpSize * 3);
    }
    // V parameters
    for (int i = 0; i < surface.paramSize; ++i)
    {
        // Get the lower and upper bounds for the current parameter
        double upper = surface.V[paraInInterval[i][1] + 1];
        double lower = surface.V[paraInInterval[i][1]];
        double alpha =
            calculate_alpha(d, x, i + surface.cpSize * 3 + surface.paramSize, upper, lower);
        d(i + surface.cpSize * 3 + surface.paramSize) =
            alpha * d(i + surface.cpSize * 3 + surface.paramSize);
    }
    return d;
}

void reassign_control_points(Bsurface &surface, const Eigen::VectorXd &x_new)
{
    std::vector<double> globV = vec_to_list(x_new);

    const int rows = surface.control_points.size();
    const int cols = surface.control_points[0].size();
    const int cpSize = rows * cols;

    for (int k = 0; k < 3; ++k)
    {
        const int base = k * surface.cpSize; // start of the k-th component block
        for (int i = 0; i < rows; ++i)
        {
            const int rowBase = i * cols; // start of the i-th row inside that block
            for (int j = 0; j < cols; ++j)
            {
                surface.control_points[i][j](k) = surface.globVars[base + rowBase + j];
            }
        }
    }
}

void getBoundingBox(const Eigen::MatrixXd &V, Eigen::Vector3d &vmin, Eigen::Vector3d &vmax)
{
    double xmin = V(0, 0);
    double xmax = V(0, 0);
    double ymin = V(0, 1);
    double ymax = V(0, 1);
    double zmin = V(0, 2);
    double zmax = V(0, 2);
    for (int i = 0; i < V.rows(); i++)
    {
        double x = V(i, 0);
        double y = V(i, 1);
        double z = V(i, 2);
        if (xmin > x)
        {
            xmin = x;
        }
        if (xmax < x)
        {
            xmax = x;
        }
        if (ymin > y)
        {
            ymin = y;
        }
        if (ymax < y)
        {
            ymax = y;
        }
        if (zmin > z)
        {
            zmin = z;
        }
        if (zmax < z)
        {
            zmax = z;
        }
    }
    vmin = Eigen::Vector3d(xmin, ymin, zmin);
    vmax = Eigen::Vector3d(xmax, ymax, zmax);
}

void reassign_parameters(const Bsurface &surface, Eigen::MatrixXd &param,
                         const Eigen::VectorXd &x_new)
{
    const int param_nbr = param.rows();
    for (int k = 0; k < 2; ++k)
    {
        for (int i = 0; i < param_nbr; ++i)
        {
            param(i, k) = x_new(i + (3 * surface.cpSize) + (k * param_nbr));
        }
    }
}

// Fairing energy doesnt consider the new control points -> needs a fix
template <typename EvalFunctionT>
double eval_energy(const Eigen::VectorXd &x_new, Bsurface &surface, PartialBasis &basis,
                   EvalFunctionT &func, const double w_fair, SparseMatrixXd &matE30)
{
    double w_fit = 1 - w_fair;
    auto f_fit = func(x_new);
    double f_fair = x_new.dot(matE30 * x_new) / 2;
    return w_fair * f_fair + w_fit * f_fit;
}

bool armijoCheck(const double f_old, const double f_new, const double alpha,
                 const Eigen::VectorXd d, const Eigen::VectorXd g_total, const double armijo_const)
{
    return f_new <= f_old + armijo_const * alpha * d.dot(g_total);
}

template <typename EvalFunctionT>
Eigen::VectorXd lineSearch(Eigen::VectorXd x, Eigen::VectorXd d, double &f_total,
                           Eigen::VectorXd g_total, Bsurface &surface, PartialBasis &basis,
                           EvalFunctionT &func, const double w_fair, SparseMatrixXd &matE30)
{
    int max_iters = 64;
    double armijo_const = 1e-4;
    double alpha = 1.0;
    double shrink = 0.8;
    double f_old = f_total;
    Eigen::VectorXd x_new;
    // Code crashes here.
    for (int i = 0; i < max_iters; ++i)
    {
        x_new = x + alpha * d;
        double f_new = eval_energy(x_new, surface, basis, func, w_fair, matE30);
        // TINYAD_ASSERT_EQ(f_new, f_new);
        if (armijoCheck(f_old, f_new, alpha, d, g_total, armijo_const))
        {
            return x_new;
        }
        else
            alpha *= shrink;
    }
    std::cout << "Line search couldn't find improvement. Gradient max norm is "
              << g_total.cwiseAbs().maxCoeff() << std::endl;
    return x;
}

// Inits the surface object
void surface_init(const std::string meshfile, const std::string tail, double delta,
                  const double per, const int target_steps, const double w_fair,
                  const bool meshInterpolation, const int nbr_pts, const int model,
                  Bsurface &surface, Eigen::MatrixXd &param,
                  std::vector<std::array<int, 2>> &paraInInterval_orig, // Add &
                  Eigen::MatrixXd &ver_orig)
{
    // Init important data
    // Bsurface gets updated by other function.
    // Mesh variables / initialization
    Eigen::MatrixXd ver;
    int nbr;
    Eigen::MatrixXi F;

    if (meshInterpolation)
    {
        mesh_parameterization(meshfile, ver, param, F);
        nbr = nbr_pts;
    }
    else
    {
        bool corners = true;
        int method = model;
        SIBSplines::examples::get_model_sample_points(nbr_pts, ver, F, param, method, corners,
                                                      meshfile);
    }
    rescale_param(param);
    const int param_nbr = param.rows();     // the number of data points
    surface.degree1 = 3;                    // degree of u direction
    surface.degree2 = 3;                    // degree of v direction
    surface.U = {{0, 0, 0, 0, 1, 1, 1, 1}}; // the initial U knot vector
    surface.V = surface.U;                  // the initial V knot vector
    bool enable_max_fix_nbr = true; // progressively update the knot vectors to make the two knot
                                    // vectors balanced in length.
    // generate knot vectors to make sure the data points can be interpolated
    surface.generate_interpolation_knot_vectors(surface.degree1, surface.degree2, surface.U,
                                                surface.V, param, delta, per, target_steps,
                                                enable_max_fix_nbr);
    std::cout << "knot vectors generated" << std::endl;
    // Solve the control points as initialization.
    PartialBasis basis(surface);
    std::cout << "Generating control points" << std::endl;
    surface.solve_control_points_for_fairing_surface(surface, param, ver, basis);
    Eigen::MatrixXd SPs_orig;
    Eigen::MatrixXi SFs_orig;
    int visual_nbr =
        200; // the discretization scale for the output surface. The mesh will be 200x200
    surface.surface_visulization(surface, visual_nbr, SPs_orig, SFs_orig);
    double precision = surface.max_interpolation_err(ver, param, surface);
    std::string path = SI_MESH_DIR;
    write_svg_pts(path + "orig_opt_pts.svg", param);
    write_points(
        path + "pts" + std::to_string(nbr_pts) + "_m_" + std::to_string(model) + "_orig.obj", ver);
    write_triangle_mesh(path + "ours_" + "p" + std::to_string(nbr_pts) + "_m_" +
                            std::to_string(model) + "_orig.obj",
                        SPs_orig, SFs_orig);
    std::cout << "Control points initialized" << std::endl;
    // Init parameter intervals for reparameterization
    std::vector<std::array<int, 2>> paraInInterval(param_nbr, {0, 0});
    // Perform binary search and return the element
    std::cout << "Init parameter interval" << std::endl;
    for (int i = 0; i < paraInInterval.size(); i++)
    {
        paraInInterval[i][0] = intervalLocator(
            surface.U, surface.degree1,
            param(i,
                  0)); // return_closest_knot_index_to_param(surface.U, param(i, 0));
        paraInInterval[i][1] = intervalLocator(
            surface.V, surface.degree2,
            param(i,
                  1)); // return_closest_knot_index_to_param(surface.V, param(i, 1));
    }

    //  Calculate the number of variables we solve for.
    //  Number of control points
    surface.paramSize = param_nbr;
    surface.cpSize =
        (surface.U.size() - 1 - surface.degree1) * (surface.V.size() - 1 - surface.degree2);
    // Number of variables = 2 * parameters (u,v) + 3 * control points (x,y,z)
    const int varSize = 2 * param_nbr + 3 * surface.cpSize;
    // Init globVars vector
    surface.cpRows = surface.control_points.size();
    surface.cpCols = surface.control_points[0].size();

    // Adding control points to globVars
    surface.globVars.resize(varSize);
    std::cout << "Setting the control points to globVars" << std::endl;
    for (int k = 0; k < 3; ++k)
    {
        for (int i = 0; i < surface.control_points.size(); ++i)
        {
            for (int j = 0; j < surface.control_points[0].size(); ++j)
            {
                // Index is k * cpSize + ith-cp
                surface.globVars[k * surface.cpSize + i * surface.control_points[0].size() + j] =
                    surface.control_points[i][j](k);
            }
        }
    }
    std::cout << "Setting parameters to globVars" << std::endl;
    for (int k = 0; k < 2; ++k)
    {
        for (int i = 0; i < param.rows(); ++i)
        {
            surface.globVars[3 * surface.cpSize + k * param.rows() + i] = param(i, k);
        }
    }
    paraInInterval_orig = std::move(paraInInterval);
    ver_orig = ver;
}

void mesh_optimization(Bsurface &surface, PartialBasis &basis, double w_fair, const int itSteps,
                       const std::vector<std::array<int, 2>> &paraInInterval,
                       Eigen::MatrixXd &param, const int method, const Eigen::MatrixXd &ver)
{
    std::vector<double> data;
    std::vector<double> d_f_fit, d_f_fair, d_f_total;
    d_f_fit.reserve(itSteps);
    d_f_fair.reserve(itSteps);
    d_f_total.reserve(itSteps);
    auto func = TinyAD::scalar_function<1>(TinyAD::range(surface.globVars.size()));
    // ##### Fitting energy ######
    // (d+1)*(d+1) cps * 3 (x,y,z) + 2 parameters (u,v).
    // For degree 3 = 50 variables
    func.add_elements<50>(
        TinyAD::range(surface.paramSize),
        [&](auto &element) -> TINYAD_SCALAR_TYPE(element)
        {
            using T = TINYAD_SCALAR_TYPE(element);
            Eigen::Index dataID = element.handle;
            // Looking for mistakes
            T parameterU =
                element.variables(variableMap(UDIR, PARAMETER, dataID, 0, 0, surface))(0, 0);
            T parameterV =
                element.variables(variableMap(VDIR, PARAMETER, dataID, 0, 0, surface))(0, 0);

            // *************************
            // This part should be fine.
            // get the uv intervals [U[uItv], U[uItv+1])
            int uItv = paraInInterval[dataID][0];
            int vItv = paraInInterval[dataID][1];
            splineBasis<T, double, T> spb; // spline basis computation
            // get the basis functions. basisU is a vector of size degree1 + 1 to represent a
            // degree1 polynomial. each element of basisU is of tinyAD type variable
            std::vector<T> basisU =
                spb.computeBasisFunctionValues(parameterU, uItv, surface.degree1, surface.U);
            std::vector<T> basisV =
                spb.computeBasisFunctionValues(parameterV, vItv, surface.degree2, surface.V);
            // **********************************

            // Looking for mistakes
            T p00 = 0, p01 = 0, p02 = 0;
            for (int i = 0; i < surface.degree1 + 1; i++)
            {
                for (int j = 0; j < surface.degree2 + 1; j++)
                {
                    // for each interval [U[uItv], U[uItv + 1]), the associated control points are
                    // P[uItv-degree],...,P[uItv]
                    int lc0 = variableMap(false, CP, uItv - surface.degree1 + i,
                                          vItv - surface.degree2 + j, 0,
                                          surface); // getTheLocationOfThe(uItv-degree1+i,
                                                    // vItv-degree2+j)-th ControlPoint_x;
                    int lc1 = variableMap(false, CP, uItv - surface.degree1 + i,
                                          vItv - surface.degree2 + j, 1,
                                          surface); // getTheLocationOfThe(uItv-degree1+i,
                                                    // vItv-degree2+j)-th ControlPoint_y;
                    int lc2 = variableMap(false, CP, uItv - surface.degree1 + i,
                                          vItv - surface.degree2 + j, 2,
                                          surface); // getTheLocationOfThe(uItv-degree1+i,
                                                    // vItv-degree2+j)-th ControlPoint_z;
                    T pt0 = 0, pt1 = 0, pt2 = 0;
                    //  get the control point
                    pt0 = element.variables(lc0)(0, 0), pt1 = element.variables(lc1)(0, 0),
                    pt2 = element.variables(lc2)(0, 0);

                    p00 += pt0 * basisU[i] * basisV[j];
                    p01 += pt1 * basisU[i] * basisV[j];
                    p02 += pt2 * basisU[i] * basisV[j];
                }
            }
            // std::cout << "DataID: " << dataID << std::endl;
            return (p00 - ver(dataID, 0)) * (p00 - ver(dataID, 0)) +
                   (p01 - ver(dataID, 1)) * (p01 - ver(dataID, 1)) +
                   (p02 - ver(dataID, 2)) * (p02 - ver(dataID, 2));
        });
    TinyAD::LinearSolver solver;
    double convergence_eps = 1e-12; // change it into 1e-6 if you want.
    // const double w_fair = 1e-6;
    const double w_fit = 1 - w_fair;
    Eigen::VectorXd x = list_to_vec(surface.globVars);

    for (int i = 0; i < itSteps; ++i)
    {
        auto [f_fit, g_fit, H_fit_proj] = func.eval_with_hessian_proj(x);
        Eigen::VectorXd ones = Eigen::VectorXd::Ones(surface.globVars.size());
        SparseMatrixXd diag = 1e-6 * SparseMatrixXd(ones.asDiagonal());
        H_fit_proj += diag;
        auto [f_fair, g_fair, H_fair] = calculate_fairing_energy(surface, basis, x);
        Eigen::VectorXd g_total = w_fit * g_fit + w_fair * g_fair;
        SparseMatrixXd H_total = H_fit_proj;
        H_total *= w_fit;
        H_total += w_fair * H_fair;
        double f_total = w_fit * f_fit + w_fair * f_fair;
        d_f_fit.push_back(f_fit);
        d_f_fair.push_back(f_fair);
        d_f_total.push_back(f_total);
        // TINYAD_DEBUG_OUT("Energy in iteration " << i << ": " << f_total);
        Eigen::VectorXd d = TinyAD::newton_direction(g_total, H_total, solver);
        d = stepBacktracker(d, paraInInterval, surface);
        // std::cout << "Steptraced d: " << d.norm() << std::endl;

        if (TinyAD::newton_decrement(d, g_total) <
            convergence_eps) // if the direction is too far from the gradient direction, break.
                             // normally this value is set as 1e-6
            break;
        w_fair = w_fair * 0.8;
        Eigen::VectorXd prev_x = x;
        x = lineSearch(x, d, f_total, g_total, surface, basis, func, w_fair, H_fair);
        if ((x - prev_x).norm() < convergence_eps) // if the step is too small, break
        {
            std::cout << "break because the line searched step is too small: "
                      << (x - prev_x).norm() << "\n";
            break;
            // w_fair = w_fair * 0.8;
        }
        // std::cout << "the dx, " << d.norm() << ", the backtraced dx " << (x - prev_x).norm()
        //<< "\n";
        // std::cout << "8" << std::endl;
        // std::cout << "--- Done with iteration: " << i << " ---" << std::endl;
    }
    std::cout << "Optimization done!" << std::endl;
    surface.globVars = vec_to_list(x);
    reassign_control_points(surface, list_to_vec(surface.globVars));
    reassign_parameters(surface, param, list_to_vec(surface.globVars));
    std::string path = SI_MESH_DIR;
    std::vector<std::string> titles = {"f_fit", "f_fair", "f_total"};
    std::vector<double> flat;
    flat.reserve(3 * d_f_fit.size()); // optional performance hint
    for (std::size_t i = 0; i < d_f_fit.size(); ++i)
    {
        flat.push_back(d_f_fit[i]);
        flat.push_back(d_f_fair[i]);
        flat.push_back(d_f_total[i]);
    }

    write_csv(path + std::to_string(method) + "energies.csv", titles, flat);

    Eigen::MatrixXd SPs;
    Eigen::MatrixXi SFs;
    int visual_nbr =
        200; // the discretization scale for the output surface. The mesh will be 200x200
    surface.surface_visulization(surface, visual_nbr, SPs, SFs);
    double precision = surface.max_interpolation_err(ver, param, surface);

    // Write to separate CSV file
    std::vector<std::string> precision_titles = {"model", "num_points", "w_fair", "iterations",
                                                 "precision"};
    std::vector<double> precision_data = {
        static_cast<double>(method), static_cast<double>(surface.paramSize), w_fair,
        static_cast<double>(d_f_fit.size()), // actual iterations performed
        precision};
    write_csv(path + "interpolation_errors.csv", precision_titles, precision_data);
    // std::cout << "maximal interpolation error "
    //           << precision << std::endl;
    if (method == -1)
        write_points(path + "pts" + std::to_string(surface.paramSize) + "_m_" + +".obj", ver);
    else
    {
        write_points(path + "pts" + std::to_string(surface.paramSize) + "_m_" +
                         std::to_string(method) + ".obj",
                     ver);
        write_triangle_mesh(path + "ours_" + "p" + std::to_string(surface.paramSize) + "_m_" +
                                std::to_string(method) + ".obj",
                            SPs, SFs);
    }
    write_svg_pts(path + "opt_pts.svg", param);
    write_svg_knot_vectors(path + "knot_vectors.svg", surface.U, surface.V);
}

void mesh_visualization(const Eigen::MatrixXd &param, const Eigen::MatrixXd &ver, Bsurface &surface,
                        std::string path, const int method)
{
    Eigen::MatrixXd SPs;
    Eigen::MatrixXi SFs;
    int visual_nbr =
        200; // the discretization scale for the output surface. The mesh will be 200x200
    surface.surface_visulization(surface, visual_nbr, SPs, SFs);
    double precision = surface.max_interpolation_err(ver, param, surface);
    std::cout << "maximal interpolation error "
              << surface.max_interpolation_err(ver, param, surface) << std::endl;
    if (method == -1)
        write_points(path + "pts" + std::to_string(surface.paramSize) + "_m_" + +".obj", ver);
    else
    {
        write_points(path + "pts" + std::to_string(surface.paramSize) + "_m_" +
                         std::to_string(method) + ".obj",
                     ver);
        write_triangle_mesh(path + "ours_" + "p" + std::to_string(surface.paramSize) + "_m_" +
                                std::to_string(method) + ".obj",
                            SPs, SFs);
    }
    write_svg_pts(path + "pts.svg", param);
    write_svg_knot_vectors(path + "knot_vectors.svg", surface.U, surface.V);
}

// Optimized for all the variables using TinyAD as a autodifferenciation.
void old(const int model, const int nbr_pts, double &per_ours, const std::string path,
         const std::string tail, const double per, const bool enable_local_energy)
{
    // Timer variables
    igl::Timer timer;
    double time_knot = 0;
    double time_solve = 0;
    double precision = 0;

    // Mesh variables / initialization
    Eigen::MatrixXd ver;
    int nbr = nbr_pts; // nbr of points
    Eigen::MatrixXi F;
    Eigen::MatrixXd param;
    int method = model;
    bool corners = true;
    SIBSplines::examples::get_model_sample_points(nbr, ver, F, param, method, corners,
                                                  path); // Inits the mesh vars
    // Vars = all control points * 2 + knots +
    std::vector<double> vars;
    auto func = TinyAD::scalar_function<2>(TinyAD::range(nbr));

    // Splines vars init
    int degree1 = 3;
    int degree2 = 3;
    std::vector<double> Uknot = {{0, 0, 0, 0, 1, 1, 1, 1}};
    std::vector<double> Vknot = Uknot;
    int perturb_itr = 0;
    int target_steps = 10;
    bool enable_max_fix_nbr = true;

    timer.start();
    std::cout << "before generating knot vectors" << std::endl;
    std::cout << "data size " << ver.rows() << std::endl;

    Bsurface surface;
    surface.generate_interpolation_knot_vectors(degree1, degree2, Uknot, Vknot, param, per_ours,
                                                per, target_steps, enable_max_fix_nbr);

    timer.stop();
    time_knot = timer.getElapsedTimeInSec();

    Eigen::MatrixXd SPs;
    Eigen::MatrixXi SFs;
    if (1)
    {
        surface.degree1 = 3;
        surface.degree2 = 3;
        surface.U = Uknot;
        surface.V = Vknot;
        std::cout << "before initialize the basis " << std::endl;
        PartialBasis basis(surface);
        std::cout << "initialize the basis done" << std::endl;
        std::cout << "before solving control points" << std::endl;
        timer.start();
        surface.solve_control_points_for_fairing_surface(surface, param, ver, basis);
        timer.stop();
        time_solve = timer.getElapsedTimeInSec();
        surface.surface_visulization(surface, 100, SPs, SFs);
        // Energy solving part
        // if (enable_local_energy)
        // {
        //     double timeitr = 0;
        //     for (int i = 0; i < 50; i++)
        //     {
        //         timer.start();
        //         Eigen::MatrixXd energy, euu, evv, euv;
        //         energy = surface.surface_energy_calculation(surface, basis, 1, euu, evv, euv);
        //         bool uorv;
        //         int which;
        //         double max_energy;
        //         surface.detect_max_energy_interval(surface, energy, euu, evv, uorv, which,
        //                                            max_energy);
        //         std::vector<double> Unew = surface.U;
        //         std::vector<double> Vnew = surface.V;
        //         if (!uorv)
        //         { // u get updated
        //             double value = (Unew[which] + Unew[which + 1]) / 2;
        //             surface.U = knot_vector_insert_one_value(Unew, value);
        //         }
        //         else
        //         {
        //             double value = (Vnew[which] + Vnew[which + 1]) / 2;
        //             surface.V = knot_vector_insert_one_value(Vnew, value);
        //         }
        //         std::cout << "knot vector get inserted" << std::endl;
        //         basis.clear();
        //         basis.init(surface);
        //         surface.solve_control_points_for_fairing_surface(surface, param, ver, basis);
        //         timer.stop();
        //         timeitr += timer.getElapsedTimeInSec();
        //         std::cout << " control points solved" << std::endl;
        //         surface.surface_visulization(surface, 100, SPs, SFs);

        //         igl::write_triangle_mesh(path + "ours_" + "p" + std::to_string(nbr) + "_refine_"
        //         +
        //                                      std::to_string(i) + "_m_" + std::to_string(method) +
        //                                      tail + ".obj",
        //                                  SPs, SFs);
        //         std::vector<std::string> titles = {{"time_knot", "time_solve", "precision", "nu",
        //                                             "nv", "cps", "time_itr", "max_energy"}};
        //         int cps = (surface.nu() + 1) * (surface.nv() + 1);
        //         precision = surface.max_interpolation_err(ver, param, surface);
        //         std::vector<double> data = {{time_knot, time_solve, precision,
        //         double(surface.nu()),
        //                                      double(surface.nv()), double(cps), timeitr,
        //                                      max_energy}};
        //         write_csv(path + "ours_" + "p" + std::to_string(nbr) + "_refine_" +
        //                       std::to_string(i) + "_m_" + std::to_string(method) + tail + ".csv",
        //                   titles, data);
        //     }
        //   return;
        //}
    }

    std::cout << "final U and V, " << surface.U.size() << " " << surface.V.size() << std::endl;
    print_vector(surface.U);
    print_vector(surface.V);
    precision = surface.max_interpolation_err(ver, param, surface);
    std::cout << "maximal interpolation error "
              << surface.max_interpolation_err(ver, param, surface) << std::endl;
    bool write_file = true;
    if (write_file)
    {
        Eigen::MatrixXi Pf;
        write_points(path + "pts" + std::to_string(nbr) + "_m_" + std::to_string(method) + ".obj",
                     ver);
        igl::write_triangle_mesh(path + "ours_" + "p" + std::to_string(nbr) + "_m_" +
                                     std::to_string(method) + tail + ".obj",
                                 SPs, SFs);

        std::vector<std::string> titles = {
            {"time_knot", "time_solve", "precision", "nu", "nv", "cps", "time"}};
        int cps = (surface.nu() + 1) * (surface.nv() + 1);
        std::vector<double> data = {{time_knot, time_solve, precision, double(surface.nu()),
                                     double(surface.nv()), double(cps), time_knot + time_solve}};
        write_csv(path + "ours_" + "p" + std::to_string(nbr) + "_m_" + std::to_string(method) +
                      tail + ".csv",
                  titles, data);
    }
    output_timing();
    if (1)
    { // write svg files
        write_svg_pts(path + "ours_" + "p" + std::to_string(nbr) + "_m_" + std::to_string(method) +
                          tail + "param.svg",
                      param);
        write_svg_knot_vectors(path + "ours_" + "p" + std::to_string(nbr) + "_m_" +
                                   std::to_string(method) + tail + "knots.svg",
                               surface.U, surface.V);
    }
    std::cout << "total time " << time_knot + time_solve << std::endl;
    polyscope::SurfaceMesh *psSurfaceMesh =
        polyscope::registerSurfaceMesh("Interpolated Surface" + std::to_string(model), SPs, SFs);
    polyscope::PointCloud *psPointCloud =
        polyscope::registerPointCloud("Old Model" + std::to_string(model), SPs);
}
} // namespace SIBSplines