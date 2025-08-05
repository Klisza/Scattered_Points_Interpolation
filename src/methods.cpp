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

#include <TinyAD/Scalar.hh>
#include <TinyAD/ScalarFunction.hh>
#include <TinyAD/Utils/LineSearch.hh>
#include <TinyAD/Utils/NewtonDecrement.hh>
#include <TinyAD/Utils/NewtonDirection.hh>

#define PARAMETER true
#define CP false

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

void write_csv(const std::string &file, const std::vector<std::string> titles,
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
// the output type, the knot vector type, and the value type.
// These are the operations to compute spline-related functions or values.
template <typename Tp, typename knotT, typename valueT> class splineBasis
{
  public:
    splineBasis() {};
    ~splineBasis() {};
    // the basis functions in each none-zero interval. in each interval there are degree + 1
    // none-zero basis functions, Each basis function is a polynomial.
    // std::vector<std::vector<std::vector<T>>> Nips;
    // the variable u is in [U_{uId}, U_{uId+1}). This function returns double(1) if uID == i,
    // otherwise double(0)
    double Ni0_func(const int i, const int uId);
    std::vector<knotT> Nip_func(const int i, const int p, const int uId,
                                const std::vector<knotT> &U, const int original_p);
    // compute the basis w.r.t. the interval, and the values of them
    std::vector<Tp> computeBasisFunctionValues(const valueT &value, const int uId, const int p,
                                               const std::vector<knotT> &U);
    // the order is the order of derivative.bvalues is the basis functions N(value) , dvalues is the
    // derivatived bvalues of order
    void computeBasisFunctionDerivativeValues(const valueT &value, const int uId, const int p,
                                              const int order, const std::vector<knotT> &U,
                                              std::vector<Tp> &bvalues, std::vector<Tp> &dvalues);
    void computeBasisFunctionDerivativeAndValues(
        const valueT &value, const int uId, const int p, const std::vector<knotT> &U,
        std::vector<std::vector<knotT>> &basisfuncs, std::vector<std::vector<knotT>> &d1funcs,
        std::vector<std::vector<knotT>> &d2funcs, std::vector<Tp> &bvalues,
        std::vector<Tp> &d1values, std::vector<Tp> &d2values);
    ply_operations<Tp, knotT, valueT> PO;
    std::vector<std::vector<std::vector<knotT>>> basisFunctions;
};
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

Eigen::VectorXd list_to_vec(std::vector<double> &vars)
{
    Eigen::VectorXd result;
    return result;
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
int return_closest_knot_index_to_param(std::vector<double> UV, double param)
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

// Gives you the varaible depending on U or V direction, variable type (control point or parameter)
int variableMap(bool Udirection, bool ParOrCp, int pos, int id)
{
    int result;
    if (ParOrCp)
    {
        result = parameterLocation(Udirection, pos, id);
    }
    else
    {
        result = controlPointLocation(Udirection, pos, id);
    }
    return result;
}

int parameterLocation(bool Udirection, int pos, int id)
{
    if (Udirection)
    {
    }
    else
    {
    }
}
int controlPointLocation(bool Udirection, int pos, int id) {}

void mesh_interpolation(std::string meshfile, double delta, double per, int target_steps)
{
    double precision = 0;
    Eigen::MatrixXd ver;
    Eigen::MatrixXi F;
    Eigen::MatrixXd param;
    std::cout << "reading mesh model: " << meshfile << std::endl;

    // mesh parametrization, and print out the parametrization result as a obj mesh.
    mesh_parameterization(meshfile, ver, param, F);
    rescale_param(param);

    // construct the surface object
    Bsurface surface;

    // set up the initial parameters.
    int param_nbr = param.rows();           // the number of data points
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
    surface.solve_control_points_for_fairing_surface(surface, param, ver, basis);
    std::cout << "Control points initialized" << std::endl;

    // Init parameter intervals for reparameterization
    std::vector<std::array<int, 2>> paraInInterval(param_nbr, {0, 0});
    // Perform binary search and return the element
    for (int i = 0; i < paraInInterval.size() - 1; i++)
    {
        paraInInterval[i][0] = return_closest_knot_index_to_param(surface.U, param(i, 0));
        paraInInterval[i][1] = return_closest_knot_index_to_param(surface.V, param(i, 1));
    }

    //  Calculate the number of variables we solve for.
    //  Number of control points
    int cpSize =
        (surface.U.size() - 1 - surface.degree1) * (surface.V.size() - 1 - surface.degree2);
    // Number of variables = 2 * parameters (u,v) + 3 * control points (x,y,z)
    int varSize = 2 * param_nbr + 3 * cpSize;
    // Init globVars vector
    int cpCols = surface.control_points.size();
    int cpRows = surface.control_points[0].size();
    // Add control points to globVars
    // Put the grid structure into a vector (i,j) directions same as (u,v) -> variableMap
    for (int k = 0; k < 3; k++)
    {
        for (int i = 0; i < cpRows; i++)
        {
            for (int j = 0; j < cpCols; j++)
            {
                surface.globVars[i * j + k * cpSize] = surface.control_points[i][j](k);
            }
        }
    }
    // Add parameters to globVars
    for (int k = 0; k < 2; k++)
    {
        for (int i = 0; i < param.size(); i++)
        {
            surface.globVars[i + k * param.size() + 3 * cpSize] = param(i, k);
        }
    }

    // Solve for the variables (parameters and control points) using the knot vectors.
    auto func = TinyAD::scalar_function<1>(TinyAD::range(varSize));
    // Fitting energy
    // (d+1)*(d+1) cps * 3 (x,y,z) + 2 parameters (u,v).
    // For degree 3 = 50 variables
    func.add_elements<50>(
        TinyAD::range(param_nbr),
        [&](auto &element) -> TINYAD_SCALAR_TYPE(element)
        {
            using T = TINYAD_SCALAR_TYPE(element);
            Eigen::Index dataID = element.handle;

            T parameterU = element.variables(variableMap())(0, 0);
            T parameterV = element.variables(variableMap())(0, 0);
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
            T p00 = 0, p01 = 0, p02 = 0;
            for (int i = 0; i < surface.degree1 + 1; i++)
            {
                for (int j = 0; j < surface.degree2 + 1; j++)
                {
                    // for each interval [U[uItv], U[uItv + 1]), the associated control points are
                    // P[uItv-degree],...,P[uItv]
                    int lc0 =
                        ; // getTheLocationOfThe(uItv-degree1+i, vItv-degree2+j)-th ControlPoint_x;
                    int lc1 =
                        ; // getTheLocationOfThe(uItv-degree1+i, vItv-degree2+j)-th ControlPoint_y;
                    int lc2 =
                        ; // getTheLocationOfThe(uItv-degree1+i, vItv-degree2+j)-th ControlPoint_z;
                    T pt0 = 0, pt1 = 0, pt2 = 0;
                    // get the control point
                    pt0 = element.variables(lc0)(0, 0), pt1 = element.variables(lc1)(0, 0),
                    pt2 = element.variables(lc2)(0, 0);

                    p00 += pt0 * basisU[i] * basisV[j];
                    p01 += pt1 * basisU[i] * basisV[j];
                    p02 += pt2 * basisU[i] * basisV[j];
                }
            }

            // std::cout<<"out add elem\n";
            return (p00 - ver(dataID, 0)) * (p00 - ver(dataID, 0)) +
                   (p01 - ver(dataID, 1)) * (p01 - ver(dataID, 1)) +
                   (p02 - ver(dataID, 2)) * (p02 - ver(dataID, 2));
        });

    // Fairing energy manual

    // solve for n iterations.
    Eigen::VectorXd x = list_to_vec(surface.globVars); // put your variables into a list x

    TinyAD::LinearSolver solver;
    double convergence_eps = 1e-12; // change it into 1e-6 if you want.
    std::cout << "check 4\n";
    for (int i = 0; i < target_steps; ++i)
    {
        auto [f, g, H_proj] =
            func.eval_with_hessian_proj(x); // compute Hessian and gradient of the fitting error
        TINYAD_DEBUG_OUT("Energy in iteration " << i << ": " << f);
        double ferror;
        s.evaluateFittingError(ferror, false); // compute the fitting error and print it out.
        std::cout << "sum of squared error " << ferror << "\n";
        // BW: assemble the fitting Hessian, fitting gradient and fairness Hessian, fairness
        // gradient together, blend them with their corresponding weights.
        Eigen::VectorXd d = TinyAD::newton_direction(g, H_proj, solver);
        // the following variable "sparsity_pattern_dirty" can use the default value "false", since
        // our sparse matrices are all in the same structure. Thus we comment the following command
        // out.
        // solver.sparsity_pattern_dirty =
        //     true; // this is crutial, since the patterns are always changing in each iteration!

        if (TinyAD::newton_decrement(d, g) <
            convergence_eps) // if the direction is too far from the gradient direction, break.
                             // normally this value is set as 1e-6
            break;
        // modify the direction vector. 2 options: 1. rescale (backtrace) the whole vector d to keep
        // the parameters staying in their regions, or 2. only re-scale the parameters that may
        // exceed their regions.
        d = s.stepBackTracer(varLevel, d, parSeparation, knotSeparation);

        // BW: write your own line search code, since the line search in TinyAD will only consider
        // about your fitting energy. we need to implement one with considering both the energies.
        x = TinyAD::line_search(x, d, f, g, func);

        if ((x - list_to_vec(s.globVars)).norm() <
            convergence_eps) // if the step is too small, break
        {
            std::cout << "break because the line searched step is too small: "
                      << (x - list_to_vec(s.globVars)).norm() << "\n";
            break;
        }
        std::cout << "the dx, " << d.norm() << ", the backtraced dx "
                  << (x - list_to_vec(s.globVars)).norm() << "\n";
    }
    /* ///////////////////////
        Data Visualization
    //////////////////////// */
    if (1)
    {
        Eigen::MatrixXd SPs;
        Eigen::MatrixXi SFs;
        int visual_nbr =
            200; // the discretization scale for the output surface. The mesh will be 200x200
        surface.surface_visulization(surface, visual_nbr, SPs, SFs);

        precision = surface.max_interpolation_err(ver, param, surface);
        std::cout << "maximal interpolation error "
                  << surface.max_interpolation_err(ver, param, surface) << std::endl;
        write_points(meshfile + "pts" + std::to_string(param_nbr) + ".obj", ver);
        write_triangle_mesh(meshfile + "_intp_" + "p" + std::to_string(param_nbr) + ".obj", SPs,
                            SFs);
        Eigen::MatrixXd verticies;
        Eigen::MatrixXi faces;
        igl::readOBJ(meshfile + "_intp_" + "p" + std::to_string(param_nbr) + ".obj", verticies,
                     faces);
        polyscope::SurfaceMesh *psSurfaceMesh =
            polyscope::registerSurfaceMesh("Interpolated Surface", verticies, faces);
    }
}

// Optimized for all the variables using TinyAD as a autodifferenciation.
void run_old_algorithm(const int model, const int nbr_pts, double &per_ours, const std::string path,
                       const std::string tail, const double per, const bool enable_local_energy)
{
    /*// Timer variables
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
    SIBSplines::examples::get_model_sample_points(nbr, ver, F, param, method, corners, path); //
    Inits the mesh vars
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

    Bsurface< surface;
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
        if (enable_local_energy)
        {
            double timeitr = 0;
            for (int i = 0; i < 50; i++)
            {
                timer.start();
                Eigen::MatrixXd energy, euu, evv, euv;
                energy = surface.surface_energy_calculation(surface, basis, 1, euu, evv, euv);
                bool uorv;
                int which;
                double max_energy;
                surface.detect_max_energy_interval(surface, energy, euu, evv, uorv, which,
    max_energy); std::vector<double> Unew = surface.U; std::vector<double> Vnew = surface.V; if
    (!uorv) { // u get updated double value = (Unew[which] + Unew[which + 1]) / 2; surface.U =
    knot_vector_insert_one_value(Unew, value);
                }
                else
                {
                    double value = (Vnew[which] + Vnew[which + 1]) / 2;
                    surface.V = knot_vector_insert_one_value(Vnew, value);
                }
                std::cout << "knot vector get inserted" << std::endl;
                basis.clear();
                basis.init(surface);
                surface.solve_control_points_for_fairing_surface(surface, param, ver, basis);
                timer.stop();
                timeitr += timer.getElapsedTimeInSec();
                std::cout << " control points solved" << std::endl;
                surface.surface_visulization(surface, 100, SPs, SFs);

                igl::write_triangle_mesh(path + "ours_" + "p" + std::to_string(nbr) + "_refine_" +
    std::to_string(i) + "_m_" + std::to_string(method) + tail + ".obj", SPs, SFs);
                std::vector<std::string> titles = {{"time_knot", "time_solve", "precision", "nu",
    "nv", "cps", "time_itr", "max_energy"}}; int cps = (surface.nu() + 1) * (surface.nv() + 1);
                precision = surface.max_interpolation_err(ver, param, surface);
                std::vector<double> data = {{time_knot, time_solve, precision,
                                             double(surface.nu()), double(surface.nv()),
    double(cps), timeitr, max_energy}}; write_csv(path + "ours_" + "p" + std::to_string(nbr) +
    "_refine_" + std::to_string(i) + "_m_" + std::to_string(method) + tail + ".csv", titles, data);
            }
            return;
        }
    }

    std::cout << "final U and V, " << surface.U.size() << " " << surface.V.size() << std::endl;
    print_vector(surface.U);
    print_vector(surface.V);
    precision = surface.max_interpolation_err(ver, param, surface);
    std::cout << "maximal interpolation error " << surface.max_interpolation_err(ver, param,
    surface) << std::endl; bool write_file = true; if (write_file)
    {
        Eigen::MatrixXi Pf;
        write_points(path + "pts" + std::to_string(nbr) + "_m_" + std::to_string(method) + ".obj",
    ver); igl::write_triangle_mesh(path + "ours_" + "p" + std::to_string(nbr) + "_m_" +
    std::to_string(method) + tail + ".obj", SPs, SFs);

        std::vector<std::string> titles = {{"time_knot", "time_solve", "precision", "nu", "nv",
    "cps", "time"}}; int cps = (surface.nu() + 1) * (surface.nv() + 1); std::vector<double> data =
    {{time_knot, time_solve, precision, double(surface.nu()), double(surface.nv()), double(cps),
                                     time_knot + time_solve}};
        write_csv(path + "ours_" + "p" + std::to_string(nbr) + "_m_" + std::to_string(method) + tail
    + ".csv", titles, data);
    }
    output_timing();
    if (1)
    { // write svg files
        write_svg_pts(path + "ours_" + "p" + std::to_string(nbr) + "_m_" + std::to_string(method) +
    tail + "param.svg", param); write_svg_knot_vectors(path + "ours_" + "p" + std::to_string(nbr) +
    "_m_" + std::to_string(method) + tail + "knots.svg", surface.U, surface.V);
    }
    std::cout << "total time " << time_knot + time_solve << std::endl;*/
}
} // namespace SIBSplines