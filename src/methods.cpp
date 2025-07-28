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

// using namespace SIBSplines;

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

void mesh_interpolation(std::string meshfile, double delta, double per, int target_steps)
{
    double precision = 0;
    Eigen::MatrixXd ver;
    Eigen::MatrixXi F;
    Eigen::MatrixXd param, paramout;
    std::cout << "reading mesh model: " << meshfile << std::endl;

    // mesh parametrization, and print out the parametrization result as a obj mesh.
    mesh_parameterization(meshfile, ver, param, F);
    paramout.resize(param.rows(), 3);
    Eigen::VectorXd param_zero = Eigen::VectorXd::Zero(param.rows());
    paramout << param, param_zero;
    write_triangle_mesh(meshfile + "_param", paramout, F);

    // construct the surface object
    Bsurface surface;

    // set up the initial parameters.
    int nbr = param.rows();                 // the number of data points
    surface.degree1 = 3;                    // degree of u direction
    surface.degree2 = 3;                    // degree of v direction
    surface.U = {{0, 0, 0, 0, 1, 1, 1, 1}}; // the initial U knot vector
    surface.V = surface.U;                  // the initial V knot vector
    bool enable_max_fix_nbr = true; // progressively update the knot vectors to make the two knot
                                    // vectors balanced in length.
    // generate knot vectors to make sure the data points can be interpolated
    // Alg 1 from the paper
    surface.generate_interpolation_knot_vectors(surface.degree1, surface.degree2, surface.U,
                                                surface.V, param, delta, per, target_steps,
                                                enable_max_fix_nbr);
    std::cout << "knot vectors generated" << std::endl;

    // Calculate the number of variables we solve for.
    int cpSize =
        (surface.U.size() - 1 - surface.degree1) * (surface.V.size() - 1 - surface.degree2);
    int varSize = 2 * nbr + 3 * cpSize;

    // Solve for the variables using the knot vectors.
    auto func = TinyAD::scalar_function<2>(TinyAD::range(varSize));
    // How many handles do I have?
    func.add_elements<2>(TinyAD::range(varSize),
                         [&](auto &element) -> TINYAD_SCALAR_TYPE(element)
                         {
                             using T = TINYAD_SCALAR_TYPE(element);
                             Eigen::Index dataID = element.handle;
                             SurfaceOpt<T, T> sOpt(surface);
                             PartialBasis<T, T> basis(surface);

                             sOpt.solve_control_points_for_fairing_surface(surface, param, ver,
                                                                           basis);
                         });
    std::cout << "surface solved" << std::endl;

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
        write_points(meshfile + "pts" + std::to_string(nbr) + ".obj", ver);
        write_triangle_mesh(meshfile + "_intp_" + "p" + std::to_string(nbr) + ".obj", SPs, SFs);
        Eigen::MatrixXd verticies;
        Eigen::MatrixXi faces;
        igl::readOBJ(meshfile + "_intp_" + "p" + std::to_string(nbr) + ".obj", verticies, faces);
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