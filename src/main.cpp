#include <iostream>
#include <string>
#include <fstream>
#include <array>

#include <polyscope/point_cloud.h>
#include <polyscope/surface_mesh.h>
#include <igl/file_dialog_open.h>
#include <igl/readOBJ.h>

#include <Eigen/Core>
#include <sparse_interp/basis.h>
#include <sparse_interp/curve.h>
#include <sparse_interp/mesh_processing.h>
#include <sparse_interp/energy.h>
#include <test.h>
#include <sparse_interp/Types.hpp>

#include <igl/opengl/glfw/Viewer.h>
#include <igl/boundary_loop.h>
#include <igl/harmonic.h>
#include <igl/write_triangle_mesh.h>


using namespace SIBSplines;
std::string root_path(SI_MESH_DIR);
std::string filename;
double sphereSize = 0.002;
bool enablePolyGUI = true;
double par = 0.9; // ours
double delta = 0.4;

using namespace SIBSplines;
std::string example_root_path(SI_MESH_DIR);
void mesh_interpolation(std::string meshfile, double per, double par, int target_steps)
{
	double precision = 0;
	Eigen::MatrixXd ver;
	Eigen::MatrixXi F;
	Eigen::MatrixXd param, paramout;
	// std::string modelname = "tiger.obj";
	// std::string meshfile = example_root_path + modelname;
	std::cout << "reading mesh model: " << meshfile << std::endl;
	// mesh parametrization, and print out the parametrization result as a obj mesh.
	mesh_parameterization(meshfile, ver, param, F);
	paramout.resize(param.rows(), 3);
	Eigen::VectorXd param_zero = Eigen::VectorXd::Zero(param.rows());
	paramout << param, param_zero;
  write_triangle_mesh(meshfile + "_param", paramout, F);
	// write_triangle_mesh(example_root_path + "param_" + modelname, paramout, F);

	// construct the surface object
	Bsurface surface;
	// set up the initial parameters.
	int nbr = param.rows();					// the number of data points
	surface.degree1 = 3;					// degree of u direction
	surface.degree2 = 3;					// degree of v direction
	surface.U = { {0, 0, 0, 0, 1, 1, 1, 1} }; // the initial U knot vector
	surface.V = surface.U;					// the initial V knot vector
	// int target_steps = 10;					// the number of iterations for constructing lists $L$.
	bool enable_max_fix_nbr = true;			// progressively update the knot vectors to make the two knot vectors balanced in length.
	// double delta = 0.4;						// the parameter to improve the solving stability
	// double per = 0.5;						// the parameter inherited from [Wen-Ke Wang et al, 2008, CAD]
	// generate knot vectors to make sure the data points can be interpolated
	surface.generate_interpolation_knot_vectors(surface.degree1, surface.degree2, surface.U, surface.V, param, delta, per, target_steps, enable_max_fix_nbr);
	std::cout << "knot vectors generated" << std::endl;

	Eigen::MatrixXd SPs;
	Eigen::MatrixXi SFs;
	int visual_nbr = 200; // the discretization scale for the output surface. The mesh will be 200x200

	// basis contains all the basis functions and their 1 and 2 order diffenrential form.
	PartialBasis basis(surface);

	// 	solve the control points to obtain the surface.
	surface.solve_control_points_for_fairing_surface(surface, param, ver, basis);
	std::cout << "surface solved" << std::endl;

	// convert B-spline surface into a triangle mesh
	surface.surface_visulization(surface, visual_nbr, SPs, SFs);

	precision = surface.max_interpolation_err(ver, param, surface);
	std::cout << "maximal interpolation error " << surface.max_interpolation_err(ver, param, surface) << std::endl;

	write_points(example_root_path + "pts" + std::to_string(nbr) + "_m_" + modelname, ver);
	write_triangle_mesh(example_root_path + "intp_" + "p" + std::to_string(nbr) + "_m_" + modelname, SPs, SFs);
}

void interpCallback() {
    if (enablePolyGUI){
      polyscope::buildPickGui();
      polyscope::buildStructureGui();
    }

    ImGui::SetNextWindowPos(ImVec2(10, 20), 0);
    // ImGui::SetNextWindowSize(ImVec2(300, 190), 0);
    ImGui::Begin("Surface Interpolation Tools",
             nullptr,
             ImGuiWindowFlags_NoMove /*| ImGuiWindowFlags_NoResize*/ | ImGuiWindowFlags_NoCollapse /*| ImGuiWindowFlags_AlwaysAutoResize*/);
  
    ImGui::Text("Mesh load options.");
    if (ImGui::Button("Import Point Cloud")) {
      filename = igl::file_dialog_open();
      if (!filename.empty()) {
        Eigen::MatrixXd V;
        Eigen::MatrixXi F;
        if (igl::readOBJ(filename, V, F)) {
          // TODO make a list of structures
          // -> array
          polyscope::PointCloud* psCloud = polyscope::registerPointCloud("mesh points", V);
          psCloud->setPointRadius(sphereSize);
          psCloud->setPointRenderMode(polyscope::PointRenderMode::Sphere);
        } else {
          std::cerr << "Failed to read OBJ." << std::endl;
        }
      }
    }
    if (ImGui::TreeNode("Interpolation")) {
      ImGui::TextWrapped("Set parameters for interpolation.");
      float f = static_cast<float>(par);
      ImGui::SliderFloat("par", &f, 0, 1);
      par = static_cast<double>(f);
      if (ImGui::Button("Interpolate Mesh")) {
        if (!filename.empty()) {
          const std::string model_filepath = filename;
          const std::string outpath = filename + "_interp";
          std::string tail = "";
          mesh_interpolation(filename, per, par);
        } else {
          std::cerr << "No mesh loaded";
        }
      }
    }
    
    ImGui::End();

    
  }

int main(int argc, char* argv[]) {
  polyscope::options::programName = "Surface Interpolation GUI";
  polyscope::options::autocenterStructures = true;
  polyscope::options::autoscaleStructures = true;
  polyscope::options::buildGui = false;
  polyscope::init();

  //run_ours();
  // Run user callback to get the mesh
  // TODO: Rewrite this to get the point cloud from the user instead of the triangle mesh
  polyscope::state::userCallback = interpCallback;

  polyscope::show();
  return 0;
}