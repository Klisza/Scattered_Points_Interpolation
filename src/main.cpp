#include <iostream>
#include <string>
#include <fstream>
#include <array>
#include <list>

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

#include <Eigen/Core>


using namespace SIBSplines;
std::string root_path(SI_MESH_DIR);
std::string filename;
double sphereSize = 0.002;
bool enablePolyGUI = true;
double user_per = 0.5; 
double user_delta = 0.9; // ours
int itSteps = 10; // default
int modelType = 0;
int nbr_of_pts = 100;
std::vector<std::pair<Eigen::MatrixXd, Eigen::MatrixXi>> meshList;

using namespace SIBSplines;
std::string example_root_path(SI_MESH_DIR);
void mesh_interpolation(std::string meshfile, double delta, double per, int target_steps)
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
	bool enable_max_fix_nbr = true;			// progressively update the knot vectors to make the two knot vectors balanced in length.
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
  write_points(meshfile + "pts" + std::to_string(nbr) + ".obj", ver);
  write_triangle_mesh(meshfile + "_intp_" + "p" + std::to_string(nbr) + ".obj", SPs, SFs);
  Eigen::MatrixXd verticies;
  Eigen::MatrixXi faces;
  igl::readOBJ(meshfile + "_intp_" + "p" + std::to_string(nbr) + ".obj", verticies, faces);
  polyscope::SurfaceMesh* psSurfaceMesh = polyscope::registerSurfaceMesh("Interpolated Surface", verticies, faces);
}

void interpCallback() {
    if (enablePolyGUI){
      //polyscope::buildPickGui();
      polyscope::buildPolyscopeGui();
      polyscope::buildStructureGui();
    }

    //  ImGui::SetNextWindowPos(ImVec2(600, 20), 0);
    // ImGui::SetNextWindowSize(ImVec2(300, 190), 0);
    ImGui::Begin("Surface Interpolation Tools",
             nullptr,
             /*ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoResize |*/ ImGuiWindowFlags_NoCollapse | ImGuiWindowFlags_AlwaysAutoResize);
  
    ImGui::Text("Mesh load options.");
    ImGui::Checkbox("Toggle Polycsope GUI", &enablePolyGUI);
    if (ImGui::Button("Import Point Cloud")) {
      filename = igl::file_dialog_open();
      if (!filename.empty()) {
        Eigen::MatrixXd V;
        Eigen::MatrixXi F;
        if (igl::readOBJ(filename, V, F)) {
          // TODO make a list of structures
          // -> array
          //for (int i = 0; i < )
          meshList.emplace_back(V, F);
          polyscope::PointCloud* psCloud = polyscope::registerPointCloud("mesh points", V);
          psCloud->setPointRadius(sphereSize);
          psCloud->setPointRenderMode(polyscope::PointRenderMode::Sphere);
        } else {
          std::cerr << "Failed to read OBJ." << std::endl;
        }
      }
    }
    ImGui::SetNextItemOpen(true, ImGuiCond_Once);
    if (ImGui::TreeNode("Interpolation")) {
      // Parameters
      ImGui::TextWrapped("Set parameters for interpolation.");
      float f1 = static_cast<float>(user_delta);
      ImGui::SliderFloat("Par", &f1, 0, 1);
      user_delta = static_cast<double>(f1);
      float f2 = static_cast<float>(user_per);
      ImGui::SliderFloat("Delta", &f2, 0, 1);
      user_per = static_cast<double>(f2);
      ImGui::SliderInt("Iteration Steps", &itSteps, 1, 20);
      // Mesh calculation
      if (ImGui::Button("Interpolate Mesh")) {
        if (!filename.empty()) {
          mesh_interpolation(filename, user_delta, user_per, itSteps);
        } else {
          std::cerr << "No mesh loaded";
        }
      }
      if (ImGui::TreeNode("Predefined Functions")) {
        ImGui::SliderInt("Function Model", &modelType, 0, 5);
        ImGui::SliderInt("Number of points", &nbr_of_pts, 10, 300);
        if (ImGui::Button("Compute Function Interpolation")) {
          run_ours(modelType, nbr_of_pts, user_delta, SI_MESH_DIR, "", user_per, false);
          Eigen::MatrixXd verticies;
          Eigen::MatrixXi faces;
          std::string prefix = "ours_p" + std::to_string(nbr_of_pts) + "_m_";
          std::string model_filename = SI_MESH_DIR + prefix + std::to_string(modelType) + ".obj";
          igl::readOBJ(model_filename, verticies, faces);
          polyscope::SurfaceMesh* psSurfaceMesh = polyscope::registerSurfaceMesh("Interpolated Surface" + std::to_string(modelType), verticies, faces);
          std::string prefix_pts = "pts" + std::to_string(nbr_of_pts) + "_m_" + std::to_string(modelType) + ".obj";
          std::string model_points = SI_MESH_DIR + prefix_pts;
          Eigen::MatrixXd verticies_pts;
          Eigen::MatrixXi faces_pts;
          igl::readOBJ(model_points, verticies_pts, faces_pts);
          polyscope::PointCloud* psPointCloud = polyscope::registerPointCloud("Model" + std::to_string(modelType), verticies_pts);
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

  polyscope::state::userCallback = interpCallback;

  polyscope::show();
  return 0;
}