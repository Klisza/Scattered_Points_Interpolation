#include <iostream>
#include <string>

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

using namespace SIBSplines;
std::string filename;

void interpCallback() {
    if (ImGui::BeginTabBar("Settings")) {

    }
    // Load and display mesh / point cloud
    if (ImGui::Button("Import Point Cloud")) {
      filename = igl::file_dialog_open();
      if (!filename.empty()) {
        Eigen::MatrixXd V;
        Eigen::MatrixXi F;
        if (igl::readOBJ(filename, V, F)) {
          polyscope::PointCloud* psCloud = polyscope::registerPointCloud("mesh points", V);
          psCloud->setPointRadius(0.002);
          psCloud->setPointRenderMode(polyscope::PointRenderMode::Sphere);
          // polyscope::registerSurfaceMesh("imported mesh", V, F);
        } else {
          std::cerr << "Failed to read OBJ." << std::endl;
        }
      }
    }
    if (ImGui::Button("Compute Mesh")) {

    }
}


int main(int argc, char* argv[]) {
  polyscope::options::programName = "Surface Interpolation GUI";
  polyscope::options::autocenterStructures = true;
  polyscope::options::autoscaleStructures = true;
  polyscope::init();

  //run_ours();
  // Run user callback to get the mesh
  // TODO: Rewrite this to get the point cloud from the user instead of the triangle mesh
  polyscope::state::userCallback = interpCallback;

  polyscope::show();
  return 0;
}