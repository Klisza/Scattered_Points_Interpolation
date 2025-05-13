#include <polyscope/polyscope.h>
#include <polyscope/surface_mesh.h>
#include <igl/readOBJ.h>
#include <igl/file_dialog_open.h>
#include <Eigen/Core>
#include <string>
#include <iostream>

int main(int argc, char* argv[]) {
  // 1) Initialize Polyscope
  polyscope::init();

  // 2) Install a per-frame ImGui callback using Polyscope's state
  polyscope::state::userCallback = []() {
    // Draw a button in the upper-left corner
    if (ImGui::Button("Import Mesh")) {

      // 3) Ask the user for an OBJ file
      std::string filename;
      // filter: “OBJ files” \0*.obj\0
      filename = igl::file_dialog_open();

      // Check if a file was selected
      if (!filename.empty()) {
        // 4) Read vertices & faces
        Eigen::MatrixXd V;
        Eigen::MatrixXi F;
        if (igl::readOBJ(filename, V, F)) {

          // 5) Register with Polyscope
          polyscope::registerSurfaceMesh("imported mesh", V, F);
        } else {
          // Log an error message to the console
          std::cerr << "Failed to read OBJ." << std::endl;
        }
      }
    }
  };

  // 6) Show the window and enter the Polyscope loop
  polyscope::show();

  return 0;
}