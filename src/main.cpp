#include <array>
#include <fstream>
#include <iostream>
#include <list>
#include <string>

#include <igl/file_dialog_open.h>
#include <igl/readOBJ.h>
#include <polyscope/point_cloud.h>
#include <polyscope/surface_mesh.h>

#include "methods.hpp"
#include <Eigen/Core>
#include <TinyAD/Scalar.hh>
#include <sparse_interp/Types.hpp>

using namespace SIBSplines;
std::string root_path(SI_MESH_DIR);
std::string filename;
double sphereSize = 0.002;
bool enablePolyGUI = true;
double user_per = 0.5;
double user_delta = 0.4; // ours
int itSteps = 10;        // default
int modelType = 0;
int nbr_of_pts = 100;
std::vector<std::pair<Eigen::MatrixXd, Eigen::MatrixXi>> meshList;

using namespace SIBSplines;
std::string example_root_path(SI_MESH_DIR);

void interpCallback()
{
    if (enablePolyGUI)
    {
        // polyscope::buildPickGui();
        polyscope::buildPolyscopeGui();
        polyscope::buildStructureGui();
    }

    //  ImGui::SetNextWindowPos(ImVec2(600, 20), 0);
    // ImGui::SetNextWindowSize(ImVec2(300, 190), 0);
    ImGui::Begin(
        "Surface Interpolation Tools", nullptr,
        /*ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoResize |*/ ImGuiWindowFlags_NoCollapse |
            ImGuiWindowFlags_AlwaysAutoResize);

    ImGui::Text("Mesh load options.");
    ImGui::Checkbox("Toggle Polycsope GUI", &enablePolyGUI);
    if (ImGui::Button("Import Point Cloud"))
    {
        filename = igl::file_dialog_open();
        if (!filename.empty())
        {
            Eigen::MatrixXd V;
            Eigen::MatrixXi F;
            if (igl::readOBJ(filename, V, F))
            {
                // TODO make a list of structures
                // -> array
                // for (int i = 0; i < )
                meshList.emplace_back(V, F);
                polyscope::PointCloud *psCloud = polyscope::registerPointCloud("mesh points", V);
                psCloud->setPointRadius(sphereSize);
                psCloud->setPointRenderMode(polyscope::PointRenderMode::Sphere);
            }
            else
            {
                std::cerr << "Failed to read OBJ." << std::endl;
            }
        }
    }
    ImGui::SetNextItemOpen(true, ImGuiCond_Once);
    if (ImGui::TreeNode("Interpolation"))
    {
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
        if (ImGui::Button("Interpolate Mesh"))
        {
            if (!filename.empty())
            {
                mesh_interpolation(filename, user_delta, user_per, itSteps);
            }
            else
            {
                std::cerr << "No mesh loaded";
            }
        }
        if (ImGui::TreeNode("Predefined Functions"))
        {
            ImGui::SliderInt("Function Model", &modelType, 0, 5);
            ImGui::SliderInt("Number of points", &nbr_of_pts, 10, 300);
            if (ImGui::Button("Compute Function Interpolation"))
            {
                run_old_algorithm(modelType, nbr_of_pts, user_delta, SI_MESH_DIR, "", user_per,
                                  true, user_delta, itSteps);
                Eigen::MatrixXd verticies;
                Eigen::MatrixXi faces;
                std::string prefix = "ours_p" + std::to_string(nbr_of_pts) + "_m_";
                std::string model_filename =
                    SI_MESH_DIR + prefix + std::to_string(modelType) + ".obj";
                igl::readOBJ(model_filename, verticies, faces);
                polyscope::SurfaceMesh *psSurfaceMesh = polyscope::registerSurfaceMesh(
                    "Interpolated Surface" + std::to_string(modelType), verticies, faces);
                std::string prefix_pts =
                    "pts" + std::to_string(nbr_of_pts) + "_m_" + std::to_string(modelType) + ".obj";
                std::string model_points = SI_MESH_DIR + prefix_pts;
                Eigen::MatrixXd verticies_pts;
                Eigen::MatrixXi faces_pts;
                igl::readOBJ(model_points, verticies_pts, faces_pts);
                polyscope::PointCloud *psPointCloud = polyscope::registerPointCloud(
                    "Model" + std::to_string(modelType), verticies_pts);
            }
            ImGui::TreePop();
        }
        ImGui::TreePop();
    }

    ImGui::End();
}

int main(int argc, char *argv[])
{
    polyscope::options::programName = "Surface Interpolation GUI";
    polyscope::options::autocenterStructures = true;
    polyscope::options::autoscaleStructures = true;
    polyscope::options::buildGui = false;
    polyscope::init();

    polyscope::state::userCallback = interpCallback;

    polyscope::show();
    return 0;
}