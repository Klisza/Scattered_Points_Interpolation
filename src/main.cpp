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

using namespace SIBSplines;
std::string example_root_path(SI_MESH_DIR);

void interpCallback()
{
    static std::string root_path(SI_MESH_DIR);
    static std::string filename;
    static double sphereSize = 0.002;
    static bool enablePolyGUI = true;
    static double user_per = 0.5;
    static double user_delta = 0.4; // ours
    static int itSteps = 50;        // default
    static int modelType = 0;
    static int nbr_of_pts = 100;
    static double w_fair = 10e-6;
    static std::vector<std::pair<Eigen::MatrixXd, Eigen::MatrixXi>> meshList;
    static Bsurface surface;
    static std::unique_ptr<PartialBasis> basis = nullptr;
    static bool surfaceInit = false;
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
                surfaceInit = false;
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
        ImGui::SliderInt("Iteration Steps", &itSteps, 0, 2000);
        ImGui::InputDouble("Weight", &w_fair);
        // Mesh calculation
        if (ImGui::Button("Interpolate Mesh"))
        {
            if (!filename.empty())
            {
                Eigen::MatrixXd param;
                std::vector<std::array<int, 2>> paraInInterval;
                Eigen::MatrixXd V;
                if (!surfaceInit)
                {
                    surface_init(filename, "", user_delta, user_per, itSteps, w_fair, true, 100, -1,
                                 surface, param, paraInInterval, V);
                    basis = std::make_unique<PartialBasis>(surface); // Initialize here
                    surfaceInit = true;
                }
                mesh_optimization(surface, *basis, w_fair, itSteps, paraInInterval, param, -1, V);
            }
        }
        if (ImGui::TreeNode("Predefined Functions"))
        {

            if (ImGui::SliderInt("Function Model", &modelType, 0, 5))
            {
                surfaceInit = false;
            }
            if (ImGui::SliderInt("Number of points", &nbr_of_pts, 10, 300))
            {
                surfaceInit = false;
            }

            if (ImGui::Button("Compute Function Interpolation"))
            {
                static Eigen::MatrixXd param;
                static std::vector<std::array<int, 2>> paraInInterval;
                static Eigen::MatrixXd V;
                if (!surfaceInit)
                {
                    surface_init("", "", user_delta, user_per, itSteps, w_fair, false, nbr_of_pts,
                                 modelType, surface, param, paraInInterval, V);
                    basis = std::make_unique<PartialBasis>(surface);
                    surfaceInit = true;
                    // -------------------------------------------
                    // Read the interpolated points into polyscope
                    std::string prefix_pts = "pts" + std::to_string(nbr_of_pts) + "_m_" +
                                             std::to_string(modelType) + "_orig.obj";
                    std::string model_points = SI_MESH_DIR + prefix_pts;
                    Eigen::MatrixXd verticies_pts;
                    Eigen::MatrixXi faces_pts;
                    igl::readOBJ(model_points, verticies_pts, faces_pts);
                    polyscope::PointCloud *psPointCloud = polyscope::registerPointCloud(
                        "Model Orig" + std::to_string(modelType), verticies_pts);
                    psPointCloud->setPointRadius(sphereSize);
                    psPointCloud->setPointRenderMode(polyscope::PointRenderMode::Sphere);
                    psPointCloud->resetTransform();
                    // -------------------------------------------
                    // Read the mesh into polyscope
                    Eigen::MatrixXd verticies;
                    Eigen::MatrixXi faces;
                    std::string prefix = "ours_p" + std::to_string(nbr_of_pts) + "_m_";
                    std::string model_filename =
                        SI_MESH_DIR + prefix + std::to_string(modelType) + "_orig.obj";
                    igl::readOBJ(model_filename, verticies, faces);
                    polyscope::SurfaceMesh *psSurfaceMesh = polyscope::registerSurfaceMesh(
                        "Interpolated Surface Orig" + std::to_string(modelType), verticies, faces);
                    psSurfaceMesh->resetTransform();
                    // -------------------------------------------
                }
                mesh_optimization(surface, *basis, w_fair, itSteps, paraInInterval, param,
                                  modelType, V);
                // -------------------------------------------
                // Read the interpolated points into polyscope
                std::string prefix_pts =
                    "pts" + std::to_string(nbr_of_pts) + "_m_" + std::to_string(modelType) + ".obj";
                std::string model_points = SI_MESH_DIR + prefix_pts;
                Eigen::MatrixXd verticies_pts;
                Eigen::MatrixXi faces_pts;
                igl::readOBJ(model_points, verticies_pts, faces_pts);
                polyscope::PointCloud *psPointCloud = polyscope::registerPointCloud(
                    "Model" + std::to_string(modelType), verticies_pts);
                psPointCloud->setPointRadius(sphereSize);
                psPointCloud->setPointRenderMode(polyscope::PointRenderMode::Sphere);
                psPointCloud->resetTransform();
                // -------------------------------------------
                // Read the mesh into polyscope
                Eigen::MatrixXd verticies;
                Eigen::MatrixXi faces;
                std::string prefix = "ours_p" + std::to_string(nbr_of_pts) + "_m_";
                std::string model_filename =
                    SI_MESH_DIR + prefix + std::to_string(modelType) + ".obj";
                igl::readOBJ(model_filename, verticies, faces);
                polyscope::SurfaceMesh *psSurfaceMesh = polyscope::registerSurfaceMesh(
                    "Interpolated Surface" + std::to_string(modelType), verticies, faces);
                psSurfaceMesh->resetTransform();
                // -------------------------------------------
            }
            if (ImGui::Button("Old Compute Function Interpolation"))
            {
                old(modelType, nbr_of_pts, user_delta, SI_MESH_DIR, "", user_per, false);
                // -------------------------------------------
                // Read the mesh into polyscope
                Eigen::MatrixXd verticies;
                Eigen::MatrixXi faces;
                std::string prefix = "ours_p" + std::to_string(nbr_of_pts) + "_m_";
                std::string model_filename =
                    SI_MESH_DIR + prefix + std::to_string(modelType) + ".obj";
                igl::readOBJ(model_filename, verticies, faces);
                polyscope::SurfaceMesh *psSurfaceMesh = polyscope::registerSurfaceMesh(
                    "Interpolated Surface" + std::to_string(modelType), verticies, faces);
                psSurfaceMesh->resetTransform();
                // -------------------------------------------
                // Read the interpolated points into polyscope
                std::string prefix_pts =
                    "pts" + std::to_string(nbr_of_pts) + "_m_" + std::to_string(modelType) + ".obj";
                std::string model_points = SI_MESH_DIR + prefix_pts;
                Eigen::MatrixXd verticies_pts;
                Eigen::MatrixXi faces_pts;
                igl::readOBJ(model_points, verticies_pts, faces_pts);
                polyscope::PointCloud *psPointCloud = polyscope::registerPointCloud(
                    "Model" + std::to_string(modelType), verticies_pts);
                // -------------------------------------------
                psPointCloud->resetTransform();
            }
            ImGui::TreePop();
        }
        ImGui::TreePop();
    }

    ImGui::End();
}

int main(int argc, char *argv[])
{
    std::cerr << "Starting main..." << std::endl;
    polyscope::options::programName = "Surface Interpolation GUI";
    polyscope::options::autocenterStructures = true;
    polyscope::options::autoscaleStructures = true;
    polyscope::options::buildGui = false;
    polyscope::init();
    std::cout << "Maybe crashes here" << std::endl;
    polyscope::state::userCallback = interpCallback;

    polyscope::show();
    return 0;
}