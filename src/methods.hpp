#include <sparse_interp/Types.hpp>

namespace SIBSplines
{
void surface_init(const std::string meshfile, const std::string tail, double delta,
                  const double per, const int target_steps, const double w_fair,
                  const bool meshInterpolation, const int nbr_pts, const int model,
                  Bsurface &surface, Eigen::MatrixXd &param,
                  std::vector<std::array<int, 2>> paraInInterval_orig,
                  const Eigen::MatrixXd &ver_orig);
void mesh_optimization(Bsurface &surface, PartialBasis &basis, const double w_fair,
                       const int itSteps, const std::vector<std::array<int, 2>> &paraInInterval,
                       Eigen::MatrixXd &param, const int method, const Eigen::MatrixXd &ver);
void mesh_visualization(const Eigen::MatrixXd &param, const Eigen::MatrixXd &ver, Bsurface &surface,
                        std::string path, const int method);
void old(const int model, const int nbr_pts, double &per_ours, const std::string path,
         const std::string tail, const double per, const bool enable_local_energy);
} // namespace SIBSplines