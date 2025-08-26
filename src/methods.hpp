#include <sparse_interp/Types.hpp>

namespace SIBSplines
{
void mesh_interpolation(std::string meshfile, double delta, double per, int target_steps);
void function_interpolation(const int model, const int nbr_pts, const std::string path,
                            const std::string tail, const double per, double &delta,
                            const int target_steps, const double w_fair);
void old(const int model, const int nbr_pts, double &per_ours, const std::string path,
         const std::string tail, const double per, const bool enable_local_energy);
} // namespace SIBSplines