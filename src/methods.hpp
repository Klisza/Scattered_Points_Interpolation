#include <sparse_interp/Types.hpp>

namespace SIBSplines {
void mesh_interpolation(std::string meshfile, double delta, double per, int target_steps);
void run_old_algorithm(const int model, const int nbr_pts, double &per_ours, const std::string path, const std::string tail,
	const double per, const bool enable_local_energy);
}