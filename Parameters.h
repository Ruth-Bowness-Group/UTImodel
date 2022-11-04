#ifndef Parameters_h
#define Parameters_h

#include <vector>
#include <string>
#include <unordered_set>
#include <random>

class Parameters
{
public:
  Parameters(std::string file_path);
  // Simulation parameters
  double time_limit;
  double timestep;
  std::string output_folder;
  double output_interval;
  const char *last_name;
  const int number_vessels = 4;
  std::vector<std::pair<int, int> > fixed_vessels = {{175, 25}, {25, 175}, {175, 175}, {25, 25}};
  std::pair<int, int> fixed_bacteria_cluster = {100, 100};


  /**
 * Diffusion Parameters
 */
  // Environment parameters
  // TODO: Grid size declared here because it needs to be known at compile time in order to allocate memory for arrays


  bool fixed_vessel_placement;
  /**
 * Initial conditions
 */

  // Fixed bacteria. true => bacteria cluster at locations in fixed_bacteria_cluster
  // false => randomly placed
  bool fixed_bacteria_placement;
  // Numbers of initial bacteria
  int initial_number_bacteria;
  int initial_number_Rmacrophages;



  /**
 * Event parameters
 */
 double spatial_step;
 double chemokine_diffusion;
 double chemokine_rmac_source;
 double chemokine_decay;

  double rmacrophages_recruitment_prob_bac;
  double probability_hmacrophage_activation;
  double activatedresident_macrophage_movement_rate;
  double rmacrophages_recruitment_prob_no_bac;
  double hmacrophages_recruitment_prob_bac;
  double hmacrophages_recruitment_prob_no_bac;
  double resting_rmacrophage_max_life_span;
  double activated_Rmacrophage_lifespan;
  double activated_Hmacrophage_lifespan;
  double activated_hmacrophage_movement_rate;
  double resting_hmacrophage_movement_rate;
  double restingresident_macrophage_movement_rate;
  double probability_rmacrophage_activation;
  double probability_neutrophil_activation;
  double probability_neutrophil_destroy_dormant_bacteria;
  double probability_neutrophil_destroy_replicating_bacteria;
  double probability_hmac_destroy_dormant_bacteria;
  double probability_hmac_destroy_replicating_bacteria;
  double activated_neutrophil_movement_rate;
  double resting_neutrophil_movement_rate;
  double resting_neutrophil_max_life_span;
  double activated_neutrophil_lifespan;
  double neutrophil_recruitment_prob_no_bac;
  double neutrophil_recruitment_prob_bac;

  double activated_Mastcell_lifespan;
  int initial_number_mastcell;
  double resting_mastcell_max_life_span;
  double restingmastcell_movement_rate;
  double mastcells_recruitment_prob_bac;
  double mastcells_recruitment_prob_no_bac;
  double activatedmastcell_movement_rate;

  double Hmacrophage_chemotactic_movement_bias;
  double resting_rmacrophage_movement_bias;
  double activated_rmacrophage_movement_bias;
  double neutrophil_chemotactic_movement_bias;

  // Max and minimum hours for replicating bacteria replication
  double bacteria_rep_hours_min_early;
  double bacteria_rep_hours_max_early;
  double bacteria_rep_hours_min_middle;
  double bacteria_rep_hours_max_middle;
  double bacteria_rep_hours_min_late;
  double bacteria_rep_hours_max_late;

  int bacteria_shedding_hours_min;
  int bacteria_shedding_hours_max;
  // Maximum distance to look for space when bacteria replicate
  int bacteria_replication_neighbourhood_max_depth;

  // Time after which bacteria can change replication state
  double time_threshold_for_bacterium_state_change;
  double prob_bacteria_penetrate_wall;

};

const int grid_size =200;
// Random number generator
static const bool output_seed = true;
// Random number seed. Replace with a fixed value to generate
// same random numbers each time.
// TODO: not certain this works
static int seed = (int)time(0);
//static const int seed = 1580155444;
// Random number generator
static auto rng = std::default_random_engine(seed);
// Uniform distributions. Call with rng to generate random number in range
// E.g. double a = distribution_0_1(rng); //a is assigned random number [0,1]
static std::uniform_real_distribution<double> distribution_0_1(0.0, 1.0);
static std::uniform_real_distribution<double> distribution_neg1_1(-1.0, 1.0);
static std::uniform_int_distribution<int> distribution_0_grid_size(0, grid_size - 1);
#endif /* Parameters_h */
