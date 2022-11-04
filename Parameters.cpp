#include <stdio.h>
#include "Parameters.h"
#include <vector>
#include <string>
#include <sstream>
#include <unordered_set>
#include <ctime>
#include "yaml-cpp/yaml.h"

using namespace std;
extern std::string file_path;


Parameters::Parameters(string file_path)
{
  time_limit = 80.0; //100.00;
  timestep = 0.005;
  output_folder = "/Volumes/Main/model_git/Outputs/Grid_codes/grid"+ to_string(std::time(0));
  last_name = output_folder.c_str();
  output_interval = 0.5;
  fixed_vessel_placement = true;
  YAML::Node config = YAML::LoadFile(file_path);
  initial_number_bacteria = config["initial_number_bacteria"].as<int>();
  initial_number_Rmacrophages = config["initial_number_Rmacrophages"].as<int>();
  fixed_bacteria_placement = true;
  bacteria_rep_hours_min_early = config["bacteria_rep_hours_min_early"].as<double>();
  bacteria_rep_hours_max_early = config["bacteria_rep_hours_max_early"].as<double>();
  bacteria_rep_hours_min_middle = config["bacteria_rep_hours_min_middle"].as<double>();
  bacteria_rep_hours_max_middle = config["bacteria_rep_hours_max_middle"].as<double>();
  bacteria_rep_hours_min_late = config["bacteria_rep_hours_min_late"].as<double>();
  bacteria_rep_hours_max_late = config["bacteria_rep_hours_max_late"].as<double>();
  probability_hmacrophage_activation = config["probability_hmacrophage_activation"].as<double>();
  probability_rmacrophage_activation = config["probability_rmacrophage_activation"].as<double>();
  resting_neutrophil_max_life_span = config["resting_neutrophil_max_life_span"].as<double>();
  activated_neutrophil_lifespan = config["activated_neutrophil_lifespan"].as<double>();
  neutrophil_recruitment_prob_no_bac = config["neutrophil_recruitment_prob_no_bac"].as<double>();
  neutrophil_recruitment_prob_bac = config["neutrophil_recruitment_prob_bac"].as<double>();
  probability_neutrophil_activation = config["probability_neutrophil_activation"].as<double>();
  probability_neutrophil_destroy_dormant_bacteria = config["probability_neutrophil_destroy_dormant_bacteria"].as<double>();
  probability_neutrophil_destroy_replicating_bacteria = config["probability_neutrophil_destroy_replicating_bacteria"].as<double>();
  probability_hmac_destroy_dormant_bacteria = config["probability_hmac_destroy_dormant_bacteria"].as<double>();
  probability_hmac_destroy_replicating_bacteria = config["probability_hmac_destroy_replicating_bacteria"].as<double>();
  activated_neutrophil_movement_rate = config["activated_neutrophil_movement_rate"].as<double>();
  resting_neutrophil_movement_rate = config["resting_neutrophil_movement_rate"].as<double>();
  activatedresident_macrophage_movement_rate = config["activatedresident_macrophage_movement_rate"].as<double>();
  bacteria_replication_neighbourhood_max_depth = config["bacteria_replication_neighbourhood_max_depth"].as<int>();
  rmacrophages_recruitment_prob_bac = config["rmacrophages_recruitment_prob_bac"].as<double>();
  hmacrophages_recruitment_prob_bac = config["hmacrophages_recruitment_prob_bac"].as<double>();
  rmacrophages_recruitment_prob_no_bac = config["rmacrophages_recruitment_prob_no_bac"].as<double>();
  hmacrophages_recruitment_prob_no_bac = config["hmacrophages_recruitment_prob_no_bac"].as<double>();
  resting_rmacrophage_max_life_span = config["resting_rmacrophage_max_life_span"].as<double>();
  activated_Rmacrophage_lifespan = config["activated_Rmacrophage_lifespan"].as<double>();
  activated_Hmacrophage_lifespan = config["activated_Hmacrophage_lifespan"].as<double>();
  activated_hmacrophage_movement_rate = config["activated_hmacrophage_movement_rate"].as<double>();
  restingresident_macrophage_movement_rate = config["restingresident_macrophage_movement_rate"].as<double>();
  resting_hmacrophage_movement_rate = config["resting_hmacrophage_movement_rate"].as<double>();
  activated_Mastcell_lifespan = config["activated_Mastcell_lifespan"].as<double>();
  initial_number_mastcell = config["initial_number_mastcell"].as<int>();
  resting_mastcell_max_life_span = config["resting_mastcell_max_life_span"].as<double>();
  restingmastcell_movement_rate = config["restingmastcell_movement_rate"].as<double>();
  mastcells_recruitment_prob_bac = config["mastcells_recruitment_prob_bac"].as<double>();
  mastcells_recruitment_prob_no_bac = config["mastcells_recruitment_prob_no_bac"].as<double>();
  activatedmastcell_movement_rate = config["activatedmastcell_movement_rate"].as<double>();
  spatial_step = config["spatial_step"].as<double>();
  chemokine_diffusion = config["chemokine_diffusion"].as<double>();
  chemokine_rmac_source = config["chemokine_rmac_source"].as<double>();
  chemokine_decay = config["chemokine_decay"].as<double>();


  Hmacrophage_chemotactic_movement_bias = config["Hmacrophage_chemotactic_movement_bias"].as<double>();
  resting_rmacrophage_movement_bias = config["resting_rmacrophage_movement_bias"].as<double>();
  activated_rmacrophage_movement_bias = config["activated_rmacrophage_movement_bias"].as<double>();
  neutrophil_chemotactic_movement_bias = config["neutrophil_chemotactic_movement_bias"].as<double>();

  prob_bacteria_penetrate_wall = config["prob_bacteria_penetrate_wall"].as<double>();
  bacteria_shedding_hours_min = config["bacteria_shedding_hours_min"].as<int>();
  bacteria_shedding_hours_max = config["bacteria_shedding_hours_max"].as<int>();


  //distribution_0_1 = std::uniform_real_distribution<double>(0.0, 1.0);
  //distribution_neg1_1 = std::uniform_real_distribution<double>(-1.0, 1.0);
  //distribution_0_grid_size = std::uniform_int_distribution<int>(0, grid_size - 1);
}
