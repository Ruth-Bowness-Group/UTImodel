
#include <iostream>
#include "Environment.h"
#include "Agent.h"
#include <stdexcept> // std::runtime_error
#include <sstream> // std::stringstream
#include <chrono>
#include <sys/types.h>
#include <sys/stat.h>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <string>
#include <unistd.h>
#include <vector>
#include <ctime>       /* time_t, struct tm, difftime, time, mktime */
#include <utility> // std::pair

#include "UTImodel_debug.h"

using namespace std;

/**
 * Neighbourhood functions:
 */

/**
* Moore Neighbourhood
*/
std::vector<std::pair<int, int>> mooreNeighbourhood(int x, int y, int level, int grid_size, std::default_random_engine rng)
{

    std::vector<std::pair<int, int>> neighbours;
    for (int i = -level; i <= level; i++)
    {
        for (int j = -level; j <= level; j++)
        {
            // Check that neighbour is not the original loc and is on the grid
            if (!(i == 0 && j == 0) and (x + i) >= 0 and (x + i) < grid_size and (y + j) >= 0 and (y + j) < grid_size)
            {
                // Create a pair of coordinates and add it to the vector
                std::pair<int, int> neighbour = std::make_pair(x + i, y + j);
                neighbours.push_back(neighbour);
            }

        }

    }
    // Shuffle the neighbours
    std::shuffle(neighbours.begin(), neighbours.end(), rng);
    return neighbours;
}

/**
 * Von Neumann Neighbourhood Function
 */

 std::vector<std::pair<int, int>> vonNeumannNeighbourhood(int x, int y, int level, int grid_size, std::default_random_engine rng)
 {
     std::vector<std::pair<int, int>> neighbours;
     for (int i = -level; i <= level; ++i)
     {
         int r_i = level - abs(i);
         for (int j = -r_i; j <= r_i; ++j)
         {
             // Check that neighbour is not the original loc and is on the grid
             if (!(i == 0 && j == 0) and (x + i) >= 0 and (x + i) < grid_size and (y + j) >= 0 and (y + j) < grid_size)
             {
                 // Create a pair of coordinates and add it to the vector
                 std::pair<int, int> neighbour = std::make_pair(x + i, y + j);
                 neighbours.push_back(neighbour);
             }
         }
     }
     // Shuffle the neighbours
     std::shuffle(neighbours.begin(), neighbours.end(), rng);
     return neighbours;
 }

 /**
  * Finite difference diffusion functions
  */
 // Calculate the rate of diffusion of a molecule in the centre of the grid (i.e not boundary) based on values at neighbours
 double diffusion_centre(double value_at_location, double value_N, double value_E, double value_S, double value_W, double diffusion_at_location, double diffusion_N, double diffusion_E, double diffusion_S, double diffusion_W, double spatial_step){
     return ((diffusion_at_location+diffusion_S)/2*(value_S-value_at_location)-(diffusion_at_location+diffusion_N)/2*(value_at_location-value_N))/(spatial_step*spatial_step)+((diffusion_at_location+diffusion_E)/2*(value_E-value_at_location)-(diffusion_at_location+diffusion_W)/2*(value_at_location-value_W))/(spatial_step*spatial_step);
 }

 // Calculate the rate of diffusion of a molecule at the boundary of the grid based on values at neighbours. 4 neighbours defined, but will be at least 1 duplicate (2 duplicates for corners)
 double diffusion_boundary(double value_at_location, double value_N1, double value_N2, double value_N3, double value_N4, double diffusion_at_location, double spatial_step){
     return diffusion_at_location*(value_N1 - 2*value_at_location + value_N2)/(spatial_step*spatial_step) + diffusion_at_location*(value_N3 - 2*value_at_location + value_N4)/(spatial_step*spatial_step);
 }



Environment::Environment(Parameters* param): p(param)
{
  // Fill each grid cell with location instatiations
  for (int x = 0; x < grid_size; ++x)
  {
    for (int y = 0; y < grid_size; ++y)
    {
      grid[x][y] = new Location(x, y, p);
      grid[x][y] -> immediate_moore = mooreNeighbourhood(x, y, 1, grid_size, rng);
    }
  }
}

Environment::~Environment()
{
  for (int x = 0; x < grid_size; x++) {
    for (int y = 0; y < grid_size; y++) {
      delete grid[x][y]; /* CHRIS: this calls the Location destructor (see Location.h and Location.cpp) */
    }
  }
  if (p != nullptr)
   p = nullptr;

 /* CHRIS: deleting the Bacterium here, and Immunecell below, calls the destructors for the Bacterium and Immunecells
  * as well as the destructor for Agent (see Agent.h and Agent.cpp)
  */
 for (Bacterium* bac : extracellular_bacteria) {
   delete bac;
 }

 for (Rmacrophage* rmac : rmacrophages) {
   delete rmac;
 }

 for (Hmacrophage* hmac : hmacrophages) {
   delete hmac;
 }
 for (Neutrophil* neut : neutrophils) {
   delete neut;
 }

 for (Mastcell* mast : mastcells) {
   delete mast;
 }
}

// Get a location based on x,y coordinates
 Location *Environment::getLocation(int x, int y){return grid[x][y];}

 std::tuple<bool, int, int> Environment::findSpace(std::vector<std::pair<int, int>> neighbours)
 {
   // TODO: may be better to return ALL values that have space and not just one
   bool space_found = false;
   std::pair<int, int> neighbour;
   // Keep searching until we find space or run out of neighbours
   while (space_found == false && !neighbours.empty())
   {
     // Pull neighbour from stack
     neighbour = neighbours.back();
     // Remove stack item
     neighbours.pop_back();
     // Choose the neighbour if its empty
     if (grid[neighbour.first][neighbour.second] -> getContents() == LocationContents::empty)
     {
       space_found = true;
     }
   }
   // Wipe the neighbour if no space found
   if (!space_found)
   {
     neighbour = std::make_pair(NULL, NULL);
   }
   // Create a tuple of space found, x coord, y coord
   return std::make_tuple(space_found, neighbour.first, neighbour.second);
 }

// Add Bacterium

Bacterium *Environment::addBacterium(int x, int y, bool is_replicating)
{
  Bacterium *bac = new Bacterium(is_replicating, p);
  //Add to grid
  try
  {
    grid[x][y] -> addAgent(bac);
  }
  catch (std::runtime_error &err)
  {
    throw std::runtime_error(std::string("addBacterium: ") + err.what());
  }
  return bac;
}

//Remove Bacterium

void Environment::eraseBacterium(Bacterium *bac)
{
  std::pair<int,int> coords = bac -> getCoordinates();
  //remove from the grid
  try
  {
    grid[coords.first][coords.second] -> removeAgent();
  }
  catch (std::runtime_error &err)
  {
    throw std::runtime_error(std::string("eraseBacterium: ") + err.what());
  }


}

//Remove neutrophil

void Environment::NeutrophilDeath(Neutrophil *neut)
{
  int x = neut->getCoordinates().first;
  int y = neut->getCoordinates().second;

  try
  {
    grid[x][y]->removeAgent();
  }
  catch (std::runtime_error &err)
  {
    throw std::runtime_error(std::string("NeutrophilDeath: ") + err.what());
  }
}
//Add Neutrophil
Neutrophil *Environment::addNeutrophil(int x, int y)
{
  //creation
  Neutrophil *neut = new Neutrophil(distribution_0_1(rng) * p->resting_neutrophil_max_life_span, p);
  //add to grid
  try
  {
    grid[x][y]->addAgent(neut);
  }
  catch (std::runtime_error &err)
  {
    throw std::runtime_error(std::string("add Neutrophil: ") + err.what());
  }
  return neut;
}

//Remove helper macrophage

void Environment::HmacrophageDeath(Hmacrophage *hmac)
{
  int x = hmac ->getCoordinates().first;
  int y = hmac->getCoordinates().second;

  try
  {
    grid[x][y]->removeAgent();
  }
  catch (std::runtime_error &err)
  {
    throw std::runtime_error(std::string("HmacrophageDeath: ") + err.what());
  }
}
//Add Helper macrophage
Hmacrophage *Environment::addHmacrophage(int x, int y)
{
  //creation
  Hmacrophage *hmac = new Hmacrophage(distribution_0_1(rng) * p->resting_rmacrophage_max_life_span, p);
  //add to grid
  try
  {
    grid[x][y]->addAgent(hmac);
  }
  catch (std::runtime_error &err)
  {
    throw std::runtime_error(std::string("add Hmacrophage: ") + err.what());
  }
  return hmac;
}


//Add Mastcell
Mastcell *Environment::addMastcell(int x, int y)
{
  //Creation (random lifespan based on parameter)
  Mastcell *mast = new Mastcell(distribution_0_1(rng) * p->resting_mastcell_max_life_span, p);
  //Add to the grid
  try
  {
    grid[x][y]->addAgent(mast);
  }
  catch (std::runtime_error &err)
  {
    throw std::runtime_error(std::string("add Mastcell: ") + err.what());
  }
  return mast;
}




//Add Resident Macrophage
Rmacrophage *Environment::addRmacrophage(int x, int y)
{
  //Creation (random lifespan based on parameter)
  Rmacrophage *rmac = new Rmacrophage(distribution_0_1(rng) * p->resting_rmacrophage_max_life_span, p);
  //Add to the grid
  try
  {
    grid[x][y]->addAgent(rmac);
  }
  catch (std::runtime_error &err)
  {
    throw std::runtime_error(std::string("add Rmacrophage: ") + err.what());
  }
  return rmac;
}


//Remove Mastcell
void Environment::MastcellDeath(Mastcell *mast)
{
  int x = mast ->getCoordinates().first;
  int y = mast->getCoordinates().second;

  try
  {
    grid[x][y]->removeAgent();
  }
  catch (std::runtime_error &err)
  {
    throw std::runtime_error(std::string("MastcellDeath: ") + err.what());
  }
}





//Remove Rmacrophage
void Environment::RmacrophageDeath(Rmacrophage *rmac)
{
  int x = rmac ->getCoordinates().first;
  int y = rmac->getCoordinates().second;

  try
  {
    grid[x][y]->removeAgent();
  }
  catch (std::runtime_error &err)
  {
    throw std::runtime_error(std::string("RmacrophageDeath: ") + err.what());
  }
}

void Environment::moveAgent(int x_from, int y_from, int x_to, int y_to)
{
  try
  {
    Agent *agent_ref = grid[x_from][y_from] -> removeAgent();
    grid[x_to][y_to] -> addAgent(agent_ref);
  }
  catch (std::runtime_error &err)
  {
    throw std::runtime_error(std::string("moveAgent: ") + err.what());
  }
}


void Environment::initialiseBacterialClusterpostShedding(int cluster_x, int cluster_y)
{
  int bacteria_placed;
  std::vector<std::pair<int,int>> cluster_locs, neighbours;
  std::tuple<bool,int,int> space_finder;
  Bacterium *bac;
  bacteria_placed = 0;

  int attempts = 0;

  while (grid[cluster_x][cluster_y] -> getContents() != LocationContents::empty && attempts < grid_size*grid_size)
  {
    neighbours = grid[cluster_x][cluster_y] -> immediate_moore;
    space_finder = findSpace(neighbours);
    if (std::get<0>(space_finder))
    {
      cluster_x = std::get<1>(space_finder);
      cluster_y = std::get<2>(space_finder);
    }
    else
    {
      // No space at the coordinates or any neighbours, so pick a random location
       cluster_x = distribution_0_grid_size(rng);
       cluster_y = distribution_0_grid_size(rng);
    }
    attempts++;
  }

  if (attempts == grid_size*grid_size)
  {
    throw std::runtime_error("no empty spaces for cluster");
  }

  bac = addBacterium(cluster_x, cluster_y, false);
  extracellular_bacteria.insert(bac);
  bacteria_placed++;
  cluster_locs.push_back(std::make_pair(cluster_x,cluster_y));

  while (bacteria_placed < bacteria_at_level_plusone)
  {
    neighbours = vonNeumannNeighbourhood(cluster_locs[0].first,cluster_locs[0].second,1,grid_size,rng);
    space_finder = findSpace(neighbours);
    if (std::get<0>(space_finder))
    {
      bac = addBacterium(std::get<1>(space_finder), std::get<2>(space_finder), false);
      bacteria_placed++;
      extracellular_bacteria.insert(bac);
      cluster_locs.push_back(std::make_pair(std::get<1>(space_finder), std::get<2>(space_finder)));
    }
    else
    {
      cluster_locs.erase(cluster_locs.begin());
    }
  }
}


void Environment::initialiseBacterialCluster(int cluster_x, int cluster_y)
{
  int bacteria_placed;
  std::vector<std::pair<int,int>> cluster_locs, neighbours;
  std::tuple<bool,int,int> space_finder;
  Bacterium *bac;
  bacteria_placed = 0;

  int attempts = 0;

  while (grid[cluster_x][cluster_y] -> getContents() != LocationContents::empty && attempts < grid_size*grid_size)
  {
    neighbours = grid[cluster_x][cluster_y] -> immediate_moore;
    space_finder = findSpace(neighbours);
    if (std::get<0>(space_finder))
    {
      cluster_x = std::get<1>(space_finder);
      cluster_y = std::get<2>(space_finder);
    }
    else
    {
      // No space at the coordinates or any neighbours, so pick a random location
       cluster_x = distribution_0_grid_size(rng);
       cluster_y = distribution_0_grid_size(rng);
    }
    attempts++;
  }

  if (attempts == grid_size*grid_size)
  {
    throw std::runtime_error("no empty spaces for cluster");
  }

  bac = addBacterium(cluster_x, cluster_y, false);
  extracellular_bacteria.insert(bac);
  bacteria_placed++;
  cluster_locs.push_back(std::make_pair(cluster_x,cluster_y));

  while (bacteria_placed < p->initial_number_bacteria)
  {
    neighbours = vonNeumannNeighbourhood(cluster_locs[0].first,cluster_locs[0].second,1,grid_size,rng);
    space_finder = findSpace(neighbours);
    if (std::get<0>(space_finder))
    {
      bac = addBacterium(std::get<1>(space_finder), std::get<2>(space_finder), false);
      bacteria_placed++;
      extracellular_bacteria.insert(bac);
      cluster_locs.push_back(std::make_pair(std::get<1>(space_finder), std::get<2>(space_finder)));
    }
    else
    {
      cluster_locs.erase(cluster_locs.begin());
    }
  }
}




void Environment::initialiseMastcell(int number_of_mastcells)
{
  int x, y;
  int mastcells_placed = 0;
  Mastcell *mast;
  while (mastcells_placed < number_of_mastcells)
  {
    // Get a random location
    x = distribution_0_grid_size(rng);
    y = distribution_0_grid_size(rng);
    if (grid[x][y] -> getContents() == LocationContents::empty)
    {
      mast = addMastcell(x,y);
      mastcells.insert(mast);
      mastcells_placed++;
    }
  }
}




void Environment::initialiseRmacrophage(int number_of_rmacrophages)
{
  int x, y;
  int rmacrophages_placed = 0;
  Rmacrophage *rmac;
  while (rmacrophages_placed < number_of_rmacrophages)
  {
    // Get a random location
    x = distribution_0_grid_size(rng);
    y = distribution_0_grid_size(rng);
    if (grid[x][y] -> getContents() == LocationContents::empty)
    {
      rmac = addRmacrophage(x,y);
      rmacrophages.insert(rmac);
      rmacrophages_placed++;
    }
  }
}

void Environment::bacterialShedding()
{
  int shedding_age,x,y;
  int cluster_x, cluster_y;
  std::vector<Bacterium *> deleted_bacteria;
  int counter_bacteria_shedded = 0;
  Bacterium *bac;
  std::uniform_int_distribution<int> shedding_dist(p->bacteria_shedding_hours_min, p->bacteria_shedding_hours_max);
  shedding_age = shedding_dist(rng) / p->timestep;

  if ((int)extracellular_bacteria.size() > 6000){
    std::cout << "||                           ||" << '\n';
    std::cout << "  ||                       ||" << '\n';
    std::cout << "    ||                   ||" << '\n';
    std::cout << "      ||               ||" << '\n';
    std::cout << "        ||           ||" << '\n';
    std::cout << "          ||SHEDDIN||" << '\n';
    for (auto iter = extracellular_bacteria.begin();iter != extracellular_bacteria.end();)
    {
      bac = *iter;

      deleted_bacteria.push_back(bac);

      counter_bacteria_shedded += 1;
      //std::cout << counter_bacteria_shedded << '\n';
      ++iter;
    }
    std::cout << counter_bacteria_shedded << '\n';
    for (auto bac:deleted_bacteria)
    {
      extracellular_bacteria.erase(bac);
      eraseBacterium(bac);
    }
    if (extracellular_bacteria.empty()){
      std::cout << "All bacteria has been successufully cleared" << '\n';
      counter_epithelial_level_shedding++;
      try {
        // Distribute bacteria in a cluster, either randomly positioned or at a set point
        if (p->fixed_bacteria_placement){
            cluster_x = p->fixed_bacteria_cluster.first;
            cluster_y = p->fixed_bacteria_cluster.second;
        } else {
            cluster_x = distribution_0_grid_size(rng);
            cluster_y = distribution_0_grid_size(rng);
        }
        initialiseBacterialClusterpostShedding(cluster_x, cluster_y);
      }
      catch (std::exception &err) {
          std::cout << "ERROR: initialisation: " << err.what() << "\n";

        }
      }

  }
}




void Environment::bacterialReplication(double time)
{
  int x,y,rep_age,neighbourhood_level, neighbour_x, neighbour_y;
  std::vector<std::pair<int,int>> neighbours;
  std::vector<Bacterium *> new_bacteria;
  // Flag to indicate a free space has been found
  bool space_found;
  std::tuple<bool, int, int> space_finder;
  Bacterium *bac;

  for (auto iter = extracellular_bacteria.begin();iter != extracellular_bacteria.end();)
  {
    bac = *iter;

    // Pull coordinates
    x = bac -> getCoordinates().first;
    y = bac -> getCoordinates().second;

    double min_rep_rate = 0.0;
    double max_rep_rate = 0.0;

    if ((int) round(time) < 10)
    {
      min_rep_rate = p->bacteria_rep_hours_min_early;
      max_rep_rate = p->bacteria_rep_hours_max_early;
    }
    else
    {
      min_rep_rate = p->bacteria_rep_hours_min_late;
      max_rep_rate = p->bacteria_rep_hours_max_late;
    }

    std::uniform_real_distribution<double> rep_dist(min_rep_rate, max_rep_rate);
    rep_age = rep_dist(rng) / p->timestep;

    if ((int)round(time / p->timestep) % rep_age == 0)
    { // Rounding needed to avoid modulo issues
      neighbourhood_level = 1;
      // Search progressive neighbourhood depths to find a free spot, stop if
      // empty neighbour found or max depth reached
      do
      {
        // Pull the (shuffled) neighbour locations based on the bacterium's neighbourhood state
        if (bac->isMooreReplicationNeighbourhood())
        {
          neighbours = mooreNeighbourhood(x,y,neighbourhood_level,grid_size,rng);
        }
        else
        {
          neighbours = vonNeumannNeighbourhood(x,y,neighbourhood_level,grid_size,rng);
        }
        // Look for space in neigbhourhood. Tuple of <space_found, neighbour_x, neighbour_y>.
        space_finder = findSpace(neighbours);

        space_found = std::get<0>(space_finder);

        // Increase the neighbourhood level
        neighbourhood_level++;
      }
      while (space_found == false && neighbourhood_level <= p->bacteria_replication_neighbourhood_max_depth);
      if (space_found)
      {
          neighbour_x = std::get<1>(space_finder);
          neighbour_y = std::get<2>(space_finder);

          try
          {
            bac = addBacterium(neighbour_x, neighbour_y, bac->isReplicating());
            new_bacteria.push_back(bac);
          }
          catch (std::runtime_error &err)
          {
            throw err;
          }
          bac->switchReplicationNeighbourhood();
      }
      else
      {
        bac->switchResting();
      }
      bac->setAge(0);
    }

    ++iter;
  }
  for (auto bac:new_bacteria)
  {
    extracellular_bacteria.insert(bac);
    //std::cout << "testing replication speed" << '\n';
  }
}


void Environment::NeutrophilRecruitment()
{
  double r;
  std::vector<std::pair<int,int>> neighbours;
  std::tuple<bool, int, int> space_finder;
  int num_bac = (int)extracellular_bacteria.size();
  Neutrophil *neut;

  for (auto vessel_coords : vessels)
  {
    r = distribution_0_1(rng);

    if ((num_bac == 0 && r < p->neutrophil_recruitment_prob_no_bac) ||
        (num_bac > 0 && r < p->neutrophil_recruitment_prob_bac))
        {
            // Get the immediate neighbours
            neighbours = grid[vessel_coords.first][vessel_coords.second]->immediate_moore;
            // Look for an empty neighbour
            space_finder = findSpace(neighbours);
            // If space has been found
            if (std::get<0>(space_finder))
            {
                // Add a new neutrophil
                try
                {
                    neut = addNeutrophil(std::get<1>(space_finder), std::get<2>(space_finder));
                    neutrophils.insert(neut);
                }
                catch (std::runtime_error const &e)
                {
                    throw e;
                } // Trycatch
            }     // if space
        }         // If recruit
    }             // for vessels
}





void Environment::HmacrophageRecruitment()
{
  double r;
  std::vector<std::pair<int,int>> neighbours;
  std::tuple<bool, int, int> space_finder;
  int num_bac = (int)extracellular_bacteria.size();
  Hmacrophage *hmac;

  for (auto vessel_coords : vessels)
  {
    r = distribution_0_1(rng);

    if ((num_bac == 0 && r < p->hmacrophages_recruitment_prob_no_bac) ||
        (num_bac > 0 && r < p->hmacrophages_recruitment_prob_bac))
        {
            // Get the immediate neighbours
            neighbours = grid[vessel_coords.first][vessel_coords.second]->immediate_moore;
            // Look for an empty neighbour
            space_finder = findSpace(neighbours);
            // If space has been found
            if (std::get<0>(space_finder))
            {
                // Add a new macrophage
                try
                {
                    hmac = addHmacrophage(std::get<1>(space_finder), std::get<2>(space_finder));
                    hmacrophages.insert(hmac);
                }
                catch (std::runtime_error const &e)
                {
                    throw e;
                } // Trycatch
            }     // if space
        }         // If recruit
    }             // for vessels
}



void Environment::MastcellRecruitment()
{
  double r;
  std::vector<std::pair<int,int>> neighbours;
  std::tuple<bool, int, int> space_finder;
  int num_bac = (int)extracellular_bacteria.size();
  Mastcell *mast;

  for (auto vessel_coords : vessels)
  {
    r = distribution_0_1(rng);

    if ((num_bac == 0 && r < p->mastcells_recruitment_prob_no_bac) ||
        (num_bac > 0 && r < p->mastcells_recruitment_prob_bac))
        {
            // Get the immediate neighbours
            neighbours = grid[vessel_coords.first][vessel_coords.second]->immediate_moore;
            // Look for an empty neighbour
            space_finder = findSpace(neighbours);
            // If space has been found
            if (std::get<0>(space_finder))
            {
              try
              {
                  mast = addMastcell(std::get<1>(space_finder), std::get<2>(space_finder));
                  mastcells.insert(mast);
              }
              catch (std::runtime_error const &e)
              {
                throw e;
              }
            }
          }
        }
}


void Environment::updateAttributes(double time){
  // Reset max values
  double max_chemokine = 0.0;
  //Temporary arrays to hold calculated values
  double temp_chemokine[grid_size][grid_size] = {{0.0}};

  double chemokine_value_at_location, chemokine_diffused_value;

  std::pair<int, int> neighbour1_coords, neighbour2_coords, neighbour3_coords, neighbour4_coords;

  Location* neighbour1;
  Location* neighbour2;
  Location* neighbour3;
  Location* neighbour4;

  for (int x = 0; x <= grid_size-1;x++){
    for (int y = 0; y <= grid_size-1;y++){
      chemokine_value_at_location = grid[x][y]->chemokine_value;

      // Diffusion
      if (x > 0 && x < grid_size-1 && y > 0 && y < grid_size-1){
        neighbour1 = grid[x][y-1];
        neighbour2 = grid[x+1][y];
        neighbour3 = grid[x][y+1];
        neighbour4 = grid[x-1][y];

        chemokine_diffused_value = diffusion_centre(chemokine_value_at_location, neighbour1->chemokine_value, neighbour2->chemokine_value, neighbour3->chemokine_value, neighbour4->chemokine_value, p->chemokine_diffusion, p->chemokine_diffusion, p->chemokine_diffusion, p->chemokine_diffusion, p->chemokine_diffusion,p->spatial_step);
      } else {
        // Location is on the boundary, so will have 2/3 neighbours
        // Left column - so only use east neighbour
        if (x==0){
          neighbour1_coords = {x+1, y};
          neighbour2_coords = {x+1, y};
        // Right column so only use West neighbour
      } else if (x == grid_size-1){
          neighbour1_coords = {x-1, y};
          neighbour2_coords = {x-1, y};
      } else {
          neighbour1_coords = {x-1, y};
          neighbour2_coords = {x+1, y};
      }

      // Top row - so only use south neighbour
      if (y == 0){
          neighbour3_coords = {x, y+1};
          neighbour4_coords = {x, y+1};
      // Bottom row - so only use north neighbour
      } else if (y == grid_size-1){
          neighbour3_coords = {x, y-1};
          neighbour4_coords = {x, y-1};
      } else {
          neighbour3_coords = {x, y-1};
          neighbour4_coords = {x, y+1};
      }

      neighbour1 = grid[neighbour1_coords.first][neighbour1_coords.second];
      neighbour2 = grid[neighbour2_coords.first][neighbour2_coords.second];
      neighbour3 = grid[neighbour3_coords.first][neighbour3_coords.second];
      neighbour4 = grid[neighbour4_coords.first][neighbour4_coords.second];

      // Calculate diffusion values using determined neighbours
      chemokine_diffused_value = diffusion_boundary(chemokine_value_at_location, neighbour1->chemokine_value, neighbour2->chemokine_value, neighbour3->chemokine_value, neighbour4->chemokine_value, p->chemokine_diffusion, p->spatial_step);
      }
      // CHEMOKINE : current value + timestep*(diffusion + source - uptake - decay)
      temp_chemokine[x][y] = chemokine_value_at_location + p->timestep * (chemokine_diffused_value + grid[x][y]->chemokine_source_at_location - p->chemokine_decay*chemokine_value_at_location);
      if (temp_chemokine[x][y] > max_chemokine){
                max_chemokine = temp_chemokine[x][y];
            }
            if (temp_chemokine[x][y] < 0){
              std::cout <<  temp_chemokine[x][y]   << '\n';
                throw std::runtime_error("Chemokine value cannot be below zero");
            }

          }
        }

      // Update the grid with new attribute values
      for (int x = 0; x <= grid_size-1; x++){
          for (int y = 0; y <= grid_size-1; y++){
              grid[x][y]->chemokine_value = temp_chemokine[x][y];
              grid[x][y]->scaled_chemokine_value = temp_chemokine[x][y] / max_chemokine * 100;

          }
      }

}




void Environment::RmacrophageRecruitment()
{
  double r;
  std::vector<std::pair<int,int>> neighbours;
  std::tuple<bool, int, int> space_finder;
  int num_bac = (int)extracellular_bacteria.size();
  Rmacrophage *rmac;

  for (auto vessel_coords : vessels)
  {
    r = distribution_0_1(rng);

    if ((num_bac == 0 && r < p->rmacrophages_recruitment_prob_no_bac) ||
        (num_bac > 0 && r < p->rmacrophages_recruitment_prob_bac))
        {
            // Get the immediate neighbours
            neighbours = grid[vessel_coords.first][vessel_coords.second]->immediate_moore;
            // Look for an empty neighbour
            space_finder = findSpace(neighbours);
            // If space has been found
            if (std::get<0>(space_finder))
            {
              try
              {
                  rmac = addRmacrophage(std::get<1>(space_finder), std::get<2>(space_finder));
                  rmacrophages.insert(rmac);
              }
              catch (std::runtime_error const &e)
              {
                throw e;
              }
            }
          }
        }
}

void Environment::setAVessel(int x, int y)
{
  grid[x][y] -> setVessel();
  vessels.push_back(std::make_pair(x,y));
}

void Environment::initialiseBloodVesselsRandom(int number_of_vessels)
{
  int vessels_placed = 0;
  while (vessels_placed < number_of_vessels)
  {
    // Get a random location
    int x = distribution_0_grid_size(rng);
    int y = distribution_0_grid_size(rng);
    // Only place if not already a vessel
    if (!(grid[x][y] -> getContents() == LocationContents::vessel))
    {
      setAVessel(x, y);
      ++vessels_placed;
    }
  }
}

// Set the initial blood vessels based on fixed positions
void Environment::initialiseBloodVesselsFixed(std::vector<std::pair<int, int>> locations)
{
  for (auto iter : locations)
  {
    setAVessel(iter.first, iter.second);
  }
}

std::pair<int, int> Environment::randomMovement(std::vector<std::pair<int, int>> neighbours)
{
  std::pair<int,int> chosen_neighbour;
  std::uniform_int_distribution<int> distribution_neighbours(0, (int)neighbours.size() - 1);
  // Chose a neighbour at random
  chosen_neighbour = neighbours[distribution_neighbours(rng)];
  return chosen_neighbour;
}

void Environment::NeutrophilAction(double time){
    int x, y, r, neighbour_count_activate;
    double movement_rate;
    std::vector<std::pair<int, int> > neighbours;
    std::pair<int, int> chosen_neighbour;
    LocationContents neighbour_contents;
    Bacterium* bac;
    Rmacrophage* rmac;
    Hmacrophage* hmac;
    Neutrophil* neut;
    Mastcell* mast;
    int counter_death = 0;

    //loop through all hmacrophages
    for(auto iter=neutrophils.begin(); iter!=neutrophils.end();){
          neut = *iter;

          if (neut->exceededLifespan()){
            ++counter_death;
            NeutrophilDeath(neut);
            iter = neutrophils.erase(iter);
            delete neut;
            continue;
          }
          switch (neut -> getState()){
            case nresting:
              movement_rate = p->resting_neutrophil_movement_rate;
              break;
            case nactivated:
              movement_rate = p->activated_neutrophil_movement_rate;
              break;
            default:
              throw std::runtime_error("Invalid Neutrophil state");
          }
          // Movement
        if ((int)round(time/p->timestep)%(int)round(movement_rate/p->timestep)
            ==0){
            // Pull macrophage co-ordinates
            //std::cout << "/* Checkpoint1 */" << '\n';
            x = neut->getCoordinates().first;
            y = neut->getCoordinates().second;
            // Get neighbours
            neighbours = grid[x][y]->immediate_moore;

            // Chose a neighbour to move to (chemotactically or randomly)
            chosen_neighbour = chemotacticMovement(neighbours, p->neutrophil_chemotactic_movement_bias);

            // Check the contents of the neighbour
            neighbour_contents = grid[chosen_neighbour.first][chosen_neighbour.second]->getContents();
            //std::cout << "###" << static_cast<int>(neighbour_contents) << std::endl;
            // Location is empty, move the neutrophil
            if (neighbour_contents == LocationContents::empty){
                //std::cout << "Macrophage moves" << std::endl;
                moveAgent(x, y, chosen_neighbour.first, chosen_neighbour.second);

            // Location has agent (i.e. a bacterium)
          } else if (neighbour_contents == LocationContents::agent){
                // For activated macrophages, bacteria may be destroyed (or otherwise left)
                //std::cout << static_cast<int>(mac->getState()) << std::endl;
                if (neut->getState()==nactivated){
                    // TODO: not keen on casting here (or below)

                    /* CHRIS: A static_cast assumes that the type which we are casting to is known at compile time (hence why it is called static).
                     * However, an Agent may be either a Bacterium or Immunecell (as they are derived from Agent). Therefore, at compile time
                     * we cannot know what grid[.][.]->getAgent() will return (it may be a Bacterium or an Immunecell). In such a case, a static_cast
                     * to a Bacterium, when the Agent is actually an Immunecell, will fail. Thus, a dynamic_cast should
                     * be used. A dynamic_cast has the additional benefit that if the casting fails then it will return a nullptr. This nullptr
                     * can be then tested for later.
                     */

                    bac = dynamic_cast<Bacterium*>(grid[chosen_neighbour.first][chosen_neighbour.second]->getAgent());

                    /* CHRIS: These are a few print statements showing when the cast to a Bacterium has succeeded or failed. If you want, you
                     * can change the dynamic cast above to a static cast and see how the output changes.
                     */
                     if (bac == nullptr) {
                       UTIMODEL_DBG_MSG("cast to bacteria has failed");
                       UTIMODEL_DBG_PRINT("\tbac: (%d, %d)\n", chosen_neighbour.first, chosen_neighbour.second);
                       Rmacrophage* rmac = dynamic_cast<Rmacrophage*>(grid[chosen_neighbour.first][chosen_neighbour.second]->getAgent());
                       if (rmac == nullptr) {
                         Neutrophil* neut = dynamic_cast<Neutrophil*>(grid[chosen_neighbour.first][chosen_neighbour.second]->getAgent());
                         UTIMODEL_DBG_MSG("cast to Neutrophil succeeded");
                         if (neut == nullptr){
                           Mastcell* mast = dynamic_cast<Mastcell*>(grid[chosen_neighbour.first][chosen_neighbour.second]->getAgent());
                           UTIMODEL_DBG_MSG("cast to Mastcell succeeded");
                           if (mast == nullptr){
                            Hmacrophage* tmp = dynamic_cast<Hmacrophage*>(grid[chosen_neighbour.first][chosen_neighbour.second]->getAgent());
                            UTIMODEL_DBG_MSG("cast to Hmacrophage succeeded");
                            if (tmp == nullptr){
                              UTIMODEL_DBG_MSG("cast to all immune cells has failed");
                            }else{
                              UTIMODEL_DBG_MSG("cast to immune cell has succeeded");
                              UTIMODEL_DBG_PRINT("\tmac: (%d, %d)\n", chosen_neighbour.first, chosen_neighbour.second);
                              //std::cout << "/* Checkpoint3 */" << '\n';
                         }
                      }
                     }
                   }
                 } else {
                       UTIMODEL_DBG_MSG("cast to bacteria has succeeded");
                       UTIMODEL_DBG_PRINT("\tbac: (%d, %d)\n", chosen_neighbour.first, chosen_neighbour.second);
                     }

                    // std::cout << "this loop is useful" << std::endl;
                    // Random number 0-1
                    r = distribution_0_1(rng);

                    /* CHRIS: As we have now used a dynamic_cast above, if the cast fails the Bacterium* bac becomes a nullptr. Therefore,
                     * we need to check for it, otherwise we would attempt to call bac->isReplicating() which would result in a seg fault
                     * due to isReplicating() not being a part of a nullptr.
                     */

                    if (bac != nullptr) { /* CHRIS: safety net in case casting fails */
                      // Determine if bacterium is destroyed based on its replication type
                      if ((bac->isReplicating() && r < p->probability_neutrophil_destroy_replicating_bacteria) ||
                          (!bac->isReplicating() && r < p->probability_neutrophil_destroy_dormant_bacteria)) {
                          // Remove from grid
                          eraseBacterium(bac);
                          // Remove from list

                          /* CHRIS: similar to Immunecell above, erasing from a Set does not clear the memory. Therefore, we need
                           * to manually delete the bac. Again, note that if we didn't have the safety net above that tested for
                           * for a nulltpr, we could potentially arrive at a seg fault
                           */

                          extracellular_bacteria.erase(bac);
                          delete bac;
                          x = neut->getCoordinates().first;
                          y = neut->getCoordinates().second;
                          // Move the neutrophil
                          moveAgent(x, y, chosen_neighbour.first, chosen_neighbour.second);
                          // Increment counter
                          counter_bacteria_destroyed_by_neutrophils++;
                      }
                    }


            }
          }

        }

        // Activation
        if (neut->getState()==nresting){
            // Reobtain co-ordinates (in case mac has moved)
            x = neut->getCoordinates().first;
            y = neut->getCoordinates().second;
            // Check immediate Moore neighbourhood
            // TODO: slight performance hit getting a neighbourhood every timestep
            neighbours = grid[x][y]->immediate_moore;
            // Count neighbours that can activate macrophages
            neighbour_count_activate = 0;
            for (auto iter:neighbours){
                if (grid[iter.first][iter.second]->getContents() == LocationContents::agent){
                    neighbour_count_activate += grid[iter.first][iter.second]->getAgent()->canActivatemacrophage();

                }
            }
            // Random number 0-1
            r = distribution_0_1(rng);

            if (r < p->probability_neutrophil_activation && neighbour_count_activate > 0){
                // Activate the macrophage
                //std::cout << "Checking for Neutrophil activation" << std::endl;
                neut->activate();
            }

          } else if (neut->getState()==nactivated){

            // Reobtain co-ordinates (in case mac has moved)
            x = neut->getCoordinates().first;
            y = neut->getCoordinates().second;
            // Check immediate Moore neighbourhood
            // TODO: slight performance hit getting a neighbourhood every timestep
            neighbours = grid[x][y]->immediate_moore;
            // Count neighbours that can activate macrophages
            neighbour_count_activate = 0;
            for (auto iter:neighbours){
                if (grid[iter.first][iter.second]->getContents() == LocationContents::agent){
                    neighbour_count_activate += grid[iter.first][iter.second]->getAgent()->canActivatemacrophage();
                }
            }
            // Random number 0-1
            r = distribution_0_1(rng);
            // Multiply the probability by the number of activating cells in neighbour
            if (neighbour_count_activate==0){
                // deactivate the macrophage
                // std::cout << "mac is getting deactivated" << std::endl;
                neut->deactivate();
            }
          }
        // Increment iterator
        ++iter;
}// immunecell loop
}


// Chemotactic movement - used for Hmacrophages and Neutrophils. Picks a neighbour to move to probabilistically based on the chemokine levels.
std::pair<int, int> Environment::chemotacticMovement(std::vector<std::pair<int, int>> neighbours, double movement_bias)
{
    double neighbour_chemokine, r;
    double total_chemokine = 0;
    double running_total = 0;
    std::pair<int, int> chosen_neighbour;

    // Check all neighbours for chemokine levels
    for (int i = 0; i < neighbours.size(); i++)
    {
        neighbour_chemokine = grid[neighbours[i].first][neighbours[i].second]->chemokine_value;
        total_chemokine += pow(neighbour_chemokine, movement_bias);
    }

    // Random movement (if all neighbours have 0 chemokine)
    if (total_chemokine == 0)
    {
        std::uniform_int_distribution<int> distribution_neighbours(0, (int)neighbours.size() - 1);
        // Chose a neighbour at random
        chosen_neighbour = neighbours[distribution_neighbours(rng)];
    }
    else
    { // Chemokine driven movement
        // Biased random walk - pick neighbour probabilistically based on chemokine levels
        r = distribution_0_1(rng) * total_chemokine;
        for (auto coords : neighbours)
        {
            running_total += pow(grid[coords.first][coords.second]->chemokine_value, movement_bias);
            if (running_total >= r)
            {
                chosen_neighbour = coords;
                break;
            }
        }
    }

    return chosen_neighbour;
}






void Environment::HmacrophageAction(double time){
    int x, y, neighbour_count_activate, empty_x, empty_y;
    double r;
    double movement_rate;
    std::vector<std::pair<int, int> > neighbours, neighbours2;
    std::tuple<bool,int,int> space_finder;
    std::pair<int, int> chosen_neighbour;
    LocationContents neighbour_contents;
    Bacterium* bac;
    Rmacrophage* rmac;
    Hmacrophage* hmac;
    Neutrophil* neut;
    int counter_death = 0;
    int neutrophils_placed = 0;
    //loop through all hmacrophages

      for(auto iter=hmacrophages.begin(); iter!=hmacrophages.end();){

            hmac = *iter;



            if (hmac->exceededLifespan()){
              ++counter_death;
              HmacrophageDeath(hmac);
              iter = hmacrophages.erase(iter);
              delete hmac;
              continue;
            }
            switch (hmac -> getState()){
              case hresting:
                movement_rate = p->resting_hmacrophage_movement_rate;
                //std::cout << p->resting_hmacrophage_movement_rate << '\n';
                break;
              case hactivated:
                movement_rate = p->activated_hmacrophage_movement_rate;
                break;
              default:
                throw std::runtime_error("Invalid Hmacrophage state");
            }
            // Movement

          if ((int)round(time/p->timestep)%(int)round(movement_rate/p->timestep)
              ==0){


              // Pull macrophage co-ordinates
              //std::cout << "/* Checkpoint1 */" << '\n';
              x = hmac->getCoordinates().first;
              y = hmac->getCoordinates().second;
              // Get neighbours
              neighbours = grid[x][y]->immediate_moore;

              // Chose a neighbour to move to (chemotactically or randomly)
              chosen_neighbour = chemotacticMovement(neighbours, p->Hmacrophage_chemotactic_movement_bias);

              // Check the contents of the neighbour
              neighbour_contents = grid[chosen_neighbour.first][chosen_neighbour.second]->getContents();
              //std::cout << "###" << static_cast<int>(neighbour_contents) << std::endl;
              // Location is empty, move the macrophage

              if (neighbour_contents == LocationContents::empty){
                  //std::cout << "Macrophage moves" << std::endl;
                  moveAgent(x, y, chosen_neighbour.first, chosen_neighbour.second);

              // Location has agent (i.e. a bacterium)
            } else if (neighbour_contents == LocationContents::agent){

                  // For activated macrophages, bacteria may be destroyed (or otherwise left)
                  //std::cout << static_cast<int>(mac->getState()) << std::endl;

                  if (hmac->getState()==hactivated){
                      // TODO: not keen on casting here (or below)

                      /* CHRIS: A static_cast assumes that the type which we are casting to is known at compile time (hence why it is called static).
                       * However, an Agent may be either a Bacterium or Immunecell (as they are derived from Agent). Therefore, at compile time
                       * we cannot know what grid[.][.]->getAgent() will return (it may be a Bacterium or an Immunecell). In such a case, a static_cast
                       * to a Bacterium, when the Agent is actually an Immunecell, will fail. Thus, a dynamic_cast should
                       * be used. A dynamic_cast has the additional benefit that if the casting fails then it will return a nullptr. This nullptr
                       * can be then tested for later.
                       */

                      bac = dynamic_cast<Bacterium*>(grid[chosen_neighbour.first][chosen_neighbour.second]->getAgent());

                      /* CHRIS: These are a few print statements showing when the cast to a Bacterium has succeeded or failed. If you want, you
                       * can change the dynamic cast above to a static cast and see how the output changes.
                       */
                       if (bac == nullptr) {
                         UTIMODEL_DBG_MSG("cast to bacteria has failed");
                         UTIMODEL_DBG_PRINT("\tbac: (%d, %d)\n", chosen_neighbour.first, chosen_neighbour.second);
                         Rmacrophage* rmac = dynamic_cast<Rmacrophage*>(grid[chosen_neighbour.first][chosen_neighbour.second]->getAgent());
                         if (rmac == nullptr) {
                           Neutrophil* neut = dynamic_cast<Neutrophil*>(grid[chosen_neighbour.first][chosen_neighbour.second]->getAgent());
                           UTIMODEL_DBG_MSG("cast to cast to rmac failed");
                          if (neut == nullptr){
                            Mastcell* mast = dynamic_cast<Mastcell*>(grid[chosen_neighbour.first][chosen_neighbour.second]->getAgent());
                             UTIMODEL_DBG_MSG("cast to cast to neut failed");
                            if (mast == nullptr){
                              Hmacrophage* tmp = dynamic_cast<Hmacrophage*>(grid[chosen_neighbour.first][chosen_neighbour.second]->getAgent());
                              UTIMODEL_DBG_MSG("cast to Hmacrophage succeeded");
                            if (tmp == nullptr){
                                UTIMODEL_DBG_MSG("cast to all immune cells has failed");
                              }else{
                                UTIMODEL_DBG_MSG("cast to immune cell has succeeded");
                                UTIMODEL_DBG_PRINT("\tmac: (%d, %d)\n", chosen_neighbour.first, chosen_neighbour.second);
                             }
                     }
                   }
                   }
                   } else {
                         UTIMODEL_DBG_MSG("cast to bacteria has succeeded");
                         UTIMODEL_DBG_PRINT("\tbac: (%d, %d)\n", chosen_neighbour.first, chosen_neighbour.second);
                       }

                      // std::cout << "this loop is useful" << std::endl;
                      // Random number 0-1
                      r = distribution_0_1(rng);

                      /* CHRIS: As we have now used a dynamic_cast above, if the cast fails the Bacterium* bac becomes a nullptr. Therefore,
                       * we need to check for it, otherwise we would attempt to call bac->isReplicating() which would result in a seg fault
                       * due to isReplicating() not being a part of a nullptr.
                       */

                      if (bac != nullptr) { /* CHRIS: safety net in case casting fails */
                        // Determine if bacterium is destroyed based on its replication type

                        //parameter for neutrophils used for macrophages for simplicity TO BE REPLACED
                        if ((bac->isReplicating() && r < p->probability_hmac_destroy_replicating_bacteria) ||
                            (!bac->isReplicating() && r < p->probability_hmac_destroy_dormant_bacteria)) {
                            //std::cout << "/* Code100- */" << '\n';
                            // Remove from grid
                            eraseBacterium(bac);
                            // Remove from list

                            /* CHRIS: similar to Immunecell above, erasing from a Set does not clear the memory. Therefore, we need
                             * to manually delete the bac. Again, note that if we didn't have the safety net above that tested for
                             * for a nulltpr, we could potentially arrive at a seg fault
                             */

                            extracellular_bacteria.erase(bac);
                            delete bac;
                            x = hmac->getCoordinates().first;
                            y = hmac->getCoordinates().second;

                            // Move the macrophage
                            moveAgent(x, y, chosen_neighbour.first, chosen_neighbour.second);
                            // Increment counter
                            counter_bacteria_destroyed_by_hmacrophages++;
                            //std::cout << "/* Code100- */" << '\n';
                        } else {
                              r = distribution_0_1(rng);
                              extracellular_bacteria.erase(bac);
                              eraseBacterium(bac);

                              //std::cout << "/* Code100- */" << '\n';
                              delete bac;
                              bacteria_at_level_plusone++;
                              //std::cout << bacteria_at_level_plusone << '\n';

                        }
                      }
              }
            }

          }
          // Activation

          if (hmac->getState()==hresting){

              // Reobtain co-ordinates (in case mac has moved)
              x = hmac->getCoordinates().first;
              y = hmac->getCoordinates().second;
              // Check immediate Moore neighbourhood
              // TODO: slight performance hit getting a neighbourhood every timestep
              neighbours = grid[x][y]->immediate_moore;
              // Count neighbours that can activate macrophages
              neighbour_count_activate = 0;
              for (auto iter:neighbours){
                  if (grid[iter.first][iter.second]->getContents() == LocationContents::agent){
                      neighbour_count_activate += grid[iter.first][iter.second]->getAgent()->canActivatemacrophage();
                  }
              }
              // Random number 0-1
              r = distribution_0_1(rng);

              if (r < p->probability_hmacrophage_activation && neighbour_count_activate > 0){
                  // Activate the macrophage
                  //std::cout << "wouldnt make sense" << std::endl;
                  hmac->activate();
                  //std::cout << mac->getState() << std::endl;
              }

            } else if (hmac->getState()==hactivated){
              // Reobtain co-ordinates (in case rmac has moved)
              x = hmac->getCoordinates().first;
              y = hmac->getCoordinates().second;
              neighbours2 = grid[x][y] -> immediate_moore;
              space_finder = findSpace(neighbours2);

              if (std::get<0>(space_finder))
              {
                empty_x = std::get<1>(space_finder);
                empty_y = std::get<2>(space_finder);
                neut = addNeutrophil(empty_x,empty_y);
                neutrophils.insert(neut);
                neutrophils_placed++;

              }


              // Reobtain co-ordinates (in case mac has moved)
              x = hmac->getCoordinates().first;
              y = hmac->getCoordinates().second;
              // Check immediate Moore neighbourhood
              // TODO: slight performance hit getting a neighbourhood every timestep
              neighbours = grid[x][y]->immediate_moore;
              // Count neighbours that can activate macrophages
              neighbour_count_activate = 0;
              for (auto iter:neighbours){
                  if (grid[iter.first][iter.second]->getContents() == LocationContents::agent){
                      neighbour_count_activate += grid[iter.first][iter.second]->getAgent()->canActivatemacrophage();
                  }
              }
              // Random number 0-1
              r = distribution_0_1(rng);
              // Multiply the probability by the number of activating cells in neighbour
              if (neighbour_count_activate==0){
                  // deactivate the macrophage
                  // std::cout << "mac is getting deactivated" << std::endl;
                  hmac->deactivate();
              }
            }
          // Increment iterator
          ++iter;
  }// immunecell loop


}

void Environment::MastcellAction(double time){
  int x, y, z, h, k, l, r, neighbour_count_activate,empty_y,empty_x;
  double movement_rate;
  std::vector<std::pair<int, int>> neighbours;
  std::pair<int, int> chosen_neighbour;
  LocationContents neighbour_contents;
  std::tuple<bool,int,int> space_finder;
  Rmacrophage* rmac;
  Hmacrophage *hmac;
  Neutrophil *neut;
  Mastcell *mast;
  int downregulation_count = 0;
  int counter_death = 0;
  int hmacrophages_placed = 0;
  int neutrophils_taken = 0;


  for (auto iter= mastcells.begin(); iter!=mastcells.end();){
    mast = *iter;
    if (mast->exceededLifespan()){
      ++counter_death;
      MastcellDeath(mast);
      iter = mastcells.erase(iter);
      delete mast;
      continue;
    }
    switch (mast->getState()) {
      case mresting:
        movement_rate = p->restingmastcell_movement_rate;
        break;
      case mactivated:
        movement_rate = p->activatedmastcell_movement_rate;
        break;
      default:
      throw std::runtime_error("invalid Mastcell state");
    }
    if ((int)round(time/p->timestep)%(int)round(movement_rate/p->timestep)==0){
      x = mast -> getCoordinates().first;
      y = mast -> getCoordinates().second;
      // Get neighbours
      neighbours = grid[x][y]->immediate_moore;

      // Chose a neighbour to move to (chemotactically or randomly)
      chosen_neighbour = randomMovement(neighbours);

      // Check the contents of the neighbour
      neighbour_contents = grid[chosen_neighbour.first][chosen_neighbour.second]->getContents();
      //std::cout << "###" << static_cast<int>(neighbour_contents) << std::endl;
      // Location is empty, move the macrophage
      if (neighbour_contents == LocationContents::empty){
          //std::cout << "Macrophage moves" << std::endl;
          moveAgent(x, y, chosen_neighbour.first, chosen_neighbour.second);
      // Location has agent (i.e. a bacterium)
    }
  }

  // Activation
  if (mast->getState()==mresting){
      // Reobtain co-ordinates (in case mac has moved)
      x = mast->getCoordinates().first;
      y = mast->getCoordinates().second;
      // Check immediate Moore neighbourhood
      // TODO: slight performance hit getting a neighbourhood every timestep
      neighbours = grid[x][y]->immediate_moore;
      // Count neighbours that can activate macrophages
      neighbour_count_activate = 0;
      for (auto iter:neighbours){
          if (grid[iter.first][iter.second]->getContents() == LocationContents::agent){
              neighbour_count_activate += grid[iter.first][iter.second]->getAgent()->canActivatemacrophage();
              downregulation_count += grid[iter.first][iter.second]->getAgent()->canActivateDRneutrophils();

          }
      }
      // Random number 0-1
      r = distribution_0_1(rng);

      if (r < p->probability_rmacrophage_activation && neighbour_count_activate > 0){
          // Activate the macrophage
          //std::cout << "wouldnt make sense" << std::endl;
          mast->activate();
          //std::cout << mac->getState() << std::endl;


      }
    } else if (mast->getState()==mactivated){
      //neighbour_count_activate = 0;
      //downregulation_count = 0;
      // Reobtain co-ordinates (in case rmac has moved)
      if (neighbour_count_activate < 2){
        x = mast->getCoordinates().first;
        y = mast->getCoordinates().second;
        neighbours = grid[x][y] -> immediate_moore;
        space_finder = findSpace(neighbours);


        if (std::get<0>(space_finder)){
          empty_x = std::get<1>(space_finder);
          empty_y = std::get<2>(space_finder);
          hmac = addHmacrophage(empty_x,empty_y);
          hmacrophages.insert(hmac);
          hmacrophages_placed++;
          //std::cout << "This part works too" << '\n';
        }
      }


      if (downregulation_count > 2){
        x = mast->getCoordinates().first;
        y = mast->getCoordinates().second;
        neighbours = grid[x][y]->immediate_moore;
        for (auto iter:neighbours){
          if (grid[iter.first][iter.second]->getContents() == LocationContents::agent && grid[iter.first][iter.second]->getAgent()->canActivateDRneutrophils()){
              neut = dynamic_cast<Neutrophil*>(grid[iter.first][iter.second]->getAgent());
              if (neut == nullptr) {
                UTIMODEL_DBG_MSG("cast to neutrophil has failed");
                UTIMODEL_DBG_PRINT("\tneut: (%d, %d)\n", chosen_neighbour.first, chosen_neighbour.second);
              }
              if (neut != nullptr){
                NeutrophilDeath(neut);
                neutrophils.erase(neut); // removes the agent from the set of neutrophils
                delete neut; // deletes the agent from memory
              }
        }

        }
      }
      // Check immediate Moore neighbourhood
      // TODO: slight performance hit getting a neighbourhood every timestep
      neighbours = grid[x][y]->immediate_moore;
      // Count neighbours that can activate macrophages
      neighbour_count_activate = 0;
      for (auto iter:neighbours){
          if (grid[iter.first][iter.second]->getContents() == LocationContents::agent){
             neighbour_count_activate += grid[iter.first][iter.second]->getAgent()->canActivatemacrophage();
          }
      }
      // Random number 0-1
      r = distribution_0_1(rng);
      // Multiply the probability by the number of activating cells in neighbour
      if (neighbour_count_activate==0){
          // deactivate the macrophage
          // std::cout << "mac is getting deactivated" << std::endl;
          mast->deactivate();
      }
    }
  // Increment iterator
  ++iter;
}// immunecell loop
}


void Environment::RmacrophageAction(double time){
  int x, y, z, h, k, l, r, neighbour_count_activate,empty_y,empty_x;
  double movement_rate, movement_bias;
  std::vector<std::pair<int, int>> neighbours, neighbours2;
  std::pair<int, int> chosen_neighbour;
  LocationContents neighbour_contents;
  std::tuple<bool,int,int> space_finder;
  std::tuple<bool,int,int> space_finder2;
  Rmacrophage* rmac;
  Hmacrophage *hmac;
  Neutrophil *neut;
  int counter_death = 0;
  int hmacrophages_placed = 0;
  int neutrophils_placed = 0;

  for (auto iter= rmacrophages.begin(); iter!=rmacrophages.end();){
    rmac = *iter;
    if (rmac->exceededLifespan()){
      ++counter_death;
      RmacrophageDeath(rmac);
      iter = rmacrophages.erase(iter);
      delete rmac;
      continue;
    }
    switch (rmac->getState()) {
      case resting:
        movement_rate = p->restingresident_macrophage_movement_rate;
        movement_bias = p->resting_rmacrophage_movement_bias;
        break;
      case activated:
        movement_rate = p->activatedresident_macrophage_movement_rate;
        movement_bias = p->activated_rmacrophage_movement_bias;
        break;
      default:
      throw std::runtime_error("invalid Macrophage state");
    }
    if ((int)round(time/p->timestep)%(int)round(movement_rate/p->timestep)==0){
      x = rmac -> getCoordinates().first;
      y = rmac -> getCoordinates().second;
      // Get neighbours
      neighbours = grid[x][y]->immediate_moore;

      // Chose a neighbour to move to (chemotactically or randomly)
      chosen_neighbour = chemotacticMovement(neighbours, movement_bias);

      // Check the contents of the neighbour
      neighbour_contents = grid[chosen_neighbour.first][chosen_neighbour.second]->getContents();
      //std::cout << "###" << static_cast<int>(neighbour_contents) << std::endl;
      // Location is empty, move the macrophage
      if (neighbour_contents == LocationContents::empty){
          //std::cout << "Macrophage moves" << std::endl;
          moveAgent(x, y, chosen_neighbour.first, chosen_neighbour.second);
      // Location has agent (i.e. a bacterium)
    }
  }

  // Activation
  if (rmac->getState()==resting){
      // Reobtain co-ordinates (in case mac has moved)
      x = rmac->getCoordinates().first;
      y = rmac->getCoordinates().second;
      // Check immediate Moore neighbourhood
      // TODO: slight performance hit getting a neighbourhood every timestep
      neighbours = grid[x][y]->immediate_moore;
      // Count neighbours that can activate macrophages
      neighbour_count_activate = 0;
      for (auto iter:neighbours){
          if (grid[iter.first][iter.second]->getContents() == LocationContents::agent){
              neighbour_count_activate += grid[iter.first][iter.second]->getAgent()->canActivatemacrophage();
          }
      }
      // Random number 0-1
      r = distribution_0_1(rng);

      if (r < p->probability_rmacrophage_activation && neighbour_count_activate > 0){
          // Activate the macrophage
          //std::cout << "wouldnt make sense" << std::endl;
          rmac->activate();
          //std::cout << mac->getState() << std::endl;


      }
    } else if (rmac->getState()==activated){
      // Reobtain co-ordinates (in case rmac has moved)
      x = rmac->getCoordinates().first;
      y = rmac->getCoordinates().second;
      neighbours2 = grid[x][y] -> immediate_moore;
      space_finder = findSpace(neighbours2);

      if (std::get<0>(space_finder))
      {
        empty_x = std::get<1>(space_finder);
        empty_y = std::get<2>(space_finder);
        hmac = addHmacrophage(empty_x,empty_y);
        hmacrophages.insert(hmac);
        hmacrophages_placed++;

        x = rmac->getCoordinates().first;
        y = rmac->getCoordinates().second;
        neighbours = grid[x][y] -> immediate_moore;
        space_finder2 = findSpace(neighbours);
        if (std::get<0>(space_finder2))
        {
          empty_x = std::get<1>(space_finder2);
          empty_y = std::get<2>(space_finder2);
          neut = addNeutrophil(empty_x,empty_y);
          neutrophils.insert(neut);
          neutrophils_placed++;
        }
      }


      // Check immediate Moore neighbourhood
      // TODO: slight performance hit getting a neighbourhood every timestep
      neighbours = grid[x][y]->immediate_moore;
      // Count neighbours that can activate macrophages
      neighbour_count_activate = 0;
      for (auto iter:neighbours){
          if (grid[iter.first][iter.second]->getContents() == LocationContents::agent){
              neighbour_count_activate += grid[iter.first][iter.second]->getAgent()->canActivatemacrophage();
          }
      }
      // Random number 0-1
      r = distribution_0_1(rng);
      // Multiply the probability by the number of activating cells in neighbour
      if (neighbour_count_activate==0){
          // deactivate the macrophage
          // std::cout << "mac is getting deactivated" << std::endl;
          rmac->deactivate();
      }
    }
  // Increment iterator
  ++iter;
}// immunecell loop
}

void Environment::bacterialStateChange(double time)
{
  (void) time;
    int x, y, neighbourhood_depth;
    std::vector<std::pair<int, int>> neighbours;
    std::tuple<bool, int, int> space_finder;
    for (auto bac : extracellular_bacteria)
    {
        x = bac->getCoordinates().first;
        y = bac->getCoordinates().second;
        // Resting state change
        if (bac->isResting())
        {
            // Start with immediate neighbours
            neighbourhood_depth = 1;
            // Keep looking at successive depths until a space found or max depth reached
            while (bac->isResting() && neighbourhood_depth <= p->bacteria_replication_neighbourhood_max_depth)
            {
                // Get neighbours
                neighbours = mooreNeighbourhood(x, y, neighbourhood_depth, grid_size, rng);
                // Look for space (only interested in presence of space, not its actual coordinates)
                space_finder = findSpace(neighbours);
                // Space found, so switch to non-resting
                if (std::get<0>(space_finder))
                {
                    bac->switchResting();
                    // No space here, so look deeper
                }
                else
                {
                    neighbourhood_depth++;
                } // If space found
            }     // while loop looking for space
            // Replication state change
        }

    }
}

void Environment::processTimestep(double time)
{
  // Update the cellular automaton attributes
    updateAttributes(time);
  for (auto iter : extracellular_bacteria)
  {
    iter->increaseAge(p->timestep);
  }
  for (auto iter : rmacrophages)
  {
    iter->increaseAge(p->timestep);
  }
  for (auto iter : hmacrophages)
  {
    iter->increaseAge(p->timestep);
  }
  for (auto iter : neutrophils)
  {
    iter->increaseAge(p->timestep);
  }
  for (auto iter : mastcells)
  {
    iter->increaseAge(p->timestep);
  }

  try
  {
    bacterialReplication(time);
    RmacrophageRecruitment();
    HmacrophageRecruitment();
    NeutrophilRecruitment();
    MastcellRecruitment();
    RmacrophageAction(time);
    //BUG
    HmacrophageAction(time);
    MastcellAction(time);

    NeutrophilAction(time);
    bacterialStateChange(time);
    bacterialShedding();
  }
  catch (const char *msg)
  {
    throw msg;
  }

}

//void Environment::outputVesselsToFile(FILE *f)
//{
//    for (int i = 0; i < grid_size; i++)
//        for (int j = 0; j < grid_size; j++)
//        {
//            if (grid[i][j]->getContents() == LocationContents::vessel)
//            {
//                fprintf(f, "%f \n", 9.0);
//            }
//            else
//            {
//                fprintf(f, "%f \n", 0.0);
//            }
//        }
//}

void Environment::recordCellNumbersToFiles(outputFiles files){
  Location* loc;
  double cell_code;
  int total_num_of_cells;
  int Type1_count = 0;
  int Type1_R_count = 0;
  int Type2_count = 0;
  int Type2_R_count = 0;
  int Type4_count = 0;
  int Type5_count = 0;
  int Type6_count = 0;
  int Type7_count = 0;
  int Type8_count = 0;
  int Type9_count = 0;
  int Type12_count = 0;
  int Type13_count = 0;
  // int Type6_count = 0;
  // int Type7_count = 0;
  // int total_intra = 0;
  // Total number of cells
  total_num_of_cells = (int)hmacrophages.size() + (int)rmacrophages.size() + (int)neutrophils.size() + (int)mastcells.size() + (int)extracellular_bacteria.size();
  // Bacteria counts
  for (auto iter: extracellular_bacteria) {
      if (iter->isReplicating()){
          Type1_count+=1;
          if (iter->isResting()){
              Type1_R_count+=1;
          }
      } else {
          Type2_count+=1;
          if (iter->isResting()){
              Type2_R_count+=1;
          }
      }
  }

  for (auto iter:mastcells){
    switch (iter->getState()){
      case mresting:
        Type12_count++;
        break;
      case mactivated:
        Type13_count++;
        break;
    }
  }

  for (auto iter:rmacrophages){
    switch (iter->getState()) {
      case resting:
        Type4_count++;
        break;
      case activated:
        Type5_count++;
        break;
    }
  }
  for (auto iter:hmacrophages){
    switch (iter->getState()) {
      case hresting:
        Type6_count++;
        break;
      case hactivated:
        Type7_count++;
        break;
    }
  }
  for (auto iter:neutrophils){
    switch (iter->getState()) {
      case nresting:
        Type8_count++;
        break;
      case nactivated:
        Type9_count++;
        break;
    }
  }
  // All bacteria (alive and destroyed
  int total_bac = (int)extracellular_bacteria.size(); //+ counter_bacteria_destroyed_by_hmacrophages + counter_bacteria_destroyed_by_neutrophils;

  //fprintf(files.f1, "%i\n", total_num_of_cells);
  //fprintf(files.f2,"%d\n", Type1_count);  // replicating
  //fprintf(files.f3,"%d\n",Type1_R_count);
  //fprintf(files.f4,"%d\n",Type2_count);  // non-replicating
  //fprintf(files.f5,"%d\n",Type2_R_count);
  //fprintf(files.fmr,"%d\n",Type4_count); // immune resting cell
  //fprintf(files.fma,"%d\n",Type5_count); // immune active cell
  //fprintf(files.fhrma,"%d\n",Type6_count);  // non-replicating
  //fprintf(files.fhama,"%d\n",Type7_count);
  //fprintf(files.frneu,"%d\n",Type8_count); // immune resting cell
  //fprintf(files.faneu,"%d\n",Type9_count); // immune active cell
  //fprintf(files.frmc,"%d\n",Type12_count); // immune active cell
  //fprintf(files.famc,"%d\n",Type13_count); // immune active cell
  //fprintf(files.fimmunekill,"%d\n",counter_bacteria_destroyed_by_hmacrophages);

  //fprintf(files.ftotalcells,"%d\n",total_bac);

  fprintf(files.analysis1,"%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %d\n",p->activatedmastcell_movement_rate,p->mastcells_recruitment_prob_no_bac,
                                                    p->mastcells_recruitment_prob_bac,
                                                    p->restingmastcell_movement_rate,p->resting_hmacrophage_movement_rate,p->restingresident_macrophage_movement_rate,
                                                    p->activated_hmacrophage_movement_rate,p->hmacrophages_recruitment_prob_no_bac,p->rmacrophages_recruitment_prob_no_bac,
                                                    p->hmacrophages_recruitment_prob_bac,p->rmacrophages_recruitment_prob_bac,p->activatedresident_macrophage_movement_rate,
                                                    p->resting_neutrophil_movement_rate,p->activated_neutrophil_movement_rate,p->probability_neutrophil_destroy_replicating_bacteria,
                                                    p->probability_neutrophil_destroy_dormant_bacteria, p->probability_hmac_destroy_replicating_bacteria,
                                                    p->probability_hmac_destroy_dormant_bacteria, p->probability_neutrophil_activation,p->neutrophil_recruitment_prob_bac,
                                                    p->neutrophil_recruitment_prob_no_bac,p->probability_rmacrophage_activation,p->probability_hmacrophage_activation,
                                                    total_bac);


  fprintf(files.analysis2,"%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %d %d\n",p->activatedmastcell_movement_rate,p->mastcells_recruitment_prob_no_bac,
                                                    p->mastcells_recruitment_prob_bac,
                                                    p->restingmastcell_movement_rate,p->resting_hmacrophage_movement_rate,p->restingresident_macrophage_movement_rate,
                                                    p->activated_hmacrophage_movement_rate,p->hmacrophages_recruitment_prob_no_bac,p->rmacrophages_recruitment_prob_no_bac,
                                                    p->hmacrophages_recruitment_prob_bac,p->rmacrophages_recruitment_prob_bac,p->activatedresident_macrophage_movement_rate,
                                                    p->resting_neutrophil_movement_rate,p->activated_neutrophil_movement_rate,p->probability_neutrophil_destroy_replicating_bacteria,
                                                    p->probability_neutrophil_destroy_dormant_bacteria, p->probability_hmac_destroy_replicating_bacteria,
                                                    p->probability_hmac_destroy_dormant_bacteria, p->probability_neutrophil_activation,p->neutrophil_recruitment_prob_bac,
                                                    p->neutrophil_recruitment_prob_no_bac,p->probability_rmacrophage_activation,p->probability_hmacrophage_activation,
                                                    counter_bacteria_destroyed_by_neutrophils,counter_bacteria_destroyed_by_hmacrophages);


  fprintf(files.analysis3,"%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %d\n",p->activatedmastcell_movement_rate,p->mastcells_recruitment_prob_no_bac,
                                                    p->mastcells_recruitment_prob_bac,
                                                    p->restingmastcell_movement_rate,p->resting_hmacrophage_movement_rate,p->restingresident_macrophage_movement_rate,
                                                    p->activated_hmacrophage_movement_rate,p->hmacrophages_recruitment_prob_no_bac,p->rmacrophages_recruitment_prob_no_bac,
                                                    p->hmacrophages_recruitment_prob_bac,p->rmacrophages_recruitment_prob_bac,p->activatedresident_macrophage_movement_rate,
                                                    p->resting_neutrophil_movement_rate,p->activated_neutrophil_movement_rate,p->probability_neutrophil_destroy_replicating_bacteria,
                                                    p->probability_neutrophil_destroy_dormant_bacteria, p->probability_hmac_destroy_replicating_bacteria,
                                                    p->probability_hmac_destroy_dormant_bacteria, p->probability_neutrophil_activation,p->neutrophil_recruitment_prob_bac,
                                                    p->neutrophil_recruitment_prob_no_bac,p->probability_rmacrophage_activation,p->probability_hmacrophage_activation,
                                                    Type1_count + Type2_count);


  fprintf(files.analysis4,"%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %d\n",p->activatedmastcell_movement_rate,p->mastcells_recruitment_prob_no_bac,
                                                    p->mastcells_recruitment_prob_bac,
                                                    p->restingmastcell_movement_rate,p->resting_hmacrophage_movement_rate,p->restingresident_macrophage_movement_rate,
                                                    p->activated_hmacrophage_movement_rate,p->hmacrophages_recruitment_prob_no_bac,p->rmacrophages_recruitment_prob_no_bac,
                                                    p->hmacrophages_recruitment_prob_bac,p->rmacrophages_recruitment_prob_bac,p->activatedresident_macrophage_movement_rate,
                                                    p->resting_neutrophil_movement_rate,p->activated_neutrophil_movement_rate,p->probability_neutrophil_destroy_replicating_bacteria,
                                                    p->probability_neutrophil_destroy_dormant_bacteria, p->probability_hmac_destroy_replicating_bacteria,
                                                    p->probability_hmac_destroy_dormant_bacteria, p->probability_neutrophil_activation,p->neutrophil_recruitment_prob_bac,
                                                    p->neutrophil_recruitment_prob_no_bac,p->probability_rmacrophage_activation,p->probability_hmacrophage_activation,
                                                    Type6_count + Type7_count);

  fprintf(files.analysis5,"%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %d\n",p->activatedmastcell_movement_rate,p->mastcells_recruitment_prob_no_bac,
                                                    p->mastcells_recruitment_prob_bac,
                                                    p->restingmastcell_movement_rate,p->resting_hmacrophage_movement_rate,p->restingresident_macrophage_movement_rate,
                                                    p->activated_hmacrophage_movement_rate,p->hmacrophages_recruitment_prob_no_bac,p->rmacrophages_recruitment_prob_no_bac,
                                                    p->hmacrophages_recruitment_prob_bac,p->rmacrophages_recruitment_prob_bac,p->activatedresident_macrophage_movement_rate,
                                                    p->resting_neutrophil_movement_rate,p->activated_neutrophil_movement_rate,p->probability_neutrophil_destroy_replicating_bacteria,
                                                    p->probability_neutrophil_destroy_dormant_bacteria, p->probability_hmac_destroy_replicating_bacteria,
                                                    p->probability_hmac_destroy_dormant_bacteria, p->probability_neutrophil_activation,p->neutrophil_recruitment_prob_bac,
                                                    p->neutrophil_recruitment_prob_no_bac,p->probability_rmacrophage_activation,p->probability_hmacrophage_activation,
                                                    Type4_count + Type5_count);


  fprintf(files.analysis6,"%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %d\n",p->activatedmastcell_movement_rate,p->mastcells_recruitment_prob_no_bac,
                                                    p->mastcells_recruitment_prob_bac,
                                                    p->restingmastcell_movement_rate,p->resting_hmacrophage_movement_rate,p->restingresident_macrophage_movement_rate,
                                                    p->activated_hmacrophage_movement_rate,p->hmacrophages_recruitment_prob_no_bac,p->rmacrophages_recruitment_prob_no_bac,
                                                    p->hmacrophages_recruitment_prob_bac,p->rmacrophages_recruitment_prob_bac,p->activatedresident_macrophage_movement_rate,
                                                    p->resting_neutrophil_movement_rate,p->activated_neutrophil_movement_rate,p->probability_neutrophil_destroy_replicating_bacteria,
                                                    p->probability_neutrophil_destroy_dormant_bacteria, p->probability_hmac_destroy_replicating_bacteria,
                                                    p->probability_hmac_destroy_dormant_bacteria, p->probability_neutrophil_activation,p->neutrophil_recruitment_prob_bac,
                                                    p->neutrophil_recruitment_prob_no_bac,p->probability_rmacrophage_activation,p->probability_hmacrophage_activation,
                                                    Type8_count + Type9_count);

  fprintf(files.analysis7,"%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %d\n",p->activatedmastcell_movement_rate,p->mastcells_recruitment_prob_no_bac,
                                                    p->mastcells_recruitment_prob_bac,
                                                    p->restingmastcell_movement_rate,p->resting_hmacrophage_movement_rate,p->restingresident_macrophage_movement_rate,
                                                    p->activated_hmacrophage_movement_rate,p->hmacrophages_recruitment_prob_no_bac,p->rmacrophages_recruitment_prob_no_bac,
                                                    p->hmacrophages_recruitment_prob_bac,p->rmacrophages_recruitment_prob_bac,p->activatedresident_macrophage_movement_rate,
                                                    p->resting_neutrophil_movement_rate,p->activated_neutrophil_movement_rate,p->probability_neutrophil_destroy_replicating_bacteria,
                                                    p->probability_neutrophil_destroy_dormant_bacteria, p->probability_hmac_destroy_replicating_bacteria,
                                                    p->probability_hmac_destroy_dormant_bacteria, p->probability_neutrophil_activation,p->neutrophil_recruitment_prob_bac,
                                                    p->neutrophil_recruitment_prob_no_bac,p->probability_rmacrophage_activation,p->probability_hmacrophage_activation,
                                                    Type12_count + Type13_count);


  for (int x=0; x<grid_size; x++){
      for (int y=0; y<grid_size; y++){
          loc = grid[x][y];
          switch (loc->getContents()) {
              case LocationContents::empty:
                  cell_code = 0.0;
                  break;
              case LocationContents::vessel:
                  cell_code = 11.0;
                  break;
              case LocationContents::agent:
                  cell_code = loc->getAgent()->getCode();
                  break;
              default:
                  break;
          }
          // Write contents
          fprintf(files.fp, "%f\n", cell_code);
          fprintf(files.analysis8,"%f\n", loc->scaled_chemokine_value);
      }
  }

}
