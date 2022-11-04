#ifndef Environment_h
#define Environment_h

#include "Parameters.h"
#include <iostream>
#include <random>
#include <set>
#include <algorithm>
#include "Location.h"
#include "Output.h"

/**
 * General functions
 */
std::vector<std::pair<int, int> > mooreNeighbourhood(int x, int y, int level, int grid_size, std::default_random_engine rng);
std::vector<std::pair<int, int> > vonNeumannNeighbourhood(int x, int y, int level, int grid_size, std::default_random_engine rng);
double diffusion_centre(double value_at_location, double value_N, double value_E, double value_S, double value_W, double diffusion_at_location, double diffusion_N, double diffusion_E, double diffusion_S, double diffusion_W, double spatial_step);
double diffusion_boundary(double value_at_location, double value_N1, double value_N2, double value_N3, double value_N4, double diffusion_at_location, double spatial_step);


class Environment
{
private:
  Parameters* p;
  Location *grid[grid_size][grid_size];
  std::vector<std::pair<int,int> > vessels;
  // Lists of agents
  std::set<Bacterium *> extracellular_bacteria;
  std::set<Rmacrophage *> rmacrophages;
  std::set<Hmacrophage *> hmacrophages;
  std::set<Neutrophil *> neutrophils;
  std::set<Mastcell *> mastcells;


  //Counts
  int counter_bacteria_destroyed_by_hmacrophages = 0;
  int counter_bacteria_destroyed_by_neutrophils = 0;
  int counter_epithelial_level_shedding = 0;
  int bacteria_at_level_plusone = 0;

  std::tuple<bool, int, int> findSpace(std::vector<std::pair<int, int> > neighbours);
  // Agent functions
  void moveAgent(int x_from, int y_from, int x_to, int y_to);
  Bacterium *addBacterium(int x, int y, bool replicating);
  void eraseBacterium(Bacterium *bac);
  Rmacrophage *addRmacrophage(int x, int y);
  Hmacrophage *addHmacrophage(int x, int y);
  Neutrophil *addNeutrophil(int x, int y);
  Mastcell *addMastcell(int x, int y);
  void MastcellDeath(Mastcell *mast);
  void HmacrophageDeath(Hmacrophage *hmac);
  void RmacrophageDeath(Rmacrophage *rmac);
  void NeutrophilDeath(Neutrophil *neut);
  // ABM events
  void bacterialReplication(double time);
  void RmacrophageRecruitment();
  void HmacrophageRecruitment();
  void NeutrophilRecruitment();
  void MastcellRecruitment();
  void RmacrophageAction(double time);
  void HmacrophageAction(double time);
  void NeutrophilAction(double time);
  void MastcellAction(double time);
  void bacterialStateChange(double time);
  void bacterialShedding();
  // Diffusion
  double diffusion(EnvironmentalAttribute attribute, int x, int y);
  std::pair<int, int> randomMovement(std::vector<std::pair<int, int> > neighbours);
  std::pair<int, int> chemotacticMovement(std::vector<std::pair<int, int>> neighbours, double movement_bias);

public:
  Environment(Parameters* param);
  ~Environment(); /* CHRIS: destructor for the Environment class */
  Location *getLocation(int x, int y);
  void initialiseBacterialCluster(int x, int y);
  void initialiseBacterialClusterpostShedding(int x, int y);
  void initialiseRmacrophage(int number_of_rmacrophages);
  void initialiseMastcell(int number_of_mastcell);
  void initialiseBloodVesselsRandom(int ves);
  void initialiseBloodVesselsFixed(std::vector<std::pair<int, int> > locations);
  void setAVessel(int x, int y);
  // ABM Events
  // Cellular automaton updates
  void updateAttributes(double time);
  void processTimestep(double time);
  void outputVesselsToFile(FILE *f);
  void recordCellNumbersToFiles(outputFiles files);
};

#endif
