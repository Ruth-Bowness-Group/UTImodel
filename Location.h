#ifndef Location_h
#define Location_h

#include <map>
#include <stdio.h>
#include "Agent.h"
#include "Parameters.h"

// Enum of the possible contents within a location
enum class LocationContents {empty, vessel, agent};
enum class EnvironmentalAttribute {chemokine};

class Location
{
private:
  Parameters* p;
  int x;
  int y;
  LocationContents contents;
  Agent* agent_ref;

public:
  Location(int x, int y, Parameters* p);
  ~Location (); /* CHRIS: Location destructor */
  // Immediate Moore neighbourhood
  std::vector<std::pair<int, int> > immediate_moore;
  LocationContents getContents();
  // Level of chemokine
  double chemokine_value = 0.0;
  double scaled_chemokine_value = 0.0;
  // Rate of chemokine source
  double chemokine_source_at_location = 0.0;
  void addAgent(Agent* agent_to_add);
  Agent* removeAgent();
  void setVessel();
  Agent* getAgent();
  std::string getAgentCode();
};

#endif
