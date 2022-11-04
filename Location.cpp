#include <map>
#include "Location.h"

/**
* Location
*/


Location::Location(int x_coord, int y_coord, Parameters* param): p(param)
{
  x = x_coord;
  y = y_coord;
  // Default: no contents, not a blood vessel, values = 0
  contents = LocationContents::empty;

  agent_ref = nullptr; /* CHRIS: generally a good idea to initialise variables, even if it is to a nullptr */
}

/* CHRIS: destructor for the Location class. This is called from the Environment destructor when the grid
 * is cleared */

Location::~Location ()
{
  if (p != nullptr)
    p = nullptr;

  if (agent_ref != nullptr)
    agent_ref = nullptr;
}

// Contents functions
LocationContents Location::getContents(){return contents;}


/**
* Set the location as containing blood vessel
*/

void Location::setVessel()
{
  if (contents != LocationContents::empty)
  {
    throw std::runtime_error("Vessel not added. Location not empty.");
  }
  // Set vessel as content of location
  contents = LocationContents::vessel;
  agent_ref = nullptr;
}
/**
* Adding an agent
*/

void Location::addAgent(Agent* agent_to_add)
{
  if (contents != LocationContents::empty)
  {
    throw std::runtime_error("addAgent: Location (" + std::to_string(x) + "," + std::to_string(y) + ") not empty");
  }
  else if (agent_to_add == nullptr)
  {
    throw std::runtime_error("addAgent: null agent");
  }
  contents = LocationContents::agent;

  agent_ref = agent_to_add;

  agent_to_add -> setCoordinates(x,y);
  // Change the chemokine source value here
  chemokine_source_at_location = agent_ref->getChemokineSource();
}

Agent* Location::removeAgent()
{
  if (contents != LocationContents::agent)
  {
    throw std::runtime_error("Cannot remove agent: no agent exists");
  }
  Agent* temp_agent = agent_ref;

  contents = LocationContents::empty;

  chemokine_source_at_location = 0.0;
  agent_ref = nullptr;

  return temp_agent;
}


/**
* Function that gets the content of the location
*/


Agent* Location::getAgent(){return agent_ref;}
