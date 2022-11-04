#ifndef Agent_h
#define Agent_h
#include <string>
#include "Parameters.h"

/**
 An agent within the model
 */
class Agent
{
private:
  double age = 0;
  int x;
  int y;
  // Agent can help activate a macrophage if in proximity
  bool can_activate_macrophage = false;
  bool can_activate_urhmacrophage = false;
  bool can_activate_drneutrophil = false;
  bool phagocytosable = false;
public:
  Parameters * p;
  int id;
  Agent(Parameters* param);
  /* CHRIS: agent destructor is virtual so that whenever a derived class, such as
   * Bacterium or Rmacrophage, is created then they must provide their own destructor.
   * Also, it ensures that whenever Bacterium or Rmacrophage is deleted then the Agent
   * destructor is also called.
   */
  virtual ~Agent ();
  double getAge();
  void setAge(double amount);
  void increaseAge(double amount);
  // Coordinates
  void setCoordinates(int x_coord, int y_coord);
  std::pair<int, int> getCoordinates();
  // Virtual code function - each subclass must specify its own code for output
  virtual double getCode();
  // Virtual code function - each subclass must specify how much chemokine it releases
  virtual double getChemokineSource();
  bool canActivateURHmacrophages();
  bool canActivateDRneutrophils();
  void setcanActivateDRneutrophils();
  bool canActivatemacrophage();
  void setcanActivateURHmacrophages();
  void setCanActivatemacrophage();
  void setPhagocytosable();
  bool isPhagocytosable();
};

/**
 * Neutrophils
 */
enum NeutrophilState
{
nresting,
nactivated
};

class Neutrophil : public Agent
{
private:
double lifespan;
// Default state is resting
NeutrophilState state = nresting;
public:
Neutrophil(double lifespan, Parameters* p);
~Neutrophil (); /* CHRIS: destructor for Hmacrophage, which is called whenever Hmacrophage is deleted */

double getCode();
NeutrophilState getState();
bool exceededLifespan();
void activate();
double getChemokineSource();
void deactivate();
};


/**
 * Mastcell
 */

// Macrophages States
enum MastcellState
{
mresting,
mactivated
};

class Mastcell : public Agent
{
private:
double lifespan;
// Default state is resting
MastcellState state = mresting;
public:
Mastcell(double lifespan, Parameters* p);
~Mastcell (); /* CHRIS: destructor for Hmacrophage, which is called whenever Hmacrophage is deleted */

double getCode();
MastcellState getState();
bool exceededLifespan();
double getChemokineSource();
void activate();
void deactivate();
};




/**
 * Hmacrophage
 */

// Macrophages States
enum HmacrophageState
{
hresting,
hactivated
};

class Hmacrophage : public Agent
{
private:
double lifespan;
// Default state is resting
HmacrophageState state = hresting;
public:
Hmacrophage(double lifespan, Parameters* p);
~Hmacrophage (); /* CHRIS: destructor for Hmacrophage, which is called whenever Hmacrophage is deleted */
double getChemokineSource();
double getCode();
HmacrophageState getState();
bool exceededLifespan();
void activate();
void deactivate();
};



                              /**
                               * Rmacrophage
                               */

// Macrophages States
enum RmacrophageState
{
resting,
activated
};

class Rmacrophage : public Agent
{
private:
double lifespan;
// Default state is resting
RmacrophageState state = resting;
public:
Rmacrophage(double lifespan, Parameters* p);
~Rmacrophage (); /* CHRIS: destructor for Rmacrophage, which is called whenever Rmacrophage is deleted */

double getCode();
RmacrophageState getState();
double getChemokineSource();
bool exceededLifespan();
void activate();
void deactivate();
};

                          /**
                           Bacterium
                           */
class Bacterium : public Agent
{
private:
   bool replicating;
   bool resting = false;
   bool moore_neighbourhood = true;

public:
   Bacterium(bool replicating, Parameters* p);
   ~Bacterium (); /* CHRIS: destructor for Bacterium, which is called whenever Bacterium is deleted */

   double getCode();
   bool isReplicating();
   void switchReplicating();
   bool isResting();
   double getChemokineSource();
   void switchResting();
   bool isMooreReplicationNeighbourhood();
   void switchReplicationNeighbourhood();

};

#endif
