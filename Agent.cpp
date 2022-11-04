#include <iostream>
#include "Parameters.h"
#include "Agent.h"

#include "UTImodel_debug.h" /* CHRIS: my header file for formatted debug printing */

// Agent ID. Incremented when an agent is created, so each agent has unique ID
static int agent_id_max = 0;

/**
 * Agent
 */

                       /**
                        * General Agent Characteristics
                        */

// Default age is zero
Agent::Agent(Parameters* param): p(param){
 age = 0;
 id = agent_id_max;
 agent_id_max++;
}

/* CHRIS: agent destructor which will be called whenever a Bacterium or Rmacrophage
 * is deleted (due to being virtual)
 */

Agent::~Agent ()
{
    if (p != nullptr)
        p = nullptr;
}


// Is Phagocytosable
void Agent::setPhagocytosable(){phagocytosable=true;}
bool Agent::isPhagocytosable(){return phagocytosable;}
//=========================================================================

// Can activate macrophages
void Agent::setCanActivatemacrophage(){can_activate_macrophage=true;}
bool Agent::canActivatemacrophage(){return can_activate_macrophage;}

void Agent::setcanActivateURHmacrophages(){can_activate_urhmacrophage=true;}
bool Agent::canActivateURHmacrophages(){return can_activate_urhmacrophage;}

void Agent::setcanActivateDRneutrophils(){can_activate_drneutrophil=true;}
bool Agent::canActivateDRneutrophils(){return can_activate_drneutrophil;}

// Set the internal coordinates
void Agent::setCoordinates(int x_coord, int y_coord){
  x = x_coord;
  y = y_coord;
}

//Retrieve Coordinates of Agent
std::pair<int, int> Agent::getCoordinates(){return std::make_pair(x, y);}

//Return Age of Agent
double Agent::getAge(){return age;}
//Set Certain Age to Agent
void Agent::setAge(double amount){age = amount;}

//Increase Agent Age by Amount
void Agent::increaseAge(double amount){age += amount;}
double Agent::getChemokineSource(){return 0.0;}

double Agent::getCode(){return -1;}


/**
                             * Neutrophils
                             */

Neutrophil::Neutrophil(double set_lifespan, Parameters* p):Agent(p){
  lifespan = set_lifespan;
  state = nresting;
  setcanActivateDRneutrophils();
}

Neutrophil::~Neutrophil ()
{

}

double Neutrophil::getCode(){
    switch (state) {
        case nresting:
            return 7.0;
        case nactivated:
            return 8.0;
        default:
            return 0.0;
  }
}

NeutrophilState Neutrophil::getState(){return state;}
double Neutrophil::getChemokineSource(){return 0.0;}

bool Neutrophil::exceededLifespan(){
    // If activated, die once exceed the set life span
    if (state==nactivated){
        if (getAge()>=p->activated_neutrophil_lifespan) {
            /* CHRIS: a few debug prints */

            UTIMODEL_DBG_MSG("exceededLifespan");
            UTIMODEL_DBG_PRINT("\tp->activated_neutrophil_lifespan = %g\n\n", p->activated_neutrophil_lifespan);
        }

        return getAge()>=p->activated_neutrophil_lifespan;
    // Otherwise, die once exceeded the internal individual lifespan
    } else {
        return getAge()>=lifespan;
    }
}

// Activation
void Neutrophil::activate(){
    if (state != nresting){
        throw std::runtime_error("Neutrophil state incompatible for activation");
    }
    state = nactivated;
}

void Neutrophil::deactivate(){
    if (state != nactivated){
        throw std::runtime_error("Neutrophil state incompatible for deactivation");
    }
    state = nresting;
}


/**
                             * Mast Cells
                             */

// Create Rmacrophage, default as resting, with lifespan specified
Mastcell::Mastcell(double set_lifespan, Parameters* p):Agent(p){
   lifespan = set_lifespan;
   state = mresting;
}

/* CHRIS: placeholder for Rmacrophage destructor. There isn't actually anything to
 * delete but is included for completeness
 */

Mastcell::~Mastcell ()
{
    /* nothing to do */
}

double Mastcell::getCode(){
    switch (state) {
        case mresting:
            return 9.0;
        case mactivated:
            return 10.0;
        default:
            return 0.0;
    }
}
// Rmacrophage state getter
MastcellState Mastcell::getState(){return state;}
double Mastcell::getChemokineSource(){return 0.0;}

// Has Rmacrophage age exceeded its lifespan, and thus should die
bool Mastcell::exceededLifespan(){
    // If activated, die once exceed the set life span
    if (state==mactivated){
        if (getAge()>=p->activated_Mastcell_lifespan) {
            /* CHRIS: a few debug prints */

            UTIMODEL_DBG_MSG("exceededLifespan");
            UTIMODEL_DBG_PRINT("\tp->activated_Mastcell_lifespan = %g\n\n", p->activated_Mastcell_lifespan);
        }

        return getAge()>=p->activated_Mastcell_lifespan;
    // Otherwise, die once exceeded the internal individual lifespan
    } else {
        return getAge()>=lifespan;
    }
}


// Activation
void Mastcell::activate(){
    if (state != mresting){
        throw std::runtime_error("Mastcell state incompatible for activation");
    }
    state = mactivated;
}

void Mastcell::deactivate(){
    if (state != mactivated){
        throw std::runtime_error("Mastcell state incompatible for deactivation");
    }
    state = mresting;
}





/**
                             * Hmacrophage
                             */

// Create Rmacrophage, default as resting, with lifespan specified
Hmacrophage::Hmacrophage(double set_lifespan, Parameters* p):Agent(p){
   lifespan = set_lifespan;
   state = hresting;
   setcanActivateURHmacrophages();
}

/* CHRIS: placeholder for Rmacrophage destructor. There isn't actually anything to
 * delete but is included for completeness
 */

Hmacrophage::~Hmacrophage ()
{
    /* nothing to do */
}

double Hmacrophage::getCode(){
    switch (state) {
        case hresting:
            return 5.0;
        case hactivated:
            return 6.0;
        default:
            return 0.0;
    }
}
// Rmacrophage state getter
HmacrophageState Hmacrophage::getState(){return state;}
double Hmacrophage::getChemokineSource(){return 0.0;}

// Has Rmacrophage age exceeded its lifespan, and thus should die
bool Hmacrophage::exceededLifespan(){
    // If activated, die once exceed the set life span
    if (state==hactivated){
        if (getAge()>=p->activated_Hmacrophage_lifespan) {
            /* CHRIS: a few debug prints */

            UTIMODEL_DBG_MSG("exceededLifespan");
            UTIMODEL_DBG_PRINT("\tp->activated_Hmacrophage_lifespan = %g\n\n", p->activated_Hmacrophage_lifespan);
        }

        return getAge()>=p->activated_Hmacrophage_lifespan;
    // Otherwise, die once exceeded the internal individual lifespan
    } else {
        return getAge()>=lifespan;
    }
}


// Activation
void Hmacrophage::activate(){
    if (state != hresting){
        throw std::runtime_error("Hmacrophage state incompatible for activation");
    }
    state = hactivated;
}

void Hmacrophage::deactivate(){
    if (state != hactivated){
        throw std::runtime_error("Hmacrophage state incompatible for deactivation");
    }
    state = hresting;
}


                            /**
                             * Rmacrophage
                             */

// Create Rmacrophage, default as resting, with lifespan specified
Rmacrophage::Rmacrophage(double set_lifespan, Parameters* p):Agent(p){
   lifespan = set_lifespan;
   state = resting;
}

/* CHRIS: placeholder for Rmacrophage destructor. There isn't actually anything to
 * delete but is included for completeness
 */

Rmacrophage::~Rmacrophage ()
{
    /* nothing to do */
}

double Rmacrophage::getCode(){
    switch (state) {
        case resting:
            return 3.0;
        case activated:
            return 4.0;
        default:
            return 0.0;
    }
}
// Rmacrophage state getter
RmacrophageState Rmacrophage::getState(){return state;}


// Determine how much chemokine the cell releases
double Rmacrophage::getChemokineSource(){
    if (state == activated){
        return p->chemokine_rmac_source;
    } else {
        return 0.0;
    }
}


// Determine how much chemokine the cell releases
//double Rmacrophage::getChemokineSource(){return 0.0;}


// Has Rmacrophage age exceeded its lifespan, and thus should die
bool Rmacrophage::exceededLifespan(){
    // If activated, die once exceed the set life span
    if (state==activated){
        if (getAge()>=p->activated_Rmacrophage_lifespan) {
            /* CHRIS: a few debug prints */

            UTIMODEL_DBG_MSG("exceededLifespan");
            UTIMODEL_DBG_PRINT("\tp->activated_Rmacrophage_lifespan = %g\n\n", p->activated_Rmacrophage_lifespan);
        }

        return getAge()>=p->activated_Rmacrophage_lifespan;
    // Otherwise, die once exceeded the internal individual lifespan
    } else {
        return getAge()>=lifespan;
    }
}


// Activation
void Rmacrophage::activate(){
    if (state != resting){
        throw std::runtime_error("Rmacrophage state incompatible for activation");
    }
    state = activated;
}

void Rmacrophage::deactivate(){
    if (state != activated){
        throw std::runtime_error("Rmacrophage state incompatible for deactivation");
    }
    state = resting;
}

                              /**
                               * Bacteria
                               */

// Agent(true) here calls Agent constructor to identify that this can be phagocytosed
Bacterium::Bacterium(bool is_replicating, Parameters* p):Agent(p){
   replicating = is_replicating;
   resting = false;
   setPhagocytosable();
   setCanActivatemacrophage();
}

/* CHRIS: placeholder for Bacterium destructor. Again, similar to Rmacrophage above,
 * there isn't anything to delete but is included for completeness
 */

Bacterium::~Bacterium ()
{
    /* nothing to do */
}



bool Bacterium::isReplicating(){return replicating;}
// Change the replication type
void Bacterium::switchReplicating(){replicating=!replicating;}
bool Bacterium::isResting(){return resting;}
// Switch the resting state
void Bacterium::switchResting(){resting=!resting;}
// Get neighbourhood type to search for space for replication
bool Bacterium::isMooreReplicationNeighbourhood(){return moore_neighbourhood;}
// Change the neighbourhood type used for replication
void Bacterium::switchReplicationNeighbourhood(){moore_neighbourhood = !moore_neighbourhood;}
double Bacterium::getChemokineSource(){return 0.0;}
// Output code (1.0/1.5/2.0/2.5)
double Bacterium::getCode(){return 1.0 + !replicating;} //+ !replicating + 0.5*resting;}
