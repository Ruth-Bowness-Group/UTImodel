// c++ -g -I/usr/local/include -L/usr/local/lib -lyaml-cpp -std=c++11 -lstdc++fs UTImodel.cpp Agent.cpp Environment.cpp Location.cpp Parameters.cpp
// ./a.out /Volumes/Mac/phd1/Model1/parameters.yml
//   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

#include <iostream>     // Input/Output stream
#include <utility>      // Implementing binary tuples using std::pair
#include <sys/types.h>      // OS types, for function prototypes
#include <sys/stat.h>       // Arguments or return values of type mode_t
#include <stdio.h>          // The header file stdio. h stands for Standard Input Output
#include <stdlib.h>         // General purpose standard library of C
#include <fstream>          // Read and write from a file
#include <string>           // For variables of the form std::string
#include <unistd.h>         // Header file that provides access to the POSIX operating system API
#include <vector>           // For variables of the form std::vector
#include <ctime>            /* time_t, struct tm, difftime, time, mktime */
#include "yaml-cpp/yaml.h"  // For YAML, parameter passing
#include <unordered_set>    // Unordered sets are containers that store unique elements in no particular order
#include <stdexcept>        // std::runtime_error
#include <sstream>          // std::stringstream
#include <chrono>           // Used for time monitoring
#include "Parameters.h"     // Model file
#include "Environment.h"    // Model file
#include "Output.h"         // Model file

using namespace std;


int main(int argc, char *argv[])        // Passing command line arguments
{
  (void) argc;

  std::string file_path = argv[1];      // File_path (command line input) points to the path where parameter files are kept
  /* Class "Parameters" is defined in both Parameters.cpp/.h files, to use it to construct
  in C++ we invoke it by passing the file_path argument to the constructor */
  Parameters p(file_path);
  // If needed to use any of the parameters or variables defined within this class we then prefix with "p."
  // p.lastname will change everytime as std::string -> Output + unique_time_tag
  mkdir(p.last_name, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
  struct stat info;


  if( stat( p.output_folder.c_str(), &info ) != 0 ){
      printf( "ERROR: Output folder \"%s\" does not exist. Please create before running. \n", p.output_folder.c_str() );
      return 9;
  }



  int cluster_x, cluster_y;

  /**
   * Output files
   * Filenames are built as strings (output folder + file name) then converted to char arrays as required by fopen using c_str()
   */
  // Struct to hold the files
  outputFiles files;



  files.fp=fopen((p.output_folder + std::string("/data_test.txt")).c_str(),"w");
  //files.fp2 = nullptr; /* CHRIS: generally good idea to initialise variables, even if it is to a nullptr */
  //files.fp4=fopen((p.output_folder + std::string("/initves.txt")).c_str(),"w");
  //files.analysis8 = fopen( (std::string("/Volumes/Main/model_git/Chemokine/ckn") + to_string(std::time(0))).c_str(),"w");

  // Types of cells
  //files.f1=fopen((p.output_folder + std::string("/Total.txt")).c_str(),"w");
  //files.f2=fopen((p.output_folder + std::string("/Type1.txt")).c_str(),"w");
  //files.f3=fopen((p.output_folder + std::string("/Type1_R.txt")).c_str(),"w");
  //files.f4=fopen((p.output_folder + std::string("/Type2.txt")).c_str(),"w");
  //files.f5=fopen((p.output_folder + std::string("/Type2_R.txt")).c_str(),"w");

  //files.fma=fopen((p.output_folder + std::string("/activemac.txt")).c_str(),"w");
  //files.fmr=fopen((p.output_folder + std::string("/restingmac.txt")).c_str(),"w");
  //files.fhrma = fopen((p.output_folder + std::string("/restinghelpermac.txt")).c_str(),"w");
  //files.fhama = fopen((p.output_folder + std::string("/activatedhelpermac.txt")).c_str(),"w");
  //files.frneu = fopen((p.output_folder + std::string("/restingneutrophil.txt")).c_str(),"w");
  //files.faneu = fopen((p.output_folder + std::string("/activatedneutrophil.txt")).c_str(),"w");

  //files.frmc = fopen((p.output_folder + std::string("/restingmastc.txt")).c_str(),"w");
  //files.famc = fopen((p.output_folder + std::string("/activatedmastc.txt")).c_str(),"w");


  //files.fimmunekill=fopen((p.output_folder + std::string("/immunekill.txt")).c_str(),"w");
  //files.ftotalcells=fopen((p.output_folder + std::string("/totalcells.txt")).c_str(),"w");

  files.analysis1 = fopen(("/Volumes/Main/model_git/Outputs/Output1/xtyt" + to_string(std::time(0)) + ".txt").c_str(),"w");
  files.analysis2 = fopen(("/Volumes/Main/model_git/Outputs/Output2/xtyt" + to_string(std::time(0)) + ".txt").c_str(),"w");
  files.analysis3 = fopen(("/Volumes/Main/model_git/Outputs/Output3/xtyt" + to_string(std::time(0)) + ".txt").c_str(),"w");
  files.analysis4 = fopen(("/Volumes/Main/model_git/Outputs/Output4/xtyt" + to_string(std::time(0)) + ".txt").c_str(),"w");
  files.analysis5 = fopen(("/Volumes/Main/model_git/Outputs/Output5/xtyt" + to_string(std::time(0)) + ".txt").c_str(),"w");
  files.analysis6 = fopen(("/Volumes/Main/model_git/Outputs/Output6/xtyt" + to_string(std::time(0)) + ".txt").c_str(),"w");
  files.analysis7 = fopen(("/Volumes/Main/model_git/Outputs/Output7/xtyt" + to_string(std::time(0)) + ".txt").c_str(),"w");
  files.analysis8 = fopen(("/Volumes/Main/model_git/Outputs/Output8/xtyt" + to_string(std::time(0)) + ".txt").c_str(),"w");

  /**
   * -------------------------------------------------------------
   * SETUP
   * -------------------------------------------------------------
   */

   clock_t simulation_start_time;

   // Create environment
   Environment env(&p);

   // Try/catch to catch any issues with placement of bac
   try {
       // Distribute blood vessels
       if (p.fixed_vessel_placement){
           env.initialiseBloodVesselsFixed(p.fixed_vessels);
       } else {
           env.initialiseBloodVesselsRandom(p.number_vessels);
         }
       // Distribute bacteria in a cluster, either randomly positioned or at a set point
       if (p.fixed_bacteria_placement){
           cluster_x = p.fixed_bacteria_cluster.first;
           cluster_y = p.fixed_bacteria_cluster.second;
       } else {
           cluster_x = distribution_0_grid_size(rng);
           cluster_y = distribution_0_grid_size(rng);
       }
       env.initialiseBacterialCluster(cluster_x, cluster_y);
   } catch (std::exception &err) {
       std::cout << "ERROR: initialisation: " << err.what() << "\n";
       return 9;
   }

   env.initialiseRmacrophage(p.initial_number_Rmacrophages);
   env.initialiseMastcell(p.initial_number_mastcell);



   /**
    * -------------------------------------------------------------
    * SIMULATION
    * -------------------------------------------------------------
    */
    simulation_start_time = clock();
    double time = 0.0;

    // Initial output
    //env.outputVesselsToFile(files.fp4);
    env.recordCellNumbersToFiles(files);


    while (time < p.time_limit){

        time += p.timestep;
        // TODO: avoid decimal point precision errors by rounding time to 3dp
        time = round(time * 1000.0) / 1000.0;

        std::cout << "Time = " << time << " hours " << std::endl;

        try {
            // Process
            env.processTimestep(time);


        } catch(std::exception &err ) {
            // Output any errors from the model
            std::cout << "ERROR: " << err.what() << "\n";
            return 9;
        }

        // Write to files on the output interval
        if ((int)round(time/p.timestep)%(int)round(p.output_interval/p.timestep)==0){
              //print_grid_to_file(env, p.output_folder, std::to_string((int)time));
            env.recordCellNumbersToFiles(files);

        }
    }

    // Simulation complete - output time taken
    double time_taken = (clock() - simulation_start_time)/CLOCKS_PER_SEC;
    printf("\nSIMULATION COMPLETE \n");
    printf("Time taken: %.2fs \t or %f hrs \n", time_taken, time_taken/3600.0);

    /* CHRIS: as far as I am aware, fclose doesn't check whether the FILE* argument is NULL or not.
     * Therefore, we should check manually.
     */

    if (files.fp != nullptr) fclose(files.fp);
    //if (files.fp2 != nullptr) fclose(files.fp2);
    //if (files.fp4 != nullptr) fclose(files.fp4);


    // Close files
    //if (files.f1 != nullptr) fclose(files.f1);
    //if (files.f2 != nullptr) fclose(files.f2);
    //if (files.f3 != nullptr) fclose(files.f3);
    //if (files.f4 != nullptr) fclose(files.f4);
    //if (files.f5 != nullptr) fclose(files.f5);
    //if (files.fma != nullptr) fclose(files.fma);
    //if (files.fmr != nullptr) fclose(files.fmr);
    //if (files.fhrma != nullptr) fclose(files.fhrma);
    //if (files.fhama != nullptr) fclose(files.fhama);
    //if (files.frneu != nullptr) fclose(files.frneu);
    //if (files.faneu != nullptr) fclose(files.faneu);
    //if (files.fimmunekill != nullptr) fclose(files.fimmunekill);
    //if (files.ftotalcells != nullptr) fclose(files.ftotalcells);
    if (files.analysis1 != nullptr) fclose(files.analysis1);
    if (files.analysis2 != nullptr) fclose(files.analysis2);
    if (files.analysis3 != nullptr) fclose(files.analysis3);
    if (files.analysis4 != nullptr) fclose(files.analysis4);
    if (files.analysis5 != nullptr) fclose(files.analysis5);
    if (files.analysis6 != nullptr) fclose(files.analysis6);
    if (files.analysis7 != nullptr) fclose(files.analysis7);
    if (files.analysis8 != nullptr) fclose(files.analysis8);
    return 0;
}
