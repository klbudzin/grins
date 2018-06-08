//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// GRINS - General Reacting Incompressible Navier-Stokes
//
// Copyright (C) 2014-2017 Paul T. Bauman, Roy H. Stogner
// Copyright (C) 2010-2013 The PECOS Development Team
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the Version 2.1 GNU Lesser General
// Public License as published by the Free Software Foundation.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc. 51 Franklin Street, Fifth Floor,
// Boston, MA  02110-1301  USA
//
//-----------------------------------------------------------------------el-


// GRINS
#include "grins_config.h"

#ifdef GRINS_HAVE_ANTIOCH

// C++
#include <iomanip>
#include <vector>
#include <fstream>

// libMesh
#include "libmesh/getpot.h"

// GRINS
#include "grins/antioch_mixture.h"
#include "grins/antioch_kinetics.h"
#include "grins/materials_parsing.h"
#include "grins/physics_naming.h"
#include "grins/antioch_evaluator.h"

#ifdef GRINS_HAVE_CANTERA
#include "grins/cantera_mixture.h"
#include "grins/cantera_kinetics.h"
#include "grins/cantera_evaluator.h"

libMesh::Real rhof(libMesh::Real P,libMesh::Real Rmix,libMesh::Real T) {
  libMesh::Real value = 0;
  value = P/(Rmix*T);
  return value;
}

int main(int argc, char* argv[])
{
   if( argc < 2 )
    {
      // TODO: Need more consistent error handling.
      std::cerr << "Error: Must specify input file." << std::endl;
      exit(1);
    }
  GetPot input( argv[1] );
  
  //Antioch Mixture and kinetics Initialization ***************************************************************
  GRINS::AntiochMixture<Antioch::CEACurveFit<libMesh::Real> >
    antioch_mixture(input,"HydrogenGas");

  const unsigned int n_species = antioch_mixture.n_species();
  
  GRINS::AntiochEvaluator<Antioch::CEACurveFit<libMesh::Real> , Antioch::IdealGasMicroThermo<Antioch::NASAEvaluator<libMesh::Real, Antioch::CEACurveFit<libMesh::Real> >, libMesh::Real> > antioch_evaluator( antioch_mixture);

  //Cantera Mixture and kinetics Initialization ***************************************************************
  GRINS::CanteraMixture cantera_mixture(input, "Hydrogen");
  GRINS::CanteraEvaluator cantera_evaluator(cantera_mixture);

  if(n_species != cantera_mixture.n_species()) 
    {
      std::cerr<< "Error: Cantera and Antioch Reading different number of species" << std::endl;
      exit(1);
    }

  //Temperature and Mass Fractions Setting ********************************************************************
  std::vector<libMesh::Real> Temperature_Dist(100,300);
  for(unsigned int i=0;i<Temperature_Dist.size();i++)
    Temperature_Dist[i] = 300 +25*i;

  std::vector<libMesh::Real> Mass_Fractions(n_species);
   for( unsigned int s = 0; s < n_species; s++ )
    {
      Mass_Fractions[s] = input( "Conditions/"+antioch_evaluator.species_name(s), 0.00);
    }

   libMesh::Real p0 = 100000; //pascals

   //Antioch Kinetic Parameters Calculations ******************************************************************
   libMesh::Real R_mix_antioch = antioch_evaluator.R_mix(Mass_Fractions);
   libMesh::Real M_mix_antioch = antioch_evaluator.M_mix(Mass_Fractions);
   std::vector<std::vector<libMesh::Real> > omega_dot_antioch(Temperature_Dist.size());
   for(unsigned int i=0;i<Temperature_Dist.size();i++)
     {
       omega_dot_antioch[i].resize(n_species);
       libMesh::Real rho = rhof(p0,R_mix_antioch,Temperature_Dist[i]);
       antioch_evaluator.omega_dot(Temperature_Dist[i],rho,Mass_Fractions,omega_dot_antioch[i]);
     }

   //Cantera Kinetic Parameters Calculations *****************************************************************
   libMesh::Real R_mix_cantera = cantera_evaluator.R_mix(Mass_Fractions);
   libMesh::Real M_mix_cantera = cantera_evaluator.M_mix(Mass_Fractions);
   std::vector<std::vector<libMesh::Real> > omega_dot_cantera(Temperature_Dist.size());
   for(unsigned int i=0; i < Temperature_Dist.size();i++)
     {
       omega_dot_cantera[i].resize(n_species,4);
       libMesh::Real rho = rhof(p0,R_mix_cantera,Temperature_Dist[i]);
       cantera_evaluator.omega_dot(Temperature_Dist[i],rho,Mass_Fractions,omega_dot_cantera[i]); //returns 0 something is horribly wrong here
     }
   
   //Printing out What we want and saving data to our files

   std::ofstream output;
   output.open("omega_dot_antioch.dat", std::ios::trunc );
   for(unsigned int s =0; s < n_species; s++)
     {
       output << antioch_mixture.species_name(s) << " ";
       for(unsigned int i = 0; i < Temperature_Dist.size();i++)
	 {
	   output << omega_dot_antioch[i][s] << " ";
	 }
       output << std::endl;
     }
   output.close();

   output.open("omega_dot_cantera.dat",std::ios::trunc );
   for(unsigned int s =0; s<n_species; s++)
     {
       output << cantera_mixture.species_name(s) << " ";
       for(unsigned int i=0; i < Temperature_Dist.size();i++)
	 {
	   output<< omega_dot_cantera[i][s] << " ";
	 }
       output << std::endl;
     }
   output.close();
     
   output.open("Cantera_Antioch_Comparison.txt");
   output << "       |Cantera         |Antioch         |" << std::endl;
   output << "R_Mix  |" <<  std::setprecision(16) << R_mix_cantera <<"|" <<  std::setprecision(16) << R_mix_antioch << "|" << std::endl;
   output << "M_Mix  |" << std::setprecision(16) << M_mix_cantera << "|" << std::setprecision(16) << M_mix_antioch << "|" << std::endl;
   //Testing to check individual species molecular weights
   for(unsigned int s=0;s<n_species;s++)
     output <<std::setprecision(7) << antioch_mixture.species_name(s) <<"|" << std::setprecision(16) << cantera_evaluator.M(s) << "|" << std::setprecision(16) << antioch_evaluator.M(s) << "|" << std::endl;

   output.close();
   
   //Testing to check individual species molecular weights
   










   return 0;
}
#endif //GRINS_HAVE_CANTERA
#endif //GRINS_HAVE_ANTIOCH
