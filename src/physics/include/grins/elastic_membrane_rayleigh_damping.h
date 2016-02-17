//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// GRINS - General Reacting Incompressible Navier-Stokes
//
// Copyright (C) 2014-2015 Paul T. Bauman, Roy H. Stogner
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

#ifndef GRINS_ELASTIC_MEMBRANE_RAYLEIGH_DAMPING_H
#define GRINS_ELASTIC_MEMBRANE_RAYLEIGH_DAMPING_H

#include "grins/elastic_membrane_base.h"

namespace GRINS
{
  template<typename StressStrainLaw>
  class ElasticMembraneRayleighDamping : public ElasticMembraneBase<StressStrainLaw>
  {
  public:

    ElasticMembraneRayleighDamping( const PhysicsName& physics_name,
                                 const GetPot& input,
                                 bool is_compressible );

    virtual ~ElasticMembraneRayleighDamping(){};

    //! Time dependent part(s) of physics for element interiors
    virtual void element_time_derivative( bool compute_jacobian,
                                          AssemblyContext& context,
                                          CachedValues& cache );

    virtual void mass_residual( bool compute_jacobian,
                                AssemblyContext& context,
                                CachedValues& /*cache*/ )
    { this->mass_residual_impl(compute_jacobian,
                               context,
                               &libMesh::FEMContext::interior_rate,
                               &libMesh::DiffContext::get_elem_solution_rate_derivative,
                               _mu_factor); }

  protected:

    libMesh::Real _lambda_factor;
    libMesh::Real _mu_factor;

  private:

    ElasticMembraneRayleighDamping();

  };

} // end namespace GRINS

#endif // GRINS_ELASTIC_MEMBRANE_RAYLEIGH_DAMPING_H
