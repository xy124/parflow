/*BHEADER**********************************************************************

  Copyright (c) 1995-2009, Lawrence Livermore National Security,
  LLC. Produced at the Lawrence Livermore National Laboratory. Written
  by the Parflow Team (see the CONTRIBUTORS file)
  <parflow@lists.llnl.gov> CODE-OCEC-08-103. All rights reserved.

  This file is part of Parflow. For details, see
  http://www.llnl.gov/casc/parflow

  Please read the COPYRIGHT file or Our Notice and the LICENSE file
  for the GNU Lesser General Public License.

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License (as published
  by the Free Software Foundation) version 2.1 dated February 1999.

  This program is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY; without even the IMPLIED WARRANTY OF
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the terms
  and conditions of the GNU General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
  USA
**********************************************************************EHEADER*/

#include "Parflow.hxx"

#include "SAMRAI/hier/VariableDatabase.h"
#include "SAMRAI/pdat/CellVariable.h"
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/hier/PatchLevel.h"
#include "SAMRAI/hier/PatchHierarchy.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/pdat/CellIndex.h"
#include "SAMRAI/mesh/StandardTagAndInitialize.h"
#include "SAMRAI/mesh/BergerRigoutsos.h"
#include "SAMRAI/mesh/TreeLoadBalancer.h"
#include "SAMRAI/mesh/GriddingAlgorithm.h"

#include "SAMRAI/tbox/MemoryDatabase.h"

#include "parflow.h"

using namespace SAMRAI;

const tbox::Dimension Parflow::d_dim[4] = { tbox::Dimension::INVALID_DIMENSION, tbox::Dimension::INVALID_DIMENSION, 
					    tbox::Dimension(3), tbox::Dimension(3)};

const std::string Parflow::VARIABLE_NAME = "variable";

const std::string Parflow::CURRENT_CONTEXT = "Current";
const std::string Parflow::SCRATCH_CONTEXT = "Scratch";

Parflow::Parflow(
      const std::string& object_name,
      tbox::Pointer<tbox::Database> input_db) :
   d_object_name(object_name),
   d_input_db(input_db)
{

}

Parflow:: ~Parflow() 
{
}

void Parflow::initializeLevelData(
   const tbox::Pointer< hier::BasePatchHierarchy > hierarchy,
   const int level_number,
   const double init_data_time,
   const bool can_be_refined,
   const bool initial_time,
   const tbox::Pointer< hier::BasePatchLevel > old_level,
   const bool allocate_data)
{
   NULL_USE(hierarchy);
   NULL_USE(level_number);
   NULL_USE(init_data_time);
   NULL_USE(can_be_refined);
   NULL_USE(initial_time);
   NULL_USE(old_level);
   NULL_USE(allocate_data);
}

void Parflow::resetHierarchyConfiguration(
   const tbox::Pointer< hier::BasePatchHierarchy > hierarchy,
   const int coarsest_level,
   const int finest_level)
{
   NULL_USE(hierarchy);
   NULL_USE(coarsest_level);
   NULL_USE(finest_level);

}

void Parflow::advanceHierarchy(
   const tbox::Pointer< hier::BasePatchHierarchy > hierarchy,
   const double loop_time, 
   const double dt) 
{
   NULL_USE(hierarchy);
   NULL_USE(loop_time);
   NULL_USE(dt);
}


void Parflow::applyGradientDetector(
   const tbox::Pointer< hier::BasePatchHierarchy > hierarchy,
   const int level_number,
   const double error_data_time,
   const int tag_index,
   const bool initial_time,
   const bool uses_richardson_extrapolation_too)
{
   NULL_USE(hierarchy);
   NULL_USE(level_number);
   NULL_USE(error_data_time);
   NULL_USE(tag_index);
   NULL_USE(initial_time);
   NULL_USE(uses_richardson_extrapolation_too);

}

tbox::Pointer<hier::PatchHierarchy > Parflow::getPatchHierarchy(int dim) const
{
   return d_patch_hierarchy[dim];
}

tbox::Pointer<mesh::GriddingAlgorithm > Parflow::getGriddingAlgorithm(int dim) const
{
   return d_patch_hierarchy[dim];
}

tbox::Array<int> Parflow::getTagBufferArray(int dim) const
{
   return d_tag_buffer_array[dim];
}

tbox::Pointer<tbox::Database> Parflow::setupGridGeometryDatabase(int dim, std::string name) 
{
   
   tbox::Pointer<tbox::Database> input_db(new tbox::MemoryDatabase(name));

   Background  *bg = GlobalsBackground;

   int lower[3];
   lower[0] = BackgroundIX(bg);
   lower[1] = BackgroundIY(bg);
   lower[2] = BackgroundIZ(bg);

   int upper[3];
   upper[0] = lower[0] +  BackgroundNX(bg) - 1;
   upper[1] = lower[1] +  BackgroundNY(bg) - 1;

   if(dim == 2) {
      upper[2] = lower[2];
   } else {
      upper[2] = lower[2] +  BackgroundNZ(bg) - 1;
   }

   tbox::DatabaseBox box(d_dim[3], lower, upper);

   tbox::Array<tbox::DatabaseBox> domain(1);
   domain[0] = box;

   input_db -> putDatabaseBoxArray("domain_boxes", domain);

   double x_lo[3];
   x_lo[0] = BackgroundXLower(bg);
   x_lo[1] = BackgroundYLower(bg);
   x_lo[2] = BackgroundZLower(bg);

   double x_up[3];
   x_up[0] = BackgroundXLower(bg) + BackgroundNX(bg) * BackgroundDX(bg);
   x_up[1] = BackgroundYLower(bg) + BackgroundNY(bg) * BackgroundDY(bg);
   x_up[2] = BackgroundZLower(bg) + BackgroundNZ(bg) * BackgroundDZ(bg);

   input_db -> putDoubleArray("x_lo", x_lo, d_dim[3]);
   input_db -> putDoubleArray("x_up", x_up, d_dim[3]);

   return input_db;
}

void Parflow::initializePatchHierarchy(double time)
{
   NULL_USE(time);

   for(int d = 2; d < 4; ++d) {

      std::string grid_geometry_name("CartesianGeometry" + tbox::Utilities::intToString(d,1) + "D");
      
      tbox::Pointer<geom::CartesianGridGeometry > grid_geometry(
	 new geom::CartesianGridGeometry(
	    d_dim[d],
	    grid_geometry_name,
	    setupGridGeometryDatabase(d, grid_geometry_name)));
      
      std::string patch_hierarchy_name("PatchHierarchy" + tbox::Utilities::intToString(d,1) + "D");
      d_patch_hierarchy[d] = new hier::PatchHierarchy(patch_hierarchy_name,
						      grid_geometry);
      
      hier::VariableDatabase *variable_database(
	 hier::VariableDatabase::getDatabase());
      
      int depth = 1;
      
      tbox::Pointer< pdat::CellVariable<double> > cell_state(
	 new pdat::CellVariable<double>(d_dim[d], VARIABLE_NAME, depth));
      
      tbox::Pointer< hier::VariableContext > current_context(
	 variable_database -> getContext(CURRENT_CONTEXT));
      
      tbox::Pointer< hier::VariableContext > scratch_context(
	 variable_database -> getContext(SCRATCH_CONTEXT));
      
      hier::IntVector ghosts(d_dim[d], 1);
      
      getFromInput(d_input_db, false);
      
      std::string standard_tag_and_initialize_name("StandardTagAndInitialize" + tbox::Utilities::intToString(d,1) + "D");
      tbox::Pointer< mesh::StandardTagAndInitialize > standard_tag_and_initialize(
	 new mesh::StandardTagAndInitialize(
	    d_dim[d],
	    standard_tag_and_initialize_name, 
	    this,
	    d_input_db -> getDatabase(standard_tag_and_initialize_name)));
      
      tbox::Pointer< mesh::BergerRigoutsos > box_generator( 
	 new mesh::BergerRigoutsos(d_dim[d]));
      
      std::string load_balancer_name("LoadBalancer" + tbox::Utilities::intToString(d,1) + "D");
      tbox::Pointer< mesh::LoadBalanceStrategy > load_balancer(
	 new mesh::TreeLoadBalancer(d_dim[d],
				    load_balancer_name,
				    d_input_db -> getDatabase(load_balancer_name)));
      
      std::string gridding_algorithm_name("GriddingAlgorithm" + tbox::Utilities::intToString(d,1) + "D");      
      d_gridding_algorithm[d] = 
	 new mesh::GriddingAlgorithm(
	    d_dim[d],
	    gridding_algorithm_name,
	    d_input_db -> getDatabase(gridding_algorithm_name),
	    standard_tag_and_initialize,
	    box_generator,
	    load_balancer);
      
      d_tag_buffer_array[d].resizeArray(d_gridding_algorithm[d]->getMaxLevels());
      for (int il = 0; il < d_gridding_algorithm[d]->getMaxLevels(); il++) {
	 d_tag_buffer_array[d][il] = 1;
      }
      
      const hier::IntVector one_vector(d_dim[d], 1);
      const hier::IntVector ratio(d_gridding_algorithm[d] -> getMaxLevels() > 1 ? 
				  d_gridding_algorithm[d] -> getRatioToCoarserLevel(1) : one_vector);

      std::vector<hier::IntVector > fine_connector_gcw;
      std::vector<hier::IntVector > peer_connector_gcw;
      d_gridding_algorithm[d] -> computeAllConnectorWidths(fine_connector_gcw,
							   peer_connector_gcw,
							   *d_patch_hierarchy[d] );
      
      d_patch_hierarchy[d]->getMappedBoxHierarchy().setMappedBoxLevelParameters(
	 0,
	 ratio,
	 fine_connector_gcw[0],
	 peer_connector_gcw[0]);
      
      tbox::Pointer< hier::MappedBoxLevel > mapped_box_level(createMappedBoxLevelFromParflowGrid());

      
      d_patch_hierarchy[d] -> makeNewPatchLevel(0, *mapped_box_level);

      tbox::Pointer<hier::PatchLevel> level(d_patch_hierarchy[d] -> getPatchLevel(0));

      for(int i = 4; i > 0; --i) {
	 level -> getMappedBoxLevel() -> getPersistentOverlapConnectors().
	    findOrCreateConnector(
	       *level -> getMappedBoxLevel(),
	       hier::IntVector(d_dim[d], i));
      }
   }
}

void Parflow::getFromInput(
   tbox::Pointer<tbox::Database> db,
   bool is_from_restart) 
{
   NULL_USE(db);
   
   if (!is_from_restart) {

   }
}

tbox::Pointer< hier::MappedBoxLevel > Parflow::createMappedBoxLevelFromParflowGrid(void)
{

   Grid *grid = CreateGrid(GlobalsUserGrid);

   // Build a box based off of Parflow grid
   const hier::Index lower( SubgridIX(GridSubgrid(grid, 0)),
			    SubgridIY(GridSubgrid(grid, 0)),
			    SubgridIZ(GridSubgrid(grid, 0)));

   const hier::Index upper ( lower[0] + SubgridNX(GridSubgrid(grid, 0)) - 1,
			     lower[1] + SubgridNY(GridSubgrid(grid, 0)) - 1,
			     lower[2] + SubgridNZ(GridSubgrid(grid, 0)) - 1);

   hier::Box box(lower, upper);

   std::cout << "In Parflow Grid create box " << box << std::endl;

   // Build a mapped box and insert into layer.  
   const int local_index = 0;
   const int my_rank = tbox::SAMRAI_MPI::getRank();

   hier::MappedBox mapped_box( box, local_index, my_rank );

   hier::MappedBoxSet mapped_box_set;
   mapped_box_set.insert(mapped_box_set.begin(), mapped_box);

   const hier::IntVector ratio(d_dim[3], 1);

   tbox::Pointer<hier::MappedBoxLevel> mapped_box_level(new hier::MappedBoxLevel(mapped_box_set,
							   ratio));

   FreeGrid(grid);

   return mapped_box_level;
}





