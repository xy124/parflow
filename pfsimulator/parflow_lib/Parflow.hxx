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

#ifndef included_parflow_Parflow
#define included_parflow_Parflow

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/mesh/StandardTagAndInitStrategy.h"
#include "SAMRAI/mesh/GriddingAlgorithm.h"
#include "SAMRAI/geom/CartesianGridGeometry.h"
#include "SAMRAI/hier/BasePatchHierarchy.h"
#include "SAMRAI/xfer/RefineAlgorithm.h"
#include "SAMRAI/xfer/CoarsenAlgorithm.h"
#include "SAMRAI/tbox/Pointer.h"
#include "SAMRAI/tbox/Database.h"

#include <string>

using namespace SAMRAI;

class Parflow:
   public mesh::StandardTagAndInitStrategy
{

public:

   Parflow(
      const std::string& object_name,
      tbox::Pointer<tbox::Database> input_db);
   
   virtual ~Parflow();

   virtual void initializeLevelData(
      const tbox::Pointer< hier::BasePatchHierarchy > hierarchy,
      const int level_number,
      const double init_data_time,
      const bool can_be_refined,
      const bool initial_time,
      const tbox::Pointer< hier::BasePatchLevel > old_level = 
      tbox::Pointer< hier::BasePatchLevel >(NULL),
      const bool allocate_data = true);

   virtual void resetHierarchyConfiguration(
      const tbox::Pointer< hier::BasePatchHierarchy > hierarchy,
      const int coarsest_level,
      const int finest_level);

   void advanceHierarchy(
      const tbox::Pointer< hier::BasePatchHierarchy > hierarchy,
      const double loop_time, 
      const double dt);

   virtual void applyGradientDetector(
      const tbox::Pointer< hier::BasePatchHierarchy > hierarchy,
      const int level_number,
      const double error_data_time,
      const int tag_index,
      const bool initial_time,
      const bool uses_richardson_extrapolation_too);

   tbox::Pointer<hier::PatchHierarchy > getPatchHierarchy(void) const;

   tbox::Pointer<mesh::GriddingAlgorithm > getGriddingAlgorithm(void) const;

   tbox::Array<int> getTagBufferArray(void) const;

   void initializePatchHierarchy(double time);

  private:

   tbox::Pointer<hier::MappedBoxLevel> createMappedBoxLevelFromParflowGrid(void);

   void setupInputDatabase(void);

   void getFromInput(
      tbox::Pointer<tbox::Database> db,
      bool is_from_restart);

   static const tbox::Dimension d_dim;

   // FIXME rename this
   static const std::string CELL_STATE;

   static const std::string CURRENT_STATE;
   static const std::string SCRATCH_STATE;


   // Not implemeted; standard C++'ism to avoid bad things
   Parflow(const Parflow&);
   void operator=(const Parflow&);
   
   tbox::Pointer<hier::PatchHierarchy > d_patch_hierarchy;

   tbox::Pointer< mesh::GriddingAlgorithm > d_gridding_algorithm;

   tbox::Pointer< xfer::RefineAlgorithm > d_boundary_fill_refine_algorithm;
   tbox::Pointer< xfer::RefineAlgorithm > d_fill_after_regrid;
   tbox::Pointer< xfer::CoarsenAlgorithm > d_coarsen_algorithm;
   
   tbox::Array< tbox::Pointer< xfer::RefineSchedule > > d_boundary_schedule_advance;
   tbox::Array< tbox::Pointer< xfer::CoarsenSchedule > > d_coarsen_schedule;

   tbox::Array<int> d_tag_buffer_array;

   int d_current_cell_state_handle;
   int d_scratch_cell_state_handle;

   std::string d_object_name;

   tbox::Pointer<tbox::Database> d_input_db;


};

#endif
