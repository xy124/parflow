#include "flowvr.h"
#include "messages.h"
#include <fca/fca.h>

#include <string.h>  // for memcpy
#include <stdlib.h>  // for malloc

FLOWVR_ACTIVE = 0;
fca_module moduleParflow;

static fca_module moduleParflowEvent;
static fca_port portIn;

void fillGridDefinition(Grid const * const grid, GridDefinition *grid_def)
{
  grid_def->nX = SubgridNX(GridBackground(grid));
  grid_def->nY = SubgridNY(GridBackground(grid));
  grid_def->nZ = SubgridNZ(GridBackground(grid));
}

void fillGridMessageMetadata(Vector const * const v, double const * const time, GridMessageMetadata *m)
{
  Grid *grid = VectorGrid(v);
  SubgridArray *subgrids = GridSubgrids(grid);
  Subgrid *subgrid;
  int g;

  ForSubgridI(g, subgrids)
  {
    subgrid = SubgridArraySubgrid(subgrids, g);
  }

  m->ix = SubgridIX(subgrid);
  m->iy = SubgridIY(subgrid);
  m->iz = SubgridIZ(subgrid);

  m->nx = SubgridNX(subgrid);
  m->ny = SubgridNY(subgrid);
  m->nz = SubgridNZ(subgrid);

  m->grid.nX = SubgridNX(GridBackground(grid));
  m->grid.nY = SubgridNY(GridBackground(grid));
  m->grid.nZ = SubgridNZ(GridBackground(grid));

  m->time = *time;
}


void NewFlowVR(void)
{
  // Refactor: shouldn't there be a GetBooleanDefault?
  NameArray switch_na = NA_NewNameArray("False True");
  char* switch_name = GetStringDefault("FlowVR", "False");

  FLOWVR_ACTIVE = NA_NameToIndex(switch_na, switch_name);
  if (FLOWVR_ACTIVE < 0)
  {
    InputError("Error: invalid print switch value <%s> for key <%s>\n",
               switch_name, "FlowVR");
    FLOWVR_ACTIVE = 0;
  }

  if (!FLOWVR_ACTIVE)
  {
    return;
  }

#ifndef HAVE_FLOWVR
  PARFLOW_ERROR("Parflow was not compiled with FlowVR but FlowVR was the input file was set to True");
  return;
#else
  D("Modname: %s, Parent: %s\n", getenv("FLOWVR_MODNAME"), getenv("FLOWVR_PARENT"));
  if (amps_size > 1)
  {
    fca_init_parallel(amps_rank, amps_size);  // TODO: amps size or amps_node_size
  }
  D("Modname: %s, Parent: %s\n", getenv("FLOWVR_MODNAME"), getenv("FLOWVR_PARENT"));
  const char* outportnamelist[] = {
    "pressure",
    "porosity",    // REM: does not really change..
    "saturation",
    "pressureSnap"
    /*"subsurf_data",         [> permeability/porosity <]*/
    /*"press",                [> pressures <]*/
    /*"slopes",               [> slopes <]*/
    /*"mannings",             [> mannings <]*/
    /*"top",                  [> top <]*/
    /*"velocities",           [> velocities <]*/
    /*"satur",                [> saturations <]*/
    /*"mask",                 [> mask <]*/
    /*"concen",               [> concentrations <]*/
    /*"wells",                [> well data <]*/
    /*"dzmult",               [> dz multiplier<]*/
    /*"evaptrans",            [> evaptrans <]*/
    /*"evaptrans_sum",        [> evaptrans_sum <]*/
    /*"overland_sum",         [> overland_sum <]*/
    /*"overland_bc_flux"      [> overland outflow boundary condition flux <]*/
    // TODO: ask BH: porosity is missing?
  };

#define n_outportnamelist (sizeof(outportnamelist) / sizeof(const char *))


  moduleParflow = fca_new_empty_module();



  for (unsigned int i = 0; i < n_outportnamelist; ++i)
  {
    fca_port port = fca_new_port(outportnamelist[i], fca_OUT, 0, NULL);
    fca_register_stamp(port, "stampTime", fca_FLOAT);
    fca_register_stamp(port, "stampFileName", fca_STRING);

    fca_append_port(moduleParflow, port);
  }
  // low: name ports? pressureOut...
  /*fca_trace trace = fca_new_trace("beginTrace", fca_trace_INT, NULL);*/
  /*if(trace != NULL) printf("Creation of trace succeded.\n"); else printf("Failed to create a trace.\n");*/

  /*fca_trace trace2 = fca_new_trace("endTrace", fca_trace_INT, NULL);*/
  /*if(trace2 != NULL) printf("Creation of trace succeded.\n"); else printf("Failed to create a trace.\n");*/

  /*fca_module modulePut = fca_new_empty_module();*/
  /*fca_append_port(modulePut, portText);*/
  /*fca_append_trace(modulePut, trace);*/
  /*fca_append_trace(modulePut, trace2);*/


  portIn = fca_new_port("in", fca_IN, 0, NULL);
  fca_append_port(moduleParflow, portIn);
  if (!fca_init_module(moduleParflow))
  {
    PARFLOW_ERROR("ERROR : init_module for moduleParflow failed!\n");
  }

  D("flowvr initialisiert.");

  /*fca_trace testTrace = fca_get_trace(modulePut,"beginTrace");*/
  /*if(testTrace == NULL) printf("ERROR : Test Trace FAIL!!\n"); else printf("Test Trace OK.\n");*/
#endif
}

#ifdef HAVE_FLOWVR

static void* translation[6];
void FlowVRinitTranslation(SimulationSnapshot *snapshot)  // TODO: macht eigentlich die uebergabe von sshot an vielen anderen stellen sinnlos!
{
  translation[VARIABLE_PRESSURE] = snapshot->pressure_out;
  translation[VARIABLE_SATURATION] = snapshot->saturation_out;
  translation[VARIABLE_KS] = NULL;  // TODO: who the fuck is KS?
  translation[VARIABLE_POROSITY] = ProblemDataPorosity(snapshot->problem_data);
  translation[VARIABLE_MANNING] = ProblemDataMannings(snapshot->problem_data);
  translation[VARIABLE_PERMEABILITY_X] = ProblemDataPermeabilityX(snapshot->problem_data);
  translation[VARIABLE_PERMEABILITY_Y] = ProblemDataPermeabilityY(snapshot->problem_data);
  translation[VARIABLE_PERMEABILITY_Z] = ProblemDataPermeabilityZ(snapshot->problem_data);
}


/// returns how much we read from buffer
size_t Steer(Variable var, Action action, const void *buffer)
{
  D("Steer");
  SteerMessageMetadata *s = (SteerMessageMetadata*)buffer;
  double *operand = (double*)(buffer + sizeof(SteerMessageMetadata));

  Vector *v = translation[var];

  Grid *grid = VectorGrid(v);
  SubgridArray *subgrids = GridSubgrids(grid);
  Subvector *subvector;
  Subgrid *subgrid;
  int g;

  ForSubgridI(g, subgrids)
  {
    subvector = VectorSubvector(v, g);
    subgrid = SubgridArraySubgrid(subgrids, g);
  }

  int nx_v = SubvectorNX(subvector);
  int ny_v = SubvectorNY(subvector);

  // TODO:speedoptimize this loop. switch to outside. maybe I can do iit with avx on whole oxes? maybe I have to change 1,1,1
  int ix = SubgridIX(subgrid);
  int iy = SubgridIY(subgrid);
  int iz = SubgridIZ(subgrid);

  int nx = SubgridNX(subgrid);
  int ny = SubgridNY(subgrid);
  int nz = SubgridNZ(subgrid);

  int i, j, k, ai = 0;
  D("operand[0]: %f", operand[0]);
  double *data;
  data = SubvectorElt(subvector, ix, iy, iz);

  size_t result = sizeof(SteerMessageMetadata) + sizeof(double) * s->nx * s->ny * s->nz;
  // TODO Check if box in this thread! only then start box loop!
  /*if (ix <= ix */

  BoxLoopI1(i, j, k, ix, iy, iz, nx, ny, nz, ai, nx_v, ny_v, nz_v, 1, 1, 1, {
    int xn = i - s->ix;
    int yn = j - s->iy;
    int zn = k - s->iz;
    // is the requested point in this chunk?
    if (xn < 0 || yn < 0 || zn < 0)
      continue;                                  // too small.
    if (xn >= s->nx || yn >= s->ny || zn >= s->nz)
      continue;                                              // too big.

    size_t index = xn + yn * s->nx + zn * s->nx * s->ny;
    switch (action)
    {
      // TODO: log steering!
      case ACTION_SET:
        if (data[ai] != operand[index])
        {
          /*D("set (%d, %d, %d) %f -> %f", i, j, k, data[ai], operand[index]);*/
        }
        data[ai] = operand[index];
        break;

      case ACTION_ADD:
        data[ai] += operand[index];
        break;

      case ACTION_MULTIPLY:
        D("diff to one: %f", operand[index] - 1.);
        data[ai] *= operand[index];
        break;

      default:
        PARFLOW_ERROR("unknown Steer Action!");
    }
  });

  // InitVectorUpdate!
  // TODO: necessary?:
  VectorUpdateCommHandle *handle;
  handle = InitVectorUpdate(v, VectorUpdateAll);
  FinalizeVectorUpdate(handle);

  D("Applied SendSteerMessage (%d) + %d + %d", sizeof(ActionMessageMetadata), sizeof(SteerMessageMetadata), sizeof(double) * s->nx * s->ny * s->nz);
  return result;
}

void SendGridDefinition(SimulationSnapshot const * const snapshot)
{
  fca_message msg = fca_new_message(moduleParflow, sizeof(GridDefinition));
  GridDefinition *g = (GridDefinition*)fca_get_write_access(msg, 0);

  fillGridDefinition(snapshot->grid, g);
  fca_put(fca_get_port(moduleParflow, "pressureSnap"), msg);
  fca_free(msg);
}


/// returns how much we read from buffer
MergeMessageParser(Interact)
{
  D("Interact %d > %d ?", size, sizeof(ActionMessageMetadata));

  if (size < sizeof(ActionMessageMetadata))
    return size;  // does not contain an action.

  SimulationSnapshot *snapshot = (SimulationSnapshot*)cbdata;
  ActionMessageMetadata *amm = (ActionMessageMetadata*)buffer;
  size_t s = sizeof(ActionMessageMetadata);

  const void *data = buffer + s;

  switch (amm->action)
  {
    case ACTION_GET_GRID_DEFINITION:
      SendGridDefinition(snapshot);
      // s+= 0
      break;

    case ACTION_TRIGGER_SNAPSHOT:
      SendSnapshot(snapshot, amm->variable);
      // s += 0;
      break;

    case ACTION_SET:
    case ACTION_ADD:
    case ACTION_MULTIPLY:
      s += Steer(amm->variable, amm->action, data);
      break;

    default:
      PARFLOW_ERROR("TODO: Unimplemented Probably somewhere adding the wrong size to s!");
  }
  /*D("processed %d / %d", s, size);*/
  return s;
}

/**
 * Executes a flowvr wait. Does the requested changes on the simulation state.
 * Returns 0 if abort was requested.
 */
int FlowVRInteract(SimulationSnapshot *snapshot)
{
  if (FLOWVR_ACTIVE)
  {
    /*D("now waiting");*/
    if (!fca_wait(moduleParflow))
      return 0;
    ParseMergedMessage(portIn, Interact, (void*)snapshot);
    // TODO: read out message on in port. do all actions that are listed there(steerings, trigger snaps...)
  }
  return 1;
}

void FreeFlowVR()
{
  if (!FLOWVR_ACTIVE)
    return;
  fca_free(moduleParflow);
}

void vectorToMessage(Vector* v, double const * const time, fca_message *result, fca_port *port)
{
  // normally really generic. low: in common with write_parflow_netcdf
  Grid *grid = VectorGrid(v);
  SubgridArray *subgrids = GridSubgrids(grid);
  Subvector *subvector;

  int g;

  ForSubgridI(g, subgrids)
  {
    subvector = VectorSubvector(v, g);
  }

  int nx_v = SubvectorNX(subvector);
  int ny_v = SubvectorNY(subvector);

  /*const fca_stamp stampMetadata = fca_get_stamp(*port, "Metadata");*/
  /*fca_write_stamp(result, stampMetadata, (void*) &stampMetadata);*/
  /*const fca_stamp stampN = fca_get_stamp(*port, "N");*/
  /*fca_write_stamp(result, stampN, 1);*/
  // write to the beginning of our memory segment
  GridMessageMetadata m;
  fillGridMessageMetadata(v, time, &m);
  size_t vector_size = sizeof(double) * m.nx * m.ny * m.nz;
  D("Sending Vector %d %d %d", m.nx, m.ny, m.nz);
  *result = fca_new_message(moduleParflow, sizeof(GridMessageMetadata) + vector_size);
  if (result == NULL)
  {
    D("Message_size: %d\n", sizeof(GridMessageMetadata) + vector_size);
    PARFLOW_ERROR("Could not create Message");
  }
  void *buffer = fca_get_write_access(*result, 0);
  D("Will write  %d bytes + %d bytes", sizeof(GridMessageMetadata), vector_size);

  memcpy(buffer, &m, sizeof(GridMessageMetadata));
  buffer += sizeof(GridMessageMetadata);


  double* buffer_double = (double*)buffer;
  double *data;
  data = SubvectorElt(subvector, m.ix, m.iy, m.iz);

  // some iterators
  int i, j, k, d = 0, ai = 0;
  BoxLoopI1(i, j, k, m.ix, m.iy, m.iz, m.nx, m.ny, m.nz, ai, nx_v, ny_v, nz_v, 1, 1, 1, { buffer_double[d] = data[ai]; d++; });
  // TODO: would be more performant if we could read the things not cell by cell I guess
}
// TODO: implement swap: do not do the memcpy but have to buffers one for read and wone for write. Change the buffers after one simulation step! (here a simulation step consists of multiple timesteps!

// Build data
typedef struct {
  const char *name;
  Vector const * const data;
} PortNameData;

void CreateAndSendMessage(SimulationSnapshot const * const snapshot, const PortNameData portnamedatas[], size_t n_portnamedata)
{
  for (unsigned int i = 0; i < n_portnamedata; ++i)
  {
    // Sometimes we do not have values for all the data...
    if (portnamedatas[i].data == NULL)
      continue;


    // Prepare the port
    fca_port port = fca_get_port(moduleParflow, portnamedatas[i].name); // TODO: maybe save all this in an array for faster acc

    // Prepare the Message
    fca_message msg;
    vectorToMessage(portnamedatas[i].data, snapshot->time, &msg, &port);


    // Create Stamps...
    const fca_stamp stampTime = fca_get_stamp(port, "stampTime");
    float time = (float)*(snapshot->time);
    fca_write_stamp(msg, stampTime, (void*)&time);

    fca_stamp stampFileName;
    if (snapshot->filename != NULL)
    {
      stampFileName = fca_get_stamp(port, "stampFileName");
      fca_write_stamp(msg, stampFileName, (void*)snapshot->filename);
    }


    // Finally send message!
    if (!fca_put(port, msg))
    {
      PARFLOW_ERROR("Could not send FlowVR-Message!");
    }
    fca_free(msg);
    D("put message!%.8f\n", *(snapshot->time));
  }
}

// REM: we are better than the nodelevel netcdf feature because during file write the other nodes are already calculating ;)
// REM: structure of nodelevel netcdf: one process per node gathers everything that has to be written and does the filesystem i/o
void DumpRichardsToFlowVR(SimulationSnapshot const * const snapshot)
{
  const PortNameData portnamedatas[] =
  {
    { "pressure", snapshot->pressure_out },
    { "porosity", snapshot->porosity_out },
    { "saturation", snapshot->saturation_out }
  };

#define n_portnamedata_ (sizeof(portnamedatas) / sizeof(const PortNameData))
  // TODO: only for those that are really connected!
  CreateAndSendMessage(snapshot, portnamedatas, n_portnamedata_);
}

void SendSnapshot(SimulationSnapshot const * const snapshot, Variable var)
{
  // TODO: extract var from snapshot
  // send snapshot!
  D("SendSnapshot");
  const PortNameData portnamedatas[] =
  {
    { "pressureSnap", snapshot->pressure_out }//,
    /*{ "porosity", porosity_out },*/
    /*{ "saturation", saturation_out }*/
  };

#define n_portnamesnapdata (sizeof(portnamedatas) / sizeof(const PortNameData))
  CreateAndSendMessage(snapshot, portnamedatas, n_portnamesnapdata);
}

void FlowVRServeFinalState(SimulationSnapshot *snapshot)
{
  NameArray switch_na = NA_NewNameArray("False True");
  char* switch_name = GetStringDefault("FlowVR.ServeFinalState", "False");


  int serve_final_state = NA_NameToIndex(switch_na, switch_name);

  if (serve_final_state < 0)
  {
    InputError("Error: invalid print switch value <%s> for key <%s>\n",
               switch_name, "FlowVR.ServeFinalState");
    serve_final_state = 0;
  }

  if (serve_final_state)
  {
    D("now serving final state.");
    while (FlowVRInteract(snapshot))
      usleep(100000);
  }
}

#endif
