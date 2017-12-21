#include "flowvr.h"
#include "messages.h"
#include <fca/fca.h>

#include <string.h>  // for memcpy
#include <stdlib.h>  // for malloc


int FLOWVR_ACTIVE;
fca_module moduleParflow;

static fca_module moduleParflowEvent;
static fca_port portIn;

void fillGridDefinition(Grid const * const grid, GridDefinition *grid_def)
{
  grid_def->nX = SubgridNX(GridBackground(grid));
  grid_def->nY = SubgridNY(GridBackground(grid));
  grid_def->nZ = SubgridNZ(GridBackground(grid));
}

void fillGridMessageMetadata(Vector const * const v, double const * const time, const Variable variable, GridMessageMetadata *m)
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
  m->variable = variable;
}

typedef struct {
  const char * port_name;
  Variable variable;
  int offset;
  int periodicity;
} Contract;

Contract *contracts;
size_t n_contracts;

/**
 * Will init out ports and create contracts from it.
 */
void InitContracts()
{
  NameArray port_names = NA_NewNameArray(GetString("FlowVR.Outports.Names"));

  n_contracts = NA_Sizeof(port_names);
  contracts = ctalloc(Contract, n_contracts);

  char key[256];
  for (size_t i = 0; i < n_contracts; ++i)
  {
    contracts[i].port_name = NA_IndexToName(port_names, i);
    sprintf(key, "FlowVR.Outports.%s.Periodicity", contracts[i].port_name);
    contracts[i].periodicity = GetInt(key);
    sprintf(key, "FlowVR.Outports.%s.Offset", contracts[i].port_name);
    contracts[i].offset = GetInt(key);
    sprintf(key, "FlowVR.Outports.%s.Variable", contracts[i].port_name);
    contracts[i].variable = NameToVariable(GetString(key));

    D("Add Contract %s executed with periodicity: %d (offset: %d)",
      contracts[i].port_name, contracts[i].periodicity, contracts[i].offset);
    // REM we do not need to free Strings obtained by GetString as they are just requestet from a database that was already cached.
  }
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
  InitContracts();
  D("Modname: %s, Parent: %s\n", getenv("FLOWVR_MODNAME"), getenv("FLOWVR_PARENT"));
  if (amps_size > 1)
  {
    fca_init_parallel(amps_rank, amps_size);  // TODO: amps size or amps_node_size
  }
  D("Modname: %s, Parent: %s\n", getenv("FLOWVR_MODNAME"), getenv("FLOWVR_PARENT"));
  /*"pressure",*/
  /*"porosity",    // REM: does not really change..*/
  /*"saturation",*/
  /*"pressureSnap"*/
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


  moduleParflow = fca_new_empty_module();

  fca_port portPressureSnap = fca_new_port("pressureSnap", fca_OUT, 0, NULL);
  fca_register_stamp(portPressureSnap, "stampTime", fca_FLOAT);
  fca_register_stamp(portPressureSnap, "stampFileName", fca_STRING);

  for (size_t i = 0; i < n_contracts; ++i)
  {
    D("Add outport %s", contracts[i].port_name);
    fca_port port = fca_new_port(contracts[i].port_name, fca_OUT, 0, NULL);
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

static void* translation[VARIABLE_LAST];
void FlowVRinitTranslation(SimulationSnapshot *snapshot)  // TODO: macht eigentlich die uebergabe von sshot an vielen anderen stellen sinnlos!
{
  translation[VARIABLE_PRESSURE] = snapshot->pressure_out;
  translation[VARIABLE_SATURATION] = snapshot->saturation_out;
  translation[VARIABLE_POROSITY] = ProblemDataPorosity(snapshot->problem_data);
  translation[VARIABLE_MANNING] = ProblemDataMannings(snapshot->problem_data);
  translation[VARIABLE_PERMEABILITY_X] = ProblemDataPermeabilityX(snapshot->problem_data);
  translation[VARIABLE_PERMEABILITY_Y] = ProblemDataPermeabilityY(snapshot->problem_data);
  translation[VARIABLE_PERMEABILITY_Z] = ProblemDataPermeabilityZ(snapshot->problem_data);
}

static inline int simple_intersect(int ix1, int nx1, int ix2, int nx2,
                                   int iy1, int ny1, int iy2, int ny2,
                                   int iz1, int nz1, int iz2, int nz2)
{
  int d;

  d = ix2 - ix1;
  if (-nx2 > d || d > nx1)
    return 0;
  d = iy2 - iy1;
  if (-ny2 > d || d > ny1)
    return 0;
  d = iz2 - iz1;
  if (-nz2 > d || d > nz1)
    return 0;
  return 1;
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

  // TODO:speedoptimize the BoxLoop! loop. switch to outside. maybe I can do it with avx on whole oxes? maybe I have to change 1,1,1
  // probably one can use IntersectSubgrids and loop only over intersction! Maybe this could influence the vectorupdate too!
  int ix = SubgridIX(subgrid);
  int iy = SubgridIY(subgrid);
  int iz = SubgridIZ(subgrid);

  int nx = SubgridNX(subgrid);
  int ny = SubgridNY(subgrid);
  int nz = SubgridNZ(subgrid);

  int i, j, k, ai = 0;
  double *data;
  data = SubvectorElt(subvector, ix, iy, iz);

  size_t read_out_size = sizeof(SteerMessageMetadata) + sizeof(double) * s->nx * s->ny * s->nz;

  // Check if box in this thread! only then start box loop!
  if (!simple_intersect(ix, nx, s->ix, s->nx,
                        iy, ny, s->iy, s->ny,
                        iz, nz, s->iz, s->nz))
    return read_out_size;

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
        /*if (data[ai] != operand[index])*/
        /*{*/
        /*D("set (%d, %d, %d) %f -> %f", i, j, k, data[ai], operand[index]);*/
        /*}*/
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
  return read_out_size;
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
  tfree(contracts);
}

void vectorToMessage(const Variable variable, double const * const time, fca_message *result, fca_port *port)
{
  Vector *v = translation[variable];
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
  fillGridMessageMetadata(v, time, variable, &m);
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


void CreateAndSendMessage(SimulationSnapshot const * const snapshot, const char * portname, Variable var)
{
  // Prepare the port
  fca_port port = fca_get_port(moduleParflow, portname); // TODO: maybe save all this in an array for faster acc

  // Prepare the Message
  fca_message msg;

  vectorToMessage(var, snapshot->time, &msg, &port);


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

// REM: we are better than the nodelevel netcdf feature because during file write the other nodes are already calculating ;)
// REM: structure of nodelevel netcdf: one process per node gathers everything that has to be written and does the filesystem i/o

void SendSnapshot(SimulationSnapshot const * const snapshot, Variable var)
{
  // TODO: extract var from snapshot
  // send snapshot!
  D("SendSnapshot");

  CreateAndSendMessage(snapshot, "pressureSnap", VARIABLE_PRESSURE);
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
// TODO: gross/kleinschreibung von function names!

/**
 * Dumps data if it has to since a contract.
 */
int FlowVRFullFillContracts(int timestep, SimulationSnapshot const * const sshot)
{
  int res = 0;

  for (size_t i = 0; i < n_contracts; ++i)
  {
    if ((abs(timestep - contracts[i].offset) % contracts[i].periodicity) == 0)
    {
      res = 1;
      CreateAndSendMessage(sshot, contracts[i].port_name, contracts[i].variable);
    }
  }
  return res;
}

#endif
