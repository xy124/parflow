from flowvrapp import *
from filters import *


#pres = FilterPreSignal("Time_PreSignal", messagetype='full')  # will be inited with one token for the beginning
pres = FilterPreSignal("Time_PreSignal", nb=1)  # will be inited with one token for the beginning

simplestarterModule = Module("simplestarter", cmdline = "python ../simplestarter/simplestarter.py")
starter_inport  = simplestarterModule.getPort("beginIt")
starter_outport = simplestarterModule.addPort("out", direction='out')

parflowModule = Module("parflow", cmdline = "tclsh ./default_richards_with_netcdf.tcl 1 1 1 --FlowVR")
#parflowModule = Module("parflow", cmdline = "tclsh ./default_single.tcl 1 1 1")
for inportname in [  # TODO: cant we use the beginit port??
    "beginPort"
    ]:
  parflowModule.addPort(inportname, direction="in");
for outportname in [
    "pressure",
    "porosity",
    "saturation" ]:
  parflowModule.addPort(outportname, direction="out");

# writer
netcdfWriterModule = Module("netcdfwriter", cmdline = "/home/xy124/MasterThesis/parflow/parflow_install/bin/netcdf-writer")
writer_inport = netcdfWriterModule.addPort("pressureIn", direction = "in")
writer_outport = netcdfWriterModule.addPort("outPort", direction = "out")


# scheme: pres -> simplestarter -(timestart, timestop)-> parflow -(pressuredata)-> netcdf-writer -(pressuredata)-> [FILE] -> e.g. visit
pres.getPort("out").link(starter_inport)
starter_outport.link(parflowModule.getPort("beginPort"))
parflowModule.getPort("pressure").link(writer_inport);
writer_outport.link(pres.getPort("in"))

app.generate_xml("parflowvr")  # is not casesensitive anyway :(
