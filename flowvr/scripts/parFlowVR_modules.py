from flowvrapp import *
from filters import *
import socket
import subprocess
import os



class FlowvrRunSlurm(Run):
    """ the id of the modules must be built from the prefix, suffixed by
    /0, /1, etc. hosts and name will not be populated automatically from
    the corresponding module.
    Optionnal parameters :
     rankfile : is used to pass a rankfile to fix the affinity of the MPI processes
    """

    def __init__(self, cmdline, prefix, hosts):
      self.cmdline = cmdline
      self.name = prefix

      self.options = '-v '

      #if bindobject != ''  and  bindnumber != '':
      #  os.environ['FLOWVR_MODULE_HWBIND'] = bindobject


      #if ownShmem:
      #  os.environ['FLOWVR_MODULE_OWNSHMEM'] = '1'

      if app.current_path: self.name = app.current_path + self.name
      hostlist = hosts.split(",")
      self.np = len(hostlist)
      self.hosts = list(set(hostlist))  # each node only once!
      self.hosts.sort()


      self.n_tasks_per_node = len(hostlist) // len(self.hosts)


    def add_environment_variable(self, variable, value = ""):
        os.environ[variable] = value

    #Note : For now we don't use the binding object here to avoid undetermined behavior
    #http://www.open-mpi.org/faq/?category=tuning#using-paffinity-v1.4
    def get_cmdline(self):
        return "srun %s -n %d --ntasks-per-node=%d --nodelist=%s --cpu-bind=rank --exclusive %s" % (self.options, self.np, self.n_tasks_per_node, ",".join(self.hosts), self.cmdline)


# helper function:
def preferLastCore(lastCore, cores):
    if lastCore and cores == '':
        import multiprocessing
        return str(multiprocessing.cpu_count()-1)
        # test by running taskset -p `pgrep netcdf-writer | head -n1`
    return cores


# Modules:
class FilterMergeGridMessages(flowvrapp.Filter):
    """A filter converting concatenated grid messages into one big grid message
    containing the full grid (nX,nY,nZ)
    of timestep t as soon as he got all the necessary data.

    WARNING: we abort if not all the data necessary was found!

    Input ports:
    -  in: Grid messages to be joined.

    Output ports:
    - out: Big grid messages.
    """

    def __init__(self, name, host = ''):
        flowvrapp.Filter.__init__(self, name, run = 'flowvr.plugins.MergeGridMessages', host = host)
        self.addPort('in', direction = 'in')
        self.addPort('out', direction = 'out')

class FilterMergeItExt(FilterWithManyInputs):
    """Merges messages received on input port into one message sent on output port if
    a message on the order port is received.

    A filter which sends ALWAYS a message (sometimes containing an empty buffer)
    if a message on the order port is received. If possible the Message consists of the
    concatenation of all messages waiting on the in-ports. The stamps of the first
    found message are hereby copied.

    Input ports:
    -  in: Messages to be filtered.
    -  order: Filtering orders  (from a  synchronizer  such as
    flowvr::plugins::GreedySynchronizer).

    Output ports:
    - out: Filtered messages.

    Use newInputPort() to get a new input port (instead of getPort()).
    """

    def __init__(self, name, host = ''):
        self.messagetype = 'full'
        self.plugin_name = 'flowvr.plugins.MergeItExt'
        FilterWithManyInputs.__init__(self, name, host = host)
        self.addPort('order', direction = 'in', messagetype = 'stamps')


class Parflow(Module):
    """Instantiate a single core parflow instance. Use outports parameter to
    set an array of out ports according to this parflow instances tcl file
    (see pfset Flowvr.Outports.Names).
    The parameter problemName is the problem's name specifying the pfidb file to load
    (must be the same as the pfidb file name without extension).

    Input ports:
    - in: receives action messages

    Output ports:
    - snap: sends snapshot data to a visualizer when requested of an incoming
      action message
    - ...: output ports according to outports parameter
    """

    def __init__(self, prefix, index=None, run=None, host=None, cmdline=None, outports=[], problemName=None):
        name = prefix + '/' + str(index) if index is not None else prefix

        if cmdline is None:
            pf_cmd = os.getenv('PARFLOW_EXE') or '$PARFLOW_DIR/bin/parflow'
            cmdline = '%s %s' % (pf_cmd, problemName)

        Module.__init__(self, name, run = run, host = host, cmdline = cmdline)

        inportnames = ['in']
        outportnames = outports + ['snap']


        for inportname in inportnames:
            self.addPort(inportname, direction = 'in')

        for outportname in outportnames:
            self.addPort(outportname, direction = 'out')

class ParflowMPI(Composite):
    """Several instances of parflow module that generate pressure, porosity... for the
    next timestep. The hosts parameter defines how many instances to launch and on which
    node. Alternatively this can be specified with a rankfile.
    If no rankfile is found and easyPinning is True, the parflow processes get pinned on
    the first cores per node. debugprefix can be used to set a debugger command that is
    prepended. MPI run args are loaded as defined during CMAKE build.
    MPIEXEC_PREFLAGS and MPIEXEC_POSTFLAGS have impact. The PARFLOW_EXE environment
    variable can be used to specify the path to use for the parflow executable"""

    def __init__(self, hosts, problemName, outports=[], debugprefix=None, rankfile='',
            easyPinning=True):
        Composite.__init__(self)

        if debugprefix is None:
            debugprefix=''
        else:
            debugprefix+=' '

        prefix = 'parflow_'+problemName
        parflowrun=None
        cmdline=None

        pf_cmd = os.getenv('PARFLOW_EXE') or '$PARFLOW_DIR/bin/parflow'

        # hosts: string with host names, separated by spaces
        if hosts.find(',') == -1:
            # only one host. start locally!
            #pf_cmd = '$PARFLOW_DIR/bin/parflow'
            cmdline = '%s %s %s' % (debugprefix, pf_cmd, problemName)
        else:
            # load runargs from file...
            command = ['bash', '-c', '. $PARFLOW_DIR/config/pf-cmake-env.sh && echo \"$MPIEXEC_PREFLAGS $MPIEXEC_POSTFLAGS\"']

            proc = subprocess.Popen(command, stdout = subprocess.PIPE)
            mpirunargs = '--report-bindings'  # -mca plm_rsh_no_tree_spawn 1 '
            if rankfile == '' and easyPinning:
                mpirunargs += ' --map-by core --bind-to core '



            for line in proc.stdout:
                mpirunargs += ' ' + line.decode('utf-8')

            proc.communicate()
            parflowrun = FlowvrRunSlurm('%s %s %s' %
                    (debugprefix, pf_cmd, problemName), hosts = hosts, prefix = prefix)


        # hosts_list: convert hosts to a list
        hosts_list = hosts.split(',')

        for i, host in enumerate(hosts_list):
            parflow = Parflow(prefix, index=i, run=parflowrun, host=host,
                    cmdline=cmdline,
                    outports=outports)

            # collect ports
            for pname in parflow.ports:
                p = parflow.ports[pname]
                if not pname in self.ports:
                    self.ports[pname] = list()
                self.ports[pname].append(p)


class VisIt(Module):
    """Module creating a visit simualation  host in ~/.visit/simulations that can be
    opened and interacted from VisIt

    Input ports:
    - in: gets grid messages and grid definition messages. Must be connected to ParFlow's
      snap port

    Output ports:
    - triggerSnap: sends action messages requesting snapshot data from ParFlow. Thus must
      be connected to ParFlow's in port"""

    def __init__(self, name, host=''):
        Module.__init__(self, name, cmdline = '$PARFLOW_DIR/bin/visit-connector-starter.sh', host=host)
        #Module.__init__(self, name, cmdline = 'xterm -e gdb $PARFLOW_DIR/bin/visit-connector')
        self.addPort('triggerSnap', direction = 'out')
        self.addPort('in', direction = 'in')

class NetCDFWriter(Module):
    """Module NetCDFWriter writes parflow output into netCDF files.
    If lastCore parameter is True it will be
    pinned on the last CPU core. Thus it is not interfering with the running ParFlow
    instances if there is a core left.
    If abortOnEnd parameter is set to true it will abort the workflow if an empty message
    is received.
    fileprefix will be the prefix of written files. This is especially useful when more
    than one NetCDFWriter is used.
    The NETCDF_WRITER_EXE environment variable can be used to specify the path to use
    for the NetCDFWriter executable

    Input ports:
    - in: receiving grid messages that will be written into a file"""
    def __init__(self, name, fileprefix='', host='', abortOnEnd=True, lastCore=True, cores=''):
# last core property is nice for pinning ;)
        # TODO: works with mpi too?!

        cores = preferLastCore(lastCore, cores)
        nw_cmd = os.getenv('NETCDF_WRITER_EXE') or '$PARFLOW_DIR/bin/netcdf-writer'
        Module.__init__(self, name, cmdline = '%s %s %s' %
                (nw_cmd, fileprefix, '--no-abort' if not abortOnEnd else ''), host=host, cores=cores)
        self.addPort('in', direction = 'in');

class NetCDFWriterMPI(Composite):
    """Module NetCDFWriterMPI writes parflow output into netCDF files.
    Used for node level fileIO. Runs on each host given in hosts a NetCDF Writer module.
    All of those share a MPI world and thus parallel netCDF is supported: multiple processes
    can open and change the same file at the same time.
    If lastCore parameter is True it will be
    pinned on the last CPU core of each node. Thus it is not interfering with the running ParFlow
    instances if there is a core left.
    If abortOnEnd parameter is set to true it will abort the workflow if an empty message
    is received.
    fileprefix will be the prefix of written files. This is especially useful when more
    than one NetCDFWriter is used.
    The NETCDF_WRITER_EXE environment variable can be used to specify the path to use
    for the NetCDFWriter executable

    Input ports:
    - in: Receiving grid messages that will be written into a file.
          One in port for each instance. getPort('in') will return a list of in ports
          with the index number corresponding to instances defined in the hosts list.
    """

    def __init__(self, prefix, hosts=[], fileprefix='', abortOnEnd=True, lastCore=True,
            rankfile=''):
        Composite.__init__(self)

        netcdfwriterrun=None
        cmdline=None

        nw_cmd = os.getenv('NETCDF_WRITER_EXE') or '$PARFLOW_DIR/bin/netcdf-writer'

        # hosts_list: convert hosts to a list
        hosts_list = hosts.split(',')

        # hosts: string with host names, separated by spaces
        if hosts.find(',') == -1:
            # only one host. start locally!
            #nw_cmd = '$PARFLOW_DIR/bin/netcdf-writer'
            cmdline = '%s %s %s' % (nw_cmd, fileprefix, '--no-abort' if not abortOnEnd else '')
        else:
            # load runargs from file...
            command = ['bash', '-c', '. $PARFLOW_DIR/config/pf-cmake-env.sh && echo \"$MPIEXEC_PREFLAGS $MPIEXEC_POSTFLAGS\"']

            proc = subprocess.Popen(command, stdout = subprocess.PIPE)
            mpirunargs = '--report-bindings'  # -mca plm_rsh_no_tree_spawn 1 '
            if rankfile == '' and lastCore:
                rankfile = os.getcwd() + '/rankfile_writer.txt'
                with open(rankfile, 'w+') as f:
                    import multiprocessing
                    for i, host in enumerate(hosts_list):
                        bind_to = multiprocessing.cpu_count()-1
                        f.write('rank %d=%s slot=%d\n' % (i, host, bind_to))

                        # bind everything to the last core ;)
                        #rank 0=aa slot=15
                        #rank 1=bb slot=15
                        #rank 2=cc slot=15

            for line in proc.stdout:
                mpirunargs += ' ' + line.decode('utf-8')

            proc.communicate()


            netcdfwriterrun = FlowvrRunSlurm('%s %s %s' %
                    (nw_cmd, fileprefix, '--no-abort' if not abortOnEnd else ''),
                    hosts=hosts,
                    prefix=prefix)




        for i, host in enumerate(hosts_list):
            name = '%s/%d' % (prefix, i)

            netcdfwriter = Module(name, cmdline=cmdline, run=netcdfwriterrun,
                    host=host)
            netcdfwriter.addPort('in', direction = 'in');

            # collect ports
            for pname in netcdfwriter.ports:
                p = netcdfwriter.ports[pname]
                if not pname in self.ports:
                    self.ports[pname] = list()
                self.ports[pname].append(p)

class Analyzer(Module):
    """Module that will analyze grid messages on the in Port and give a Steer Proposition
    on the out port. The log port can be used with the logger module to log and plot
    floating point outputs in real time. If lastCore parameter is True it will be
    pinned on the last CPU core. Thus it is not interfering with the running ParFlow
    instances if there is a core left.

    Input ports:
    - in: receiving grid messages, typically from one of ParFlow's out ports

    Output ports:
    - out: sends action messages (steer messages) to typically ParFlow's in port"""

    def __init__(self, name, cmdline, host='', lastCore=True, cores=''):
        cores = preferLastCore(lastCore, cores)
        Module.__init__(self, name, cmdline = cmdline, host=host, cores=cores)
        self.addPort('in', direction = 'in')
        p = self.addPort('out', direction = 'out')#, blockstate='nonblocking') # send out thaat I need a trigger ;)
        #p.blockstate='non blocking'
        self.addPort('log', direction = 'out')

class Logger(Module):
    """Module that will log the given stamp value in a graph. Used mainly for debugging.
    stampNames: space separated line of stamp names. Must be float stamps.
    The float content of those stamps will be logged and also plotted if showWindows is
    True.

    Input ports:
    - in: receives stamp messages whose float stamps are logged according to the
      stampNames given on module initialization
    """
    def __init__(self, name, stampNames, showWindows=True, host=''):
        Module.__init__(self, name, cmdline = 'python $PARFLOW_DIR/bin/parflowvr/logger.py %s %s' % (stampNames, '--show-windows' if showWindows else ''), host=host)
        self.addPort('in', direction = 'in')  #, messagetype='stamps')

class AbortOnEmpty(Module):
    """Module that will abort the FlowVR application when receiving an empty message.

    Input ports:
    - in: receives messages that when empty will lead to an abort of the FlowVR application
    """
    def __init__(self, name, host=''):
        Module.__init__(self, name,
                cmdline = 'python $PARFLOW_DIR/bin/parflowvr/abort-on-empty.py',
                host=host)
        self.addPort('in', direction = 'in')

class Ticker(Module):
    """Module that sends an empty message every T seconds. Used for debugging and module
    testing.

    Output ports:
    - out: port that sends the message
    """
    def __init__(self, name, size=0, T=0.5, host=''):
        Module.__init__(self, name, cmdline = 'python $PARFLOW_DIR/bin/parflowvr/ticker.py %d %f' % (size, T), host=host)
        #self.addPort('beginIt', direction = 'in')
        self.addPort('out', direction = 'out')
