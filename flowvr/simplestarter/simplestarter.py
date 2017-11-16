import flowvr
import time
import sys

outport = flowvr.OutputPort("out")

stampStartTime = outport.addStamp("stampStartTime", float)
stampStopTime = outport.addStamp("stampStopTime", float)
inport = flowvr.InputPort("beginIt")  # really strange that we have to init it here...

v = flowvr.vectorPort()
v.append(inport)
v.append(outport)


simplestarterModule = flowvr.initModule(v)


print ("-------simple starter now waiting!")
while simplestarterModule.wait() :
    print ("-------got beginIt")
    #rm = inport.get()
    #rm.clear()


    # I do not know how to make an own stamplist atm...
    m = flowvr.MessageWrite(outport.stamps)


    startTime = float(sys.argv[1])
    stopTime = float(sys.argv[2])
    print("----setting stamp %f - %f" % (startTime, stopTime))
    m.setStamp("stampStartTime", startTime)
    m.setStamp("stampStopTime",  stopTime)

    # we need to initizlize data also for stamps messages ;)
    m.data = simplestarterModule.alloc(0)
    outport.put(m)

    m.clear()

    print ("-------putting message")
    time.sleep(1)
    break



print ("-----quit simple starter")

simplestarterModule.close()
#simplestarterModule.abort()



