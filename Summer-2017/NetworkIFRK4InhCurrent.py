from pylab import *
import numpy as np
from math import *
from NeuronClass import *
import random as random
from matplotlib import *
from FunctionsIFInhCurrent import *


def record(aggregateData,parameters):
    for dataVector in aggregateData:
        text_file = open("Data_W/"+ str(dataVector[0])+".txt", "w")
        for item in dataVector[1]:
          text_file.write(str(item) + "\n")
        text_file.close()

    text_file = open("Data_W/Parameters.txt","w")
    for item in parameters:
        text_file.write(str(item) + "\n")
    text_file.close()

def recordRaster(rasterList):
    text_file = open("Data_W/2D/RasterPlot.txt", "w")
    for list in rasterList:
        for item in list:
            text_file.write(str(item) + " ")
        text_file.write("\n")
    text_file.close()

def main():
  #Parameters
  N_E=4000
  N_I=1000
  J_EE=1.0
  J_IE=1.0
  J_EI=-2.0
  J_II=-1.8
  # Number of connections is made to match up with the overall number of neurons
  connectivity=0.02*(N_E+N_I)
  #connectivity=2
  Jmatrix=setupJmatrix(N_E,N_I,J_EE,J_EI,J_IE,J_II,connectivity)
  #print Jmatrix
  #JmatrixNP = np.array(Jmatrix)
  #for i in range(N_E+N_I):
#      print JmatrixNP[:,i]

  # Linear IF properties
  Rm=1
  Cm=1
  tau_m=Rm*Cm
  Vth=1
  Vrest=0
  t_0=0
  injectCurrent = 1.5

  T = 200
  binSize = 20
  dt=0.01

  #2D data


  # Ratio of inputs for a randomly chosen neuron fixed throughout time.
  #EIRatio = []

  # Inputs for the same randomly chosen neuron.
  TotalInput =[]
  ExcitatoryInput = []
  InhibitoryInput = []

  #3D data; all data points collected are time-averaged after each run; mean activity annotated as example.

  # Adaptation variables

  # Decay variable
  lambdaValueList=[]
  # Jum variable
  phiValue=[]

  # Mean population activities; averaged over time after having run an entire run for a 2-tuple of adaptation variables.
  MeanInhActivity3D=[]
  MeanInhActivity3DSD=[]
  MeanExcActivity3D=[]
  MeanExcActivity3DSD=[]
  # Input ratios.
  EIRatio3D=[]
  EIRatio3DSD=[]
  # Inputs
  MeanExcInput3D=[]
  MeanExcInput3DSD=[]
  MeanInhInput3D=[]
  MeanInhInput3DSD=[]

  # Mean threshold data accoording to populations
  networkMeanWExc3D=[]
  networkMeanWExc3DSD=[]
  networkMeanWInh3D=[]
  networkMeanWInh3DSD=[]


  # Setting up and recording adaptation parameters:

  #for j in range(10,69,6):
  for j in range(1,2):
      #phi=float(j)/100.0
      phi=j*0.3
      #for i in range(30,405,41):
      for i in range (1,2):
        lambdaValue= i*0.03
        lambdaValueList.append(lambdaValue)
        phiValue.append(phi)

        #Setting up network structure: appending the list of neurons.
        allNeurons = []
        Vm=[]
        for i in range(N_E):
          Vm= random.random()
          neuron=Neuron('Exc', Rm,Cm,Vth,Vm,Vrest,i)
          allNeurons.append(neuron)
        for i in range(N_I):
          Vm = random.random()*0.7
          neuron=Neuron('Inh',Rm,Cm,Vth*0.7,Vm,Vrest,i+N_E)
          allNeurons.append(neuron)
        #allNeurons=ExcitNeurons+InhibNeurons


        # Start recording the passage of time; starts with t_0 = 0
        t_history=[]
        t_history.append(t_0)
        # Grab a random neuron and record needed parameters
        randIndex=random.random()*len(allNeurons)
        #randNeuron=allNeurons[randIndex]

        #manually setup randNeuron
        randNeuron = allNeurons[0]
        randThresholdHistory=[]

        # Rastor plot related variables.
        RastorPlotList=[]
        RasExc=[]
        RasInh=[]
        RasterList = []

        # Mean Activity variables
        listActExcPercent=[]
        listActInhPercent=[]

        # Gain curve
        gainCurveData=[]

        # The time at which we record spikes
        spikeHistory=[]

        # Mean Network Threshold: this is a 2D data list
        MeanWExc=[]
        MeanWInh=[]

        MeanWExc.append(0)
        MeanWInh.append(0)

        # Average Active Population according to update times.
        MeanExcActivity=[]
        MeanInhActivity=[]

        EIRatio = []

        #time
        Time = []
        Time_Single = []

        #Inputs
        MeanExcInput=[]
        MeanInhInput=[]

        # Record the passage of time
        timeElapsed=0
        # Voltage history for the selected neuron
        singleVoltageHistory=[]
        singleWHistory=[]

        # Count the number of rounds passed
        roundCounter=0

        #Binning. Bin size controls the amount of steps in one bin.
        binIndex = 0
        binExcWChange =0
        binInhWChange =0
        binDataCount = [0]*4


        # Set up a stop time and reinitiate values
        while timeElapsed <T:
            binIndex +=1
            excInput=0
            inhInput=0
            roundCounter+=1

            excWChange = 0
            inhWChange = 0

            timeElapsed+=dt
            # voltageSum = 0

            ######FORMAT: [excCount, inhCount, excInput, inhInput]
            dataCount=[0]*4

            firedList = []


            print "The current system time is at %f."%timeElapsed

            # Main update mechanism: run through each neuron
            for i, neuron in enumerate(allNeurons):

                #neuorn.updateVoltageAndThreshold(dt,injectCurrent,lambdaValue)
                if neuron.population == "Excitatory":
                    excWChange -= neuron.W
                    neuron.updateVoltageAndInhCurrent(dt,injectCurrent,lambdaValue)
                    excWChange += neuron.W
                else:
                    inhWChange -= neuron.W
                    neuron.updateVoltageAndInhCurrent(dt,injectCurrent,lambdaValue)
                    inhWChange += neuron.W

                if i==1:
                    Time_Single.append(timeElapsed)
                    singleVoltageHistory.append(neuron.V_m)
                    singleWHistory.append(neuron.W)
                #Trigger adaptation once the neuron under examination fires
                if neuron.V_m >= neuron.threshold:
                    #reset the voltage
                    neuron.V_m = neuron.Vrest
                    firedList.append(neuron)
                    neuron.W += phi
                    if neuron.population == "Excitatory":
                        excWChange+=phi
                    else:
                        inhWChange+= phi
                    #print "the %d th neuron spiked" %(i)
                    # Examine neighboring neurons
                    data = neighborSpike(neuron,randNeuron,allNeurons,Jmatrix,phi,firedList,injectCurrent,dt)
                    #print "We have our current data[0] at %d." %(data[0])
                    #print "We have our current data[1] at %d." %(data[1])
                    binDataCount[0] += data[0]
                    binDataCount[1] += data[1]
                    binDataCount[2] += data[2]
                    binDataCount[3] += data[3]

            for item in firedList:
                RasterList.append([item.index,timeElapsed])

            binExcWChange += excWChange
            binInhWChange += inhWChange

            if (binIndex == binSize):
                excInput=binDataCount[2]
                inhInput=binDataCount[3]
                excInput += injectCurrent*dt*(binSize)
                binExcWChange += binDataCount[0]*phi
                binInhWChange += binDataCount[1]*phi

                ExcitatoryInput.append(excInput)
                InhibitoryInput.append(inhInput)
                TotalInput.append(excInput+inhInput)

                Time.append(timeElapsed)

                #Mean Activity & Mean Threshold
                MeanExcActivity.append(binDataCount[0]/float(N_E))
                MeanInhActivity.append(binDataCount[1]/float(N_I))
                MeanWExc.append((MeanWExc[-1]*float(N_E) + binExcWChange)/float(N_E))
                MeanWInh.append((MeanWInh[-1]*float(N_I) + binInhWChange)/float(N_I))
                #reset bin
                binIndex = 0
                binExcWChange = 0
                binInhWChange = 0
                binDataCount = [0]*4


        #3D Data collection, after the update loop is finished
        for neuron in allNeurons:
            #print "exc input sum: ", neuron.getExcInputSum()
            #print "inh input sum: ", neuron.getInhInputSum()
            if (neuron.population == "Inhibitory"):
                ratio = float(neuron.getExcInputSum()+(injectCurrent*0.8*timeElapsed)) / neuron.getInhInputSum()
            else:
                ratio = float(neuron.getExcInputSum()+(injectCurrent*timeElapsed)) / neuron.getInhInputSum()
            #ratio = float(neuron.getExcInputSum()/timeElapsed+injectCurrent)/(neuron.getInhInputSum()/timeElapsed)
            EIRatio.append(ratio)

        #Averaged EI ratio for the entire pop.
        avgEIRatio = sum(EIRatio)/len(EIRatio)

        print "Averaged EI ratio: ", avgEIRatio



        appendMean(MeanExcActivity3D,MeanExcActivity)
        appendStandardDeviation(MeanExcActivity3DSD,MeanExcActivity)
        appendMean(MeanInhActivity3D,MeanInhActivity)
        appendStandardDeviation(MeanInhActivity3DSD,MeanInhActivity)

        appendMean(EIRatio3D,EIRatio)
        appendStandardDeviation(EIRatio3DSD,EIRatio)

        appendMean(networkMeanWExc3D,MeanWExc)
        appendStandardDeviation(networkMeanWExc3DSD,MeanWExc)
        appendMean(networkMeanWInh3D,MeanWInh)
        appendStandardDeviation(networkMeanWInh3DSD,MeanWInh)

        # appendMean(MeanExcInput3D,ExcitatoryInput)
        # appendStandardDeviation(MeanExcInput3DSD,ExcitatoryInput)
        # appendMean(MeanInhInput3D,InhibitoryInput)
        # appendStandardDeviation(MeanInhInput3DSD,InhibitoryInput)
        #appendMean(MeanExcInput3D,MeanExcInput)
        #appendStandardDeviation(MeanExcInput3DSD,MeanExcInput)
        #appendMean(MeanInhInput3D,MeanInhInput)
        #appendStandardDeviation(MeanExcInput3DSD,MeaInhInput)




  aggregateData = []
  #2D
  aggregateData.append(["2D/EI_Ratios",EIRatio])
  aggregateData.append(["2D/Exc_Mean_Atv" ,MeanExcActivity])
  aggregateData.append(["2D/Inh_Mean_Atv",MeanInhActivity])
  aggregateData.append(["2D/Exc_Mean_W",MeanWExc[1:]])
  aggregateData.append(["2D/Inh_Mean_W",MeanWInh[1:]])
  aggregateData.append(["2D/Single_Total_Input",TotalInput])
  aggregateData.append(["2D/Single_Exc_Input",ExcitatoryInput])
  aggregateData.append(["2D/Single_Inh_Input",InhibitoryInput])
  aggregateData.append(["2D/Single_W",singleWHistory])
  aggregateData.append(["2D/Single_Voltage",singleVoltageHistory])
  aggregateData.append(["2D/Time_Single",Time_Single])
  aggregateData.append(["2D/Time",Time])

  #3D
  # aggregateData.append(["3D/Phi",phiValue])
  # aggregateData.append(["3D/Lambda",lambdaValueList])
  # aggregateData.append(["3D/EI_Ratios_Mean",EIRatio3D])
  # aggregateData.append(["3D/EI_Ratios_SD",EIRatio3DSD])
  # aggregateData.append(["3D/Exc_EM_activity",MeanExcActivity3D])
  # aggregateData.append(["3D/Inh_EM_activity",MeanInhActivity3D])
  # aggregateData.append(["3D/Exc_EM_SD",MeanExcActivity3DSD])
  # aggregateData.append(["3D/Inh_EM_SD",MeanInhActivity3DSD])
  # aggregateData.append(["3D/Inh_MeanW",networkMeanWInh3D])
  # aggregateData.append(["3D/Exc_MeanW",networkMeanWExc3D])
  # aggregateData.append(["3D/Inh_WSD",networkMeanWInh3DSD])
  # aggregateData.append(["3D/Exc_WSD",networkMeanWExc3DSD])

  #parameters
  parameters = []
  parameters.append("N_E: " + str(N_E))
  parameters.append("N_I: " + str(N_I))
  parameters.append("K: " + str(connectivity))
  parameters.append("External Drive: " + str(injectCurrent))
  parameters.append("dt: " + str(dt))
  parameters.append("T: " + str(T))
  parameters.append("Bin size: " + str(binSize))


  record(aggregateData,parameters)
  recordRaster(RasterList)



main()
