from pylab import *
from numpy import *
from math import *
from NeuronClass import *
from random import *
from matplotlib import *
from FunctionsIFWithAdaptation import *
def main():
  #Setting up parameters
  N_E=800
  N_I=200
  J_EE=1.0
  J_IE=1.0
  J_EI=-2.0
  J_II=-2.0
  connectivity=0.02*(N_E+N_I)
  #connectivity=2
  Jmatrix=setupJmatrix(N_E,N_I,J_EE,J_EI,J_IE,J_II,connectivity)

  #LIF properties
  Rm=1
  Cm=1
  tau_m=Rm*Cm
  Vth=1
  Vrest=0
  t_0=0
  decayConstant=0.2
  phi=0
  #Initialize neurons
  ExcitNeurons=[]
  InhibNeurons=[]
  Vm=[]
  for i in range(N_E):
    Vm=(randrange(10000000))/10000000.00
    neuron=Neuron('Exc', Rm,Cm,Vth,Vm,Vrest,i)
    ExcitNeurons.append(neuron)
  for i in range(N_I):
    Vm=(randrange(7000000))/10000000.00
    neuron=Neuron('Inh',Rm*0.9,Cm,Vth*0.7,Vm,Vrest,i)
    InhibNeurons.append(neuron)

  allNeurons=ExcitNeurons+InhibNeurons

  t_history=[]
  t_history.append(t_0)



  updateTimes=1000
  randIndex=randrange(len(allNeurons))
  randNeuron=allNeurons[randIndex]
  randThresholdHistory=[]
  RastorPlotList=[]
  ISI=[]
  listActExcPercent=[]
  listActInhPercent=[]
  gainCurveData=[]
  gainCurveDataNoBin=[]
  spikeHistory=[]
  currentHistory=[]
  RasExc=[]
  RasExcNoBin=[]
  RasInh=[]
  RasInhNoBin=[]
  for i in range(2,3):
    currentHistory.append(i)
    ratioActive=[]
    #injectCurrent=1*sqrt(connectivity)*i/2 ##Where 1 is a constant
    injectCurrent=1.1

    t_0=0
    resetBin=0
    binWidth=0.3
    nextBin=resetBin+binWidth
    ExcInputData=[]
    InhInputData=[]
    TotalInputData=[]

    ExcMeanActivity=[]
    ExcMeanActivityNoBin=[]
    InhMeanActivity=[]
    InhMeanActivityNoBin=[]
    thresholdHistory=[]
    randNeuron.set_OnceExc(0)
    randNeuron.set_OnceInh(0)
    ExcMeanActivityCounter=0
    InhMeanActivityCounter=0
    binCount=0

    for i in range(updateTimes):
      originalVm=randNeuron.getVm()
      print "o"*100
      print "In main. injectCurrent=%f" %injectCurrent
      #t_spike,index=NextUpdateNeuronNewton(t_0,allNeurons,injectCurrent,decayConstant,Vth)
      t_spike,index=NextUpdateNeuron(t_0,allNeurons,injectCurrent)
      timeDifference=t_spike-t_0
      updateNeurons(allNeurons, timeDifference, injectCurrent,Vth,decayConstant,phi)
      neuronFiredList=[]


      data=neighborSpike(Jmatrix,index,allNeurons,connectivity,injectCurrent,neuronFiredList,randNeuron,ExcMeanActivityCounter,InhMeanActivityCounter,decayConstant,phi,Vth)
      threshold=randNeuron.getThreshold()
      ExcMeanActivityCounter=data[0]
      InhMeanActivityCounter=data[1]
      RasThisRound=neuronFiredList
      ExcCounter=0
      InhCounter=0
      for neuron in neuronFiredList:
          if neuron.getPop()=="Excitatory":
              ExcCounter+=1
          elif neuron.getPop()=="Inhibitory":
              InhCounter+=1
      ExcMeanActivityNoBin.append(ExcCounter)
      InhMeanActivityNoBin.append(InhCounter)

      if t_spike > nextBin:
          ExcThisBin=randNeuron.get_OnceExc()/binWidth+injectCurrent
          InhThisBin=randNeuron.get_OnceInh()/binWidth

          ExcInputData.append(ExcThisBin)
          InhInputData.append(InhThisBin)

          TotalInputData.append(ExcThisBin+InhThisBin)
          resetBin=nextBin
          nextBin=nextBin+binWidth
          ExcMeanActivity.append(float(ExcMeanActivityCounter)/float(N_E))
          InhMeanActivity.append(float(InhMeanActivityCounter)/float(N_I))


          randNeuron.set_OnceExc(0)
          randNeuron.set_OnceInh(0)

          ExcMeanActivityCounter=0
          InhMeanActivityCounter=0

          binCount+=1
      '''
      for neuron in RasThisRound:
          if neuron.getPop()=="Excitatory":
              ExcMeanActivityCounter+=1
          elif neuron.getPop()=="Inhibitory":
              InhMeanActivityCounter +=1
      '''
      t_0=t_spike
      print t_0
      spikeHistory.append(t_0)
      thresholdHistory.append(threshold)
      print "This is round %d." %i

      ###Stuff to do with collecting excitatory input and inhibitory input by different time bins

      #EI ratio stuff
      ExcSumList = []
      InhSumList = []

      for neuron in allNeurons:
        ExcSumList.append(neuron.getExcInputSum()/t_0+injectCurrent)
        InhSumList.append(neuron.getInhInputSum()/t_0)

      EIRatioList = []
      for i in range(len(ExcSumList)):
        if InhSumList[i]!=0:
          EIRatioList.append(float(ExcSumList[i])/InhSumList[i])



    sumRatio=0
    for item in ExcMeanActivity:
      sumRatio+=item
    averageOverTimes=sumRatio/len(ExcMeanActivity)/t_0
    RasExc.append(averageOverTimes)

    sumRatio=0
    for item in InhMeanActivity:
        sumRatio+=item
    averageOverTimes=sumRatio/len(InhMeanActivity)/t_0
    RasInh.append(averageOverTimes)

    sumRatio=0
    for item in ExcMeanActivityNoBin:
        sumRatio+=item
    averageOverTimes=sumRatio/len(ExcMeanActivityNoBin)/t_0
    RasExcNoBin.append(averageOverTimes)
    sumRatio=0
    for item in InhMeanActivityNoBin:
        sumRatio+=item
    averageOverTimes=sumRatio/len(InhMeanActivityNoBin)/t_0
    RasInhNoBin.append(averageOverTimes)
  print RasExc
  print RasInh

  plot(spikeHistory,thresholdHistory)
  title("Threshold against Update Times")
  ylabel("Threshold (mV)")
  xlabel("Update Time (mSec)")
  show()


  plot(ExcInputData)
  title("Exc Input History")
  ylabel("Input Strength (mV)")
  xlabel("Update Bins")

  plot(InhInputData)
  title("Inh Input History")
  ylabel("Input Strength (mV)")
  xlabel("Update Bins")
  show()

  plot(RasExcNoBin)
  title("Exc Mean Activity In Different Current Strengths")
  ylabel("Mean Activity Ratio")
  xlabel("Current Strength")
  plot(RasInhNoBin)
  title("Inh Mean Activity In Different Current Strengths")
  ylabel("Mean Activity Ratio")
  xlabel("Current Strength")
  show()


  plot(RasExc)
  title("Mean Activity In Different Current Strengths")
  ylabel("Mean Activity Ratio")
  xlabel("Current Strength")

  plot(RasInh)
  title("Mean Activity In Different Current Strengths")
  ylabel("Mean Activity Ratio")
  xlabel("Curent Strength")
  show()

  """

    #Editing the format of the print for the Rastor Plots as in description:
  RasDescriptList=RastorPlotList
  for i in range(len(RastorPlotList)):
    for j in range(len(RastorPlotList[i])):
      if RastorPlotList[i][j]>= N_E:
        RasDescriptList[i][j]="Inh: %d" %(RastorPlotList[i][j] -N_E+1)
      elif RastorPlotList[i][j]<N_E:
        RasDescriptList[i][j]="Exc: %d" %(RastorPlotList[i][j]+1)
  #print RasDescriptList
    #Draw actual Raster Plots
  print "--------------------->"
  #print RastorPlotList
  numExcFired=0
  numInhFired=0
  for i in range(len(RastorPlotList)):
    for j in range(len(RastorPlotList[i])):
      if RastorPlotList[i][j][:3]=="Exc":
        numExcFired+=1
      elif RastorPlotList[i][j][:3]=="Inh":
        numInhFired+=1
  print "%d excitatory neurons fired; %d inhibitory neurons fired." %(numExcFired,numInhFired)

  plot(RastorPlotList)
  title("Rastor Plot")
  ylabel("Indices of Firing Neurons")
  xlabel("Spike Times")
  show()"""

  plot(ExcMeanActivity)
  title("Percentage of Active Excitatory Neurons Per Update")
  ylabel("Percentage")
  xlabel("Update Bins")
  show()

  plot(InhMeanActivity)
  title("Percentage of Active Inhibitory Neurons Per Update")
  ylabel("Percentage")
  xlabel("Update Times")
  show()
  connectionCounter=0

  scatter(ExcSumList,InhSumList)
  title("EI Ratios")
  ylabel("I Input")
  xlabel("E Input")
  show()

  savetxt("ExcSumList",ExcSumList)
  savetxt("InhSumList",InhSumList)
  savetxt("EIRatioList", EIRatioList)


  """for i in range(N_E+N_I):
      for j in range(N_E+N_I):
          if Jmatrix[i][j]!=0:
              connectionCounter+=1
  print "We count %d number of connections."%connectionCounter

    #print "We have the Rastor Plot information:"
    #print RasDescriptList
    #listConnected=IFNetwork.getConnected(randNeuron)
    #print "TotalExcitatoryInput:"
    #print randNeuron.getTotalExcitatoryInput()
    #print "TotalInhibitoryInput"
    #print randNeuron.getTotalInhibitoryInput()


    totalExcitatoryInput=randNeuron.getTotalExcitatoryInput()
    totalInhibitoryInput=randNeuron.getTotalInhibitoryInput()
    EIRatios=[]
    for i in range(len(totalExcitatoryInput)):
      if totalInhibitoryInput[i]==0:
        ratio="NaN"
      else:
        ratio=totalExcitatoryInput[i]/totalInhibitoryInput[i]
      EIRatios.append(ratio)

  print "-------------------------> \n Gain Curve Data"
  print gainCurveData
  #print "EI Ratios: \n"
  #print EIRatios
  #print "InterSpikeIntervals: \n"
  #print ISI
  #print "Excitatory Mean Activity:\n"
  #print listActExcPercent
  #print "Inhibitory Mean Activity:\n"
  #print listActInhPercent
  #print RastorPlotList
  #plot(randNeuronVoltageHistory)
  #title('Membrane potential history for a randomly chosen neuron')
  #ylabel('Membrane Potential (V)')
  #xlabel('Time (msec)')
  #show()
  #Plotting TotalInput for this random Neuron:
  plot(randNeuronTotalInputHistory)
  title('Total Input for a randomly chosen neuron')
  ylabel('Membrane Potential (V)')
  xlabel('Time (msec)')
  ylim([-1,1])
  show()
    # print listInhConnected
  plot(randNeuron.getTotalExcitatoryInput())
  title('Total Excitatory Input for a randomly chosen neuron')
  ylabel('Membrane Potential (V)')
  xlabel('Time (msec)')
  ylim([-1,1])
  show()
"""
"""

      #t_history.append(t_0)
  plot(randNeuron.getTotalInhibitoryInput())
  title('Total Inhibitory Input for a randomly chosen neuron')
  ylabel('Membrane Potential (V)')
  xlabel('Time (msec)')
  ylim([-1,1])
  show()

  plot(ISI)
  title('Interspike Intervals')
  ylabel('Interval Length (msec)')
  xlabel('Experiment Time (msec)')
  ylim([-1,1])
  show()

  plot(listActExcPercent)
  title('Mean Activity Per Time')
  ylabel('Interval Length (msec)')
  xlabel('Experiment Time (msec)')
  ylim([0,0.05])

  plot(listActInhPercent)
  title('Mean Activity Per Time')
  ylabel('Interval Length (msec)')
  xlabel('Experiment Time (msec)')
  ylim([0,0.05])

  show()

  #print "This run: \n Neurons: %d \n Update Times: %d \n I=%f.\n Update History:"%(N_E+N_I,updateTimes,injectCurrent)
  #print t_history
  #print "Random Neuron Voltage History:"
  #print randNeuronVoltageHistory


  plt.semilogy(, P, 'b')
  plt.semilogy(t, Pb, 'r')
  plt.axvline(10, color='k', linestyle='solid')
  plt.xlabel('Time [s]')
  plt.ylabel('index')
  plt.grid()
  plt.show()
"""
main()
