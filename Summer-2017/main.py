from pylab import *
import numpy as np
from math import *
from NeuronClass import *
import random as random
from matplotlib import *
from FunctionsIFCorrect import *

def main():
    #Parameters
    N_E=5000
    N_I=1000
    J_EE=1.0
    J_IE=1.0
    J_EI=-2.0
    J_II=-1.8
    # Number of connections is made to match up with the overall number of neurons
    K=0.02*(N_E+N_I)

    # Linear IF properties
    Rm=1
    Cm=1
    tau_m=Rm*Cm
    Vth=1
    Vrest=0
    t_0=0
    injectCurrent = 1.5

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
    networkMeanThresholdExc3D=[]
    networkMeanThresholdExc3DSD=[]
    networkMeanThresholdInh3D=[]
    networkMeanThresholdInh3DSD=[]


  # Setting up and recording adaptation parameters:

  #for j in range(10,69,6):
  for j in range(10,11):
      #phi=float(j)/100.0
      phi=0
      #for i in range(30,405,41):
      for i in range (30,31,17):
        lambdaValue=float(i)/1000.0
        lambdaValueList.append(lambdaValue)
        phiValue.append(phi)


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


        '''
        appendMean(MeanExcActivity3D,MeanExcActivity)
        appendStandardDeviation(MeanExcActivity3DSD,MeanExcActivity)
        appendMean(MeanInhActivity3D,MeanInhActivity)
        appendStandardDeviation(MeanInhActivity3DSD,MeanInhActivity)
        appendMean(EIRatio3D,EIRatio)
        appendStandardDeviation(EIRatio3DSD,EIRatio)
        #appendMean(networkMeanThresholdExc3D,networkMeanThresholdExc)
        #appendStandardDeviation(networkMeanThresholdExc3DSD,networkMeanThresholdExc)
        #appendMean(networkMeanThresholdInh3D,networkMeanThresholdInh)
        #appendStandardDeviation(networkMeanThresholdInh3DSD,networkMeanThresholdInh)
        appendMean(MeanExcInput3D,ExcitatoryInput)
        appendStandardDeviation(MeanExcInput3DSD,ExcitatoryInput)
        appendMean(MeanInhInput3D,InhibitoryInput)
        appendStandardDeviation(MeanInhInput3DSD,InhibitoryInput)
        #appendMean(MeanExcInput3D,MeanExcInput)
        #appendStandardDeviation(MeanExcInput3DSD,MeanExcInput)
        #appendMean(MeanInhInput3D,MeanInhInput)
        #appendStandardDeviation(MeanExcInput3DSD,MeaInhInput)
        '''


  MeanExcInput3D=[]
  MeanExcInput3dSD=[]
  MeanInhInput3D=[]
  MeanInhInput3DSD=[]


  plot(EIRatio)
  title("2DEIRatio")
  ylabel("Ratio")
  xlabel("Time")
  show()

  plot(MeanExcActivity)
  title("Mean Excitatory Activity")
  ylabel("Percentage")
  xlabel("Time")

  plot(MeanInhActivity)
  title("Mean Inhibitory Activity")
  ylabel("Percentage")
  xlabel("Time")
  show()

  plot(ExcitatoryInput)
  title("Excitatory Input")
  ylabel("Signal Strength")
  xlabel("Time")

  plot(InhibitoryInput)
  title("Inhibitory Input")
  ylabel("Signal Strength")
  xlabel("Time")
  show()

  plot(singleVoltageHistory)
  title("singleVoltageHistory")
  ylabel("Voltage")
  xlabel("Time")
  show()

main()
