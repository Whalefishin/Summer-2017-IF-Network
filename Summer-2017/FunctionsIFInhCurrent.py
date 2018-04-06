from pylab import *
from numpy import *
from math import *
from NeuronClass import *
import random as random
from matplotlib import *

import collections

"""
Function: find the next update Neuron
Algorithm:
  1. Set up an arbitrarily large t_spike
  2. Compute the next spiking time for each neuron
  3. For every spike time, compare with the t_spike
  4. If spike time is less than t_spike, set t_spike to be this spike time
  5. This should return a smallest spike time associated with a neuron of a certain index; return this index to locate the neuron
"""


###########NEWTON#########

def NextUpdateNeuronNewton(t_0,allNeurons,injectCurrent,lbda,originalThresh):
  deltaT=10000000000
  index=0
  for i in range(len(allNeurons)):
    currentNeuron=allNeurons[i]
    Vth=currentNeuron.getThreshold()
    Vrest=currentNeuron.getVrest()
    V_0=currentNeuron.getVm()
    Rm=currentNeuron.getRm()
    tau_m=Rm*(currentNeuron.getCm())
    population=currentNeuron.getPop()
    if population=="Inhibitory":
        I=injectCurrent*0.8
        originalT=originalThresh*0.7
    elif population=="Excitatory":
        I=injectCurrent
        originalT=originalThresh
    #t_sp=t_0+log(float(Vrest+I*Rm-originalT)/( (Vth-originalT)*exp(-lbda) -(V_0-Vrest-I*Rm)*exp(-1.000/tau_m    )))
    deltaTNewton= newton(20,3,lbda,tau_m, Vth, originalT,V_0,Vrest,injectCurrent,Rm)
    if deltaT>=deltaTNewton:
      deltaT=deltaTNewton
      index=i
    t_spike=t_0+deltaT
  print "Our computed spike time is %f." %t_spike
    #instVoltage=(Vrest+injectCurrent*Rm)+(V_0-Vrest-injectCurrent*Rm)*exp(-(t_spike-t_0)/tau_m)
    #currentNeuron.setV_m(instVoltage)
    #currentNeuron.setV_m(instVoltage)
  descript=(t_spike,index)

  return descript


#############ANALYTICAL###############
def NextUpdateNeuron(t_0,allNeurons,injectCurrent):
  t_spike=100000
  index=0
  for i in range(len(allNeurons)):
    currentNeuron=allNeurons[i]
    Vth=currentNeuron.getThreshold()
    Vrest=currentNeuron.getVrest()
    V_0=currentNeuron.getVm()
    Rm=currentNeuron.getRm()
    tau_m=Rm*(currentNeuron.getCm())
    population=currentNeuron.getPop()
    omegaValue=currentNeuron.get_omegaValue()
    print "Current omega value is %f." %(omegaValue)
    if population=="Inhibitory":
        I=injectCurrent*0.8
    elif population=="Excitatory":
        I=injectCurrent
    t_sp=t_0-tau_m*log(float(Vth-Vrest-(I-omegaValue)*Rm)/float(V_0-Vrest-(I-omegaValue)*Rm))

    #newton's method
    #t_sp = newton(100,t_0+0.1,)
    if t_spike>=t_sp:
      t_spike=t_sp
      index=i
    #instVoltage=(Vrest+injectCurrent*Rm)+(V_0-Vrest-injectCurrent*Rm)*exp(-(t_spike-t_0)/tau_m)
    #currentNeuron.setV_m(instVoltage)
    #currentNeuron.setV_m(instVoltage)
  descript=(t_spike,index)

  return descript



def setupJmatrix(N_E,N_I,J_EE,J_EI,J_IE,J_II, K):
    Jmatrix=[]
    for i in range (N_E+N_I):
      Jmatrix.append([])
    for i in range(N_E+N_I):
      for j in range(N_E+N_I):
        Jmatrix[i].append(0)

    for i in range(N_E+N_I):
      for j in range(N_E+N_I):
        if i != j:
          indicValue = random.random()
          #print indicValue
          if i<N_E and j < N_E:
            if indicValue <= float(K)/N_E:
              Jmatrix[i][j]=J_EE/sqrt(K)
          elif i>=N_E and j< N_E:
            if indicValue <= float(K)/N_E:
              Jmatrix[i][j]=J_IE/sqrt(K)
          elif i<N_E and j>=N_E:
            if indicValue <= float (K)/N_I:
              Jmatrix[i][j]=J_EI/sqrt(K)
          elif i>=N_E and j >= N_E:
            if indicValue<= float(K)/N_I:
              Jmatrix[i][j]=J_II/sqrt(K)
    return Jmatrix


def updateNeurons(allNeurons,timeDifference,injectCurrent,Vth,lbda,phi,decayConstant):
  counter=0
  listNegative=[]
  for neuron in allNeurons:
    population=neuron.getPop()
    if population=="Inhibitory":
        I=injectCurrent*0.8
        originalThresh=Vth*0.7
    elif population=="Excitatory":
        I=injectCurrent
        originalThresh=Vth
    omegaValue=neuron.get_omegaValue()
    V_0=neuron.getVm()
    V_rest=neuron.getVrest()
    Rm=neuron.getRm()
    tau=Rm*neuron.getCm()
    neuron.setV_m(V_rest+(I-omegaValue)*Rm+(V_0-V_rest-(I-omegaValue)*Rm)*exp(-timeDifference/decayConstant))
    newOmegaValue=omegaValue*exp(-timeDifference/tau)
    neuron.set_omegaValue(newOmegaValue)
    V_oldT=neuron.getThreshold()
    newThreshold=(V_oldT -originalThresh)*exp(-lbda*timeDifference)+originalThresh
    neuron.setThreshold(newThreshold)
    #print "This neuron originally has threshold at %f; currently has threshold %f." %(V_oldT,newThreshold)


"""
def neighborSpike(Jmatrix,index,allNeurons,K,injectCurrent,listFiredNeurons,randNeuron,ExcMeanActivityCounter,InhMeanActivityCounter,lbda,phi,Vth):
    ourNeuron=allNeurons[index]
    print "Neuron %d has fired!" %index
    population=ourNeuron.getPop()
    if population=="Excitatory":
        originalThresh=Vth
    elif population=="Inhibitory":
        originalThresh=Vth*0.7
    V_oldT=ourNeuron.getThreshold()
    newThreshold=originalThresh+(V_oldT-originalThresh+phi*originalThresh)
    ourNeuron.setThreshold(newThreshold)
    if population=="Excitatory":
        ExcMeanActivityCounter+=1
    else:
        InhMeanActivityCounter+=1
    listFiredNeurons.append(ourNeuron)
    ourNeuron.setV_m(0)
    for i in range(len(allNeurons)):
      if Jmatrix[i][index]!=0:
        neuron=allNeurons[i]
        V_0=neuron.getVm()
        V_rest=neuron.getVrest()
        Rm=neuron.getRm()
        tau=Rm*neuron.getCm()

        if neuron not in listFiredNeurons:
          neuron.setV_m(neuron.getVm()+Jmatrix[i][index])
          #EI ratio stuff

          if Jmatrix[i][index] >=0:
            neuron.addExcInputSum(Jmatrix[i][index])
          else:
            neuron.addInhInputSum(Jmatrix[i][index])

        if neuron==randNeuron:
          if Jmatrix[i][index] >=0:
            neuron.set_OnceExc(neuron.get_OnceExc()+Jmatrix[i][index])
            print "*"*50
            print "Neuron %d is sending selected neuron signal." %index
          else:
            neuron.set_OnceInh(neuron.get_OnceInh()+Jmatrix[i][index])
            print "*"*50
            print"Neuron %d is sending selected neuron signal." %index

        if neuron.getVm()>=neuron.getThreshold():
          index=allNeurons.index(neuron)
          neighborSpike(Jmatrix,index,allNeurons,K,injectCurrent, listFiredNeurons,randNeuron,ExcMeanActivityCounter,InhMeanActivityCounter,lbda,phi,Vth)
    data=[]
    data.append(ExcMeanActivityCounter)
    data.append(InhMeanActivityCounter)
    return data

"""

def reset(allNeurons):
  for neuron in allNeurons:
    if neuron.getVm()>=neuron.getThreshold():
      neuron.setV_m(0)


def restNeuronUpdate(injectCurrent,firedList,updatedList,allList,timeDifference):
  for neuron in allList:
    if (neuron not in firedList) and (neuron not in updatedList):
      V_rest=neuron.getVrest()
      V_0=neuron.getVm()
      Rm=neuron.getRm()
      tau=Rm*neuron.getCm()
      neuron.setV_m((V_rest+injectCurrent*Rm)+(V_0-V_rest-injectCurrent*Rm)*exp(-timeDifference/tau))


def neighborSpike(startNeuron,neuronToRecord,allNeurons,Jmatrix,phi,firedList,injectCurrent,dt):
    excCount, inhCount, excInput, inhInput = 0,0,0,0
    current = startNeuron
    #queue is a special kind of list
    queue = collections.deque()
    queue.append(current)

    #while queue is not empty
    while queue:
        #set the first neuron in the queue as current
        current = queue.popleft()
        #get the current neuron's neighbors, then iterate through each neighbor
        neighbors, signals = getNeighbors(current,allNeurons,Jmatrix)
        for index, neigh in enumerate(neighbors):
            signal = signals[index]
            #print "sent signal to %d th neuron" %(neigh.index)
            #print "signal is %f" %(signal)
            if neigh not in firedList:
                #For EI ratio
                if signal > 0:
                    #neigh.addExcInputSum(signal/dt + injectCurrent)
                    neigh.addExcInputSum(signal)
                else:
                    neigh.addInhInputSum(signal)
                #For total input of one neuron
                if neigh == neuronToRecord:
                    if signal >0:
                        excInput += signal
                        #excInput += signal/dt + injectCurrent
                    else:
                        inhInput += signal

                fired = update(neigh,signal,phi)
                if fired:
                    #print "the %d th neuron fired in neighborSpike" %(neigh.index)
                    if neigh.population == "Excitatory":
                        excCount+=1
                    else:
                        inhCount+=1
                    firedList.append(neigh)
                    queue.append(neigh)

    #return all the data collected
    data = [excCount,inhCount,excInput,inhInput]
    return data


def getNeighbors(neuron,allNeurons,Jmatrix):
    neighbors,signals = [],[]
    startIndex = allNeurons.index(neuron)
    for i,neighbor in enumerate(allNeurons):
        if Jmatrix[i][startIndex]!=0:
            neighbors.append(neighbor)
            signal = Jmatrix[i][startIndex]
            signals.append(signal)

    return neighbors, signals


def update(neuron,signal,phi):
    fired = False
    #if neuron.V_m +signal >=0:
    neuron.V_m += signal
    if neuron.V_m >= neuron.threshold:
        neuron.V_m = neuron.Vrest
        neuron.W += phi
        fired = True

    return fired



def computeMean(vector):
    return sum(vector)/len(vector)

    
def appendMean(summary,data):
    Mean=computeMean(data)
    summary.append(Mean)

def computeStandardDeviation(vector):
    summation=0
    mean = computeMean(vector)
    size = float(len(vector))
    for item in vector:
        summation+= (item-mean)**2
    returnValue=(summation/size)**.5
    return returnValue

def appendStandardDeviation(summary,data):
    StandardDeviation=computeStandardDeviation(data)
    summary.append(StandardDeviation)



"""
def newton(n,guess,lbda,tau, V_oldT, originalThreshold,V_0,Vrest,injectCurrent,Rm):
    x_0 = guess
    for i in range(n):
        x_0 = x_0 - float((V_oldT-originalThreshold)*exp(-lbda*x_0)-(V_0-Vrest-injectCurrent*Rm)*exp(-float(x_0)/tau)-Vrest-injectCurrent*Rm+originalThreshold)/(-lbda*(V_oldT-originalThreshold)*exp(-lbda*x_0)+(1.0/tau)*(V_0-Vrest-injectCurrent*Rm)*exp(-float(x_0)/tau))
    return x_0
"""
