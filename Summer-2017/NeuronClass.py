"""
This file defines the class of neurons.

Sida Li
"""

class Neuron:
    def __init__(self,descript,Rm,Cm,Vth,V_m,Vrest,\
        index):
        if descript=='Exc':
            self.population="Excitatory"
        elif descript=='Inh':
            self.population="Inhibitory"
        self.Rm=Rm
        self.Cm=Cm
        self.threshold=Vth
        # self.spike=V_spike
        # self.spikeTo=V_T
        self.totalInput=[]
        self.totalExcitatoryInput=[]
        self.totalInhibitoryInput=[]
        self.tau=Rm*Cm
        self.timeToBeUpdated=[]
        self.V_m=V_m
        self.Vrest=Vrest
        self.index=index
        self.OnceInhInput=None
        self.OnceExcInput=None
        self.originalThresh=Vth
        #EI ratio plot
        self.ExcInputSum=0
        self.InhInputSum=0
        #self.omegaValue=omegaCoefficient*Vth
        self.W = 0

    def getIndex(self):
        return self.index

    def getTau(self):
        return self.tau
    def get_originalThresh(self):
        return self. originalThresh
    def set_OnceInh(self,onceInh):
        self.OnceInhInput=onceInh

    def get_omegaValue(self):
      return self.omegaValue

    def set_omegaValue(self,omegaValue):
      self.omegaValue=omegaValue

    def set_OnceExc(self,onceExc):
        self.OnceExcInput=onceExc

    def get_OnceInh(self):
        return self.OnceInhInput

    def get_OnceExc(self):
        return self.OnceExcInput

    def setV_m(self,Vm):
        self.V_m=Vm

    def getPop(self):
        return self.population

    def getThreshold(self):
        return self.threshold

    def setThreshold(self,threshold):
        self.threshold=threshold

    def getVrest(self):
        return self.Vrest

    def getVm(self):
        return self.self.V_m

    def getRm(self):
        return self.Rm

    def getCm(self):
        return self.Cm

    def addTotalInput(self, totalInput):
        (self.totalInput).append(totalInput)

    def addTotalExcitatoryInput(self,totalExcitatoryInput):
        (self.totalExcitatoryInput).append(totalExcitatoryInput)

    def addTotalInhibitoryInput(self,totalInhibitoryInput):
        (self.totalInhibitoryInput).append(totalInhibitoryInput)

    def setTimeToBeUpdated(self,timeToBeUpdated):
        self.timeToBeUpdated.append(timeToBeUpdated)

    def getTotalExcitatoryInput(self):
        return self.totalExcitatoryInput

    def getTotalInhibitoryInput(self):
        return self.totalInhibitoryInput

    def getTotalInput(self):
        return self.totalInput

    def getExcInputSum(self):
        return self.ExcInputSum

    def getInhInputSum(self):
        return self.InhInputSum

    def addInhInputSum(self, new):
        self.InhInputSum += new

    def addExcInputSum(self,new):
        self.ExcInputSum += new


    def updateVoltageAndThreshold(self,dt,injectCurrentExternal,lambdaValue):
        injectCurrent = injectCurrentExternal
        if self.population == "Inhibitory":
            injectCurrent = injectCurrentExternal*0.8
        k_1=(-(self.V_m -self.Vrest)+injectCurrent)/self.tau
        k_2=(-(self.V_m+k_1/2.0*dt -self.Vrest)+injectCurrent)/self.tau
        k_3=(-(self.V_m+k_2/2.0*dt-self.Vrest)+injectCurrent)/self.tau
        k_4=(-(self.V_m + k_3*dt-self.Vrest)+injectCurrent)/self.tau
        self.V_m=self.V_m+dt/6.0*(k_1+2.0*k_2+2.0*k_3+k_4)

        k_1_t=dt*lambdaValue*(self.originalThresh-self.threshold)
        k_2_t=dt*lambdaValue*(self.originalThresh-(self.threshold + k_1_t/2.0))
        k_3_t=dt*lambdaValue*(self.originalThresh-(self.threshold + k_2_t/2.0))
        k_4_t=dt*lambdaValue*(self.originalThresh-(self.threshold+k_3_t))
        self.threshold=self.threshold+1.0/6*(k_1_t+2*k_2_t+2*k_3_t+k_4_t)


    def updateVoltageAndInhCurrent(self,dt,injectCurrentExternal,inverseLambda):
        lambdaValue = 1.0/inverseLambda
        injectCurrent = injectCurrentExternal
        if self.population == "Inhibitory":
            injectCurrent = injectCurrentExternal*0.8
        k_1=(-(self.V_m -self.Vrest)+injectCurrent-self.W)/self.tau
        k_2=(-(self.V_m+k_1/2.0*dt -self.Vrest)+injectCurrent-self.W)/self.tau
        k_3=(-(self.V_m+k_2/2.0*dt-self.Vrest)+injectCurrent-self.W)/self.tau
        k_4=(-(self.V_m + k_3*dt-self.Vrest)+injectCurrent-self.W)/self.tau
        self.V_m=self.V_m+dt/6.0*(k_1+2.0*k_2+2.0*k_3+k_4)

        k_1_w=dt*(-1*self.W / lambdaValue)
        k_2_w=dt*(-1*(self.W+k_1_w/2.0)/lambdaValue)
        k_3_w=dt*(-1*(self.W+k_2_w/2.0)/lambdaValue)
        k_4_w=dt*(-1*(self.W+k_3_w)/lambdaValue)
        self.W=self.W+1.0/6*(k_1_w+2.0*k_2_w+2.0*k_3_w+k_4_w)
        #print "We compute the updated omega value to be %f." %omega_update



    def __str__(self):
        # ourstr="Neuron:\n Population: %s \n total Input: %l \n total ExcitatoryInput: %l \n total Inhibitory Input: %l \n timeToBeUpdated: %l \n delta:%f" %(self.population,self.totalInput, self.totalExcitatoryInput,self.totalInhibitoryInput,self.timeToBeUpdated)
        ourstr="This Neuron:"
        ourstr=ourstr+self.population +"%d"%(self.index+1) +"\nTotal Excitatory Input:"+ ', '.join(str(e) for e in (self.totalExcitatoryInput)) +"\nTotal InhibitoryInput" ', '.join(str(e) for e in (self.totalInhibitoryInput)) +"\nTimes to be Updated:" ', '.join(str(e) for e in (self.timeToBeUpdated))#+"\nVoltage History:"+', '.join(str(e) for e in (self.V_history))

        return ourstr
