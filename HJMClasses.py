# -*- coding: utf-8 -*-
"""
Created on Tue Oct 16 21:06:53 2018

@author: Akhil
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from Numerics import curvatures, evalSpline, trapInteg

'''takes the Fwd Curve F(t,T) from BOE as input and returns 
continuous Interpolated F(t,T) function as output using Cubic Spline between From and To points
Inputs: f = F(t,T) from BOE at discrete points; RatesTenor = Tenors from BOE curve


'''
def getInstantFwdCurve(f, RatesTenor, DateFrom, DateTo, Plot = 0, step = 0.01):
    k = curvatures(RatesTenor,f)
    InstFwds = []
    plottenors = np.arange(DateFrom, DateTo + step, step)
    for i in np.arange(DateFrom, DateTo + step, step):
        InstFwds.append(evalSpline(RatesTenor,  f,k,i))
        
      ##Test the Interpoltion  -the curve should go through interpolated points  
    if Plot ==1:
        plt.plot(plottenors, InstFwds, 'red')
        plt.plot(RatesTenor, f, 'go')
        plt.legend(fontsize=10)
        plt.xlabel('Maturities $ T $')
        plt.ylabel('$F(t,T)$')
        plt.show()
        
    return InstFwds

    '''
Function that returns discretised Fhat (as defined in Glasserman) from initial instantaneous continuous Forward Curve
Input: DiscreteMats are the discretised maturities and tenors which we are using for simulating the full curve
f, RatesTenor - f = F(t,T) from BOE at discrete points; RatesTenor = Tenors from BOE curve
'''
        
def getDiscreteFwdCurves(DiscreteMats,f,RatesTenor, dTau, step=  0.01): #Plot this in a Graph, Also Show in a table along with Initial F(t,T)

    Fhat = np.zeros(len(DiscreteMats))
    
    for l, t in enumerate(DiscreteMats):
        y = getInstantFwdCurve(f, RatesTenor, t, t+dTau, 0, step)
        Fhat[l] = trapInteg(y,step)*(1/dTau)
    return Fhat


class Product(object):

    def __init__(self):
        self.Notional = 1

    def cashflow(self, time):
        return self.Notional


class ZeroCouponBond(Product):

    def __init__(self, maturity):
        super(ZeroCouponBond, self).__init__()
        self.maturity = maturity

    def cashflow(self, time):
        if self.maturity == time:
            return self.Notional
        else:
            return 0.0

    def get_pay_schedule(self, AsofDate = 0.0 ):
        return [self.maturity]

    def __repr__(self):
        return "%sy ZCB" % self.maturity


class IRS(Product):

    def __init__(self, maturity, fixedrate, freq):
        super(IRS, self).__init__()
        self.maturity = maturity
        self.fixedrate = fixedrate
        self.freq = freq
        self.paymentdates = []
        for i in np.arange(1, int(np.ceil(self.maturity / self.freq)) + 1):
            self.paymentdates.append(i * self.freq)

    def get_pay_schedule(self,AsOfDate ):#needed because for CVA residual maturity of IRS decreases
        self.paymentdates = np.arange(AsOfDate + self.freq, self.maturity + self.freq, self.freq)
        return self.paymentdates

    def cashflow(self, libor, time):
        if time in self.paymentdates:
            return (libor - self.fixedrate) * self.Notional * self.freq

    def __repr__(self):
        return "%sy %s%% Pay Fix IRS" % (self.maturity, self.fixedrate * 100)

def calcPriceAnalytic(DFCurve, dt, product): #Returns both Par Swap rate based on Initial Curve as well as PV of Swap for given fixed rate
    fixedpv = 0
    floatingpv = 0
    pay_schedule = product.get_pay_schedule(0.0)
    for idx, item in enumerate(pay_schedule):
        df = DFCurve.loc[item]
        fixedpv += df 
    if hasattr(product, 'freq'):
        floatingpv = 1 - DFCurve.loc[item] #item will be the final payment date
        pv  =  (floatingpv - (fixedpv*product.fixedrate*product.Notional*product.freq))
    else:
        pv = fixedpv*product.Notional
    return pv


def calcPrice(Curve, AsOfDate, dt, product):
    pv = 0
    pay_schedule = product.get_pay_schedule(AsOfDate)
    for idx, item in enumerate(pay_schedule):
        df = Curve.getDiscountFactor(dt, AsOfDate, item) #enter as off date here
        if hasattr(product, 'freq'):
            libor = Curve.getFwdLiborRate(dt, AsOfDate, item)
            pv += df * product.cashflow(libor, item)
        else:
            pv += df * product.cashflow(item)
    return pv



class Curve(object):

    def __init__(self, method,Tenors, input_curve):
        self.method = method
        self.starting_curve = input_curve
        self.drift = None
        self.Tenors = Tenors
        self.forwardCurve = None

 
    
    def createForwardCurve(self, time):
        return self.forwardCurve
    
    def getFwdLiborRate(self, time):
        return print ("error getFwdLiborRate")

    def getDiscountFactor(self, time):
        return print ("error getDiscountFactor")



class continuousHJM(Curve):
    

    
    ''' Function to calculate the Forward Curve under the continuous HJM Model given the Drift and Volatility vectors
    Num factors is the number of independent factors used from PCA. Currently this is 3 but was written to accomodate more factors in future
    '''
    def calcForwardCurve(self, Drift, V1, V2, V3, dfTenor, numfactors = 3):
        dt = 0.01
        dtau = 0.5
        Maturity = 5
        dates = np.arange(0, Maturity+dt, dt)
        numrow = len(dates)
        numcol = len(dfTenor)
        curve = pd.DataFrame(data = np.ones((numrow, numcol)), columns = dfTenor)
        #InitialCurve = inputforwardsBOE.iloc[-1]*.01 #THIS WILL BE REMOVED AS INPUT CURVE IS SEPARATE
        curve.loc[0] = np.array(self.starting_curve)
        for i in range(1,len(dates)):
            # dfbar =  m(t)*dt+SUM(Vol_i*phi*SQRT(dt))+dF/dtau*dt
            f = curve.loc[i-1]
            Phi = np.random.standard_normal(numfactors)
            deltaf = np.ones(len(dfTenor))
            for k, t in enumerate(dfTenor):
                if k == len(dfTenor)-1: deltaf[k] = f.iloc[k] - f.iloc[k-1] #for terminal we do reverse delta
                else: deltaf[k] = f.iloc[k+1] - f.iloc[k]
            dfbar = Drift*dt + (V1* Phi[0] +V2* Phi[1] +V3* Phi[2])*np.sqrt(dt) + deltaf/dtau*dt 
            curverow = curve.loc[i-1] + dfbar
            
            curverow = curverow.clip(lower = 0.0000001)
            #print curverow
            curve.loc[i] = curverow
        curve = curve.set_index(dates)
        self.forwardCurve = curve
        return     #Forward Curve created in the curve object 
    
    def sampleCurvePlots(self):
        ###ALSO Save the Plots 
        '''Cut this into a function to reduce clutter on Main function'''
        self.forwardCurve.plot(kind = 'box')
        plt.xlabel('$Maturities$')
        plt.ylabel('$F(t,T)$')
        plt.title('Distribution of the forward curve for 1 simulation', fontsize=15)
        plt.show()
        
        maturitiesToPlot = [0,1,2,10,15, 20, 25, 30, 40, 50] 
        self.forwardCurve.iloc[:, maturitiesToPlot].plot()
        plt.xlabel('Evolution of time $t$')
        plt.ylabel('$F(t,T)$')
        plt.title('Evolution of forward rates with time for a simulation', fontsize=15)
        plt.show()
        
        for i in np.arange(0.0, 5.5, 0.5):
            plt.plot(self.Tenors, self.forwardCurve.loc[i, :], label = i)
        plt.legend(fontsize=10)
        plt.xlabel('Maturities $ T $')
        plt.ylabel('$F(t,T)$')
        plt.title('Forward Curve shapes with evolution of time from Simulation', fontsize=15)
        plt.show()

    def saveCurve(self):
        self.forwardCurve.to_csv('ContinuousHJMFwds.csv')
    
    def getFwdLiborRate(self, tenor, asofDay, PaymentDate):
        tenor = 0.5   #Means we are looking for 6 month libor rate for every payment
        try:
            row_idx = asofDay +  (PaymentDate - tenor)
            col_idx =  tenor
            instfwdrate = self.forwardCurve.loc[row_idx,col_idx]   #we have only 1 rate over 6 months. we can chose to average over Ti, Ti+1
            fwdLibor = 2*(np.exp(0.5*instfwdrate) -1) 
            return  fwdLibor
        except:
            print ("error getFwdLiborRate") 
            return 

    def getDiscountFactor(self, tenor, asofDay,  PaymentDate, method = 1): #Code for DF(0,0)
        tenor = 0.5
        try:
            if method == 1: #We take DF through the entire forward rate curve Z(0,T) = exp(-sum(f(0,ti)*deltaT)
               row_idx = asofDay
               col_start = 0
               col_end = PaymentDate
               rate =  self.forwardCurve.loc[row_idx, col_start:col_end].sum()*tenor
               DiscountFactor = np.exp(-rate)
               return DiscountFactor
            else: #here we can implement discounting by the short rate
                return -1
        except:
            print ("error getDiscountFactor") 
            return 


class discretisedHJM(Curve):
    
    def __init__(self,  method, Tenors, input_curve,maxMaturity ,DiscreteMats):
        super(discretisedHJM, self).__init__(method,Tenors, input_curve)
        self.LOIS = None
        self.maxMaturity = maxMaturity
        self.DiscreteMats = DiscreteMats

    def saveCurve(self):
        self.forwardCurve.to_csv('DiscreteHJMFwds.csv')
 
    def sampleCurvePlots(self):
        for i in self.Tenors:
            plt.plot(self.forwardCurve.loc[i, 0.0:4.5-i],'o-' ,label = i)
        plt.legend(fontsize=10)
        plt.xlabel('Maturities $ T $')
        plt.ylabel('$F^(ti,Ti)$')
        plt.title('Forward Curve shapes with evolution of time from Simulation - using Discrete HJM formulation', fontsize=15)   
        
    

    
        '''
    Returns one simulated value of the discretised Forward Curve as a Data Frame based on Glasserman's method
    It - 
     - Gets initial Forward Curve discretised at T=0
     - calculates discrete drift and generates simulated values of Fhat for each tenor uptil maturity
     - As tenor increases, number of maturities reduce. Rest are all written -1
     
    Inputs: Calibrated Vol functions V1, V2, V3; dfTenor = Global at which the Tenors were FIT 
    '''
    #Drift not being used in this function
    def calcForwardCurve(self, Drift, V1, V2, V3, dfTenor, numfactors = 3):
        #dt = 0.01
        dtau = 0.5 #dtau could also be 1m, 3m; hence keep it flexible
        
        #DiscreteMats = np.arange(0, self.maxMaturity, dtau) 
        #dates = arange(0, Maturity+dt, dt)
        numrow = len(self.DiscreteMats)
        numcol = len(self.DiscreteMats)
        curve = pd.DataFrame(data =np.ones((numrow, numcol))*-1, columns = self.DiscreteMats)
        curve = curve.set_index(self.DiscreteMats)
        curve.loc[0] = np.array(self.starting_curve)
        s = pd.DataFrame(data = np.column_stack((V1, V2, V3)), columns = range(3), index = dfTenor)  #this is the Vol Matrix
        
        def getDiscreteDrift(Ti): #s is the Data Frame of Vols, 50 x3 ; Col = 0,1,2 Index = Mats ; Ti = 0.5, 1, 1.5...4.5; i = 0 is the initial curve
            remainMats = np.arange(0, self.maxMaturity - Ti , dtau) #0, 0.5, uptil M-i
            Anext = np.zeros(numfactors)
            Aprev = np.zeros(numfactors)
            Bprev = 0.0
            M = np.zeros(len(remainMats))
            for ind,j in enumerate(remainMats):
                Bnext = 0
                for k in range(numfactors):
                    Anext[k] = Aprev[k] + s.loc[j,k]*dtau  
                    Bnext = Bnext + Anext[k]*Anext[k] 
                    Aprev[k] = Anext[k]
                M[ind] = (Bnext - Bprev)/(2*dtau)
                Bprev = Bnext
            return M
        
            
        for i in range(1,len(self.DiscreteMats)):
            Ti = self.DiscreteMats[i] #next discretised time i
            remainMats = np.arange(0, self.maxMaturity - Ti , dtau)
            Phi = np.random.standard_normal(numfactors)
            m = getDiscreteDrift(Ti)
            for ind, j in enumerate(remainMats):
                S = 0        
                S = (s.loc[j,0]* Phi[0] +s.loc[j,1]* Phi[1] +s.loc[j,2]* Phi[2])*np.sqrt(dtau)
                curve.loc[Ti,j] = curve.loc[Ti - dtau,j+dtau] + m[ind]*dtau + S  #KEY Step , relative maturity reduces by dtau, Hence F(mat  = 0) depends on previous F(mat = 0.5)
        self.forwardCurve = curve
        return     #Forward Curve created in the curve object 
        
    ###(self, dt, asofDay, PaymentDate)
    def getFwdLiborRate(self, dtau, asofDay, PaymentDate):
        #dtau = 0.5   #Means we are looking for 6 month libor rate for every payment; 
        #As of Date not really needed as Payment Date is always abolute b/w 0.5 - 5.0
        #if the Fwd rates were discretised monthly/quarterly - this will need to be adjusted
        try:
            row_idx =  (PaymentDate - dtau)
            col_idx =  0.0
            discretefwdrate = self.forwardCurve.loc[row_idx,col_idx]   #
            fwdLibor = (1/dtau)*(np.exp(dtau*discretefwdrate) -1) 
            return  fwdLibor
        except:
            print ("error getFwdLiborRate") 
            return 
    
    ###(self, dt, asofDay,  PaymentDate, method = 1)
    def getDiscountFactor(self, dtau, asofDay,  PaymentDate, method = 2): #Code for DF(0,0)
        #dtau = 0.5 d
        try:
            if method == 2: #here we can implement discounting by the short rate   = exp(-sum(r(t,t)*deltaT) = exp(-sum(F(ti,ti)*deltaT)
               col_idx = 0
               row_start = asofDay
               row_end = PaymentDate - dtau
               rate =  self.forwardCurve.loc[row_start:row_end, col_idx].sum()*dtau
               DiscountFactor = np.exp(-rate)
               return DiscountFactor
            elif method == 3: #Implement OIS Discounting by adjusting LOIS = exp(-sum(rOIS(t,t)*deltaT) = exp(-sum(F(ti,ti)*deltaT) = exp(-sum{(F(ti,ti) - LOIS(dtau)}*deltaT) 
               col_idx = 0    #Short rate picked from first column
               row_start = asofDay
               row_end = PaymentDate - dtau   ##All Indices may need to be rounded in order to check if we are doing this for monthly tenors
               periods = len(self.forwardCurve.loc[row_start:row_end, col_idx])
               rate =  self.forwardCurve.loc[row_start:row_end, col_idx].sum()
               LOISindex = round(dtau, 5)  #Index is checked on rounded value
               LOISAdjustedrate = rate - periods*self.LOIS[LOISindex] 
               #DiscountFactor = np.exp(-rate*dtau)
               DiscountFactor = np.exp(-LOISAdjustedrate*dtau)
               return DiscountFactor
            else:  #OIS Discounting but on Todays curve
                return -1
     
        except:
            print ("error getDiscountFactor") 
            return -1 


   
    def getTodaysDiscountFactorCurve(self, dtau, PaymentDates): #Code for DF(0,0)
        #dtau = 0.5 dtau = 1/12.0
        
        try:
            row_idx = 0.0
            col_start = 0
            DFcurve = []
            Tenors = []
            DFcurve.append(1.0)
            Tenors.append(0.0)
            for tenor in PaymentDates:
                col_end = tenor - dtau
                rate =  self.forwardCurve.loc[row_idx, col_start:col_end].sum()*dtau
                DiscountFactor = np.exp(-rate)
                DFcurve.append(DiscountFactor)
                Tenors.append(tenor)
            DFcurve =  pd.Series(data = np.array(DFcurve), index = Tenors)
            plt.plot(DFcurve)
            
            return DFcurve
    
        except:
            print ("error getDiscountFactor") 
            return -1 

