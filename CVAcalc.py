
import numpy as np
import matplotlib.pyplot as plt
'''
Function to Bootstrap Survival probabilities from given CDS spreads
Use simple linear interpolation at semi-annual points 
'''
def funcBootstrapCDS(DiscountFactors, deltaT = 0.5):
    maturity = 5.0
    InputCDS = [0, 141.76,  165.36, 188.56, 207.32, 218.38] #Input CDS in Basis Points
    InputCDS  = 0.0001*np.array(InputCDS) 
        

   #t = np.array(range(6))
    t = np.array(np.arange(0, maturity +deltaT, deltaT))

    S = np.interp(t, range(6), InputCDS) #this step linearly interpolates CDS rates

    
    SumPL = 0.0
    SumDL = 0.0
    RR = 0.4
    
    DF = np.array(DiscountFactors) #pick from DF matrix Interpolation
    P = np.ones(len(t))
    
    for i in range(1,len(t)):
        P[i] = (SumDL - SumPL + ((1-RR)*P[i-1]*DF[i]) )/(S[i]*DF[i]*deltaT + (1-RR)*DF[i]  )
        SumPL = 0.0
        SumDL = 0.0
        if i < 10 :
            for j in range(1,i+1):
                SumPL = SumPL + (S[i+1]*P[j]*DF[j]*deltaT)
                SumDL = SumDL + (1-RR)*(P[j-1] - P[j])*DF[j]
        print("T = %r Survival Probability = %8.7f " %(i , P[i]))
    plt.plot(t, P, 'go-')
    plt.xlabel('Time in years')
    plt.ylabel('Survival Probability of DB')
    plt.title('Plot DB survival probability with time, bootstrapped from CDS')
    plt.show()

#Test - validate that the Survival Probabilities are correct. Both Premium and Default Legs should be equal using calculated Ps for 5 yr CDS
    SumPL = 0.0
    SumDL = 0.0
    for k in range(1, 11):
        SumPL = SumPL + (S[10]*P[k]*DF[k]*deltaT)
        SumDL = SumDL + (1-RR)*(P[k-1] - P[k])*DF[k]
        #print(k, SumPL, SumDL)
    print ("Test: Premium Leg Value = %8.10f and Default Leg Value = %8.10f for 5 years CDS are equal" %(SumPL, SumDL))
    return P




##### Use P to get Lamda = Term structure of Hazard rates

def getHazardRates(P):
    Lamda = np.ones(10)
    deltaT = 0.25
    for i in range(1,11):
        Lamda[i-1] = -np.log(P[i]/P[i-1])/deltaT
        print (i,Lamda[i-1])
    return Lamda

def piecewiseHazard(x, Lamda):
    
    return Lamda[int(x)]

def plotHazardRates():
    P = funcBootstrapCDS()
    Lamda = getHazardRates(P)
    
    time = np.arange(0.025,10., 0.025)

    y = []
    for i in range(len(time)):
        y.append(piecewiseHazard(time[i], Lamda))

    plt.plot(time,y,c='red', ls='', ms=5, marker='.')
    ax = plt.gca()
    ax.set_xlim([0, 10])
    ax.set_ylim([0,0.05])
    plt.xlabel('Time in years')
    plt.ylabel('Hazard Rates of DB')
    plt.title('Term Structure of Hazard Rates for DB - piecewise constant Lambda')

    plt.show()
    return

'''
Calculate CVA by Reimann sum of the theoretical Integral
Hence calculates average of Expected exposure at the 

'''


def calcCVA(ExpectedExposure,SurvivalProb, DiscountFactors, exposureTenors, RR, Tenor = 0.5):

    if not len(ExpectedExposure) == len(exposureTenors):
        print ('Indexes are not of same shape')
    elif not len(DiscountFactors) == len(exposureTenors):
        print ('Indexes are not of same shape')
    elif not len(SurvivalProb) == len(exposureTenors):
        print ('Indexes are not of same shape')
    
    LogDF = np.log(np.array(DiscountFactors)) #this step used for Log Linear Interpolation of Discount Factors


    CVA_period = np.zeros(shape=(len(exposureTenors)-1),dtype=np.float)  #Calculate CVA per period
    for i in range(len(CVA_period)):
        MidExposure = (ExpectedExposure.iloc[i+1] + ExpectedExposure.iloc[i]) * 0.5
        DefaultProbability = SurvivalProb[i] - SurvivalProb[i+1]
        MidTenor = (exposureTenors[i+1] + exposureTenors[i]) * 0.5   #assume that default happens at mid-point of the period on an average
        LogDFMid = np.interp(MidTenor, exposureTenors, LogDF)
        DF = np.exp(LogDFMid)
        CVA_period[i] = MidExposure * DF * (1-RR) * DefaultProbability
        #print i, MidExposure, DefaultProbability, LogDFMid, DF
    print ("CVA for the 5 Year IRS  notional = 1.0 = %8.8f " % sum(CVA_period))
    return CVA_period