import pandas as pd
import numpy as np
import os
from Numerics import calculate_pca
import matplotlib.pyplot as plt
import timeit
from Numerics import trapInteg
from CVAcalc import funcBootstrapCDS, calcCVA
import HJMClasses

'''
Function that predicts the values of the Fitted Vol function using the regression coefficients
'''

def funcVolfunction(intercept, coefficients, T):
    V = intercept
    for i in range(len(coefficients)):
        V = V+ coefficients[i]*pow(T,i+1)
    return V

'''
Function that  checks whether Vol1, Vol2, Vol3 need to be calculated using the fitted regression functions
and returns a value by calling the previous function 
'''    

def V123(T, n):
    if n ==1 :
        try:   #try catch so that code caters to both float inputs as well as list inputs
            V1 = np.ones(len(T))
            V1 = V1*Vol1
        except: 
            V1 = Vol1
        return V1
    elif n ==2:
        return funcVolfunction(Intercept2,Coeff2, T)
    elif n ==3:
        return funcVolfunction(Intercept3,Coeff3, T)
    else: 
        print ("Unknown value of n")
        return

def funcM(Tau):
    dt = 0.01
    steps = int((Tau/dt)) +1
    T = np.linspace(0, Tau, steps)
    
    V1 = V123(T, 1)
    V2 = V123(T, 2)
    V3 = V123(T, 3)
    
         
    intV1 = trapInteg(V1, dt)
    intV2 = trapInteg(V2, dt)
    intV3 = trapInteg(V3, dt)
    
    M1 = V123(Tau,1)*intV1
    M2 = V123(Tau,2)*intV2
    M3 = V123(Tau,3)*intV3
    
    M = M1 + M2 +M3
    return M
    

'''
    function that returns the Risk Neutral drift Vector for the entire Vector of Maturities
'''    
def funcMVector(Tenors):
    M = []
    for i, t in enumerate(Tenors):
        M.append(funcM(t))
    return np.array(M)

# function to generate the Price over multiple MC simulations. 

def getPricefromMC(NSims, product, inputforwards, HJMMethodology = 'continuous' , asofDay = 0, maxMaturity = 5.0,DiscreteMats = None):
     prices = []
     tstart = timeit.default_timer()
     for i in range(0,NSims):
         #print (i)
         if HJMMethodology == 'continuous':
             forwardcurve = HJMClasses.continuousHJM("continuous", dfTenor, inputforwards)
         elif HJMMethodology == 'discrete':
             forwardcurve = HJMClasses.discretisedHJM("discrete", dfTenor, inputforwards, maxMaturity,DiscreteMats)
             
         forwardcurve.calcForwardCurve(Drift, V1, V2, V3, dfTenor, 3)
         price = HJMClasses.calcPrice(forwardcurve, 0.0, dtau, product) # Price at at As of Date T0 = 0 
         if price != -1:
             prices.append(price)
     Price_series = pd.Series(prices)
     timeend = timeit.default_timer()
     SimTime = timeend  - tstart
         
     print("Simulation completed using %r Paths, time taken %3.5f" % (NSims, SimTime) )
 
     return Price_series

##Function to show convergence of MC price series
def MTMConvergence(NumberofSimulations,PriceSeries):
    prices = []
    Simsteps = int(NumberofSimulations/10)
    Sims = range(Simsteps, NumberofSimulations +Simsteps , Simsteps)
    Variance = []
    VarSim = []
    for i in Sims:
        price = PriceSeries[0:i].mean()
        prices.append(price)
        if i > NumberofSimulations/2:
            std = np.std(prices)
            Variance.append(std)
            VarSim.append(i)
           
    plt.plot(Sims, prices, label = 'IRS price by MC')
    plt.xlabel('Number of HJM Simulation')
    plt.ylabel('IRS MTM')
    plt.title('Convergence of IRS Price with Number of Simulations')
    plt.legend(fontsize=10)
    plt.text(1200,0.026,'Analytic Value of IRS is $0.0317$')
    plt.show()
    
    plt.plot(VarSim, Variance, label = 'Standard deviation of MC price')
    plt.xlabel('Number of HJM Simulation')
    plt.title('Std deviation of MC price with Number of Simulations')
    plt.show()
    return 

'''Function to calculate IRS MTM values for multiple number of simulations, say 10,000
At each Simulation, it first simulates the forward Curve and calls the IRSExposures function
returns a Data Frame with each Simulated MTM profile as a row with Columns as the tenors
'''    
def ExposureProfile(NSim, product, inputforwards): #Build Exposure profile on basis of MC simulations. Each simulation = row of Data Frame

    ##Function to calculate the IRS MTM at various payment dates from T = 0 to T= 4.5 for 1 simulated value of the FWD Curve;

    def IRSExposures(forwardcurve,product, dtau = 0.5):
         exposureDates = np.arange(0, product.maturity, product.freq) #from T = 0 to T= 4.5; T = 5 - IRS  MTM = 0
         ExposureValues = []
         for i, asofDate in enumerate(exposureDates):
             #MTM = priceIRS2(FwdCurve, FixedRate, asofDate,Maturity)
             MTM = HJMClasses.calcPrice(forwardcurve, asofDate, dtau, product) # Price at at As of Date T0 = 0 
             ExposureValues.append(MTM)
         Exposure_series = pd.Series(ExposureValues)
         return Exposure_series

    ExposureProfile = pd.DataFrame()
    tstart = timeit.default_timer()
    for i in np.arange(0,NSim):
        FwdCurve = HJMClasses.discretisedHJM("discrete", dfTenor, inputforwards, maxMaturity,DiscreteMats)
        FwdCurve.calcForwardCurve(Drift, V1, V2, V3, dfTenor, 3)
        Exposures = IRSExposures(FwdCurve,product)
        ExposureProfile = ExposureProfile.append([Exposures],ignore_index=True) 
    timeend = timeit.default_timer()
    SimTime = timeend  - tstart
         
    print("Simulation completed using %r Paths, time taken %3.5f" % (NSim, SimTime) )
         
    return ExposureProfile
 
def ConvergenceDiscreteHJM(NumberofSimulations):
    '''Convergence of IRS under Discrete HJM'''
    
    IRSSeries = MTMDataFrame.iloc[:,0 ]
    pricesIRS = []
    Simsteps = int(NumberofSimulations/10)
    Sims = range(Simsteps, NumberofSimulations +Simsteps , Simsteps)
    for i in Sims:
        price = IRSSeries[0:i].mean()
        pricesIRS.append(price)
          
    plt.plot(Sims, pricesIRS, label = 'IRS price by MC')
    plt.xlabel('Number of HJM Simulation')
    plt.ylabel('IRS MTM')
    plt.title('Discrete HJM - Convergence of IRS Price with Number of Simulations')
    plt.legend(fontsize=10)
    analyticPriceIRS = HJMClasses.calcPriceAnalytic(DiscountFactors, 0.5, IRStoPrice)
    SimulatedPriceIRS = IRSSeries.mean()
    print ('Analytic Price of IRS is %8.6F and Simulated Price over 10,000 simulations is %8.6F.' %(analyticPriceIRS,SimulatedPriceIRS))
    
    
def PlotPCAs():
    #Plot the first 3 Eigen vectors (principal componenets to check if they are consistent with our interpretation as level, slope and curvature
    '''Cut this into a function to reduce clutter on Main function'''
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_position('zero')
    plt.plot(dfTenor,b.iloc[:, 0], label='PC1 = Level')
    plt.plot(dfTenor,b.iloc[:, 1], label="PC2 = Slope")
    plt.plot(dfTenor,b.iloc[:, 2], label="PC3 = Curvature")
    plt.legend(fontsize=10)
    plt.title('Principal Components of Forward Curve - Weekly Diff', fontsize=15)
    plt.show() 

    
def PlotFittedVolFunctions():
    v1 = funcVolfunction( Intercept1, Coeff1, X['Tau'])
    plt.plot( X['Tau'] , VolFunc.iloc[:,0], 'o' ,label='Vol1')
    plt.plot( X['Tau'] , v1,label='First Vol function fit using $quadratic Polynomial$ - REJECTED' )
    plt.legend(fontsize=10)
    plt.text(15,.001,'$RSquare = 0.541$')
    plt.title('First Volatility Function fit using quadratic - Rejected', fontsize=12)
    '''Cut this into a function to reduce clutter on Main function'''
    v2 = funcVolfunction( Intercept2, Coeff2, X['Tau'])
    plt.plot( X['Tau'] , VolFunc.iloc[:,1], 'o' ,label='Vol2')
    plt.plot( X['Tau'] , v2, label  = 'Vol2 function fit using Cubic Spline')
    
    v3 = funcVolfunction( Intercept3, Coeff3, X['Tau'])
    plt.plot( X['Tau'] , VolFunc.iloc[:,2], 'o', label = 'Vol3')
    plt.plot( X['Tau'] , v3, label =  'Vol3 function fit using Cubic Spline')
    plt.xlabel('$Tenors$')
    plt.legend(fontsize=10)
    plt.title('Volatility Functions Fit using Cubic Splines', fontsize=15)
    plt.show() 

def PlotDrift():
    plt.plot( dfTenor, Drift, 'bo-')
    plt.xlabel('$Maturity$')
    plt.ylabel('$Risk Neutral Drift$')
    #plt.legend(fontsize=10)
    plt.title('Risk Neutral Drift from continuous HJM', fontsize=15)
    plt.show() 

'''
Main Project thread
'''

os.getcwd()
'''Add the Path to the directory all Data Files are stored'''

path="C:\\Users\Akhil\Documents\CQF\Project"
os.chdir(path)
os.getcwd()


#### Load Input Forward Rates from BOE
inputforwardsBOE = pd.read_excel('InputRatesAll.xlsx', 'fwd curve')  # this file contains last 4 and half years of data from January 2014
inputforwardsBOE = inputforwardsBOE.drop(['Tenor'], axis = 1)
inputforwardsBOE.dropna(inplace=True)
inputforwardsBOE = inputforwardsBOE.reset_index(drop = True)


dfTenor = np.array(inputforwardsBOE.columns, dtype = float)
dfTenor[0] = 0.0  #Normalise T = 0.08333 = 0 to serve as short rate

'''calculate first 3 Principal Components using weekly differences (since daily differences was not giving inappropriate PCs
a = EigenValues
b = EigenVectors
c = Cumulative Sum of EigenValues to see how much variance can be explained by first 3 components
'''
a,b,c = calculate_pca(inputforwardsBOE, 5, 3)

print ("First PC explains  %8.2f , Second PC explains %8.2f  and Third PC explains = %8.2f percent of the variance " %(c.iloc[0]*100, c.iloc[1]*100, c.iloc[2]*100))
Lambda = a.iloc[0:3]
PrinComps = b.iloc[:,0:3]  #the vector of Principal components that we chose to fit the volatility
VolFunc = np.sqrt(Lambda)*PrinComps  #Array Multiplication. This is the Final Volatility Function that is the result of Calibration

PlotPCAs()

'''Next Step - Fitting of Vol Functions
'''

from sklearn import linear_model

X = pd.DataFrame(dfTenor, columns = ['Tau'])
X['Tau2'] = dfTenor**2

Y1 = VolFunc.iloc[:,0]

lm1 = linear_model.LinearRegression()
model = lm1.fit(X,Y1)
Rsq1 = lm1.score(X,Y1)
Coeff1 = lm1.coef_
Intercept1 = lm1.intercept_

'''Trying to fit the second and third volatility functions using Cubic Splines
'''
X = pd.DataFrame(dfTenor, columns = ['Tau'])
X['Tau2'] = dfTenor**2
X['Tau3'] = dfTenor**3

Y2 = VolFunc.iloc[:,1]
lm2 = linear_model.LinearRegression()
model = lm2.fit(X,Y2)
Rsq2 = lm2.score(X,Y2)
Coeff2 = lm2.coef_
Intercept2 = lm2.intercept_



Y3 = VolFunc.iloc[:,2]
lm3 = linear_model.LinearRegression()
model = lm3.fit(X,Y3)
Rsq3 = lm3.score(X,Y3)
Coeff3 = lm3.coef_
Intercept3 = lm3.intercept_

#This gives us the fitted Vol Functions
#plots for fitted functions
PlotFittedVolFunctions()

print ('Coefficients of Fitting the Volatility 2 Function' ,Intercept2, Coeff2 , Rsq2)
print ('Coefficients of Fitting the Volatility 3 Function', Intercept3, Coeff3 , Rsq3)


'''Next Step to calculate the functions V1, V2, V3 and No Arbitrage Drift that are used globaly as part of the SDE
to calculate the dF(t,T) = the change in forward rate
'''
Vol1 = np.median(VolFunc.iloc[:,0])
V1 =   V123(dfTenor, 1)
V2 = V123(dfTenor, 2)
V3 = V123(dfTenor, 3)

Drift = funcMVector(dfTenor)  
PlotDrift()

'''Calculate and Plot a sample Simulated Forward Curve'''

CurveA = HJMClasses.continuousHJM("continuous", dfTenor, inputforwardsBOE.iloc[-1]*.01)
CurveA.calcForwardCurve(Drift, V1, V2, V3, dfTenor, 3)
CurveA.sampleCurvePlots()

#Price a sample IRS
IRStoPrice = HJMClasses.IRS(5.0, 0.02, 0.5)
HJMClasses.calcPrice(CurveA, 0.0, 0.01, IRStoPrice) # Price at at As of Date T0 = 0 

'''Note - it takes 25 minutes for 2000 Simulations!!!'''
NumberofSimulations = 20
HJMmethod = 'Discrete'
if HJMmethod == 'Continuous':
    PriceSeries1 = getPricefromMC(NumberofSimulations, IRStoPrice, inputforwardsBOE.iloc[-1]*.01, 'continuous')
    print (PriceSeries1)
    MTMConvergence(NumberofSimulations,PriceSeries1)

'''Simulation using the discrete HJM formulation in Glasserman'''
#Converting to arrays to use in Spline Interpolation functions
RatesTenor = np.array(inputforwardsBOE.columns, dtype = float)
f = np.array(inputforwardsBOE.iloc[-1:]*.01)

maxMaturity = IRStoPrice.maturity
dtau = IRStoPrice.freq
DiscreteMats = np.arange(0, maxMaturity, dtau) 
InitialFHatsCurve = HJMClasses.getDiscreteFwdCurves(DiscreteMats,f[0],RatesTenor, dtau,  0.01) #f is the instataneous forward curve at time 0

''' Simulate 1 forward Curve and save as CSV'''
CurveB = HJMClasses.discretisedHJM("discrete", dfTenor, InitialFHatsCurve, maxMaturity,DiscreteMats)
CurveB.calcForwardCurve(Drift, V1, V2, V3, dfTenor, 3)
#CurveB.sampleCurvePlots()
CurveB.saveCurve()

HJMClasses.calcPrice(CurveB, 0.0, dtau, IRStoPrice) # Price at at As of Date T0 = 0 
NumberofSimulations = 1000
PriceSeries2 = getPricefromMC(NumberofSimulations, IRStoPrice, InitialFHatsCurve, 'discrete', 0.0, maxMaturity,DiscreteMats)
MTMConvergence(NumberofSimulations,PriceSeries2)

#Discount Factor curve at T = 0, used to get analytic value of IRS, ZCB
DiscountFactors = CurveB.getTodaysDiscountFactorCurve(IRStoPrice.freq, IRStoPrice.get_pay_schedule(0.0))              

'''Test for Zero Coupon Bond Convergence and hence testing No-Arbitrage of the HJM discrete model'''

'''
Calculate and Plots exposures. Also calculates the expected exposure as E[Max(MTM, 0)]
'''
MTMDataFrame = ExposureProfile(NumberofSimulations, IRStoPrice, InitialFHatsCurve)
MTMDataFrame.to_csv('IRSMTMsDiscreteHJM.csv')

exposureDates = np.arange(0, IRStoPrice.maturity, IRStoPrice.freq) 
MTMDataFrame.columns = exposureDates
PositiveMTM = MTMDataFrame
PositiveMTM = PositiveMTM.clip(lower=0)
ExpectedExposure = PositiveMTM.mean()

'''Plots
'''     
MTMDataFrame.hist( bins=50)  # MTM Distribution- Histograms for all Tenors

PositiveMTM.plot(kind = 'box')  #Box Plot for exposures

plt.plot(exposureDates, ExpectedExposure, 'o-')  ##Plot in Bold
plt.xlabel('SwapTenor')
plt.ylabel('IRS Expected Exposure')
plt.title('IRS Expected exposure acros swap tenor under discrete HJM ', fontsize=12)   
plt.show()
    
exposureOutPercentiles = 0.01*np.array([1, 2.5, 5, 10, 25, 50, 75, 90, 95, 97.5, 99])
expPercent = PositiveMTM.quantile(exposureOutPercentiles)
expMeans = PositiveMTM.mean()
expStd = PositiveMTM.std()
    
print ('Exposure Percentiles', expPercent)
print ('Expected Exposure', expMeans)
print ('Standard Deviation  = ',expStd)



ExpectedExposure['5.0'] = 0.0
SurvivalProb = funcBootstrapCDS(np.array(DiscountFactors))
calcCVA(ExpectedExposure,SurvivalProb, DiscountFactors, DiscountFactors.index, 0.4, Tenor = 0.5)
#Show convergence of IRS exposure to Analytic IRS price
ConvergenceDiscreteHJM(NumberofSimulations)

