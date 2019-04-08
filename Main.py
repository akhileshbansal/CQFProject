# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import os
import pandas as pd
import numpy as np
from scipy.optimize import minimize, LinearConstraint, Bounds
from scipy.integrate import quad


def getDFfromForwards():
    DFcurve = []
    df = 1.0
    for i in Tenors:
        df = df*(1/(1 + inputForwardVol.loc[i]['Forward']))
        DFcurve.append(df)
    DFcurve =  pd.Series(data = np.array(DFcurve), index = Tenors)
    inputForwardVol['DF'] = DFcurve
    return DFcurve
    
'''
Alpha and Beta are in Absolute times from T0
Hence you add the option expiry and swap maturity to get the Beta

***TEST that sum of weights  = 1 for each Series
'''
def getWeights(DFcurve):
    Weights = {}
    for swaptionexpiry in np.arange(1, 11.0):
        alpha = swaptionexpiry
        Weights[alpha] = {}
        for swapmaturity in np.arange(1, 11.0):
            beta = swaptionexpiry + swapmaturity
            Weights[alpha][beta] = {}
        
            dfdcf= 0.0
            for i in np.arange(alpha+1, beta+1):
                dfdcf += DFcurve.loc[i]
            for i in np.arange(alpha+1, beta+1):
                Wi = DFcurve.loc[i]/dfdcf
                Weights[alpha][beta][i] = Wi
        
    return Weights

'''Alpha and Beta are in Absolute times from T0

***TEST that the Sx1 forward swap rates are same as the Forwards rates that are input 
'''
    
def getForwardSwapRates(Weights, Forwards):
    SwapRate = {}
    WeightAlphaBeta = {}
    for swaptionexpiry in np.arange(1, 11.0):
        alpha = swaptionexpiry
        SwapRate[alpha] = {}
        for swapmaturity in np.arange(1, 11.0):
            beta = swaptionexpiry + swapmaturity
            WeightAlphaBeta = Weights[alpha][beta]
            sumWiFi = 0.0
            for i in np.arange(alpha+1, beta+1):
                sumWiFi += Forwards.loc[i]*WeightAlphaBeta[i]
            SwapRate[alpha][beta] = sumWiFi
    return SwapRate

    #Integral for parametric volatility formulation using a,b,c,d
    #     
def getfuncI(t, Params, T):
    a,b,c,d = Params
    f = (((a*(T-t) +d)*np.exp(-b*(T-t)))+c)**2
    return f
 
def getfProdI(t, Params, Ti, Tj):
    a,b,c,d = Params
    f = (((a*(Ti-t) +d)*np.exp(-b*(Ti-t)))+c)*(((a*(Tj-t) +d)*np.exp(-b*(Tj-t)))+c)
    return f
    

'''
PSI, Vol, PHI will all be indexed 0-18 (total 19 values). Ideally would index them from 1-19
Vol[0] means Vol1 = caplet expiring at 1 year i.e. F (1x2) i.e. Vol of F[t,1,2]

i.e. Sigma[0,0] = Vol of F2, obviously only for 1 period from 0-1

Sigma[18,t] = Vol of F20 that expires in 19 years and applicable for t-t+1
'''
#PSIs are parameters for Volatility
# represents psi values for non-parametric formulation and a,b,c,d for parametric formulation
def getPHIs(MarketCapletVols, PSIs, Method = 'Non-Parametric'):
    tenorindex = np.arange(19.0)
    PHI = []
    tempsum = 0.0
    if Method == 'Parametric':
        
        for i, vol in enumerate(MarketCapletVols):
            I = quad(getfuncI, 0, i+1, args = (PSIs,i+1))  #would include a tau if dcf was not 1
            phi = np.sqrt((i+1)*vol**2/I[0]) #Ti*Volmktcaplet^2
            PHI.append(phi)
    else:

        for i, vol in enumerate(MarketCapletVols):
            tempsum += PSIs[i]**2  #would include a tau if dcf was not 1
            phi = np.sqrt((i+1)*vol**2/tempsum) #Ti*Volmktcaplet^2
            PHI.append(phi)
    PHI = pd.Series(PHI, index = tenorindex)    
    return PHI

'''
indexed to Ints - Ideally make it indexed to Tenors - implicitly done as you index it to PHIs
'''        
def getInstantVols(PHIs, PSIs):
    
        numrow = len(PHIs)
        numcol = len(PHIs)
        InstVolMatrix = pd.DataFrame(data =np.ones((numrow, numcol))*-1, columns = PHIs.index)
        InstVolMatrix = InstVolMatrix.set_index(PHIs.index)
        for i in range(numrow):
            for j in range(i+1):
                InstVolMatrix.loc[i,j] = PHIs[i]*PSIs[i-j] #if INdex 0-18, then this is i-j; if Index = 1-19, then i-j+1
        return InstVolMatrix
                
'''
indexed to Tenors - currently hard coded to 19.0 x 19.0
'''
def getInstantCorr(Thetas):
        tenorindex = np.arange(19.0)
        numrow = len(Thetas)
        numcol = len(Thetas)
        InstCorrMatrix = pd.DataFrame(data =np.ones((numrow, numcol))*-1, columns = tenorindex)
        InstCorrMatrix = InstCorrMatrix.set_index(tenorindex)
        for i in range(numrow):
            for j in range(numcol):
                InstCorrMatrix.loc[i,j] = np.cos(Thetas[i] - Thetas[j]) #if INdex 0-18, then this is i-j; if Index = 1-19, then i-j+1
        return InstCorrMatrix

'''Remember to check this Swaption Vol is Variance or Vol that can be compared to market
Same check needed for Caplet Vols. Also the Time factor needs to be considered'''
    
def RebonatoSwaption(PSIs, Thetas, MarketCapletVols, SwapRates, Forwards, Weights,Method = 'Non-Parametric' ):
    InstCorrMatrix = getInstantCorr(Thetas) 
    PHIs = getPHIs(MarketCapletVols, PSIs, Method)
    if Method =='Non-Parametric':    InstVolMatrix = getInstantVols(PHIs, PSIs)
    
    dtau = 1.0
    tenors = np.arange(1, 11.0)
    numrow = len(tenors)
    numcol = len(tenors)

    SwaptionLFM = pd.DataFrame(data = np.ones((numrow, numcol)), columns = tenors)
    SwaptionLFM = SwaptionLFM.set_index(tenors)
    for swaptionexpiry in np.arange(1, 11.0): # change it to swaptionexpiry in tenors
        if (swaptionexpiry== 6.0) or (swaptionexpiry== 8.0) or (swaptionexpiry== 9.0):continue
        alpha = swaptionexpiry
        for swapmaturity in np.arange(1, 11.0):
            beta = swaptionexpiry + swapmaturity
            WeightAlphaBeta = Weights[alpha][beta]
            swaptionVol = 0.0
            for i in np.arange(alpha +1, beta+1):
                for j in np.arange(alpha+1, beta+1):
                    sumVol = 0.0
                    if Method =='Parametric': sumVol = quad(getfProdI, 0, alpha, args = (PSIs,i-1, j-1))[0]*PHIs.loc[i-2]*PHIs.loc[j-2] #Integral(Phi)
                    else:
                        for k in np.arange(alpha): 
                            sumVol += InstVolMatrix.loc[i-2,k]*InstVolMatrix.loc[j-2,k]*dtau
                    swaptionVol += (sumVol*WeightAlphaBeta[i]*WeightAlphaBeta[j]*Forwards[i]*Forwards[j]*InstCorrMatrix.loc[i-2,j-2])/(SwapRates[alpha][beta]**2 * alpha)
            SwaptionLFM.loc[swaptionexpiry, swapmaturity] = np.sqrt(swaptionVol)
    return SwaptionLFM
            
                        
'''
Function takes Market Vols as Inputs
Calculates LFM Vols using built in function
For the given market vols - calculates the sum of squares of differences between market vols and LFM
Can be used in Optimization routine
'''
def minimizationfunc(PSIs, Thetas, MarketCapletVols, SwapRates, Forwards, Weights, MarketSwaptions):
    expiries = MarketSwaptions.index
    swapTenors = MarketSwaptions.columns
    swaptionLFM = RebonatoSwaption(PSIs, Thetas, MarketCapletVols, SwapRates, Forwards, Weights)
    
    SumSquare = 0.0
    for i in expiries:
        for j in swapTenors:
            SumSquare +=  (MarketSwaptions.loc[i,j]*0.01)**2 - swaptionLFM.loc[i,j]**2
    return SumSquare
 
    
class VolCalibration(object):

    def __init__(self, method):
        self.method = method
        self.Sol = None
        self.optimizedPSI = None
        self.optimizedTheta = None
        self.X0  = None

     
    
    def objectivefunc(self, x):
        return print ("error getFwdLiborRate")
    
class nonParametricVol(VolCalibration):
    
    def setInitValues(self):
        PSIs0 = np.ones(len(MarketCapletVols))
        Thetas0 = np.ones(len(MarketCapletVols))*0.5*np.pi        
        self.X0 = np.hstack((PSIs0, Thetas0))        
    
    def objectivefunc(self, x):
        PSIs = x[0:19]
        Thetas = x[19:]
    
        #Use this section to add constraints and Bounds
        thetaBound = 0.25*np.pi
        for psi in PSIs:
            if (psi < 0.0): return 1000
        for i in range(len(Thetas)-1):
            if abs(Thetas[i] - Thetas[i+1]) > thetaBound: return 1000
        
        
        
        expiries = MktSwaptionVol.index
        swapTenors = MktSwaptionVol.columns
        swaptionLFM = RebonatoSwaption(PSIs, Thetas, MarketCapletVols, ForwardSwapRate0, forwards, Weights0)
        
        SumSquare = 0.0
        for i in expiries:
            for j in swapTenors:
                if j == 1: continue #Skipping Sx1 vols, as we are using caplets vols
                SumSquare +=  (MktSwaptionVol.loc[i,j]*0.01 - swaptionLFM.loc[i,j])**2
        return SumSquare
    
    def outputParams(self):
        self.optimizedPSI = self.Sol.x[0:19]
        self.optimizedTheta = self.Sol.x[19:]
        return


class ParametricVol(VolCalibration):
    
    def setInitValues(self):
        PSIs0 = np.ones(4)
        Thetas0 = np.ones(len(MarketCapletVols))*0.5*np.pi
        Testa = 0.0285
        Testb = 0.2004
        Testc = 0.11
        Testd = 0.0570
        PSIs0 = Testa, Testb, Testc, Testd
        
        self.X0 = np.hstack((PSIs0, Thetas0))  
        return

    def objectivefunc(self, x):
        PSIs = x[0:4]
        Thetas = x[4:]
    
        #Use this section to add constraints and Bounds
        thetaBound = 0.333*np.pi
        PhiLBound = 0.7
        PhiUBound = 1.1
        PHIs = getPHIs(MarketCapletVols, PSIs, self.method)
        
        for phi in PHIs:
            if (phi < PhiLBound) or (phi > PhiUBound) : return 1000
        for i in range(len(Thetas)-1):
            if abs(Thetas[i] - Thetas[i+1]) > thetaBound: return 1000
        
        
        
        expiries = MktSwaptionVol.index
        swapTenors = MktSwaptionVol.columns
        swaptionLFM = RebonatoSwaption(PSIs, Thetas, MarketCapletVols, ForwardSwapRate0, forwards, Weights0, self.method)
        
        SumSquare = 0.0
        for i in expiries:
            for j in swapTenors:
                if j == 1: continue #Skipping Sx1 vols, as we are using caplets vols
                SumSquare +=  (MktSwaptionVol.loc[i,j]*0.01 - swaptionLFM.loc[i,j])**2
        return SumSquare
    
    def outputParams(self):
        self.optimizedPSI = self.Sol.x[0:4]
        self.optimizedTheta = self.Sol.x[4:]
        return




def objectivefunc(x):
    PSIs = x[0:19]
    Thetas = x[19:]

    #Use this section to add constraints and Bounds
    thetaBound = 0.25*np.pi
    for psi in PSIs:
        if (psi < 0.0): return 1000
    for i in range(len(Thetas)-1):
        if abs(Thetas[i] - Thetas[i+1]) > thetaBound: return 1000
    
    
    
    expiries = MktSwaptionVol.index
    swapTenors = MktSwaptionVol.columns
    swaptionLFM = RebonatoSwaption(PSIs, Thetas, MarketCapletVols, ForwardSwapRate0, forwards, Weights0)
    
    SumSquare = 0.0
    for i in expiries:
        for j in swapTenors:
            if j == 1: continue #Skipping Sx1 vols, as we are using caplets vols
            SumSquare +=  (MktSwaptionVol.loc[i,j]*0.01 - swaptionLFM.loc[i,j])**2
    return SumSquare

'''So the issue is with Alphas or option expiries.. Forward swap rates are correct..coz we checked for Sx1
 issue should be with integral part on the rebonato.. either integral or the vols...
 psis and phis have been input.. so not there
 

'''
#corr(Fi-Talpha)(Fj-Talpha) = rho i,j * integral(0-Talpha)Sigmai(t)*Sigmaj(t)/(integral(0-Talpha)Sigmai(t))*integral(0-Talpha)SigmaJ(t)
def calcTerminalCorrelation(Talpha, Thetas,PSIs):
    InstCorrMatrix = getInstantCorr(Thetas) 
    PHIs = getPHIs(MarketCapletVols, PSIs)
    InstVolMatrix = getInstantVols(PHIs, PSIs)
    sumCross = 0.0
    sumii = 0.0
    sumjj = 0.0
    tenorindex = np.arange(19.0)
    numrow = len(tenorindex)
    numcol = len(tenorindex)
    TermCorrMatrix = pd.DataFrame(data =np.ones((numrow, numcol))*-1, columns = tenorindex)
    TermCorrMatrix = InstCorrMatrix.set_index(tenorindex)
    for ii in tenorindex:
        for jj in tenorindex:
            for k in np.arange(Talpha):
                sumCross += InstVolMatrix.loc[ii,k]*InstVolMatrix.loc[jj,k]
                sumii += InstVolMatrix.loc[ii,k]
                sumjj += InstVolMatrix.loc[jj,k]
        termCorr = sumCross * InstCorrMatrix.loc[ii,jj] /(np.sqrt(sumii)*np.sqrt(sumjj))
        TermCorrMatrix.loc[ii,jj] = termCorr
    return TermCorrMatrix

def setBoundsConstraints():
#CONSTRAINED OPTIMIZATION
#Setting constraints

    A = np.zeros((18,38))
    for index in range(18):
        A[index, index] = 1
        A[index, index+1] = -1
    A = np.matrix(A)
    
    ConstLB = -0.5*np.pi#*np.ones(18)
    ConstUB = 0.5*np.pi#*np.ones(18)
    
    const = LinearConstraint(A,ConstLB, ConstUB )
    
    
    BndLB = -1*np.ones(38)*np.inf
    BndUB = np.ones(38)*np.inf
    
    for index in range(19): BndLB[index] = 0.0
    bnds = Bounds(BndLB,BndUB)
    
    return (const, bnds)


    

'''
Main Project thread
'''

os.getcwd()
'''Add the Path to the directory all Data Files are stored'''

path="C:\\Users\Akhil\Documents\CQF\LMM"
os.chdir(path)
os.getcwd()



inputForwardVol = pd.read_excel('InputData.xlsx', 'Rates')  # contains Annual forwards and caplet vols
inputForwardVol = inputForwardVol.set_index('Maturity')
Tenors = inputForwardVol.index
forwards = pd.Series(inputForwardVol['Forward'])
DF = getDFfromForwards()
Weights0 = getWeights(DF)
ForwardSwapRate0 = getForwardSwapRates(Weights0,forwards)
MarketCapletVols = inputForwardVol.iloc[1:]['Caplet Vol'].reset_index(drop = True)
PSIs0 = np.ones(len(MarketCapletVols))
Thetas0 = np.ones(len(MarketCapletVols))*0.5*np.pi
MktSwaptionVol = pd.read_excel('InputData.xlsx', 'SwaptionVol')  # contains Annual forwards and caplet vols
MktSwaptionVol.set_index('Expiry', inplace = True)
#sumin = minimizationfunc(PSIs, Thetas, MarketCapletVols, ForwardSwapRate0, forwards, Weights0, MktSwaptionVol)

X0 = np.hstack((PSIs0, Thetas0))

'''
Comment out while testing
'''

#UNCONSTRAINED OPTIMIZATION
#Sol = minimize(objectivefunc, X0, method = 'Nelder-Mead', tol = 1e-5, options={'maxiter': 5000, 'disp': True} )

optim = ParametricVol('Parametric')
optim.setInitValues()
optim.Sol = minimize(optim.objectivefunc, optim.X0, method = 'Nelder-Mead', tol = 1e-5, options={'maxiter': 100, 'disp': True} )
optim.outputParams()
optimizedVols = RebonatoSwaption(optim.optimizedPSI, optim.optimizedTheta, MarketCapletVols, ForwardSwapRate0, forwards, Weights0, optim.method)
Diffs = (MktSwaptionVol - optimizedVols*100)*100/ MktSwaptionVol

#CONSTRAINED OPTIMIZATION
#consts, bounds = setBoundsConstraints
#Sol = minimize(objectivefunc, X0, method = 'trust-constr', bounds = bnds, tol = 1e-5, options={'maxiter': 100, 'disp': True} )







'''
Test
TestParams = pd.read_excel('InputData.xlsx', 'TestParams')  # contains Annual forwards and caplet vols
TestPsi = np.array(TestParams['Psi'])
TestTheta = np.array(TestParams['Theta'])

optimizedVols = RebonatoSwaption(TestPsi, TestTheta, MarketCapletVols, ForwardSwapRate0, forwards, Weights0)
'''
'''
TestParams = pd.read_excel('InputData.xlsx', 'TestParams2')  # contains Annual forwards and caplet vols
TestTheta = np.array(TestParams['Theta'])
TestPsi = np.ones(4)
Testa = 0.29342753
Testb = 1.25080230
Testc = 0.13145869
Testd = 0.0
TestPsi = Testa, Testb, Testc, Testd
optimizedVols = RebonatoSwaption(TestPsi, TestTheta, MarketCapletVols, ForwardSwapRate0, forwards, Weights0, 'Parametric')
'''


