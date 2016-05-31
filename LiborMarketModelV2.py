'''
Created on Jan 21, 2016

port of following C code:
http://www.quantcode.com/modules/mydownloads/singlefile.php?lid=63


[1] Linking Caplet and Swaption Volatilities in a BGM/J Framework: Approximate Solutions by Peter Jackel, Riccardo Rebonato
@author: David Dallaire

Purpose is to learn pricing swaptions using Libor Market Model with Monte Carlo
'''
import numpy as np
import math as math
import sys
import random
import csv


def GetForwardRatesAndTerms(P,tau):
    theRates = np.zeros(P.shape[0]-1)
    theTerms = np.zeros(P.shape[0]-1)
    #lose one observation at end
    for i in range(0,P.shape[0]-1):
        #print(i)
        val = (P[i,]/P[i+1,]-1)/tau
        theRates[i,] = val
        theTerms[i,] = i*tau
        #print(val)
    return theTerms,theRates

def CreateCorrelation(beta,terms):
    numberOfTerms= len(terms)
    corrMatrix = np.zeros((numberOfTerms,numberOfTerms))
    countX = 0
    for i in terms:
        countY = 0
        for j in terms:
            val = math.exp(-beta * abs(i - j))
            #print(str(i) + " " + str(j) + " " + str(val))
            corrMatrix[countX,countY] = val
            countY += 1
        countX += 1
    return corrMatrix

#equation 19 in [1]
def corr(Tj,Tk,beta):
    return math.exp(-beta * abs(Tj - Tk))


#see equation 18 in [1]
def vol(Tj,T_o):
    a = -0.05
    b = 0.5
    c = 1.5
    d = 0.15
    return ((a + b * (Tj - T_o)) * math.exp(-c * (Tj - T_o)) + d)

def GetSwapRate(F,alpha,beta):
    alpha = int(alpha)
    beta = int(beta)
    #print(alpha)
    #print(beta)
    tmp = 1
    tmpSum = 0
    swapRate = 1
    for j in range(alpha,beta):
        tmp *= (1/(1 + tau * F[j,]))
    #Floating Leg
    SR = 1 - tmp
    for i in range(alpha,beta):
        tmp = 1
        for j in range(alpha,i+1):
            tmp *= (1/(1+tau*F[j,]))
        tmpSum +=  tau * tmp
    SR = SR/tmpSum
    return SR   


def MonteCarloLiborMarketModel(beta,tau,end,swapPrices,numberOfScenarios):

    terms,forwardRates = GetForwardRatesAndTerms(swapPrices,tau)
    
    swapRate = GetSwapRate(forwardRates,g_alpha,g_beta)
    terms = np.arange(0,end,tau)
 
    corrMatrix = CreateCorrelation(beta,terms)
    cholMatrix = (np.linalg.cholesky(corrMatrix)).T
    
    #anti Thetic Flag starts at false
    antiTheticFlag = False

    timePoints = ((g_alpha * tau) / dt)
    
    rand_nos_mat=np.zeros((timePoints,g_beta))  #store cholRand for use in anti
    
    #initialize arrays
    libor_simulations = np.zeros((numberOfScenarios,g_beta))
    finalFVec = np.zeros((1,swapPrices.shape[0]))
    discountCurve = np.zeros((1,swapPrices.shape[0]))
    
    sqrDt = math.sqrt(dt)
    payoff_sum = 0
    
    
    useRandom = True #This is used for debuging. When false, random numbers are ones.

    np.random.seed(100000)
    
    tmpCount = 0
    
    for scenario in range(0,numberOfScenarios): 
        currLog = np.log(forwardRates)
        currLog_t = np.log(forwardRates)
    
        if antiTheticFlag == True:
            antiTheticFlag = False
        else:
            antiTheticFlag = True
            
        mOnes = np.ones((1,14))
        
        t = 0
        count = 0
       
        while(t < (g_alpha * tau)):
            if(antiTheticFlag):
                randNumbers = np.random.randn(g_beta)
               
                if (useRandom):
                    randChol = np.dot(randNumbers,cholMatrix)
                else:
                    randChol = np.dot(mOnes,cholMatrix)
                rand_nos_mat[count,:] = randChol
            else:
                randChol = -1 * rand_nos_mat[count,:]
            
            nextResetIdx = math.floor(t/tau) + 1

            for k in range(int(nextResetIdx),int(g_beta)):
                driftSum = 0
                for j in range(nextResetIdx,k + 1):
                    tmp = corr(terms[k],terms[j],beta) * tau * vol(terms[j],t) * math.exp(currLog_t[j])
                    tmp = tmp / ( 1 + tau * math.exp(currLog_t[j]))
                    driftSum += tmp           
                dLogF=0
                vol_Tk_t = vol(terms[k],t)
                dLogF += vol_Tk_t*driftSum*dt
                dLogF -= 0.5*vol_Tk_t*vol_Tk_t * dt
                
                if (useRandom):
                    dLogF = dLogF + vol_Tk_t * randChol[k] * sqrDt
                else:
                    dLogF += vol_Tk_t*mOnes[0,k] * sqrDt
                currLog[k] += dLogF
                
            for i in range(0,currLog_t.shape[0]):
                currLog_t[i] = currLog[i]
            count += 1     
            t = t + dt
            
        for i in range(0,int(g_beta)):
            libor_simulations[scenario,i] = math.exp(currLog_t[i])
        
        for i in range(0,int(g_beta)):
            finalFVec[0,i] =  math.exp(currLog_t[i])
        discountCurve[0,0] = 1 / (1 + tau * finalFVec[0,0])
        
        for i in range(1,int(g_beta)):
            discountCurve[0,i] = discountCurve[0,i-1] / (1 + tau * finalFVec[0,i])
        
        payoff = 0
        
        for i in range(int(g_alpha),int(g_beta)):
            payoff += (swapRate - finalFVec[0,i]) * tau * discountCurve[0,i]
        
        payoff_sum += max(payoff,0)
       
    return 100*payoff_sum/numberOfScenarios 

        
if __name__=='__main__':
    numberOfScenarios = 1000
 
        
    start = 2.5 #in years
    end = 7 #in years
    
    print("Swap start = " + str(start) + " year")
    print("Swap end = " + str(end) + " year")   
    beta = 0.1
    tau = 0.5 #term
    
    g_alpha = start / tau
    g_beta = end / tau
    dt = 0.1
    
    swapPrices = np.array([1.000000,0.969514,0.939441,0.909913,0.881024,0.852807,0.825482,0.799100,
                  0.773438,0.749042,0.725408,0.702527,0.680361,0.659402,0.639171,0.619580,
                  0.600668,0.582455,0.564873,0.547888,0.531492,0.515651,0.500360,0.485543,
                  0.471240,0.457861,0.444977,0.432554,0.420575,0.409019,0.397888,0.387341,
                  0.377196,0.367435,0.358056,0.348978,0.340292,0.331614,0.323265,0.315460,
                  0.307945,0.300321])
    #print(P.shape)
    
    #flatCurve = np.array([1.000000,0.966736,0.934579,0.903492,0.873439,0.844385,0.816298,
    #                     0.789145,0.762895,0.737519,0.712986,0.689270,0.666342,0.644177,
    #                     0.622750,0.602035,0.582009,0.562649,0.543934,0.525841,0.508349,
    #                     0.491440,0.475093,0.459290,0.444012,0.429243,0.414964,0.401161,
    #                     0.387817,0.374917,0.362446,0.350390,0.338735,0.327467,0.316574,
    #                     0.306044,0.295864,0.286022,0.276508,0.267311,0.258419,0.249823])
    
    #print(flatCurve.shape)
    swaptionPrice = MonteCarloLiborMarketModel(beta,tau,end,swapPrices,numberOfScenarios)
    
    print("swaption price = " + str(swaptionPrice))
