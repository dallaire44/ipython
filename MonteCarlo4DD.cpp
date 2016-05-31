
// Monte Carlo simulations to find the Price and Greeks of a simple
// Vanilla Call Option.
// 
// David Dallaire - Example from FADBAD++ but modified
// 			- to include dividend yield
//                      - use C++14 random generator that default to Mersenne Twister 
//                  Answer trivial but coding not so
/*

	Answer from xterm

	BadSimpleMonteCarlo

	Price=9.11906
	dPriceSpot=0.579618
	dPriceVol=37.4675
	dPriceR=48.8427
	dDividend=-57.9618

	FADBADSimpleMonteCarlo

	Price=9.11906
	dPriceSpot=0.579618
	dPriceVol=37.4675
	dPriceR=48.8427
	dDividend=-57.9618
	dPrice2Spot=0
*/

#include <cstdlib>
#include <cmath>

#include "badiff.h"
#include "fadiff.h"
#include <iostream>
#include <random>

//compile with following
//g++ -I/media/ddu/HD1/MyMain/usr/FADBAD++ -std=c++14  -O3 -o MonteCarlo4DD MonteCarlo4DD.cpp -lm

using namespace std;
using namespace fadbad;

template <typename ArithT>
ArithT SimpleMonteCarlo(const double Expiry,
	const double Strike,
	const ArithT& Spot,
	const ArithT& Dividend,
	const ArithT& Vol,
	const ArithT& r,
	unsigned long NumberOfPaths)
{
	random_device rd;
  	//mt19937 gen(rd());
	//set constant see
	mt19937 gen(1);
  	normal_distribution<> d(0,1);
  
	ArithT variance = Vol*Vol*Expiry;
	ArithT rootVariance = sqrt(variance);
	ArithT itoCorrection = -0.5*variance;
	ArithT movedSpot = Spot*exp((r-Dividend)*Expiry +itoCorrection);
	ArithT thisSpot;
	ArithT runningSum=0;
	for (unsigned long i=0; i < NumberOfPaths; i++)
	{
		thisSpot = movedSpot*exp( rootVariance*d(gen));
		ArithT thisPayoff = thisSpot - Strike;
		thisPayoff = thisPayoff >0 ? thisPayoff : 0;
		runningSum += thisPayoff;
	}
	ArithT mean = runningSum / double(NumberOfPaths);
	mean *= exp(-r*Expiry);
	return mean;
}



template <typename ArithT>
void BADSimpleMonteCarlo(const double Expiry,
	const double Strike,
	const ArithT& Spot,
	const ArithT& Dividend,
	const ArithT& Vol,
	const ArithT& r,
	unsigned long NumberOfPaths,
	ArithT& Price,
	ArithT& dPriceSpot,
	ArithT& dPriceVol,
	ArithT& dPriceR,
	ArithT& dDividend)
{
	B< ArithT > BSpot(Spot);
	B< ArithT > BVol(Vol);
	B< ArithT > Br(r);
        B< ArithT > BDividend(Dividend);
	B< ArithT > BPrice=SimpleMonteCarlo(Expiry,Strike,BSpot,BDividend,BVol,Br,NumberOfPaths);
	BPrice.diff(0,1);
	Price=BPrice.x();
	dPriceSpot=BSpot.d(0);
	dPriceVol=BVol.d(0);
	dPriceR=Br.d(0);
	dDividend = BDividend.d(0);
}

template <typename ArithT>
void FADBADSimpleMonteCarlo(const double Expiry,
	const double Strike,
	const ArithT& Spot,
	const ArithT& Dividend,
	const ArithT& Vol,
	const ArithT& r,
	unsigned long NumberOfPaths,
	ArithT& Price,
	ArithT& dPriceSpot,
	ArithT& dPriceVol,
	ArithT& dPriceR,
	ArithT& dDividend,
	ArithT& dPrice2Spot)
{

	F< ArithT > FSpot(Spot);
	F< ArithT > FVol(Vol);
	F< ArithT > Fr(r);
	F< ArithT > FDividend(Dividend);
	FSpot.diff(0,1);
	F< ArithT > FPrice;
	F< ArithT > FdPriceSpot;
	F< ArithT > FdPriceVol;
	F< ArithT > FdPriceR;

	BADSimpleMonteCarlo(Expiry,Strike,FSpot,FDividend,FVol,Fr,NumberOfPaths,
		FPrice, FdPriceSpot, FdPriceVol, FdPriceR, FDividend);

	Price=FPrice.x();
	dPriceSpot=FdPriceSpot.x();
	dPriceVol=FdPriceVol.x();
	dPriceR=FdPriceR.x();
	dDividend = FDividend.x();
	dPrice2Spot=FdPriceSpot.d(0);
}


int main()
{
	double Strike=100;
	double Expiry=1.0;
        double Dividend = 0.03;
	double Spot=100;
	double Vol=.2;
	double r=0.06;
	double Price;
	double dPriceSpot;
	double dPriceVol;
	double dPriceR;	
	double dPrice2Spot;
        double dDividend;

	srand(0);
	
	BADSimpleMonteCarlo(Expiry,Strike,Spot,Dividend,Vol,r,100000,
		Price,dPriceSpot,dPriceVol,dPriceR,dDividend);
	cout << "\nBadSimpleMonteCarlo" << endl;
	cout << "\nPrice=" << Price << endl;
	cout << "dPriceSpot=" << dPriceSpot << endl;
	cout << "dPriceVol=" << dPriceVol << endl;
	cout << "dPriceR=" << dPriceR << endl;
	cout << "dDividend=" << dDividend << endl;

	srand(0);
	
	// NOTICE: that the gamma-value is wrong since the 
	// derivative of the intrinsic value as a function of 
	// spot is not Lipschitz.

	srand(0);

	FADBADSimpleMonteCarlo(Expiry,Strike,Spot,Dividend,Vol,r,100000,
		Price,dPriceSpot,dPriceVol,dPriceR,dDividend,dPrice2Spot);

	cout << "\nFADBADSimpleMonteCarlo" << endl;
	cout << "\nPrice=" << Price << endl;
	cout << "dPriceSpot=" << dPriceSpot << endl;
	cout << "dPriceVol=" << dPriceVol << endl;
	cout << "dPriceR=" << dPriceR << endl;
	cout << "dDividend=" << dDividend << endl;
	cout << "dPrice2Spot=" << dPrice2Spot << endl;
	
	return 0;

}
