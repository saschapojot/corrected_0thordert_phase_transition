//
// Created by polya on 5/10/24.
//

#ifndef T_PHASE_TRANSITION_VERSION1LJPOT2ATOM_HPP
#define T_PHASE_TRANSITION_VERSION1LJPOT2ATOM_HPP

#include <algorithm>
#include <armadillo>
#include <array>

#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/filesystem.hpp>
#include <boost/json.hpp>
#include <boost/python.hpp>
#include <boost/serialization/complex.hpp>
#include <boost/serialization/vector.hpp>

#include <cmath>
#include <chrono>
#include <cstdlib>
#include <cxxabi.h>

#include <fstream>
#include <iostream>
#include <memory>
#include <msgpack.hpp>
#include <random>
#include <regex>
#include <sstream>
#include <string>
#include <typeinfo>
#include <vector>

namespace fs = boost::filesystem;
//this subroutine computes the mc evolution for a 1d system, 2-atom, Lennard-Jones+quartic potential +PBC

class potentialFunction {
    //base class for potential function
public:
    potentialFunction (const double &alpha1Val, const double &alpha2Val, const double &beta1Val,
                       const double &beta2Val, const double &p1Val, const double &p2Val, const double &q1Val, const double &q2Val, const double &r0) {

        this->alpha1 = alpha1Val;
        this->alpha2 = alpha2Val;
        this->beta1 = beta1Val;
        this->beta2 = beta2Val;
        this->p1 = p1Val;
        this->p2 = p2Val;
        this->q1 = q1Val;
        this->q2 = q2Val;
        this->r0=r0;


        std::cout<<"alpha1="<<alpha1<<", beta1="<<beta1
                 <<", p1="<<p1<<", q1="<<q1
                 <<", alpha2="<<alpha2<<", beta2="<<beta2
                 <<", p2="<<p2<<", q2="<<q2<<std::endl;

    }//end of constructor
public:
    virtual double operator()(const arma::dcolvec &xA, const arma::dcolvec &xB, const double &L) const = 0;
    virtual double dVEst(const double &r, const unsigned long long &N)const = 0;

    virtual ~ potentialFunction() {};


public:
    double alpha1;
    double alpha2;
    double beta1;
    double beta2;
    double p1;
    double p2 ;
    double q1 ;
    double q2 ;
    double r0;// eq distance

};

class LJPotPBC : public potentialFunction {
public:
    LJPotPBC(const double &alpha1Val, const double &alpha2Val, const double &beta1Val,
             const double &beta2Val, const double &p1Val, const double &p2Val, const double &q1Val, const double &q2Val, const double &r0Val):potentialFunction(alpha1Val, alpha2Val, beta1Val, beta2Val, p1Val, p2Val, q1Val, q2Val, r0Val)  {
        this->alpha1 = alpha1Val;
        this->alpha2 = alpha2Val;
        this->beta1 = beta1Val;
        this->beta2 = beta2Val;
        this->p1 = p1Val;
        this->p2 = p2Val;
        this->q1 = q1Val;
        this->q2 = q2Val;
        this->r0=r0Val;

    }//end of constructor

public:
    ///
    /// @param xA positions of atom A
    /// @param xB positions of atom B
    /// @return potential energy
    double operator()(const arma::dcolvec &xA, const arma::dcolvec &xB, const double& L) const override {
        return V1Total(xA, xB) + V2Total(xA, xB,L);

    }


    ///
    /// @param xA positions of atom A
    /// @param xB positions of atom B
    /// @return the sum of all V1 energy
    double V1Total(const arma::dcolvec &xA, const arma::dcolvec &xB) const {
        arma::dcolvec rVec = arma::abs(xA - xB);


        arma::dcolvec vecPart1 = alpha1 * arma::pow(rVec, -p1);
        arma::dcolvec vecPart2 = -beta1 * arma::pow(rVec, -q1);
        arma::dcolvec vecPart3 = arma::pow(rVec, 4);
        std::cout<<"vecPart1="<<vecPart1<<", vecPart2="<<vecPart2<<", vecPart3="<<vecPart3<<std::endl;
        double val = arma::sum(vecPart1) + arma::sum(vecPart2) + arma::sum(vecPart3);
//        std::cout<<"V1Total="<<val<<std::endl;
        return val;


    }

    ///
    /// @param xA positions of atom A
    /// @param xB positions of atom B
    /// @param L length of the PBC loop
    /// @return the sum of all V2 energy under OBC
    double V2Total(const arma::dcolvec &xA, const arma::dcolvec &xB, const double &L) const {
        int N = xB.size();
        if (N <= 1) {
            return 0;
        }
        arma::dcolvec sliceA = xA.subvec(1, N - 1);
//    std::cout<<"sliceA="<<sliceA<<std::endl;
        arma::dcolvec sliceB = xB.subvec(0, N - 2);
        arma::dcolvec rVec = arma::abs(sliceA - sliceB);
//        std::cout<<"sliceB="<<sliceB<<std::endl;

//        std::cout<<"V2Total: rVec="<<rVec<<std::endl;

        arma::dcolvec vecPart1 = alpha2 * arma::pow(rVec, -p2);
        arma::dcolvec vecPart2 = -beta2 * arma::pow(rVec, -q2);
        arma::dcolvec vecPart3 = arma::pow(rVec, 4);

        double rLastAbs=std::abs(L-(xB(N-1)-xA(0)));
        double valBoundary=alpha2*std::pow(rLastAbs,-p2)-beta2*std::pow(rLastAbs,-q2)+std::pow(rLastAbs,4);

        double val = arma::sum(vecPart1) + arma::sum(vecPart2) + arma::sum(vecPart3)+valBoundary;
//        std::cout<<"V2Total="<<val<<std::endl;

        return val;

    }

    double dV1(const double &r)const{
        return -alpha1*p1*std::pow(r,(-p1-1))+beta1*q1*std::pow(r,(-q1-1))
               +4*std::pow(r,3);

    }

    double dV2(const double &r)const{


        return -alpha2*p2*std::pow(r,(-p2-1))+beta2*q2*std::pow(r,(-q2-1))
               +4*std::pow(r,3);
    }

    double dVEst(const double &r, const unsigned long long &N)const{
        double val=static_cast<double>(N)*(dV1(r)+ dV2(r));
        return val;

    }


public:
    double alpha1;
    double alpha2;
    double beta1;
    double beta2 ;
    double p1 ;
    double p2 ;
    double q1 ;
    double q2 ;
    double r0;// eq distance


};


class version1dLJPot2Atom {
public:
    version1dLJPot2Atom(int rowNum, double temperature, unsigned long long cellNum,
                        const std::shared_ptr<potentialFunction> &funcPtr) {
        this->rowNum = rowNum;
        this->T = temperature;
        this->beta = 1 / T;
//        this->h = stepSize;
        this->potFuncPtr = funcPtr;
        double rEst=funcPtr->r0;

        std::cout<<"rEst="<<rEst<<std::endl;
        double dValEst=2;

        double stepSize=dValEst*T/(std::abs(funcPtr->dVEst(rEst,cellNum)));
        if (stepSize>0.005){
            stepSize=0.005;
        }
        this->h=0.005;//stepSize;
        std::cout<<"h="<<h<<std::endl;

//        this->diag=isDiag;
        this->N = cellNum;
        this->stddev =h;// std::sqrt(2.0 * h);
    }


public:

    ///
    /// @param xA positions of atom A
    /// @param xB positions of atom B
    /// @param L
    /// @return beta*potential
    double f(const arma::dcolvec &xA, const arma::dcolvec &xB, const double & L);

    ///
    /// @param xACurr positions of atom A
    /// @param xBCurr positions of atom B
    /// @param LCurr
    /// @param zANext proposed positions of atom A
    /// @param zBNext proposed positions of atom B
    /// @param LNext
    void proposal(const arma::dcolvec &xACurr, const arma::dcolvec &xBCurr, const double &LCurr,
                  arma::dcolvec &zANext, arma::dcolvec &zBNext, double &LNext);


    ///
    /// @param filename xml file name of vec
    ///@param vec vector to be saved
    static void saveVecToXML(const std::string &filename, const std::vector<double> &vec);

    ///
    /// @param filename bin file name of vec
    /// @param vec vector to be saved
    static void saveVecToBin(const std::string &filename, const std::vector<double> &vec);

    ///
    /// @param cmd python execution string
    /// @return signal from the python
    static std::string execPython(const char *cmd);


    ///
    /// @param xA current positions of atom A
    /// @param xB current positions of atom B
    /// @param LCurr
    /// @param zA proposed positions of atom A
    /// @param zB proposed positions of atom B
    /// @param LNext
    /// @return
    double acceptanceRatio(const arma::dcolvec &xA, const arma::dcolvec &xB,const double& LCurr,
                           const arma::dcolvec &zA, const arma::dcolvec &zB, const double &LNext);

    ///
    /// @param xAInit initial positions of A
    /// @param xBInit initial positions of B
    /// @param LInit
    void initPositionsEquiDistance(arma::dcolvec &xAInit, arma::dcolvec &xBInit, double &LInit);

    ///
    /// @param filename  xml file name of vecvec
    /// @param vecvec vector<vector> to be saved
    static void saveVecVecToXML(const std::string &filename, const std::vector<std::vector<double>> &vecvec);

    ///
    /// @param filename bin file name of vecvec
    /// @param vecvec vector<vector> to be saved
    static void saveVecVecToBin(const std::string &filename, const std::vector<std::vector<double>> &vecvec);

    ///
    /// @param lag decorrelation length
    /// @param loopTotal total mc steps
    /// @param equilibrium whether equilibrium has reached
    /// @param same whether all values of potential are the same
    /// @param xALast last positions of atom A
    /// @param xBLast last positions of atom B
    /// @param LLast
    void readEqMc(int &lag, int &loopTotal, bool &equilibrium, bool &same, std::vector<double> &xALast,
                  std::vector<double> &xBLast, double &LLast);

    ///
    /// @param lag decorrelation length
    /// @param loopEq total loop numbers in reaching equilibrium
    /// @param xA_init xA from readEqMc
    /// @param xB_init xB from readEqMc
    /// @param LInit
    void executionMCAfterEq(const int &lag, const int &loopEq, const std::vector<double> &xA_init,
                            const std::vector<double> &xB_init, const double &LInit);

    std::string demangle(const char *name) {
        int status = -1;
        char *demangled = abi::__cxa_demangle(name, NULL, NULL, &status);
        std::string result(name);
        if (status == 0) {
            result = demangled;
        }
        std::free(demangled);
        return result;
    }

    ///
    /// @param rowNum row number
    static void parseCSV(const int &rowNum, double &alpha1, double &beta1, double &p1, double &q1,
                         double &alpha2, double &beta2, double &p2, double &q2,double &r0);

    double T;// temperature
    double beta;
    int moveNumInOneFlush = 3000;// flush the results to python every moveNumInOneFlush iterations
    int flushMaxNum = 7000;
    int dataNumTotal = 8000;
    double h;// step size
//    double a=5.3;//stiffness
//    bool diag=true;// whether the quadratic form of energy is diagonal
    int N;//number of unit cells

    double lastFileNum = 0;
    std::shared_ptr<potentialFunction> potFuncPtr;
    double stddev;
    int rowNum;
    int startingFileInd=0;

};


void save_array_to_pickle(double* ptr, std::size_t size, const std::string& filename);
///to msgpack bin file
void save_to_bin_file(double* data, unsigned long long size, const std::string& filename);

#endif //T_PHASE_TRANSITION_VERSION1LJPOT2ATOM_HPP
