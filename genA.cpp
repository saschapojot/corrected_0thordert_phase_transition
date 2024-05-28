#include "./version1/LJPotPBC/version1LJPotPBC2Atom.hpp"
#include <boost/python.hpp>
#include <fstream>

void printVec(const std::vector<double >& vec, std::ostream& os){
    for (size_t i = 0; i < vec.size(); ++i) {
        os << vec[i];
        if (i != vec.size() - 1) {
            os << ",";
        }
    }
}




int main(int argc, char *argv[]) {
    if (argc != 3) {
        std::cout << "wrong arguments" << std::endl;
        std::exit(2);
    }
    double T = std::stod(argv[1]);
    int rowNum=std::stoi(argv[2]);

    double  alpha1;
    double alpha2;
    double beta1;
    double beta2;

    double p1;
    double q1;

    double p2;
    double  q2;
    double r0=-1;
    version1dLJPot2Atom::parseCSV(rowNum,alpha1,beta1,p1,q1,alpha2,beta2,p2,q2,r0);

//    std::cout<<"alpha1="<<alpha1<<", beta1="<<beta1
//    <<", p1="<<p1<<", q1="<<q1
//    <<", alpha2="<<alpha2<<", beta2="<<beta2
//    <<", p2="<<p2<<", q2="<<q2<<std::endl;
    std::cout.precision(11);
    auto LJFunc=LJPotPBC(alpha1,alpha2,beta1,beta2,p1,p2,q1,q2,r0);

//    double  h=0.005;
    int cellNum = 10;
    auto v1Obj=version1dLJPot2Atom(rowNum,T,cellNum,std::make_shared<LJPotPBC>(alpha1,alpha2,beta1,beta2,p1,p2,q1,q2,r0));



}