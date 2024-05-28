//
// Created by polya on 5/13/24.
//

#include "./version1/LJPotPBC/parseXML.hpp"
std::vector<std::string> scanFiles(const int& rowNum){
    std::string searchPath="./version1Data/1d/funcLJPotPBC/row"+std::to_string(rowNum)+"/";
    std::vector<std::string> TDirs;
    if(fs::exists(searchPath) && fs::is_directory(searchPath)){
        for(const auto &entry:fs::directory_iterator(searchPath)){
            if(entry.path().filename().string()[0]=='T'){
                TDirs.push_back(entry.path().filename().string());
            }

        }
    }
//    for(const auto&s:TDirs)
//    {
//        std::cout<<s<<std::endl;
//    }

    return TDirs;

}



int main(int argc, char *argv[]) {
    if (argc != 2) {
        std::cerr << "wrong number of arguments" << std::endl;
        exit(1);
    }
    int rowNum = std::stoi(argv[1]);

    std::vector<std::string> TDirs = scanFiles(rowNum);
    for (const auto &s: TDirs) {
        std::cout<<"file is "<<s<<std::endl;
//        std::regex TPattern("T([+-]?\\d*(\\.\\d+)?)");
//        std::smatch T_match;
//        std::regex_search(s,T_match,TPattern);
//        if (std::stod(T_match.str(1))>0.45 or std::stod(T_match.str(1))<0.35){
//            continue;
//        }

        const auto tCStart{std::chrono::steady_clock::now()};
        auto rd = reader(rowNum, s);
        rd.searchFiles();
        rd.sortFiles();

        rd.parseSummary();
//        std::string smrAfterEq = rd.searchSummaryAfterEq();
//        rd.parseSummaryAfterEq(smrAfterEq);
        rd.UAndxFilesSelected();
        rd.parseUFiles();
        rd.parsexAxB();
        rd.data2json();
        rd.colmeans();
        rd.computeGAA();
        rd.computeGAB();
        rd.computeGBB();

    }

//    arma::dmat x{0,1,2,3,4,5,6,7,8,9,10,11};
//    x.reshape(4,3);
//    std::cout<<"x="<<x.t()<<std::endl;

}
