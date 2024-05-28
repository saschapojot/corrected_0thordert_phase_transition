//
// Created by polya on 5/13/24.
//
#include "parseXML.hpp"

///UAll, xA_All, xB_All folder's files
void reader::searchFiles() {

//    this->UPath = this->TDir + "/UAll/";
    this->UPath = this->TDir + "/UAllBin/";
//    this->xAPath = this->TDir + "/xA_All/";
    this->xAPath = this->TDir + "/xA_AllBin/";
//    this->xBPath = this->TDir + "/xB_All/";
    this->xBPath = this->TDir + "/xB_AllBin/";
    std::cout<<"UPath="<<UPath<<std::endl;
//    std::cout<<xAPath<<std::endl;
//    std::cout<<xAPath<<std::endl;
    for (const auto &entry: fs::directory_iterator(UPath)) {
        this->UFilesAll.push_back(entry.path().string());
    }
//    printVec(UFilesAll);

    for (const auto &entry: fs::directory_iterator(xAPath)) {
        this->xAFilesAll.push_back(entry.path().string());
    }

    for (const auto &entry: fs::directory_iterator(xBPath)) {
        this->xBFilesAll.push_back(entry.path().string());
    }


}


///
/// @param path the path containing xml files
/// @return sorted xml files by starting loop
std::vector<std::string> reader::sortOneDir(const std::vector<std::string> &allFiles) {
    std::vector<int> startingLoopsAll;
    for (const std::string &name: allFiles) {
        std::regex startPattern("loopStart(\\d+)loopEnd");
        std::smatch matchPattern;
        if (std::regex_search(name, matchPattern, startPattern)) {
            startingLoopsAll.push_back(std::stoi(matchPattern.str(1)));
        }

    }

    std::vector<size_t> inds = this->argsort<int>(startingLoopsAll);;

    std::vector<std::string> sortedFiles;
    for (const auto &i: inds) {
        sortedFiles.push_back(allFiles[i]);
    }

    std::cout<<"sortedFiles[0]="<<sortedFiles[0]<<std::endl;
    std::cout<<"sortedFiles[1]="<<sortedFiles[1]<<std::endl;

    return sortedFiles;


}


///sort files by starting loop
void reader::sortFiles() {
    this->sorted_UFilesAll = this->sortOneDir(this->UFilesAll);
    this->sorted_xAFilesAll = this->sortOneDir(this->xAFilesAll);
    this->sorted_xBFilesAll = this->sortOneDir(this->xBFilesAll);
    std::cout<<"sorted_UFilesAll[0]="<<sorted_UFilesAll[0]<<std::endl;
    std::cout<<"sorted_UFilesAll[1]="<<sorted_UFilesAll[1]<<std::endl;

}


void reader::parseSummary() {
    std::cout<<"entering parseSummary()"<<std::endl;
    std::string smrPath = TDir + "/summary.txt";
    std::regex lagPattern("lag=([+-]?\\d+)");
    std::regex lastFilesNumPattern("lastFileNum=(\\d+)");

    std::smatch matchLag;
    std::smatch matchFileNum;


    std::ifstream smrIn(smrPath);
    for (std::string line; std::getline(smrIn, line);) {
        //extract lag value
        if (std::regex_search(line, matchLag, lagPattern)) {
            this->lag = std::stoi(matchLag.str(1));
            std::cout << "lag=" << lag << std::endl;
        }
        //extract lastFilesNum
        if (std::regex_search(line, matchFileNum, lastFilesNumPattern)) {

            this->lastFileNum = std::stoi(matchFileNum.str(1));

            std::cout<<"lastFilesNum="<<lastFileNum<<std::endl;
        }

    }//end readline for


std::cout<<"leaving parseSummary()"<<std::endl;
}

std::string reader::searchSummaryAfterEq() {

    std::regex afterEqPattern("summaryAfterEq");
    std::smatch matchAfter;

    for (const auto &entry: fs::directory_iterator(this->TDir)) {
        std::string fileName = entry.path().string();
        if (std::regex_search(fileName, matchAfter, afterEqPattern)) {
            return fileName;

        }


    }//end for

    return "";


}

//void reader::parseSummaryAfterEq(const std::string &afterEqPath) {
//
//    std::regex loopPattern("total loop number:\\s*(\\d+)");
//    std::smatch matchLoop;
//    if (afterEqPath.size() > 0) {
//        std::ifstream afterIn(afterEqPath);
//        for (std::string line; std::getline(afterIn, line);) {
//            if (std::regex_search(line, matchLoop, loopPattern)) {
//                this->loopNumAfterEq = std::stoi(matchLoop.str(1));
//                std::cout<<"loopNumAfterEq="<<loopNumAfterEq<<std::endl;
//                break;
//
//            }//end if
//
//        }//end for
//
//
//    }//end if
//
//}


void reader::UAndxFilesSelected() {
//    loopNumToInclude = moveNumInOneFlush * lastFileNum + loopNumAfterEq;
//    std::cout<<"loopNumToInclude="<<loopNumToInclude<<std::endl;
//    double loopNumToIncludeDB = static_cast<double >(loopNumToInclude);
//    double moveNumInOneFlushDB = static_cast<double >(moveNumInOneFlush);

//    fileNumSelected = static_cast<int>(std::ceil(loopNumToIncludeDB / moveNumInOneFlushDB));
    int numfilePart0=2000;
    int numFilePart0Selected=numfilePart0-lastFileNum;
    fileNumSelected=sorted_UFilesAll.size()-numfilePart0+numFilePart0Selected;

    std::cout<<"fileNumSelected="<<fileNumSelected<<std::endl;


    for (int i = sorted_UFilesAll.size() - fileNumSelected; i < sorted_UFilesAll.size(); i++) {
        this->UFilesSelected.push_back(sorted_UFilesAll[i]);
    }
// std::cout<<"fileNumSelected="<<fileNumSelected<<std::endl;
// std::cout<<"len(UFilesSelected)="<<UFilesSelected.size()<<std::endl;

    for (int i = sorted_xAFilesAll.size() - fileNumSelected; i < sorted_xAFilesAll.size(); i++) {
        this->xAFilesSelected.push_back(sorted_xAFilesAll[i]);
    }

    for (int i = sorted_xBFilesAll.size() - fileNumSelected; i < sorted_xBFilesAll.size(); i++) {
        this->xBFilesSelected.push_back(sorted_xBFilesAll[i]);
    }


}


void reader::parseUFiles() {
    const auto tUStart{std::chrono::steady_clock::now()};
    UIn.reserve(this->UFilesSelected.size() * moveNumInOneFlush);
    for (const std::string &oneUFile: this->UFilesSelected) {
        std::vector<double> vecInOneFile = readMsgBinVec(oneUFile);
//        std::ifstream ifs(oneUFile);
//        if (!ifs.is_open()) {
//            std::cerr << "cannot open " << oneUFile << std::endl;
//            return;
//        }
//        boost::archive::xml_iarchive ia(ifs);
//        ia >> BOOST_SERIALIZATION_NVP(vecInOneFile);

        UIn.insert(UIn.end(), vecInOneFile.begin(), vecInOneFile.end());


    }

    const auto tUEnd{std::chrono::steady_clock::now()};
    const std::chrono::duration<double> elapsed_secondsAll{tUEnd - tUStart};
    std::cout << "parse U time: " << elapsed_secondsAll.count() << " s" << std::endl;

    for (int i = 0; i < UIn.size(); i += lag + 1) {
        USelected.push_back(UIn[i]);
    }
//    std::cout<<"lag="<<lag<<std::endl;
//    std::cout<<"len(UIn)="<<UIn.size()<<std::endl;
//    std::cout<<"len(USelected)="<<USelected.size()<<std::endl;
//    armaU= arma::dcolvec (USelected);

//    std::cout<<armaU<<std::endl;


}


void reader::parsexAxB() {

    //A
    const auto tAStart{std::chrono::steady_clock::now()};
    //reserve lengths
    std::vector<std::vector<double>> zerothVecVec = readMsgBinVecVec(xAFilesSelected[0]);

    cellNum = zerothVecVec[0].size();
    std::cout << "cellNum=" << cellNum << std::endl;
    xAIn.reserve(moveNumInOneFlush * xAFilesSelected.size());
    std::vector<double> initVecVal(cellNum, 0);
    for (int i = 0; i < moveNumInOneFlush * xAFilesSelected.size(); i++) {
        xAIn.push_back(initVecVal);
    }
    xBIn.reserve(moveNumInOneFlush * xAFilesSelected.size());
    for (int i = 0; i < moveNumInOneFlush * xAFilesSelected.size(); i++) {
        xBIn.push_back(initVecVal);
    }
    int AStart = 0;

    for (const std::string &onexAFile: xAFilesSelected) {
        std::vector<std::vector<double>> vecVecInOneFile = readMsgBinVecVec(onexAFile);
//    std::ifstream ifs(onexAFile);
//        if (!ifs.is_open()) {
//            std::cerr << "cannot open "<<onexAFile << std::endl;
//            return;
//        }
//        boost::archive::xml_iarchive ia(ifs);
//        ia >> BOOST_SERIALIZATION_NVP(vecVecInOneFile);
//

//        xAIn.insert(xAIn.end(),vecVecInOneFile.begin(),vecVecInOneFile.end());
        int lengthTmp = vecVecInOneFile.size();
        for (int i = AStart; i < AStart + lengthTmp; i++) {
            xAIn[i] = vecVecInOneFile[i - AStart];
        }
        AStart += lengthTmp;

    }

//    std::cout<<"cellNum="<<cellNum<<std::endl;

    const auto tAEnd{std::chrono::steady_clock::now()};
    const std::chrono::duration<double> A_elapsed_secondsAll{tAEnd - tAStart};
    std::cout << "parse A time: " << A_elapsed_secondsAll.count() << " s" << std::endl;

    int counterA = 0;
    for (int i = 0; i < xAIn.size(); i += lag + 1) {
//        xASelected.push_back(xAIn[i]);
        xASelectedFlat.insert(xASelectedFlat.end(), xAIn[i].begin(), xAIn[i].end());
        counterA++;
    }
    arma_xA = ((arma::dmat(xASelectedFlat)).reshape(cellNum, counterA)).t();


//    std::cout<<"T="<<std::to_string(T)<<std::endl;
//    std::cout<<"A size=("<<arma_xA.n_rows<<", "<<arma_xA.n_cols<<")"<<std::endl;




    //B
    int BStart = 0;
    const auto tBStart{std::chrono::steady_clock::now()};
    for (const std::string &onexBFile: xBFilesSelected) {
        std::vector<std::vector<double>> vecVecInOneFile = readMsgBinVecVec(onexBFile);
//        xBIn.insert(xBIn.end(),vecVecInOneFile.begin(),vecVecInOneFile.end());
        int lengthTmp = vecVecInOneFile.size();
        for (int i = BStart; i < BStart + lengthTmp; i++) {
            xBIn[i] = vecVecInOneFile[i - BStart];
        }
        BStart += lengthTmp;

    }

    const auto tBEnd{std::chrono::steady_clock::now()};
    const std::chrono::duration<double> B_elapsed_secondsAll{tBEnd - tBStart};
    std::cout << "parse B time: " << B_elapsed_secondsAll.count() << " s" << std::endl;

    for (int i = 0; i < xBIn.size(); i += lag + 1) {
        xBSelectedFlat.insert(xBSelectedFlat.end(), xBIn[i].begin(), xBIn[i].end());
    }

    arma_xB = ((arma::dmat(xBSelectedFlat)).reshape(cellNum, counterA)).t();


//    arma::drowvec meanB=arma::mean(arma_xB,0);
//    std::cout<<"B size=("<<arma_xB.n_rows<<", "<<arma_xB.n_cols<<")"<<std::endl;
//    std::cout<<"-------------------------"<<std::endl;
}


///
/// @param filename bin file name
/// @return vector
std::vector<double> reader::readMsgBinVec(const std::string &filename) {
    std::ifstream infile(filename, std::ios::binary);
    if (!infile.is_open()) {
        throw std::runtime_error("Unable to open file for reading");
    }

    std::stringstream buffer;
    buffer << infile.rdbuf();
    std::string serialized_data = buffer.str();
    infile.close();

    // Step 2: Deserialize the binary data using MessagePack
    msgpack::object_handle oh = msgpack::unpack(serialized_data.data(), serialized_data.size());
    msgpack::object obj = oh.get();

    std::vector<double> vec;
    obj.convert(vec);

    return vec;


}

///
/// @param filename bin file name
/// @return vector<vector<double>>
std::vector<std::vector<double>> reader::readMsgBinVecVec(const std::string &filename) {

// Step 1: Read the binary data from the file
    std::ifstream infile(filename, std::ios::binary);
    if (!infile.is_open()) {
        throw std::runtime_error("Unable to open file for reading");
    }

    std::stringstream buffer;
    buffer << infile.rdbuf();
    std::string serialized_data = buffer.str();
    infile.close();

    // Step 2: Deserialize the binary data using MessagePack
    msgpack::object_handle oh = msgpack::unpack(serialized_data.data(), serialized_data.size());
    msgpack::object obj = oh.get();

    std::vector<std::vector<double>> nested_vec;
    obj.convert(nested_vec);

    return nested_vec;


}


///data to json, json as input to plot
void reader::data2json() {

    std::string jsonPath = this->TDir + "/jsonData/";

    //write U
    std::string UJsonPath = jsonPath + "/jsonU/";
    if (!fs::is_directory(UJsonPath) || !fs::exists(UJsonPath)) {
        fs::create_directories(UJsonPath);
    }

    std::string UJsonFile = UJsonPath + "/UData.json";

    boost::json::object objU;
    boost::json::array arrU;
    for (const auto &val: USelected) {
        arrU.push_back(val);
    }
    objU["U"] = arrU;
    std::ofstream ofsU(UJsonFile);
    std::string UStr = boost::json::serialize(objU);
    ofsU << UStr << std::endl;
    ofsU.close();


    //write xA, xB
    for (int i = 0; i < cellNum; i++) {
        std::string cellPathTmp = jsonPath + "jsonUnitCell" + std::to_string(i) + "/";
        if (!fs::is_directory(cellPathTmp) || !fs::exists(cellPathTmp)) {
            fs::create_directories(cellPathTmp);
        }
        boost::json::object obj_xAxB;

        std::string cellJsonFile = cellPathTmp + "xAxBData.json";

        arma::dcolvec colA = arma_xA.col(i);
        boost::json::array oneCol_xA;
        for (const auto &val: colA) {
            oneCol_xA.push_back(val);
        }
        arma::dcolvec colB = arma_xB.col(i);
        boost::json::array oneCol_xB;
        for (const auto &val: colB) {
            oneCol_xB.push_back(val);
        }
        obj_xAxB["xA"] = oneCol_xA;

        obj_xAxB["xB"] = oneCol_xB;

        std::ofstream ofsxAxB(cellJsonFile);
        std::string xAxBStr = boost::json::serialize(obj_xAxB);
        ofsxAxB << xAxBStr << std::endl;
        ofsxAxB.close();


    }


}



///compute the column means of arma_xA, arma_xB
void reader::colmeans(){
    this->E_xARow=arma::mean(arma_xA,0);
    this->E_xBRow=arma::mean(arma_xB,0);

    this->E_xACol=E_xARow.t();
    this->E_xBCol=E_xBRow.t();

    this->E_xA2=E_xACol*E_xARow;
    this->E_xB2=E_xBCol*E_xBRow;
//    std::cout<<"cellNum="<<cellNum<<std::endl;

//    E_xARow.print("mean xA:");
//    E_xBRow.print("mean xB");
//
//E_xA2.print("EA2:");
//    E_xB2.print("EB2:");


}

///compute correlation functions GAA
void reader::computeGAA() {

    arma::dmat YA = arma::zeros(cellNum, cellNum);

    int Q = arma_xA.n_rows;

    for (int q = 0; q < Q; q++) {
        arma::drowvec rowTmp = arma_xA.row(q);
        arma::dcolvec colTmp = rowTmp.t();
        YA += colTmp * rowTmp;

    }

    double QDB = static_cast<double>(Q);

    YA /= QDB;

//    YA.print("YA:");

    arma::dmat GAA = YA - E_xA2;

    std::string outGAA=TDir+"/GAA.csv";
    std::ofstream ofs(outGAA);

    printMat(GAA,ofs);
    ofs.close();


}


///compute correlation functions GAB
void reader::computeGAB(){

    arma::dmat YAB = arma::zeros(cellNum, cellNum);
    int Q = arma_xA.n_rows;

    for(int q=0;q<Q;q++){
        arma::drowvec rowATmp = arma_xA.row(q);
        arma::dcolvec colATmp=rowATmp.t();

        arma::drowvec rowBTmp=arma_xB.row(q);

        YAB+=colATmp*rowBTmp;

    }

    double QDB = static_cast<double>(Q);

    YAB/=QDB;

    arma::dmat GAB=YAB-E_xACol*E_xBRow;
    std::string outGAB=TDir+"/GAB.csv";
    std::ofstream ofs(outGAB);

    printMat(GAB,ofs);
    ofs.close();


}


///compute correlation functions GBB
void reader::computeGBB(){
    arma::dmat YB = arma::zeros(cellNum, cellNum);

    int Q = arma_xB.n_rows;

    for(int q=0;q<Q;q++){
        arma::drowvec rowTmp = arma_xB.row(q);
        arma::dcolvec colTmp = rowTmp.t();
        YB += colTmp * rowTmp;


    }
    double QDB = static_cast<double>(Q);
    YB/=QDB;

    arma::dmat GBB=YB-E_xB2;
    std::string outGBB=TDir+"/GBB.csv";
    std::ofstream ofs(outGBB);

    printMat(GBB,ofs);
    ofs.close();

}