#include <iostream>
#include <map>
#include <string>
#include <vector>
#include <fstream>
#include <cstdlib>
#include <algorithm>
#include <tuple>
#include <regex>
#include <unistd.h>

using namespace std;

struct Lib {
    string name;
    int width;
    int height;
};

struct Tech {
    string name;
    int numLib;
    int id;
    vector<Lib> libCells;
    map<string, Lib> libMap;
};

struct DieSize {
    unsigned long long width;
    unsigned long long height;
};

struct Die {
    string name;
    string techname;
    int util;
    int id;
};

struct Net;

struct Cell {
    string name;
    string libname;
    int id;
    vector<Net> netlist;
};

struct Net {
    string name;
    int numCell;
    int id;
    vector<Cell> netcells;
    int dist[2]={0,0};
};
//依面積由大到小排序 相同面積net數少的排前面
bool cmp_area(const tuple<int, int,int>& ca, const tuple<int, int,int>& cb) {
    int ca2 = get<1>(ca);
    int cb2 = get<1>(cb);
    if(ca2==cb2){
        int ca3 = get<2>(ca);
        int cb3 = get<2>(cb);
        return ca3 < cb3;
    }
    return ca2 > cb2; 
}
//依cell連接的net數由小到大排序 相同net數面積大的排前面
bool cmp_nets(const tuple<int, int,int>& ca, const tuple<int, int,int>& cb) {  //0:id 1:area 2:netnum
    int ca2 = get<2>(ca);
    int cb2 = get<2>(cb);
    if(ca2==cb2){
        int ca1 = get<1>(ca);
        int cb1 = get<1>(cb);
        return ca1 > cb1;
    }
    return ca2 < cb2; 
}
//依cell id由小到大排序
bool cmp_id(const tuple<int, int,int>& ca, const tuple<int, int,int>& cb) {
    int ca1 = get<0>(ca);
    int cb1 = get<0>(cb);
    return ca1 < cb1; 
}

int main(int argc, char *argv[]){
    bool flag = false;
    int count = 0;
    ifstream inputFile(argv[1]);

    if (!inputFile) {
        cerr << "Failed to open the input file." << endl;
        return 1;
    }

    string s;
    int numTechs;
    inputFile >> s >> numTechs;
    
    //紀錄每個tech裡面的lib資訊
    map<string, Tech> techMap;

    for (int i = 0; i < numTechs; i++) {
        Tech tech;
        inputFile >> s >> tech.name >> tech.numLib;
        tech.id = i;
        for (int j = 0; j < tech.numLib; j++) {
            Lib libCell;
            inputFile >> s >> libCell.name >> libCell.width >> libCell.height;
            tech.libCells.push_back(libCell);
            tech.libMap[libCell.name] = libCell;
        }
        techMap[tech.name] = tech;
    }
   
    string line;
    getline(inputFile, line); 
    if (line.empty()){} // Skip blank lines
        
    DieSize diesize;
    inputFile >> s >> diesize.width >> diesize.height;

    //紀錄每個die使用哪個tech
    map<int, Die> dieMap;

    for (int i = 0; i < 2; i++) {
        Die die;
        inputFile >> die.name >> die.techname >> die.util;
        die.id = i;
        dieMap[die.id] = die;
    }

    getline(inputFile, line); 
    if (line.empty()) {} // Skip blank lines

    //紀錄cell屬於哪個lib
    int numCells;
    inputFile >> s >> numCells;
    map<string, Cell> cellMap; //可用cell name查詢
    map<int, Cell> cellidMap; //可用cell id查詢
    vector<Cell> cells;

    for (int i = 1; i <=numCells; i++) {
        Cell cell;
        inputFile >> s >> cell.name >> cell.libname;
        cell.id = i;
        //cout<<"cell_id= " <<cell.id <<endl;
        cells.push_back(cell);
        cellMap[cell.name] = cell;
        cellidMap[cell.id] = cell;
    }
    
    getline(inputFile, line); 
    if (line.empty()) {} // Skip blank lines

    //紀錄每個net包含哪些cell
    int numNets;
    inputFile >> s >> numNets;
    map<string, Net> netMap;
    vector<Net> nets;
    for (int i = 0; i < numNets; i++) {
        Net net;
        inputFile >> s >> net.name >> net.numCell;
        net.id = i;
        for (int j = 0; j < net.numCell; j++) {
            Cell cell;
            inputFile >> s >> cell.name;
            Cell& cellInfo = cellMap.at(cell.name);
            cellInfo.netlist.push_back(net);
            int cid = cellInfo.id;
            //cout<<"cid="<<cid<<endl;
            Cell& cellId = cellidMap.at(cid);
            cellId.netlist.push_back(net);
            net.netcells.push_back(cell);
            //cout<<"cellname="<<cellId.name<<endl;
        }
        netMap[net.name] = net;
        nets.push_back(net);
    }
    inputFile.close();
    
    //construct vwgts 用DieA的tech對應的libcell長寬計算
    unsigned long long  die_area = diesize.width*diesize.height;
    int *vwgts = new int[numCells];
    Tech tech;
    Lib lib;
    tech = techMap[dieMap[0].techname];
    for(int i=0; i<numCells; i++){
        lib = tech.libMap[cells[i].libname];
        vwgts[i] = lib.width * lib.height;
        if(vwgts[i]>= die_area *0.02){   //看有沒有單一cell面積大於die面積的2% 會影響後面sorting選擇
            flag = true;
            count +=1 ;
        }
    }

    //準備shmetis需要的inputfile hyperedge graph
    ofstream HGraphFile("hyper.hgr");
    
    HGraphFile << numNets << " " << numCells << " "<<"10" <<endl;
    for (const auto& net : nets) {
        for (const Cell& cells : net.netcells) {
            if(cellMap.find(cells.name) != cellMap.end()){
                const Cell& cellInfo = cellMap.at(cells.name);
                //cout << "Cell libame: " << cellInfo.libname << endl;
                HGraphFile << cellInfo.id << " ";
            }
        }
        HGraphFile<< endl;
    }
    for(int i=0;i<numCells; i++){
        HGraphFile<< vwgts[i]*0.1 << endl;
    }
    HGraphFile.close();
    
    /*FILE *pipe = popen("../src/shmetis hyper.hgr 2 1", "r");
    if (!pipe) {
        cerr << "Error opening pipe to hmetis" << endl;
        return 1;
    }*/
    char cwd[1024];
    int hyperedgeCut = -1;
    int part0 = -1;
    int part1 = -1;
    string command;
    if (getcwd(cwd, sizeof(cwd)) != NULL) {
        string hgrfile = string(cwd) + "/hyper.hgr"; 
        const char* shmetisEct = "../src/shmetis";  
        //string shmetispath = std::string(cwd) + "/" + shmetisEct;

        command = string(shmetisEct) + " " + hgrfile + " 2 1";
    cout<<command<<endl;
        FILE* pipe = popen(command.c_str(), "r");
        if (!pipe) {
            cerr << "Error opening pipe to shmetis" << endl;
            return 1;
        }
    

    char buffer[128];
    regex pattern("Hyperedge Cut:\\s+(\\d+)");

    while (fgets(buffer, sizeof(buffer), pipe) != nullptr) {
        string line(buffer);
        smatch hmatch;
        if (regex_search(line, hmatch, pattern)) {
            hyperedgeCut = stoi(hmatch[1]);
        }else if (line.find("Partition Sizes & External Degrees:") != string::npos) {
            if (fgets(buffer, sizeof(buffer), pipe) != nullptr) {
                sscanf(buffer, " %d[%*d] %d[%*d]", &part0, &part1);
            }
        }
    }


    /*cout << "Hyperedge Cut: " << hyperedgeCut << endl;
    cout << "Partition Size 1: " << part0 << endl;
    cout << "Partition Size 2: " << part1 << endl;*/

    // Close the pipe
    pclose(pipe);
    }

    string part2FilePath = string(cwd) + "/hyper.hgr.part.2";
    ifstream part2File("hyper.hgr.part.2");

    //紀錄每個cell被分到哪一邊
    int partition[numCells];
    for (int i = 0; i <numCells; i++) {
        part2File >> partition[i];
    }
    
    //util高的die就對應shmetis分出來面積大那一邊
    string part0name, part1name;
    if(part0>=part1) {
        part0name = dieMap[0].util>=dieMap[1].util? dieMap[0].name:dieMap[1].name;
        part1name = dieMap[0].util<dieMap[1].util? dieMap[0].name:dieMap[1].name;
    }else{
        part0name = dieMap[0].util<dieMap[1].util? dieMap[0].name:dieMap[1].name;
        part1name = dieMap[0].util>=dieMap[1].util? dieMap[0].name:dieMap[1].name;
    }
    //cout<< part0name << part1name;
    //k代表dieA那一邊
    int k = (part0name==dieMap[0].name)? 0:1;
    //cout<<"k= "<<k<<endl;
    unsigned long long sum_areaA=0; //加總被分到part0的vertex weight
    unsigned long long sum_areaB=0;
    
     //紀錄每個die上面的cell id
    vector<tuple<int,int,int>> DiecellsA;  //<id,weight,# of nets>
    vector<tuple<int,int,int>> DiecellsB;
    int vwgtsb[numCells];
    int weight = 0;

    //根據cell所在的die計算它的面積 以及dieA dieB目前總面積
    for(int i=0; i<numCells; i++){
        if(partition[i]==k){
            sum_areaA += vwgts[i];
            Cell& cell = cellidMap.at(i+1);
            int netnum = cell.netlist.size();
            DiecellsA.push_back(make_tuple(i+1,vwgts[i], netnum));
        }else{
            if(numTechs==1){weight=vwgts[i];}   //如果tech只有一種就不會重算
            else{
            Tech tech;
            Lib lib;
            tech = techMap[dieMap[1].techname];
            lib = tech.libMap[cells[i].libname];
            weight = lib.width * lib.height;
            }
            sum_areaB += weight;
            Cell& cell = cellidMap.at(i+1);
            int netnum = cell.netlist.size();
            DiecellsB.push_back(make_tuple(i+1,weight, netnum));
        }
    }

    //計算每個cell在dieB的面積
    for(int i=0; i<numCells; i++){
        if(numTechs==1){vwgtsb[i]=vwgts[i];}
        else{
            Tech tech;
            Lib lib;
            tech = techMap[dieMap[1].techname];
            lib = tech.libMap[cells[i].libname];
            weight = lib.width * lib.height;
            vwgtsb[i]= weight;
        }
    }
  
    //計算每個net在兩個die的distribution
    for (int i = 0; i < numCells; i++) {
    Cell &cell = cellidMap[i + 1]; 
    int cellp = partition[i]; 
    //cout<<"cellid= "<<i+1<<",partition= "<<cellp<<endl;
        for (auto &net : cell.netlist) {
            Net& netInfo = netMap.at(net.name);
            if (cellp == k) {
                netInfo.dist[0] += 1;
            }else {
                netInfo.dist[1] += 1;
            }
        //cout<<"cellid= "<<i+1<<"netid= "<<net.id<<"in a:"<<netInfo.dist[0]<<" in b: "<<netInfo.dist[1]<<endl;
        }
    }
    int num_DieA = DiecellsA.size();
    int num_DieB = DiecellsB.size();
  
    //計算兩個die目前的util
    double utilA = sum_areaA/(double)die_area;
    double utilB = sum_areaB/(double)die_area;
    //cout<<utilA << utilB<<endl;

    //如果有超過10個cell的單一面積大於整個die的2% 就依numnet由小到大排序
    //相同net數中面積較大的排前面
    if(flag==true && count>=10){
        sort(DiecellsA.begin(), DiecellsA.end(), cmp_nets);
        sort(DiecellsB.begin(), DiecellsB.end(), cmp_nets);
    }else{
        sort(DiecellsA.begin(), DiecellsA.end(), cmp_area); //一般情況依面積由大到小排序
        sort(DiecellsB.begin(), DiecellsB.end(), cmp_area);
    }
    
    //讓兩邊符合util規定
    while(utilA > dieMap[0].util*0.01 || utilB > dieMap[1].util*0.01){

        while(utilB > dieMap[1].util*0.01){ //取面積最大但net數量較少的換到另一邊
            int gains = 0;
            tuple<int,int,int> chtup = DiecellsB[0];
            int cid = get<0>(chtup);
            //cout<<"cid= "<<cid<<endl;
            Cell& cellInfo = cellidMap.at(cid);
            for(auto& net: cellInfo.netlist){  //<Net>
                Net& netInfo = netMap.at(net.name);
                int gain = 0;
                int ba = netInfo.dist[0];
                int bb = netInfo.dist[1];
                if(ba==0) gain=-1;
                if(bb==1) gain=1;
                netInfo.dist[0]+=1;
                netInfo.dist[1]-=1;
                partition[cid-1] ^=1;
                gains += gain;
            }
            DiecellsA.push_back(chtup);
            DiecellsB.erase(DiecellsB.begin());
            hyperedgeCut -= gains; 
            //update area and cutsize
            sum_areaA += vwgts[cid-1];
            sum_areaB -= vwgtsb[cid-1];
            utilA = sum_areaA/(double)die_area;
            utilB = sum_areaB/(double)die_area;
        }
        while(utilA > dieMap[0].util*0.01){ //取面積最大但net數量較少的換到另一邊
            int gains = 0;
            tuple<int,int,int> chtup = DiecellsA[0];
            int cid = get<0>(chtup);
            Cell& cellInfo = cellidMap.at(cid);
            for(auto& net: cellInfo.netlist){  //<Net>
                Net& netInfo = netMap.at(net.name);
                int gain = 0;
                int ba = netInfo.dist[0];
                int bb = netInfo.dist[1];
                if(bb==0) gain=-1;
                if(ba==1) gain=1;
                netInfo.dist[1]+=1;
                netInfo.dist[0]-=1;
                partition[cid-1] ^=1;
                gains += gain;
            }
            DiecellsB.push_back(chtup);
            DiecellsA.erase(DiecellsA.begin());
            hyperedgeCut -= gains; 
            //update area and cutsize
            sum_areaA -= vwgts[cid-1];
            sum_areaB += vwgtsb[cid-1];
            utilA = sum_areaA/(double)die_area;
            utilB = sum_areaB/(double)die_area;
        }
    }
    cout<<"utilA= "<< utilA <<"utilB= "<< utilB<<endl; 
    cout<<"cutsize= "<<hyperedgeCut<<endl;

    //計算DieA DieB各有多少個cells
    num_DieA = DiecellsA.size();
    num_DieB = DiecellsB.size();
   
    // 依據cell id由小到大排序準備寫檔
    sort(DiecellsA.begin(), DiecellsA.end(), cmp_id);
    sort(DiecellsB.begin(), DiecellsB.end(), cmp_id);

    //output to the outputfile
    ofstream outputFile(argv[2]);
    outputFile<<"CutSize "<< hyperedgeCut<<endl;
    outputFile<<"DieA "<<num_DieA<<endl;
    //outputFile<<command<<endl; 
    //out put cells name
    for(auto&cells: DiecellsA){
        int cid = get<0>(cells);
        Cell cell;
        cell = cellidMap.at(cid);
        outputFile<<cell.name<<endl;
    }
    outputFile<<"DieB "<<num_DieB<<endl;
    for(auto&cells: DiecellsB){
        int cid = get<0>(cells);
        Cell cell;
        cell = cellidMap.at(cid);
        outputFile<<cell.name<<endl;
    }
    outputFile.close();
    return 0;
}
