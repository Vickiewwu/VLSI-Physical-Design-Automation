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
#include "function.h"
#include "struct.h"

using namespace std;
unsigned int chipwidth;
unsigned int chipheight;
map<string,unsigned int> mixedmap;
map<string,unsigned int> smap;
vector<unsigned int>hcontour(100000,0);
vector<unsigned int>lcontour(100000,0);
vector<unsigned int>tcontour(100000,0);
void simulated_annealing(BST* bst, Chip* chip, std::vector<SoftModule*>smodvec,std::vector<FixedModule*>fmodvec,std::vector<Net*> nets,char* name);
void writeoutput(char* outputfilename,BST bst);
vector<Net*> nets;
vector<FixedModule*> fmods;

int main(int argc, char *argv[]){

    ifstream inputFile(argv[1]);
    if (!inputFile) {
        cerr << "Failed to open the input file." << endl;
        return 1;
    }

    string s;
    Chip* chip = new Chip; 
    inputFile >> s >> chip->width >> chip->height ;
    chipwidth = chip->width;
    chipheight = chip->height;

    string line;
    getline(inputFile, line); 
    if (line.empty()){} // Skip blank lines

    int numSMod;
    inputFile >> s >> numSMod;
    vector<SoftModule*> smods;
    
    for(int i=1; i<=numSMod; i++){
        SoftModule* smod = new SoftModule;
        inputFile >> s >> smod->name >> smod->min_area;
        smod->id = i;
        smod->x1 = smod->y1 = 0;
        smod->area = 0;
        smod->getwh();
        smod->getx2y2();
        smod->getmiddle();
        smods.push_back(smod);
        smap[smod->name]=i;
        mixedmap[smod->name]=i;
    }
    
    getline(inputFile, line); 
    if (line.empty()){} // Skip blank lines

    int numFMod;
    inputFile >> s >> numFMod;
   
    for(int i=1; i<=numFMod; i++){
        FixedModule* fmod = new FixedModule;
        inputFile >> s >> fmod->name >> fmod->x1 >>fmod->y1 >>fmod->width>>fmod->height;
        fmod->x2 = fmod->x1+fmod->width;
        fmod->y2 = fmod->y1+fmod->height;
        fmod->area = fmod->width*fmod->height;
        fmod->id = i+numSMod;
        fmod->getmiddle();
        fmods.push_back(fmod);
        mixedmap[fmod->name]=i+numSMod;
    }

    
    getline(inputFile, line); 
    if (line.empty()){} // Skip blank lines

    int numNet;
    inputFile >> s >> numNet;
    
    for(int i=0; i<numNet; i++){
        Net* net = new Net;
        net->id = i;
        inputFile >> s >> net->m1 >>net->m2 >>net->weight;
        nets.push_back(net);
    }
    inputFile.close();
        
    BST bst;
    bst.init(smods,fmods,chip->width,chip->height,numSMod); 
    cout<<"Finishing init bstree.\n";
    bst.save();
    
    srand(1);
    simulated_annealing(&bst,chip,smods,fmods,nets,argv[2]);
    cout<<"Finishing sa.\n";

    return 0;
}
