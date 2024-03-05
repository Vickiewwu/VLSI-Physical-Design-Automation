#include <iostream>
#include <fstream>
#include <map>
#include <string>
#include <vector>
#include <cstdlib>
#include <algorithm>
#include <tuple>
#include <regex>
#include <unistd.h>
#include <random>
#include <iterator> 
//#include <ctime>
#include "function.h"
#include "struct.h"

using namespace std;
extern unsigned int chipwidth;
extern unsigned int chipheight;
extern std::vector<Net*> nets;
extern std::vector<FixedModule*> fmods;
extern map<string,unsigned int> smap;
//extern std::map<std::string,unsigned int> mixedmap;
void writeoutput(char* outputfilename,BST* bst);
//std::vector<unsigned int>contour;

/*bool cmpx(const FixedModule& a, const FixedModule& b) {
    if(a.x==b.x){
        return a.y<b.y;
    }
    return a.x < b.x;
}
bool cmpy(const FixedModule& a, const FixedModule& b) {
    if(a.y==b.y){
        return a.x<b.x;
    }
    return a.y < b.y;
}*/
/*void readfile(const char* inputfile)
{
    ifstream inputFile(inputfile);

    if (!inputFile) {
        cerr << "Failed to open the input file." << endl;
    }

    string s;
    Chip chip; 
    inputFile >> s >> chip.width >> chip.height ;

    string line;
    getline(inputFile, line); 
    if (line.empty()){} // Skip blank lines

    int numSMod;
    inputFile >> s >> numSMod;
    vector<string>smods;
    vector<int>smodsid;
    map<string,SoftModule>smodmap;
    map<int,SoftModule>smodidmap;
    for(int i=0; i<numSMod; i++){
        SoftModule smod;
        inputFile >> s >> smod.name >> smod.min_area;
        smod.id = i;
        smods.push_back(smod);
        smodsid.push_back(i);
        smodmap[smod.name] = smod;
        smodidmap[i] = smod;
    }

    
    getline(inputFile, line); 
    if (line.empty()){} // Skip blank lines

    int numFMod;
    inputFile >> s >> numFMod;
    vector<string>fmods;
    map<string,FixedModule>fmodmap;
    map<int,FixedModule>fmodidmap;
    for(int i=0; i<numFMod; i++){
        FixedModule fmod;
        inputFile >> s >> fmod.name >> fmod.x1 >>fmod.y1 >>fmod.width>>fmod.height;
        fmod.id = i;
        fmods.push_back(fmod);
        fmodmap[fmod.name] = fmod;
        fmodidmap[i] = fmod;
    }
    sort(fmods.begin(), fmods.end(), cmpx);
    sort(fmods.begin(), fmods.end(), cmpy);


    
    getline(inputFile, line); 
    if (line.empty()){} // Skip blank lines

    int numNet;
    inputFile >> s >> numNet;
    for(int i=0; i<numNet; i++){
        Net net;
        string m;
        inputFile >> s >> net.m1 >>net.m2;
    }
    inputFile.close();
}*/

/*floorplan(char* inputfile){

    //readfile(inputfile);

    simulated_annealing();

    writefile(outputfile);
}*/
void BST::init(std::vector<SoftModule*>smodvec,std::vector<FixedModule*>fmodvec,unsigned int chipw,unsigned chiph,unsigned int numsmod){
    this->numSMod = numsmod;
    
    this->bstree.push_back(Node()); //head
    this->smodlist.push_back(NULL);
    unsigned int size = smodvec.size();
    //cout<<"size="<<size<<endl;
    unsigned int index = 1;
    for(unsigned int i=0;i<size;i++){
        SoftModule* cur = smodvec[i];
        this->smodlist.push_back(cur);
        unsigned int parent = index/2;
        if(index==parent*2){  //left child
            this->bstree[parent].left = index;
        }else{ //root is the right child of head
            this->bstree[parent].right = index;
        }
        this->bstree.push_back(Node(parent,0,0,cur->width,cur->height));
        this->idsmap.insert(std::pair<int,int>(cur->id,index));
        index++;
    }
    //cout<<"smodlist.size="<<this->smodlist.size()<<endl;
    
    unsigned int bottombound = 0;
    unsigned int leftbound = 0;
    unsigned int topbound = chiph;
    unsigned int rightbound = chipw;
    //cout<<"leftbound="<<leftbound<<"rightbound="<<rightbound<<endl;
    //cout<<"topbound="<<topbound<<"bottombound="<<bottombound<<endl;
    index = 0;
    for(unsigned int i=0;i<fmodvec.size();i++){
        FixedModule* cur = fmodvec[i];
        this->fmodlist.push_back(cur);
        //this->bstree.pushback(Node(0,0,0,cur->width,cur->height));
        this->idfmap.insert(std::pair<int,int>(cur->id,index));
        if(cur->y1==0){
            if(cur->height>bottombound){
                bottombound = cur->height;
            }
        }
        if(cur->x1==0){
            if(cur->width>leftbound){
                leftbound = cur->width;
            }
        }
        else if(cur->x2==chipw){
            if(cur->x1 < rightbound){
                rightbound = cur->x1;
            }
        }
        else if(cur->y2==chiph){
            if(cur->y1 < topbound){
                topbound = cur->y1;
            }
        }
        index++;
    }
    //cout<<"after leftbound="<<leftbound<<"rightbound="<<rightbound<<endl;
    //cout<<"topbound="<<topbound<<"bottombound="<<bottombound<<endl;
    this->origx = leftbound;
    this->origy = bottombound;
    this->newcwidth = rightbound-leftbound;
    this->newcheight = topbound-bottombound;
    this->topbound = topbound;
    this->currmaxx = 0;
    this->currmaxy = 0;
    this->fwidth = 0;
    this->fheight = 0;
    //cout<<"init fwidth="<<this->fwidth<<"fheight="<<this->fheight<<endl;
    this->placed = 0;
    this->tofloorplan();
    
}
void BST::save(){
    savedbst.resize(bstree.size());
    auto first = bstree.begin();
    auto last = bstree.end();
    std::copy(first,last,savedbst.begin());
} 

void BST::savebest(){
    bestbst.resize(bstree.size());
    auto first = bstree.begin();
    auto last = bstree.end();
    std::copy(first,last,bestbst.begin());
} 

void BST::loadprev(){
    bstree.resize(savedbst.size());
    auto first = savedbst.begin();
    auto last = savedbst.end();
    std::copy(first,last,bstree.begin());
    this->placed = false;
} 

void BST::loadbest(){
    bstree.resize(bestbst.size());
    auto first = bestbst.begin();
    auto last = bestbst.end();
    std::copy(first,last,bstree.begin());
    this->placed = false;
} 

void BST::rotateNode(SoftModule* smod){
    unsigned int index = idsmap[smod->id];
    bstree[index].rotate();
    smodlist[index]->getx2y2();
    
    this->placed = false;
}

void BST::moveNode(SoftModule* smod1,SoftModule* smod2){
    
    /*unsigned int moveindex = idsmap[smod1->id];
    unsigned int destindex = idsmap[smod2->id];
    if(moveindex==destindex) return;
    deleteNode(moveindex);
    cout<<"here a"<<endl;
    unsigned int direction = rand()%2;
    cout<<"here b"<<endl;
    insertNode(moveindex,destindex,direction);
    cout<<"here c"<<endl;*/

    unsigned int moveid = smod1->id;
    unsigned int destid = smod2->id;
    //cout<<"moveID"<<moveid<<"destid="<<destid<<endl;
    if(moveid==destid) return;
    //cout<<"move 1"<<endl;
    this->deleteNode(moveid);
    //cout<<"move 2"<<endl;
    unsigned int direction = rand()%2;
    //cout<<"move 3"<<endl;
    this->insertNode(moveid,destid,direction);
    //cout<<"move 4"<<endl;
    this->placed = false;
}

void BST::deleteNode(unsigned int id){
    //no child
    unsigned int index1, index2;
    
    if(bstree[id].left==0 && bstree[id].right ==0){
        unsigned int parent = bstree[id].parent;
        if(bstree[parent].left == id){
            bstree[parent].left = 0;
        }
        else if(bstree[parent].right == id){
            bstree[parent].right = 0;
        }
        
    }
    //2 children
    else if(bstree[id].left!=0 && bstree[id].right !=0){
        bool direction;
        while(bstree[id].left!=0 || bstree[id].right!=0 ){
            if(bstree[id].left!=0 && bstree[id].right !=0){
                direction = rand()%2; //0:left,1:right
            }
            else if(bstree[id].left>0){
                direction = false;
            }
            else{
                direction = true;
            }
            index1 = this->idsmap[id];
            if(direction==0){
              
                index2 = this->idsmap[bstree[id].left];
                this->swapNode(this->smodlist[index1],this->smodlist[index2]);
                
            }else{
                index2 = this->idsmap[bstree[id].right];
                this->swapNode(this->smodlist[index1], this->smodlist[index2]);
                
            }
            
        }
        unsigned int parent = bstree[id].parent;
        if(bstree[parent].left == id){
            bstree[parent].left = 0;
        }
        else if(bstree[parent].right == id){
            bstree[parent].right = 0;
        }
        
    }
    //1 child
    else{
        
        unsigned int child;
        if(bstree[id].left >0){
            child = bstree[id].left;
        }
        else {
            child = bstree[id].right;
        }
        
        unsigned int parent = bstree[id].parent;
        if(parent>0){
            if(bstree[parent].left == id){
            bstree[parent].left = child;
            }
            else if(bstree[parent].right == id){
                bstree[parent].right = child;
            }
        }
        bstree[child].parent = parent;
        /*if(index==root){
            root = child;
        }*/
        
    }
}
void BST::insertNode(unsigned int moveid, unsigned int destid,bool direction){
    //left
    if(direction==0){
        //cout<<"destid.left="<<bstree[destid].left<<endl;
        if(bstree[destid].left >0 && bstree[destid].left <=numSMod){
            unsigned int leftchild = bstree[destid].left;
            bstree[leftchild].parent=moveid;
            bstree[moveid].left = bstree[destid].left;
        }
        else{
            bstree[moveid].left = 0;
        }
        //cout<<"move.left="<<bstree[moveid].left<<endl;
        bstree[moveid].right = 0;
        bstree[moveid].parent = destid;
        bstree[destid].left = moveid;
    }
    else{
        //cout<<"destid.right="<<bstree[destid].right<<endl;
        if(bstree[destid].right >0 && bstree[destid].left <=numSMod){
            unsigned int rightchild = bstree[destid].right;
            bstree[rightchild].parent=moveid;
            bstree[moveid].right = bstree[destid].right;
        }
        else{
            bstree[moveid].right = 0;
        }
        //cout<<"move.right="<<bstree[moveid].right<<endl;
        bstree[moveid].left = 0;
        bstree[moveid].parent = destid;
        bstree[destid].right = moveid;
    }
   
}

void BST::swapNode(SoftModule* smod1,SoftModule* smod2){
    //unsigned int index1 = idsmap[smod1->id];
    //nsigned int index2 = idsmap[smod2->id];
    unsigned int id1 = smod1->id;
    unsigned int id2 = smod2->id;
    unsigned int right1 = bstree[id1].right;
    unsigned int right2 = bstree[id2].right;
    unsigned int left1 = bstree[id1].left;
    unsigned int left2 = bstree[id2].left;
    unsigned int parent1 = bstree[id1].parent;
    unsigned int parent2 = bstree[id2].parent;

    if(bstree[parent1].left == id1){
        bstree[parent1].left = id2;
    }
    else{
        bstree[parent1].right = id2;
    }
    if(bstree[parent2].left == id2){
        bstree[parent2].left = id1;
    }
    else{
        bstree[parent2].right = id1;
    }

    if(left1 >0){
        bstree[left1].parent = id2;
    }
    if(left2 >0){
        bstree[left2].parent = id1;
    }
    if(right1 >0){
        bstree[right1].parent = id2;
    }
    if(right2 >0){
        bstree[right2].parent = id1;
    }

    bstree[id1].parent = parent2;
    bstree[id1].right = right2;
    bstree[id1].left = left2;
    bstree[id2].parent = parent1;
    bstree[id2].right = right1;
    bstree[id2].left = left1;

    if(bstree[id1].parent==id1){
        bstree[id1].parent = id2;
    }
    else if(bstree[id1].left==id1){
        bstree[id1].left= id2;
    }
    else if(bstree[id1].right==id1){
        bstree[id1].right= id2;
    }
    if(bstree[id2].parent==id2){
        bstree[id2].parent = id1;
    }
    else if(bstree[id2].left==id2){
        bstree[id2].left= id1;
    }
    else if(bstree[id2].right==id2){
        bstree[id2].right= id1;
    }
    this->placed = false;
}
void BST::resizeNode(SoftModule* smodt){
    unsigned int index = idsmap[smodt->id];
    double w,h,ratio;
    if(index==0) return;
    w = ceil(smodlist[index]->width*1.00001);
    h = ceil(smodlist[index]->min_area/(double)w);
    ratio = (double)h/smodlist[index]->width;
    if(ratio< 0.5 || ratio >=1.98){
        if(w/2+1>h){
            smodlist[index]->width = w;
            smodlist[index]->height = w/2+1;
        }
    }else{
        smodlist[index]->width = w;
        smodlist[index]->height = h;
    }
    
    //cout<<"inresize, w="<<smodlist[index]->width<<", h="<<smodlist[index]->height<<endl;
    smodlist[index]->getx2y2();
    if(smodlist[index]->x2 > this->currmaxx){
        this->currmaxx = smodlist[index]->x2;
        this->fwidth = this->currmaxx - this->origx;
    }
     if(smodlist[index]->y2 > this->currmaxy){
        this->currmaxy = smodlist[index]->y2;
        this->fheight = this->currmaxy - this->origy;
    }
    this->placed = false;
}
/*void BST::resizeNode1(SoftModule* smodt){
    unsigned int index = idsmap[smodt->id];
    //Node& cur = bstree[smodt->id];
    unsigned int ex2 = this->newcwidth;
    unsigned int dx2 = 0;
    unsigned int ay2 = this->newcheight;
    unsigned int cy2 = 0;
    cout << "smodt->id="<<smodlist[index]->id<<endl;
    for(auto&smod:this->smodlist){
        if(smod ==NULL) continue;
        cout << "smod->id="<<smod->id<<" smod->x2= "<< smod->x2<<endl;
        if (smod->x2>smodlist[index]->x2 && smod->x2<ex2){
            cout << "smodt->x2= "<< smodlist[index]->x2<< "ex2="<<ex2<<endl;
            ex2 = smod->x2;
            
        }
        if (smod->x2<smodlist[index]->x2 && smod->x2>dx2){
            cout << "smodt->x2= "<< smodlist[index]->x2<< "dx2="<<dx2<<endl;
            dx2 = smod->x2;
        }
        if (smod->y2>smodlist[index]->y2 && smod->y2<ay2){
            cout << "smodt->y2= "<< smodlist[index]->y2<< "ay2="<< ay2<<endl;
            ay2 = smod->y2;
        }
        if (smod->y2<smodlist[index]->y2 && smod->y2>cy2){
            cout << "smodt->y2= "<< smodlist[index]->y2<< "cy2="<<cy2<<endl;
            cy2 = smod->y2;
        }
    }
    cout << "ex2= "<< ex2<< "dx2 ="<<dx2 << "ay2="<<ay2<< "cy2="<<cy2<<endl;
    bool Rflag = 1;
    bool Lflag = 1;
    bool Bflag = 1;
    bool Tflag = 1;
    unsigned int R = 0;
    if(ex2>(smodlist[index]->x1)){
        R = ex2-(smodlist[index]->x1);
        Rflag = 0;
    }
    unsigned int L = 0;
    if(dx2>(smodlist[index]->x1)){
        L = dx2-(smodlist[index]->x1);
        Lflag = 0;
    }
    
    unsigned int T = 0;
    if(ay2>(smodlist[index]->y1)){
        T = ay2-(smodlist[index]->y1);
        Tflag = 0;
    }
    unsigned int B = 0;
    if(cy2>(smodlist[index]->y1)){
        B = cy2-(smodlist[index]->y1);
        Bflag = 0;
    }
    unsigned int h,w;
    double ratio;
    cout << "R= "<< R<< "L="<<L<< "T="<<T<< "B="<<B<<endl;  
    srand(time(NULL)); 
    unsigned int which = rand()%4;
    switch(which){
    case 0:
        if(!Rflag){
            smodlist[index]->width = R;
            h = ceil(smodlist[index]->min_area/(double)R);
            ratio = R/(double)h;
            if(ratio>=0.5 && ratio<2){
                smodlist[index]->height = h;
            }
            else if(R>=h){
                smodlist[index]->height = R/2;
            }
            else{
                smodlist[index]->height = R*2;
            }
        }
        break;
    case 1:
        if(!Lflag){
            smodlist[index]->width = L;
            h = ceil(smodlist[index]->min_area/(double)L);
            ratio = L/(double)h;
            if(ratio>=0.5 && ratio<2){
                smodlist[index]->height = h;
            }
            else if(L>=h){
                smodlist[index]->height = L/2;
            }
            else{
                smodlist[index]->height = L*2;
            }
        }
        break;
    case 2:
        if(!Tflag){
            smodlist[index]->height = T;
            w = ceil(smodlist[index]->min_area/(double)T);
            ratio = T/(double)w;
            if(ratio>=0.5 && ratio<2){
                smodlist[index]->width = w;
            }
            else if(T>=w){
                smodlist[index]->width = T/2;
            }
            else{
                smodlist[index]->width = T*2;
            }
        }
        break;
    case 3:
        if(!Bflag){
            smodlist[index]->height = B;
            w = ceil(smodlist[index]->min_area/(double)B);
            ratio = B/(double)w;
            if(ratio>=0.5 && ratio<2){
                smodlist[index]->width = w;
            }
            else if(B>=w){
                smodlist[index]->width = B/2;
            }
            else{
                smodlist[index]->width = B*2;
            }
        }
        break;
    }
    smodlist[index]->getx2y2();
    if(smodlist[index]->x2 > this->currmaxx){
        this->currmaxx = smodlist[index]->x2;
        this->fwidth = this->currmaxx - this->origx;
    }
     if(smodlist[index]->y2 > this->currmaxy){
        this->currmaxy = smodlist[index]->y2;
        this->fheight = this->currmaxy - this->origy;
    }
 
}*/
unsigned int getwirelength(BST* bst){
    //cout << "here j "<< endl;
    unsigned int wl =0;
    unsigned int mx1, mx2, my1, my2;  //middle coordinate of the module
    for(unsigned int i=1;i<=bst->numSMod;i++){
        bst->smodlist[i]->getmiddle();
        //cout<<" i= "<<i<<endl;
    }
    for(auto&net:nets){
        //cout << "netid= "<< net->id<< endl;
        unsigned int id1 = mixedmap[net->m1];
        unsigned int id2 = mixedmap[net->m2];
        //cout << "id1= "<< id1<<"id2="<<id2<< endl;
        if(id1>=(bst->numSMod+1)){
            unsigned int index1 = bst->idfmap[id1];
            mx1 = bst->fmodlist[index1]->middle_x;
            //cout << "mx1= "<< mx1<< endl;
            my1 = bst->fmodlist[index1]->middle_y;
            //cout << "my1= "<< my1<< endl;
        }else{
            unsigned int index1 = bst->idsmap[id1];
            mx1 = bst->smodlist[index1]->middle_x;
            //cout << "mx1= "<< mx1<< endl;
            my1 = bst->smodlist[index1]->middle_y;
            //cout << "my1= "<< my1<< endl;
        }
        if(id2>=(bst->numSMod+1)){
            unsigned int index2= bst->idfmap[id2];
            mx2 = bst->fmodlist[index2]->middle_x;
            //cout << "mx2= "<< mx2<< endl;
            my2 = bst->fmodlist[index2]->middle_y;
            //cout << "my2= "<< my2<< endl;
        }else{
            unsigned int index2 = bst->idsmap[id2];
            mx2 = bst->smodlist[index2]->middle_x;
            //cout << "mx2= "<< mx2<< endl;
            my2 = bst->smodlist[index2]->middle_y;
            //cout << "my2= "<< my2<< endl;
        }
        
        unsigned int d = 0;
        d += std::abs(static_cast<int>(mx1-mx2))+std::abs(static_cast<int>(my1-my2));
        //cout << "d= "<<d<< endl;
        wl += d*net->weight;  //weighted wirelength
        //cout << "wl= "<<wl<< endl;
    }
    return wl;
}
void BST::tofloorplan(){
    if(this->placed) return;
    //std::list<std::pair<unsigned int, unsigned int>> contour;
    //contour.push_back(std::pair<unsigned int, unsigned int>(origx,origy));
    //contour.push_back(std::pair<unsigned int, unsigned int>(chip.width,origy));
    //preorder(root,origx,origy,contour);
    
    std::fill(hcontour.begin(), hcontour.end(), 0);
    std::fill(lcontour.begin(), lcontour.end(), 0);
    std::fill(tcontour.begin(), tcontour.end(), this->topbound);
    //cout<<"newcwidth="<<newcwidth<<endl;
  
    for(unsigned int i=0;i<fmods.size();i++){
        //cout<<"size= "<<fmods.size()<< endl;
        FixedModule* cur = fmods[i];
        //cout<<"x1="<<cur->x1<<"y1="<<cur->y1<<"width="<<cur->width<<"height="<<cur->height<<endl;
        if(cur->y1==0){
            unsigned int startx = cur->x1;
            unsigned int endx = startx + cur->width;
            for(unsigned int i=startx; i<endx;i++){
                hcontour[i] = cur->height;    
                //cout<<"h[i]="<<hcontour[i]<<endl;
            }
        }
        if(cur->x1==0){
            unsigned int starty = cur->y1;
            unsigned int endy = starty + cur->height;
            for(unsigned int i=starty; i<endy;i++){
                lcontour[i] = cur->width; 
                //cout<<"l[i]="<<lcontour[i]<<endl;   
            }
        }    
        /*if(cur->x1+cur->width==chip->width){
            unsigned int starty = cur->y1;
            unsigned int endy = starty + cur->height;
            for(unsigned int i=starty; i<endy;i++){
                rcontour[i] = cur->x1; 
                //cout<<"l[i]="<<lcontour[i]<<endl;   
            }
        } */
        if(cur->y1+cur->height==chipheight){
            unsigned int startx = cur->x1;
            unsigned int endx = startx + cur->width;
            for(unsigned int i=startx; i<endx;i++){
                tcontour[i] = cur->y1; 
                //cout<<"l[i]="<<lcontour[i]<<endl;   
            }
        }    
    }
    //cout<<"here"<<endl;
    //root
    unsigned int id = bstree[0].right;
    unsigned int index = this->idsmap[id];
    unsigned int maxx = 0;
    for(unsigned int i=0;i<this->smodlist[index]->height;i++){
        if (lcontour[i] > maxx){
            maxx = lcontour[i];
        }
         //cout<<"here1"<<endl;
    }
     //cout<<"here2"<<endl;
    this->smodlist[index]->x1 = maxx;
    unsigned int maxy = 0;
    for(unsigned int i=maxx;i<maxx+this->smodlist[index]->width;i++){
        if (hcontour[i] > maxy){
            maxy = hcontour[i];
        }
        hcontour[i]+=this->smodlist[index]->width;
    }
     //cout<<"here3"<<endl;
    this->smodlist[index]->y1 = maxy;
    //cout<<"x1"<<this->smodlist[index]->x1<<"y1="<<this->smodlist[index]->y1<<endl;
    //this->smodlist[index]->x1 = this->origx;
    //this->smodlist[index]->y1 = this->origy;
    this->smodlist[index]->getx2y2();
    
    //cout<<"index="<<index<<endl;*/
    for(unsigned int i=maxy;i<this->smodlist[index]->y2;i++){
        lcontour[i] += smodlist[index]->width;
    }
    this->currmaxx = this->smodlist[index]->x2;
    this->currmaxy = this->smodlist[index]->y2;
    this->fwidth = this->currmaxx - this->origx;
    this->fheight = this->currmaxy - this->origy;

    /*if(this->smodlist[index]->x1!=0 && this->smodlist[index]->y1!=0){ //root 不在0 0
        this->handle = 1;
    }*/
    //cout<<"origy4"<<endl;
    if(bstree[id].left>0){
        preorder(bstree[id].left,false);
    }
    //cout<<"origy5"<<endl;
    if(bstree[id].right>0){
        preorder(bstree[id].right,true);
    }
    //cout<<"origy6"<<endl;
    this->placed = 1 ;
    
}

void BST::preorder(unsigned int id, bool direction){//std::vector<unsigned int>&contour
    //Softmodule* curmd = this->smodlist[index];
    //Node& curnd = this->bstree[index];
    //unsigned int width = curnd.width;
    //unsigned int height = curnd.height;
    //cout<<"smodsize="<<smodlist.size()<<endl;
    /*if (index >= smodlist.size()) {
        // Add error handling for an invalid index
        cout << "Invalid index in smodlist: " << index << endl;
        return;
    }
    if (bstree[id].parent >= smodlist.size()) {
        // Add error handling for an invalid parent index
        cout << "Invalid parent index in bstree: " << bstree[index].parent << endl;
        return;
    }*/

    unsigned int parentid = bstree[id].parent;
    unsigned int index = this->idsmap[id];
    unsigned int indexp = this->idsmap[parentid];
    //cout<<"id="<<id<<", index="<<index<<", indexp="<<indexp<<endl;
    //if(index<=0) return;
   
    if(direction==0 ){ //left
        this->smodlist[index]->x1 = this->smodlist[indexp]->x1 + this->smodlist[indexp]->width;
        unsigned int startx = this->smodlist[index]->x1;
        unsigned int endx = startx + this->smodlist[index]->width;
        unsigned int maxy = 0;
        for(unsigned int i=startx; i<endx;i++){
            if(hcontour[i]>maxy){
                maxy = hcontour[i];
                //cout<<"i="<<i<<"maxy="<<maxy<<endl;
            }
        }
        this->smodlist[index]->y1 = maxy;
        this->smodlist[index]->getx2y2();
        //if(this->smodlist[index]->y2 > this->topbound){
            /*if(this->smodlist[index]->y2 - this->topbound < this->smodlist[index]->height/2){
                resizeNode(this->smodlist[index]);
                this->smodlist[index]->getx2y2();
                if(this->smodlist[index]->y2 > this->topbound){
                    moveNode(this->smodlist[index],this->smodlist[numSMod]);
                }
            }
            else{
                moveNode(this->smodlist[index],this->smodlist[numSMod]);
            }*/

       // }
        /*else{
            maxy += this->smodlist[index]->height;
            for(unsigned int i=startx; i<endx;i++){
                hcontour[i] = maxy;
            }
            if( this->smodlist[index]->x2 > this->currmaxx){
                this->currmaxx =  this->smodlist[index]->x2;
            }
            if( this->smodlist[index]->y2 > this->currmaxy){
                this->currmaxy =  this->smodlist[index]->y2;
            }
            this->fwidth = this->currmaxx - this->origx;
            this->fheight = this->currmaxy - this->origy;
        }*/
        /*maxy += this->smodlist[index]->height;
        for(unsigned int i=startx; i<endx;i++){
            hcontour[i] = maxy;
        }*/
        maxy += this->smodlist[index]->height;
        for(unsigned int i=startx; i<endx;i++){
            hcontour[i] = maxy;
        }
        /*for(unsigned int i=maxy;i<this->smodlist[index]->y2;i++){
            lcontour[i] += smodlist[index]->width;
        }*/
    } 
    else {  //right
        this->smodlist[index]->x1 = this->smodlist[indexp]->x1;
        unsigned int startx = this->smodlist[index]->x1;
        unsigned int endx = startx + this->smodlist[index]->width;
        unsigned int maxy = 0;
        for(unsigned int i=startx; i<endx;i++){
            if(hcontour[i]>maxy){
                maxy = hcontour[i];
                //cout<<"i="<<i<<"maxy="<<maxy<<endl;
            }
        }
        this->smodlist[index]->y1 = maxy;
        this->smodlist[index]->getx2y2();
        /*if(this->smodlist[index]->y2 > this->topbound){
            /*if(this->smodlist[index]->y2 - this->topbound < this->smodlist[index]->height/2){
                resizeNode(this->smodlist[index]);
                if(this->smodlist[index]->y2 - this->topbound < this->smodlist[index]->height/2){
                    moveNode(this->smodlist[index],this->smodlist[numSMod]);
                }
            }
            else{
                moveNode(this->smodlist[index],this->smodlist[numSMod]);
            }*/
            
            //resizeNode(this->smodlist[index]);
            
            //this->smodlist[index]->getx2y2();
            /*if(smodlist[index]->x1 >= chipwidth/2){
                unsigned int maxx = 0;
                for(unsigned int i=0;i<this->smodlist[index]->y1;i++){
                    if (lcontour[i] > maxx){
                        maxx = lcontour[i];
                    }
                    //cout<<"here1"<<endl;
                }
                this->smodlist[index]->x1 = maxx;
            }
            else{
                //this->smodlist[index]->x1 = this->smodlist[indexp]->x1 + 500;
            }*/
            //unsigned int startx = this->smodlist[index]->x1;
           // unsigned int endx = startx + this->smodlist[index]->width;
            /*unsigned int maxy = 0;
            //maxy <= topbound-height
            for(unsigned int i=startx; i<endx;i++){
                if(hcontour[i]>maxy){
                    maxy = hcontour[i];
                    //cout<<"i="<<i<<"maxy="<<maxy<<endl;
                }
            }
            this->smodlist[index]->y1 = maxy;
            this->smodlist[index]->getx2y2();
            startx = this->smodlist[index]->x1;
            endx = startx + this->smodlist[index]->width;*/
            maxy += this->smodlist[index]->height;
            for(unsigned int i=startx; i<endx;i++){
                hcontour[i] = maxy;
            }
            for(unsigned int i=maxy;i<this->smodlist[index]->y2;i++){
                lcontour[i] += smodlist[index]->width;
            }
        }
        // else if(direction==0 && handle == 1){ //left && handle != 1
        // this->smodlist[index]->x1 = this->smodlist[indexp]->x1 + this->smodlist[indexp]->width;
        // unsigned int startx = this->smodlist[index]->x1;
        // unsigned int endx = startx + this->smodlist[index]->width;
        // unsigned int maxy = 0;
        // for(unsigned int i=startx; i<endx;i++){
        //     if(hcontour[i]>maxy){
        //         maxy = hcontour[i];
        //         //cout<<"i="<<i<<"maxy="<<maxy<<endl;
        //     }
        // }
        // this->smodlist[index]->y1 = maxy;
        // this->smodlist[index]->getx2y2();

        // while(this->smodlist[index]->y1+this->smodlist[index]->height > this->topbound){   //放太高 往右往下走  //考慮右邊也要檢查
        //             unsigned int maxy = 0;
        //             unsigned int newx = this->smodlist[index]->x1+50;
        //             for(unsigned int i=newx;i< newx + this->smodlist[index]->width;i++){
        //                 if(hcontour[i]>maxy){
        //                     maxy = hcontour[i];
        //                     newx = i;
        //                     //cout<<"i="<<i<<"maxy="<<maxy<<endl;
        //                 }   
        //             }
        //             this->smodlist[index]->y1 = maxy;
        //             this->smodlist[index]->x1 = newx;
        //             this->smodlist[index]->getx2y2();
        // }

        // if(this->smodlist[index]->y2 > this->topbound){
        //     if(this->smodlist[index]->y2 - this->topbound < this->smodlist[index]->height/2){
        //         resizeNode(this->smodlist[index]);
        //         this->smodlist[index]->getx2y2();
        //         if(this->smodlist[index]->y2 > this->topbound){
        //             moveNode(this->smodlist[index],this->smodlist[numSMod]);
        //         }
        //     }
        //     else{
        //         moveNode(this->smodlist[index],this->smodlist[numSMod]);
        //     }

        // }
        /*else{
            maxy += this->smodlist[index]->height;
            for(unsigned int i=startx; i<endx;i++){
                hcontour[i] = maxy;
            }
            if( this->smodlist[index]->x2 > this->currmaxx){
                this->currmaxx =  this->smodlist[index]->x2;
            }
            if( this->smodlist[index]->y2 > this->currmaxy){
                this->currmaxy =  this->smodlist[index]->y2;
            }
            this->fwidth = this->currmaxx - this->origx;
            this->fheight = this->currmaxy - this->origy;
        }*/
        /*maxy += this->smodlist[index]->height;
        for(unsigned int i=startx; i<endx;i++){
            hcontour[i] = maxy;
        }*/
        //maxy += this->smodlist[index]->height;
    //     startx = this->smodlist[index]->x1;
    //     endx = startx + this->smodlist[index]->width;
    //     for(unsigned int i=startx; i<endx;i++){
    //         hcontour[i] += smodlist[index]->height;
    //     }
    //     for(unsigned int i=this->smodlist[index]->y1;i<this->smodlist[index]->y2;i++){
    //         lcontour[i] += smodlist[index]->width;
    //     }
    // } 
    //     else if(direction==1 && this->handle==1){ //原本要放parent上面
    //         this->smodlist[index]->x1 = this->smodlist[indexp]->x1;
    //         unsigned int startx = this->smodlist[index]->x1;
    //         unsigned int endx = startx + this->smodlist[index]->width;
    //         unsigned int maxy = 0;
    //         for(unsigned int i=startx; i<endx;i++){
    //             if(hcontour[i]>maxy){
    //             maxy = hcontour[i];
    //             //cout<<"i="<<i<<"maxy="<<maxy<<endl;
    //             }
    //         }
    //         this->smodlist[index]->y1 = maxy;

    //         unsigned int maxx = 0;
    //         for(unsigned int i=maxy;i<this->smodlist[index]->height;i++){
    //             if (lcontour[i] > maxx){
    //                 maxx = lcontour[i];
    //             }
    //         }
    //         this->smodlist[index]->x1 = maxx;
    //         this->smodlist[index]->getx2y2();

    //         // maxy += this->smodlist[index]->height;
    //         // for(unsigned int i=startx; i<endx;i++){
    //         //     hcontour[i] = maxy;
    //         // }
    //         // for(unsigned int i=maxy;i<this->smodlist[index]->y2;i++){
    //         //     lcontour[i] += smodlist[index]->width;
    //         // }
            
    //         while(this->smodlist[index]->y1+this->smodlist[index]->height > this->topbound){   //放太高 往右往下走
    //                 unsigned int maxy = 0;
    //                 unsigned int newx = this->smodlist[index]->x1+50;
    //                 for(unsigned int i=newx;i< newx + this->smodlist[index]->width;i++){
    //                     if(hcontour[i]>maxy){
    //                         maxy = hcontour[i];
    //                         newx = i;
    //                         //cout<<"i="<<i<<"maxy="<<maxy<<endl;
    //                     }   
    //                 }
    //                 this->smodlist[index]->y1 = maxy;
    //                 this->smodlist[index]->x1 = newx;
    //                 this->smodlist[index]->getx2y2();
    //         }

    //         startx = this->smodlist[index]->x1;
    //         endx = startx + this->smodlist[index]->width;
    //         for(unsigned int i=startx; i<endx;i++){
    //             hcontour[i] += smodlist[index]->height;
    //         }
    //         for(unsigned int i=this->smodlist[index]->y1;i<this->smodlist[index]->y2;i++){
    //             lcontour[i] += smodlist[index]->width;
    //         }
    //     }

            // else{
            //     //this->smodlist[index]->x1 = this->smodlist[indexp]->x1 + 500;
            // }
            //unsigned int startx = this->smodlist[index]->x1;
           // unsigned int endx = startx + this->smodlist[index]->width;
            /*unsigned int maxy = 0;
            //maxy <= topbound-height
            for(unsigned int i=startx; i<endx;i++){
                if(hcontour[i]>maxy){
                    maxy = hcontour[i];
                    //cout<<"i="<<i<<"maxy="<<maxy<<endl;
                }
            }
            this->smodlist[index]->y1 = maxy;
            this->smodlist[index]->getx2y2();
            startx = this->smodlist[index]->x1;
            endx = startx + this->smodlist[index]->width;*/
            /*maxy += this->smodlist[index]->height;
            for(unsigned int i=startx; i<endx;i++){
                hcontour[i] = maxy;
            }
            for(unsigned int i=maxy;i<this->smodlist[index]->y2;i++){
                lcontour[i] += smodlist[index]->width;
            }
        }*/
        

       
        /*maxy += this->smodlist[index]->height;
        for(unsigned int i=startx; i<endx;i++){
            hcontour[i] = maxy;
        }
        for(unsigned int i=maxy;i<this->smodlist[index]->y2;i++){
            lcontour[i] += smodlist[index]->width;
        }*/
       /* else{
            maxy += this->smodlist[index]->height;
            for(unsigned int i=startx; i<endx;i++){
                hcontour[i] = maxy;
            }
            if( this->smodlist[index]->x2 > this->currmaxx){
                this->currmaxx =  this->smodlist[index]->x2;
            }
            if( this->smodlist[index]->y2 > this->currmaxy){
                this->currmaxy =  this->smodlist[index]->y2;
            }
            this->fwidth = this->currmaxx - this->origx;
            this->fheight = this->currmaxy - this->origy;
        }*/
        
        /*unsigned int starty = maxy;
        unsigned int endy = starty + this->smodlist[index]->height;
        unsigned int maxx = x1;
        for(unsigned int i=starty; i<endy;i++){
            if(lcontour[i]>maxx){
                maxx = lcontour[i];
                lcontour[i]+=this->smodlist[index]->width;
                //cout<<"i="<<i<<"maxy="<<maxy<<endl;
            }
        }
        this->smodlist[index]->x1 = maxx;*/
    //}
    //cout<<"x1="<<this->smodlist[index]->x1<<endl;

    /*unsigned int startx = this->smodlist[index]->x1;
    unsigned int endx = startx + this->smodlist[index]->width;
    //cout<<"id="<<this->smodlist[index]->id<<"width="<<this->smodlist[index]->width<<endl;
    //cout<<"startx="<<startx<<"endx="<<endx<<endl;
    unsigned int maxy = 0;
    for(unsigned int i=startx; i<endx;i++){
        if(hcontour[i]>maxy){
            maxy = hcontour[i];
            //cout<<"i="<<i<<"maxy="<<maxy<<endl;
        }
    }*/
    
    //this->smodlist[index]->y1 = maxy;
    //cout<<"y1="<<this->smodlist[index]->y1<<endl;
    //this->smodlist[index]->getx2y2();
    //cout<<"x2="<<this->smodlist[index]->x2<<endl;
    //cout<<"y2="<<this->smodlist[index]->y2<<endl;

    /*if( this->smodlist[index]->x2 > this->currmaxx){
        this->currmaxx =  this->smodlist[index]->x2;
    }
    if( this->smodlist[index]->y2 > this->currmaxy){
        this->currmaxy =  this->smodlist[index]->y2;
    }
    this->fwidth = this->currmaxx - this->origx;
    this->fheight = this->currmaxy - this->origy;*/
    //cout<<"in tofloorplan currmaxx="<<this->currmaxx<<"currmaxy="<<this->currmaxy<<endl;
    //cout<<"fwidth="<<this->fwidth<<"fheight="<<this->fheight<<endl;

    /*maxy += this->smodlist[index]->height;
    for(unsigned int i=startx; i<endx;i++){
        hcontour[i] = maxy;
    }*/
    
    if( this->smodlist[index]->x2 > this->currmaxx){
        this->currmaxx =  this->smodlist[index]->x2;
    }
    if( this->smodlist[index]->y2 > this->currmaxy){
        this->currmaxy =  this->smodlist[index]->y2;
    }
    this->fwidth = this->currmaxx - this->origx;
    this->fheight = this->currmaxy - this->origy;
    
    if(bstree[id].left>0){
        //cout<<"id="<<id<<"left="<<bstree[id].left<<endl;
        preorder(bstree[id].left,false);
        
    }
    
    if(bstree[id].right>0){
        //cout<<"id="<<id<<"right="<<bstree[id].right<<endl;
        preorder(bstree[id].right,true);
    }
    //cout<<"end preorder"<<endl;

}
unsigned int BST:: verify()
{
    for(unsigned int i=1;i<=this->numSMod;i++){
        SoftModule* smod = this->smodlist[i];
        unsigned int x1 = smod->x1;
        unsigned int y1 = smod->y1;
        unsigned int xr1 = smod->x1+smod->width;
        unsigned int yr1 = smod->y1+smod->height;
        for(unsigned int j=1;j<=this->numSMod;j++){
        if(i==j) continue;
            SoftModule* smod = this->smodlist[j];
            unsigned int x2= smod->x1;
            unsigned int y2 = smod->y1;
            unsigned int xr2 = smod->x1+smod->width;
            unsigned int yr2 = smod->y1+smod->height;

            if(!(xr1<=x2 || x1>=xr2 ||yr1 <= y2 ||y1>=yr2)){
                return 1;
            }
        }
        for (auto fmod : fmods) {
            unsigned int fx2 = fmod->x1;
            unsigned int fy2 = fmod->y1;
            unsigned int fxr2 = fmod->x1 + fmod->width;
            unsigned int fyr2 = fmod->y1 + fmod->height;

            if (!(xr1 <= fx2 || x1 >= fxr2 || yr1 <= fy2 || y1 >= fyr2)) {
                return 1;  // Overlapping SoftModule with FixedModule
            }
        }

    }
    return 0;
}
double getprob(double deltaC, double t) {
        double power = -deltaC / t;
        double answer = exp(power);
        return ( answer < 1 ) ? answer : 1;
    }   
    bool acceptbad(double prob) {
        const int g = 100000000;
        int r = rand() % g;
        const int bound = g * prob;

        if ( r < bound ) return true;
        else return false;
    }   
void simulated_annealing(BST* bst, Chip* chip, std::vector<SoftModule*>smodvec,std::vector<FixedModule*>fmodvec,std::vector<Net*> nets,char* outputfilename){
    
    //cost=aA+bL+(1-a-b)(RS-R)^2
    //double a_base = 0.5;
    //unsigned int nf,n;
    //double a=a_base+(1-nf)/n;
    double a = 0.3;
    double b=0.4;
    //unsigned int A = compute_area(fmods);
    //unsigned int L=0;
    double RS = chip->height/(double)chip->width;
    //double RS = bst->newcheight/(double)bst->newcwidth;
    //double R=0;
    double initP=0.99;
    unsigned int perturb_num = 7;
    
    unsigned int init_area = bst->newcwidth * bst->newcheight;
   
    unsigned int init_wl = getwirelength(bst);
    //cout<<"hi"<<endl;
    double init_R = bst->newcheight/(double)bst->newcwidth;
    //cout<<"hi1"<<endl;
    double init_ratio = (RS>init_R)? RS-init_R:init_R-RS;
    //cout<<"hi2"<<endl;
    double init_c = a*init_area+b*init_wl+(1-a-b)*init_ratio;
   //cout<<"hi3"<<endl;
    double T1 = 10000;
    double T = T1;
    int c = 100;
    int k = 7;
    int r = 2;
    unsigned int stage2rd = 2400;
    unsigned int stage3rd = 4400;

    double delta_upavg = 0;
    double delta_cavg = 0;
    double area_norm = init_area;
    double wl_norm = init_wl;
    double ratio_norm = 0;

    unsigned int uphillcnt = 0;
    double costn = 1;
    double newcostn = 0;
    double mincostn = 1;
    double delta_cn = 0;
    double cost = init_c;
    double newcost = 0;
    double mincost = 0;
    double delta_c = 0;
    double Tmin=10;
    unsigned int area = init_area;
    unsigned int wl = init_wl;
    double R = init_R;
    double ratio = init_ratio;
    double prob;
    bool accept;
    bool best = 0;
    double overlap;
    unsigned int min_wl = init_wl*100;
    cout<<"min_wl="<<min_wl<<endl;
    /*unsigned int init_rounds = k * numSMod;
    for(unsigned int i=0; i<init_rounds; i++){
        unsigned int success = 0;
        while(success==0){
            //unsigned int init_dice_roll = rolldice();
            unsigned int M = rand()%4;
            switch(M){
            case 0:  //rotate
                unsigned int t = rand() % numSMod;
                success = bst.rotateNode(t);
                break;
            case 1:  //move
                unsigned int t1 = rand() % numSMod;
                do{
                    unsigned int t2 = rand() % numSMod;
                }while(t2==t1);
                success = bst.move(t1,t2);
                break;
            case 2:  //swap
                unsigned int t1 = rand() % numSMod;
                do{
                    unsigned int t2 = rand() % numSMod;
                }while(t2==t1);
                success = bst.swap(t1,t2);
                break;
            case 3:  //resize
                unsigned int t = rand() % numSMod;
                success = bst.resize(t);
                break;
            }
        }
        bst.floorplan();
        init_wirelength = getwirelength(nets);*/

        //sa
        //cout << "here "<< endl;
        unsigned int t1, t2,index1,index2;
        for(unsigned int n=0;n<=stage2rd+stage3rd;n++){
            double delta_cavg=0;
           //cout << "here1 "<< endl;
            for(unsigned int i=0;i<k*bst->numSMod;i++){
                //srand(static_cast<unsigned int>(time(nullptr)));
                unsigned int M = rand()%6;  
                //cout << "M= "<< M<<endl;
                switch(M){  //要改回M
                case 0:
                 //rotate
                    //cout << "here2 "<< endl;
                    t1 = 1+rand() % (bst->numSMod-1); //不能選到0 
                    index1=bst->idsmap[t1];
                    bst->rotateNode(bst->smodlist[index1]);
                    //cout << "here3 "<< endl;
                    break;
                case 1:     
                case 2:
                case 3:
                
                //move
                    t1 = 1+rand() % (bst->numSMod-1);
                    index1=bst->idsmap[t1];
                    do{
                    t2 = 1+rand() % (bst->numSMod-1);
                    }while(t2==t1);
                    index2=bst->idsmap[t2];
                    //cout << "here4 "<< endl;
                    bst->moveNode(bst->smodlist[index1],bst->smodlist[index2]);
                    //cout << "here5 "<< endl;
                    break;
                //swap
                
                case 4: 
                case 5:
                    t1 = 1+rand() % (bst->numSMod-1);
                    index1=bst->idsmap[t1];
                    do{
                    t2 = 1+rand() % (bst->numSMod-1);
                    }while(t2==t1);
                    index2=bst->idsmap[t2];
                    //cout << "here6 "<< endl;
                    bst->swapNode(bst->smodlist[index1],bst->smodlist[index2]);
                    //cout << "here7 "<< endl;
                    break;
                case 8:  //resize
                    t1 = 1+rand() % (bst->numSMod-1);  //t1:id
                    //cout << "here08 "<< endl;
                    index1 = bst->idsmap[t1];
                    //cout << "here008"<< endl;
                    bst->resizeNode(bst->smodlist[index1]);
                    //cout << "here09 "<< endl;
                    break;
                }
                //cout << "here8 "<< endl;
                bst->tofloorplan();
                //cout << "here9 "<< endl;
                area = bst->fwidth*bst->fheight;
                //cout << "area= "<< area<< endl;
                wl = getwirelength(bst);
                //cout << "wl "<< wl<<endl;
                R = bst->fheight/(double)bst->fwidth;
                ratio = (RS>R)? RS-R:R-RS;
                overlap = bst->verify();
                costn = newcostn;
                newcostn = a*area/area_norm+b*wl/wl_norm+(1-a-b)*overlap;
                //cout << "newcostn "<< newcostn<<endl;
                delta_cn = newcostn-costn;
                delta_cavg += delta_cn;
               
                cost = newcost;
                //newcost = a*area+b*wl+(1-a-b)*ratio;
                newcost = a*area+b*wl+(1-a-b)*overlap*100;
                if(mincost==0) mincost=newcost;
                delta_c = newcost-cost;

                prob = getprob(delta_c, T);
                accept = acceptbad(prob);
                if(overlap==0 && bst->currmaxx<=chip->width && bst->currmaxy<= chip->height){
                        //mincost = newcost;
                        cout<<"updated"<<endl;
                        cout<< "bst->currmaxx="<<bst->currmaxx<<"chip->width"<<chip->width<<endl;
                        cout<< "bst->currmaxy="<<bst->currmaxy<<"chip->height"<<chip->height<<endl;
                        if(wl<=min_wl){
                        min_wl = wl;
                        cout<<"updated1"<<endl;
                        bst->savebest();
                        }
                        best = 1;
                    }
                if(newcost<cost){
                    bst->save();
                    // if(overlap==0/*newcost<=mincost*/ && bst->currmaxx<=chip->width && bst->currmaxy<= chip->height){
                    //     //mincost = newcost;
                    //     cout<<"updated"<<endl;
                    //     cout<< "bst->currmaxx="<<bst->currmaxx<<"chip->width"<<chip->width<<endl;
                    //     cout<< "bst->currmaxy="<<bst->currmaxy<<"chip->height"<<chip->height<<endl;
                    //     if(wl<=min_wl){
                    //     cout<<"updated1"<<endl;
                    //     bst->savebest();
                    //     }
                    //     best = 1;
                    // }
                }
                else{
                    /*double random = ((double)rand())/RAND_MAX;
                    if(random < exp(-delta_c*10/T)){
                        uphillcnt++;
                        bst->save();
                    }*/
                    if(accept==true){
                        bst->save();
                    }
                    else{
                        bst->loadprev();
                        //raw_cost = raw_last_cost;
                        newcostn = costn;
                    }
                }
                
            }

            delta_cavg /= bst->numSMod*k;

            if(n<=stage2rd){
                T = T1*delta_cavg/((n+1)*c);
            }
            else{
                T = T1*delta_cavg/(n+1);
                //a = 0.2;
            }
        }
        //area = bst->fwidth*bst->fheight;
        if(best==1){
            bst->loadbest();
        }else{
            bst->loadprev();
        }
        
        //bst->loadbest();
        bst->tofloorplan();
        unsigned int illegal = bst->verify();
        bst->wirelength = getwirelength(bst);
        std::cout<<"bestwl="<<bst->wirelength<<endl;
        cout<< "bst->currmaxx="<<bst->currmaxx<<"chip->width"<<chip->width<<endl;
        cout<< "bst->currmaxy="<<bst->currmaxy<<"chip->height"<<chip->height<<endl;
        //ratio = bst->fheight/(double)bst->fwidth;

        //if(bst->currmaxx<=chip->width && bst->currmaxy<= chip->height){
        if(illegal==0){
            std::cout<<"legal flooplan found!"<<endl;
            writeoutput(outputfilename,bst);
            //printfloorplan();
            //break;
        }else{
            std::cout<<"illegal flooplan!"<<endl;
            writeoutput(outputfilename,bst);
        }

    }
    
    /*delta_c /= init_rounds;
    delta_upavg /= uphillcnt;
    


    BST::BST bst1;
    BST::BST best;
    bst1.root = BST::BST(smods[0],smods.size());
    float T=2000;
    //adaptive fast-SA
    //stage 1, r=1


    while(T<Tmin && r<k){
        BST* newBST = perturb(bst1);
        place(newBST);
        update_contour();
        newcost = compute_cost();
        modify_weight();
        update_T();
        delta_c = abs(newcost-cost);
        T = T*delta_c/n*c;
        cost = newcost;
        ++r;
    }
    while(T<Tmin){
        perturb();
        update_contour();
        compute_cost();
        modify_weight();
        update_T();
        T = T*delta_c/n;
    }
    return best;
}*/

void writeoutput(char* outputfilename,BST* bst){
    ofstream outputFile(outputfilename);
    outputFile<< "Wirelength " << bst->wirelength<<endl;
    outputFile<<endl;
    outputFile<< "NumSoftModules " << bst->numSMod<<endl;
    for(unsigned int i=1;i<=bst->numSMod;i++){
        SoftModule* smod = bst->smodlist[i];
        outputFile<< smod->name <<" "<<smod->x1<<" "<<smod->y1<<" "<<
            smod->width<<" "<<smod->height<<endl;
    }
    /*for(auto&smod:bst.smodlist){
        outputFile<< smod->name <<" "<<smod->x1<<" "<<smod->y1<<" "<<
            smod->width<<" "<<smod->height<<endl;
    }*/
    outputFile.close();
}

/*float compute_cost(a,b,RS,BST* bst){
    float cost = 0;
    //cost function=aA+BL+(1-a-b)(RS-R)^2
    unsigned int A = compute_area();
    unsigned int L = compute_HPWL();
    unsigned int R = compute_ratio();
    cost = a*A+b*L+(1-a-b)(RS-R)^2;
    return cost;
}

int compute_area(vector<string>vec){
    packing();
    int area = max_x*max_y;
    return area;
}
//每個NET的wirelength是連接的兩個module的中心點的x座標差+y座標差
//程式最後要輸出weight*HPWL
int compute_HPWL(vector<string>vec){
    int HPWL = 0;
    for (auto&net:vec){
        HPWL = abs(net.m1.middle_x-net.m2.middle_x)+
                abs(net.m1.middle_y-net.m2.middle_y);
    }
    
    return HPWL;
}

perturb(vector<string>vec){
    move_block();
    swap_block();
    resize();
}

resize(){
    srand(static_cast<unsigned int>(time(nullptr)));
    int option = rand() % 4;
}

place(BST& bst){
    Node* root = bst.root;
    sort(fmods.begin(), fmods.end(), cmpx);
    FixedModule fm1 = fmods.front();
    fx = fm1.x;
    sort(fmods.begin(), fmods.end(), cmpy);
    FixedModule fm2 = fmods.front();
    fy = fm2.y;
    if(fx!=0 && fy!=0 && root.area<=fx*fy){
        root->x=0;
        root->y=0;
    }
    else{
        root->x=fx+fm1.width;
        root->y=0;
    }
    Node* node;
    postorder(root){
        if(node!=nullptr){
            postorder(root->left);  //左子放右邊 x=x1+w  y=y1
            postorder(root->right); //右子放上面 x=x1  y=y1+h  
            //check x and y coordinate of fixedmodules
            if(lchild==1){  //用bit紀錄是左子還是右子
                node.x = parent.x+parent.width;
            }
            for(auto&fmod:fmods){
                if(node.x<=fmod.x && node.x+node.width>=fmod.x){
                    node.x = fmod.x+fmod.width;
                    if(node.x+node.width>=chip.width){

                    }
                }
            }
        }
    }
}*/
