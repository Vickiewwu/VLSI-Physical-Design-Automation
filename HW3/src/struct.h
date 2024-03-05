#ifndef STRUCT_H
#define STRUCT_H
#include <string>
#include <vector>
#include <map>
#include <cmath>
#include <iostream>

struct Chip {
    unsigned int width;
    unsigned int height;
};

struct SoftModule {
    std::string name;
    //左下座標
    unsigned int x1;
    unsigned int y1;
    //右上座標
    unsigned int x2;
    unsigned int y2;
    unsigned int width;
    unsigned int height;
    unsigned int min_area;
    unsigned int area;
    unsigned int id;
    //module 中心點x.y座標 要取floor
    unsigned int middle_x;
    unsigned int middle_y;
    SoftModule(): x1(0),y1(0),x2(0),y2(0),min_area(0),area(0),id(0),middle_x(0),middle_y(0) {}
    SoftModule(std::string s, unsigned int a, unsigned int i): name(s), x1(0),y1(0),x2(0),y2(0),min_area(a),area(0),id(i),middle_x(0),middle_y(0) {}
    void getmiddle(){
        middle_x = std::floor((x1+x2)/2);
        middle_y = std::floor((y1+y2)/2);
    }
    void getx2y2(){
        x2 = x1+width;
        y2 = y1+height;
    }
    void getwh(){
        int d = 1;
        for(unsigned int i = 2; i<min_area; i++){
            if(i*i==min_area) width = i;
            if(i*i<min_area && min_area%i==0) d = i; 
            if(i*i>min_area) width = d;
        }
        height = min_area/width;
        if(height/width<=0.5 || height/width>2){
            std::cout<<"regetwh"<<std::endl;
            for(unsigned int i = 2; i<min_area; i++){
            if(i*i==min_area) width = i;
            if(i*i<min_area ) d = i; 
            if(i*i>min_area) width = d+1;
            height = min_area/width +1;
        }

        }
    }
   
};

struct FixedModule {
    std::string name;
    //左下座標
    unsigned int x1;
    unsigned int y1;
    //右上座標
    unsigned int x2;
    unsigned int y2;
    unsigned int width;
    unsigned int height;
    unsigned int area;
    unsigned int id;
    //module 中心點x.y座標
    unsigned int middle_x;
    unsigned int middle_y;
    FixedModule(): x1(0),y1(0),x2(0),y2(0),area(0),id(0),middle_x(0),middle_y(0) {}
    FixedModule(std::string s, unsigned int x, unsigned int y,unsigned int w, unsigned int h,unsigned int i): name(s), x1(x),y1(y),x2(x+w),y2(y+h),area(w*h),id(i),middle_x(0),middle_y(0) {}
    void getmiddle(){
        middle_x = floor((x1+x2)/2);
        middle_y = floor((y1+y2)/2);
    }
};

struct Net {
    unsigned int id;
    unsigned int weight;
    std::string m1;
    std::string m2;
};

struct Node{
    unsigned int parent;
    unsigned int left;
    unsigned int right;
    unsigned int width;
    unsigned int height;
    Node(): parent(0),left(0),right(0),width(0),height(0){}
    Node(unsigned int p,unsigned int l,unsigned int r,unsigned int w,unsigned int h) :
    parent(p),left(l),right(r),width(w),height(h){}
    void rotate(){
        unsigned int tmp;
        tmp = height;
        height = width;
        width = tmp;
    }
    void resize(unsigned int w,unsigned int h){
        width = w;
        height = h;
    }
};

class BST{
public:
    void init(std::vector<SoftModule*>smodvec,std::vector<FixedModule*>fmodvec,unsigned int chipw,unsigned chiph,unsigned int numsmod);
    void rotateNode(SoftModule*);
    void resizeNode(SoftModule*);
    void moveNode(SoftModule*,SoftModule*);
    void swapNode(SoftModule*,SoftModule*);
    void savebest();
    void loadbest();
    void save();
    void loadprev();
    void tofloorplan();
    void preorder(unsigned int,bool);
    std::map<unsigned int,unsigned int> idsmap; 
    std::map<unsigned int,unsigned int> idfmap;
    std::vector<Node> bstree;
    std::vector<Node> bestbst;
    std::vector<Node> savedbst;
    std::vector<SoftModule*> smodlist;
    std::vector<FixedModule*> fmodlist;

    unsigned int numSMod;
    unsigned int topbound;
    unsigned int origx;
    unsigned int origy;
    unsigned int newcwidth;
    unsigned int newcheight;
    unsigned int fwidth;
    unsigned int fheight;
    unsigned int currmaxx;
    unsigned int currmaxy;
    unsigned int wirelength;
    bool placed;
    bool handle;
    void insertNode(unsigned int,unsigned int,bool);
    void deleteNode(unsigned int);
    unsigned int verify();
};
#endif 
