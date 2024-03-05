#include "GlobalPlacer.h"
#include "ExampleFunction.h"
#include "NumericalOptimizer.h"
#include <cmath>


ExampleFunction::ExampleFunction(wrapper::Placement &placement)
: _placement(placement)
{
    num = _placement.numModules();

    double Corewidth = _placement.boundryRight() - _placement.boundryLeft();
    double Coreheight = _placement.boundryTop() - _placement.boundryBottom();
    binwidth = Corewidth / binnum;
    binheight = Coreheight / binnum;

    a = 0.01*Corewidth;
    beta = 0;
}


double ExampleFunction::getTb(const double width, const double height)
{
    core_area = width * height;

    double area = 0;
    for (unsigned int i = 0; i < num; i++)
        area += _placement.module(i).area();

    return area / core_area;
}

double ExampleFunction::getWL(const vector<double> &x)
{
    double wl = 0;

    vector<double> exp_x(num* 4);
    for (unsigned int i = 0; i < num; i++) {
        exp_x[i * 4] = exp(x[i * 2] / a);    
        exp_x[i * 4 + 1] = exp(-x[i * 2] / a);    
        exp_x[i * 4 + 2] = exp(x[i * 2 + 1] / a);    
        exp_x[i * 4 + 3] = exp(-x[i * 2 + 1] / a);    
    }

    for (unsigned int net_id = 0; net_id < _placement.numNets(); net_id++) {
        
        double xsum = 0;
        double xnsum = 0;
        double ysum = 0;
        double ynsum = 0;
        const unsigned int num_pins = _placement.net(net_id).numPins();
        for (unsigned int index = 0; index < num_pins; index++) {
            const unsigned int module_id = _placement.net(net_id).pin(index).moduleId();

            xsum += exp_x[module_id * 4];
            xnsum += exp_x[module_id * 4 + 1];
            ysum += exp_x[module_id * 4 + 2];
            ynsum += exp_x[module_id * 4 + 3];
        }

        wl += log(xsum) + log(xnsum) + log(ysum) + log(ynsum);
    }

    return a * wl;
}


double ExampleFunction::grad_WL(const vector<double> &x, vector<double> &g)
{
    double wl = 0;

    vector<double> exp_x(num * 4);
    for (unsigned int i = 0; i < num; i++) {
        exp_x[i * 4] = exp(x[i * 2] / a);
        exp_x[i * 4 + 1] = exp(-x[i * 2] / a);
        exp_x[i * 4 + 2] = exp(x[i * 2 + 1] / a);
        exp_x[i * 4 + 3] = exp(-x[i * 2 + 1] / a);
    }

    for (unsigned int j = 0; j < _placement.numNets(); j++) {
       
        double xsum = 0;
        double xnsum = 0;
        double ysum = 0;
        double ynsum = 0;
       
        for (unsigned int pid = 0; pid < _placement.net(j).numPins(); pid++) {
            unsigned int id = _placement.net(j).pin(pid).moduleId();

            xsum += exp_x[id * 4];
            xnsum += exp_x[id * 4 + 1];
            ysum += exp_x[id * 4 + 2];
            ynsum += exp_x[id * 4 + 3];
        }

        wl += log(xsum) + log(xnsum) + log(ysum) + log(ynsum);

        for (unsigned int pid = 0; pid < _placement.net(j).numPins(); pid++) {
            unsigned int id = _placement.net(j).pin(pid).moduleId();

            g[id * 2] += exp_x[id * 4] / xsum;
            g[id * 2] -= exp_x[id * 4 + 1] / xnsum;
            g[id * 2 + 1] += exp_x[id * 4 + 2] / ysum;
            g[id * 2 + 1] -= exp_x[id * 4 + 3] / ynsum;
        }
    }

    return a * wl;
}



double ExampleFunction::getBinDensity(const vector<double> &x)
{
    double D = 0;
    
    for (unsigned int i = 0; i < binnum; i++) {
        for (unsigned int j = 0; j < binnum; j++) {
            double bx = _placement.boundryLeft() + (0.5 + j) * binwidth;
            double by = _placement.boundryBottom() + (0.5 + i) * binheight;
            double d = 0;
            for (unsigned int k = 0; k < num; k++) {
                const double wk = _placement.module(k).width();
                const double jk = _placement.module(k).height();

                double xov = 0;
                double yov = 0;
                double dx = abs(x[k * 2] - bx);
                double dy = abs(x[k * 2 + 1] - by);

                if (dx <= binwidth / 2 + wk / 2) {
                    double ax = 4 / ((binwidth + wk) * (2 * binwidth + wk));
                    xov = 1 - ax * dx * dx;
                }
                else if (dx <= binwidth + wk / 2) {
                    double bx = 4 / (binwidth * (2 * binwidth + wk));
                    xov = bx * (dx - binwidth - wk / 2) * (dx - binwidth - wk / 2);
                }
                else {
                    xov = 0;
                }

                if (dy <= binheight / 2 + jk / 2) {
                    double ay = 4 / ((binheight + jk) * (2 * binheight + jk));
                    yov = 1 - ay * dy * dy;
                }
                else if (dy <= binheight + jk / 2) {
                    double by = 4 / (binheight * (2 * binheight + jk));
                    yov = by * (dy - binheight - jk / 2) * (dy - binheight - jk / 2);
                }
                else {
                    yov = 0;
                }

                d += xov * yov;
            }

            // if(d>Tb) d=d-Tb;
            D += d * d;  
        }
    }

    return beta * D;
}


double ExampleFunction::grad_BinDensity(const vector<double> &x, vector<double> &g)
{
    double D = 0;
    
    for (unsigned int i = 0; i < binnum; i++) {
        for (unsigned int j = 0; j < binnum; j++) {
            vector<double> g_bin(g.size(), 0);
            double bx = _placement.boundryLeft() + (0.5 + j) * binwidth;
            double by = _placement.boundryBottom() + (0.5 + i) * binheight;
            double d = 0;
            for (unsigned int k = 0; k < num; k++) {
                const double wk = _placement.module(k).width();
                const double hk = _placement.module(k).height();

                double xov = 0;
                double yov = 0;
                double g_xov = 0;
                double g_yov = 0;
                double dx = x[k * 2] - bx;
                double dy = x[k * 2 + 1] - by;
                double adx = abs(dx);
                double ady = abs(dy);
                double ax = 4 / ((binwidth + wk) * (2 * binwidth + wk));
                double bx = 4 / (binwidth * (2 * binwidth + wk));
                double ay = 4 / ((binheight + hk) * (2 * binheight + hk));
                double by = 4 / (binheight * (2 * binheight + hk));

                if (adx <= binwidth / 2 + wk / 2) {
                    xov = 1 - ax * adx * adx;
                    g_xov = (-2) * ax * dx;
                }
                else if (adx <= binwidth + wk / 2) {
                    xov = bx * (adx - binwidth - wk / 2) * (adx - binwidth - wk / 2);
                    if (dx > 0)
                        g_xov = 2 * bx * (dx - binwidth - wk / 2);
                    else
                        g_xov = 2 * bx * (-1)*(dx - binwidth - wk / 2);
                }
                else {
                    xov = 0;
                    g_xov = 0;
                }

                if (ady <= binheight / 2 + hk / 2) {
                    yov = 1 - ay * ady * ady;
                    g_yov = -2 * ay * dy;
                }
                else if (ady <= binheight + hk / 2) {
                    yov = by * (ady - binheight - hk / 2) * (ady - binheight - hk / 2);
                    if (dy > 0)
                        g_yov = 2 * by * (dy - binheight - hk / 2);
                    else
                        g_yov = 2 * by * (-1)*(dy - binheight - hk / 2);
                }
                else {
                    yov = 0;
                    g_yov = 0;
                }

                d += xov * yov;

                g_bin[k * 2] = g_xov * yov;
                g_bin[k * 2 + 1] = xov * g_yov;
            }

        //    if(d>Tb) d=d-Tb;
            D += d * d;

            for (unsigned int m = 0; m < num; m++) {
                g[m * 2] += beta * 2 * d * g_bin[m * 2];
                g[m * 2 + 1] += beta * 2 * d * g_bin[m * 2 + 1];
            }
        }
    }

    return beta * D;
}

void ExampleFunction::evaluateFG(const vector<double> &x, double &f, vector<double> &g)
{
    g = vector<double>(g.size(), 0);
    f = grad_WL(x, g);
    f += grad_BinDensity(x, g);
}

void ExampleFunction::evaluateF(const vector<double> &x, double &f)
{
    f = getWL(x) + getBinDensity(x);
}

unsigned ExampleFunction::dimension()
{
    return num * 2; 
}



// // double ExampleFunction::WL(){
// //     //printf("e1\n");  
// //     double wl=0.0;
// //     double a = 10;
// //     double xsum,xnsum,ysum,ynsum;
// //     //xsum=xnsum=ysum=ynsum=0.0;
// //     for(int i=0;i<_placement.numNets(); ++i){
// //         xsum=xnsum=ysum=ynsum=0.0;
// //         for(int j=0;j<_placement.net(i).numPins();++j){
// //             double x = _placement.net(i).pin(j).x();
// //             double y = _placement.net(i).pin(j).y();
// //             //printf("x=%lf,y=%lf\n", x,y);
// //             xsum += exp(x/a);
// //             xnsum += exp((x*(-1))/a);
// //             ysum += exp(y/a);
// //             ynsum += exp((y*(-1))/a);
// //             // unsigned  expX = static_cast<unsigned >(exp(x / a));
// //             // unsigned  expXn = static_cast<unsigned >(exp((-x) / a));
// //             // unsigned  expY = static_cast<unsigned >(exp(y / a));
// //             // unsigned  expYn = static_cast<unsigned >(exp((-y) / a));

// //             // xsum += expX;
// //             // xnsum += expXn;
// //             // ysum += expY;
// //             // ynsum += expYn;
// //             //printf("xsum=%lf,xnsum=%lf,ysum=%lf,ynsum=%lf\n", xsum,xnsum,ysum,ynsum);

// //         }
// //         //printf("lxsum=%llu,lxnsum=%llu,lysum=%llu,lynsum=%llu\n",log(xsum),log(xnsum),log(ysum),log(ynsum));
// //         wl += a*(log(xsum)+log(xnsum)+log(ysum)+log(ynsum));
// //         //printf("WL=%lf\n", wl);
// //     }
// //     //wl = a*(log(xsum)+log(xnsum)+log(ysum)+log(ynsum));
// //     printf("totalWL=%lf\n", wl);
// //     return wl;
// // }
// // void ExampleFunction::WL(vector<double> &g){
// //     //printf("e1\n");  
// //     f=0;
// //     double coreWidth = _placement.boundryRight() - _placement.boundryLeft();
// //     double a = 0.01*coreWidth;
// //     double xsum,xnsum,ysum,ynsum,cx,cy;
// //     //xsum=xnsum=ysum=ynsum=0.0;
// //     for(int i=0;i<_placement.numNets(); ++i){
// //         xsum=xnsum=ysum=ynsum=0.0;
// //         for(int j=0;j<_placement.net(i).numPins();++j){
// //             double x = _placement.module(_placement.net(i).pin(j).moduleId()).centerX();
// //             double y = _placement.module(_placement.net(i).pin(j).moduleId()).centerY();
// //             //printf("x=%lf,y=%lf\n", x,y);
            
// //             xsum += exp(x/a);
// //             xnsum += exp((x*(-1))/a);
// //             ysum += exp(y/a);
// //             ynsum += exp((y*(-1))/a);
           
// //             //printf("xsum=%lf,xnsum=%lf,ysum=%lf,ynsum=%lf\n", xsum,xnsum,ysum,ynsum);
// //         }
// //         //printf("lxsum=%llu,lxnsum=%llu,lysum=%llu,lynsum=%llu\n",log(xsum),log(xnsum),log(ysum),log(ynsum));
// //         f += log(xsum)+log(xnsum)+log(ysum)+log(ynsum);
// //         //printf("WL=%lf\n", wl);
// //         for(int j=0;j<_placement.net(i).numPins();++j){
// //             int moduleid = _placement.net(i).pin(j).moduleId();
// //             cx = _placement.module( moduleid).centerX();
// //             cy = _placement.module( moduleid).centerY();
// //             //printf("x=%lf,y=%lf\n", x,y);
// //             if(!_placement.module(moduleid).isFixed()){
// //                 g[moduleid*2] += exp(cx/a)/(xsum*a);
// //                 g[moduleid*2] -= exp(cx/a)/(xnsum*a);
// //                 g[moduleid*2+1] += exp(cx/a)/(ysum*a);
// //                 g[moduleid*2+1] -= exp(cx/a)/(ynsum*a);
// //             }
// //         }
// //     }
// //     //wl = a*(log(xsum)+log(xnsum)+log(ysum)+log(ynsum));
// //     //printf("totalWL=%lf\n", wl);
  
// // }
// // void ExampleFunction::grad_WL( vector<double> &g){
// //     //printf("w2\n");
// //     double wl=0.0;
// //     double coreWidth = _placement.boundryRight() - _placement.boundryLeft();
// //     double a = 0.01*coreWidth;
// //     double xsum,xnsum,xdir,xndir,ysum,ynsum,ydir,yndir,gx,gnx,gy,gny;
// //     xsum=xnsum=xdir=xndir=ysum=ynsum=ydir=yndir=gx=gnx=gy=gny=0.0;
// //     int netid=0;
// //     double x,y,xk,yk;
// //     for(int k=0;k<_placement.numModules(); ++k){
// //         //printf("w8\n");
// //         if (_placement.module(k).isFixed()){
// //             g[k*2] = 0;
// //             g[k*2+1] = 0;
// //             continue;
// //         }
// //         xk = _placement.module(k).centerX();
// //         yk = _placement.module(k).centerY();
       
                
// //         //printf("w3\n");
// //         for(int i=0;i<_placement.module(k).numPins(); ++i){
// //             //printf("w4\n");   //numnets
// //             netid = _placement.module(k).pin(i).netId();
// //             //printf("w5\n");
// //             xsum=xnsum=xdir=xndir=ysum=ynsum=ydir=yndir=gx=gnx=gy=gny=0.0;
// //             for(int j=0;j<_placement.net(netid).numPins();++j){
// //                 //printf("w6\n");
// //                 x = _placement.module(_placement.net(netid).pin(j).moduleId()).centerX();
// //                 y = _placement.module(_placement.net(netid).pin(j).moduleId()).centerY();
// //                 //printf("w7\n");
// //                 xsum += exp(x/a);
// //                 xnsum += exp(-x/a);
// //                 ysum += exp(y/a);
// //                 ynsum += exp(-y/a);
// //                 //printf("xsum=%lf,xnsum=%lf,ysum=%lf,ynsum=%lf\n", xsum,xnsum,ysum,ynsum);
// //             }
// //             xdir = exp(xk/a)/xsum;
// //             xndir = exp(-xk/a)/xnsum;
// //             ydir = exp(yk/a)/ysum;
// //             yndir =  exp(-yk/a)/ynsum;
// //             gx += a*(xdir-xndir);
// //             gy += a*(ydir-yndir);
// //             //printf("xdir=%lf,xndir=%lf,ydir=%lf,yndir=%lf\n", xdir,xndir,ydir,yndir);
// //             //printf("gx=%lf,gy=%lf\n", gx,gy);
// //         }
// //         g[k*2] += gx*0.4;
// //         g[k*2+1] += gy*0.4;
// //     }
// //     grad_BinDensity(g);
// // }
// // void ExampleFunction::grad_WL( vector<double> &g){
// //        //printf("e2\n");
// //     double wl=0.0;
// //     double a = 10;
// //     double xsum,xnsum,xdir,xndir,ysum,ynsum,ydir,yndir,gx,gnx,gy,gny;
// //     xsum=xnsum=xdir=xndir=ysum=ynsum=ydir=yndir=gx=gnx=gy=gny=0.0;
// //     int netid=0;
// //     double x,y,xk,yk;
// //     for(int k=0;k<_placement.numModules(); ++k){
// //         xk = _placement.module(k).centerX();
// //         yk = _placement.module(k).centerY();
// //         for(int i=0;i<_placement.numPins(); ++i){   //numnets
// //             netid = _placement.module(k).pin(i).netId();
// //             xsum=xnsum=xdir=xndir=ysum=ynsum=ydir=yndir=gx=gnx=gy=gny=0.0;
// //             for(int j=0;j<_placement.net(netid).numPins();++j){
// //                 x = _placement.module(_placement.net(i).pin(j).moduleId()).centerX();
// //                 y = _placement.module(_placement.net(i).pin(j).moduleId()).centerY();
// //                 xsum += exp(x/a);
// //                 xnsum += exp(-x/a);
// //                 ysum += exp(y/a);
// //                 ynsum += exp(-y/a);
// //             }
// //             xdir = exp(xk/a)/xsum;
// //             xndir = exp(-xk/a)/xnsum;
// //             ydir = exp(yk/a)/ysum;
// //             yndir =  exp(-yk/a)/ynsum;
// //             gx += xdir-xndir;
// //             gy += ydir-yndir;
// //         }
// //         g[k*2] += gx;
// //         g[k*2+1] += gy;
// //     }
// //     for(int i=0;i<_placement.numNets(); ++i){
// //         for(int j=0;j<_placement.net(i).numPins();++j){
// //             xdir = exp(_placement.net(i).pin(j).x()/a)/xsum;
// //             xndir = exp((-1)*(_placement.net(i).pin(j).x()/a))/xnsum;
// //             ydir = exp(_placement.net(i).pin(j).y()/a)/ysum;
// //             yndir = exp((-1)*(_placement.net(i).pin(j).y()/a))/ynsum;
// //             g[2*_placement.net(i).pin(j).moduleId()] = a*(xdir-xndir);
// //             g[1+2*_placement.net(i).pin(j).moduleId()] = a*(ydir-yndir);
// //         }
// //     }
// //     for(int i=0;i<_placement.numModules(); ++i){
// //         if (!_placement.module(i).isFixed()){
// //             g[i*2] /= _placement.module(i).numPins();
// //             g[i*2+1] /= _placement.module(i).numPins();
// //         }else{
// //             g[i*2] = 0;
// //             g[i*2+1] = 0;
// //         }
// //     }
// // }
// double ExampleFunction::getTb(){
//     double coreWidth = _placement.boundryRight() - _placement.boundryLeft();
//     double coreHeight = _placement.boundryTop() - _placement.boundryBottom();
//     unsigned binnum = ceil(sqrt(_placement.numModules()));
//     double coreArea = coreWidth*coreHeight;
//     double area= 0.0;
    
//     for(int k=0;k<_placement.numModules(); ++k){
//         area += _placement.module(k).area();
//     }
//     return area/coreArea+0.1;
// }
// double ExampleFunction::BinDensity(){
//     //printf("e3\n");
//     double xsum,xnsum,ysum,ynsum,d,D,c,ax,bx,ay,by;
//     double xb,yb,wb,hb, xj,yj,dx,dy,xov,yov,wj,hj;
//     xsum=xnsum=ysum=ynsum=d=D=c=ax=bx=ay=by=0.0;
//     xb=yb=wb=hb= xj=yj=dx=dy=xov=yov=wj=hj=0.0;
//     double coreWidth = _placement.boundryRight() - _placement.boundryLeft();
//     double coreHeight = _placement.boundryTop() - _placement.boundryBottom();
//     unsigned binnum = ceil(sqrt(_placement.numModules())/4);
//     //printf("binnum =%lu\n",binnum);
//     double binwidth = coreWidth/binnum;  
//     //printf("binwidth=%lf\n",binwidth);
//     double binheight = coreHeight/binnum;
//     //printf("binheight=%lf\n",binheight);
//     double tb = getTb();
//     double * bindensity;
//     bindensity = new double[binnum*binnum];
//     //printf("tb=%lf\n", tb);
//     //c = 2;
//     for(int i = 0; i < binnum*binnum; ++i){
// 		    bindensity[i] = 0;   
// 	}
//     for(int i=0;i<binnum; ++i){
//         for(int j=0;j<binnum; ++j){
//             wb = binwidth;
//             hb = binheight;
//             xb = _placement.boundryLeft() + i*wb + 0.5*wb;
//             yb = _placement.boundryBottom() + j*hb + 0.5*hb;
//             //printf("wb=%lf,hb=%lf\n", wb,hb);
//             //printf("xb=%lf,yb=%lf\n", xb,yb);
//             d=0;
//             for(int k=0;k<_placement.numModules(); ++k){
//                 //printf("id=%d\n", k);
//                 //printf("wb=%lf,hb=%lf\n", wb,hb);
//                 //printf("xb=%lf,yb=%lf\n", xb,yb);
//                 xj = _placement.module(k).centerX();
//                 yj = _placement.module(k).centerY();
//                 wj = _placement.module(k).width();
//                 hj = _placement.module(k).height();
//                 c = _placement.module(k).area()/(binwidth*binheight);
//                 //printf("xj=%lf,yj=%lf\n", xj,yj);
//                 //printf("wj=%lf,hj=%lf\n", wj,hj);
//                 dx = abs(xb-xj);
//                 dy = abs(yb-yj);
//                 //printf("dx=%lf,dy=%lf\n", dx,dy);
//                 ax = 4/((4*wb+wj)*(2*wb+wj));
//                 bx = 4/(wb*(4*wb+wj));
//                 ay = 4/((4*hb+hj)*(2*hb+hj));
//                 by = 4/(hb*(4*hb+hj));
//                 //printf("ax=%lf,bx=%lf,ay=%lf,by=%lf\n",ax,bx,ay,by);

//                 if(0<=dx && dx<=wb+wj/2){
//                     xov = 1-ax*(dx*dx);
//                 }else if(wb+wj/2<=dx && dx<=2*wb+wj/2){
//                     xov = bx*((dx-2*wb-wj/2)*(dx-2*wb-wj/2));
//                 }else if(2*wb+wj/2<=dx){
//                     xov = 0;
//                 }
//                 if(0<=dy && dy<=hb+hj/2){
//                     yov = 1-ay*(dy*dy);
//                 }else if(hb+hj/2<=dy && dy<=2*hb+hj/2){
//                     yov = by*((dy-2*hb-hj/2)*(dy-2*hb-hj/2));
//                 }else if(2*hb+hj/2<=dy){
//                     yov = 0;
//                 }
//                 //printf("xov=%lf,yov=%lf\n", xov,yov);
//                 bindensity[binnum*j+i] += xov*yov*c;  //一個bin的density
//                 //printf("d=%lf\n", d);
//                 f+=b*(bindensity[binnum*j+i]-tb)*(bindensity[binnum*j+i]-tb);

//                 for(int m = 0; m < _placement.numModules(); ++m){
// 				g[2 * m    ] += b * 2 * (binDensity[binnum*j+i] - tb) * g_temp[2 * m    ];
// 				g[2 * m + 1] += b * 2 * (binDensity[binnum*j+i] - tb) * g_temp[2 * m + 1];
// 			}
//             }
//             // if(d>tb){
//             //     //printf("d=%lf\n", d);
//             //     D += (d-tb)*(d-tb);  //加總bins den
//             // }
//             //printf("D=%lf\n", D);
//         }

//     } 
//     //printf("Dsum=%lf\n", D);
//     return D;
// }
// void ExampleFunction::grad_BinDensity(vector<double> &g){
//     for(int m=0;m<_placement.numModules(); ++m){
//         //printf("gx%lf,gy=%lf\n",g[m*2],g[m*2+1]);
//     }
//     //printf("e4\n");
//     double xsum,xnsum,ysum,ynsum,d,D,c,ax,bx,ay,by,axi,bxi,sign,db,Db;
//     double xb,yb,wb,hb, xj,yj,dx,dy,xov,yov,wj,hj,xi,yi,wi,hi,dxi;
//     double d1,dyi,xov1,yov1,ayi,byi,sign1;
//     d1=dyi=xov1=yov1=ayi=byi=sign1=0.0;
//     xsum=xnsum=ysum=ynsum=d=db=D=Db=c=ax=bx=ay=by=axi=bxi=sign=0.0;
//     xb=yb=wb=hb= xj=yj=dx=dy=xov=yov=wj=hj=xi=yi=wi=hi=dxi=0.0;
//     double coreWidth = _placement.boundryRight() - _placement.boundryLeft();
//     double coreHeight = _placement.boundryTop() - _placement.boundryBottom();
//     unsigned binnum = ceil(sqrt(_placement.numModules())/4);
//     //printf("binnum =%lu\n",binnum);
//     double binwidth = coreWidth/binnum;  
//     //printf("binwidth=%lf\n",binwidth);
//      double binheight = coreHeight/binnum;
//     //printf("binheight=%lf\n",binheight);
//     double result =0;
//     double tb = getTb();
//     //c = 2;
//     for(int i=0;i<binnum; ++i){
//         for(int j=0;j<binnum; ++j){
//             wb = binwidth;
//             hb = binheight;
//             xb = _placement.boundryLeft() + i*wb + 0.5*wb;
//             yb = _placement.boundryBottom() + j*hb + 0.5*hb;
//             //printf("xb=%lf,yb=%lf\n",xb,yb);
//             db=0;
//             for(int k=0;k<_placement.numModules(); ++k){
//                 xj = _placement.module(k).centerX();
//                 yj = _placement.module(k).centerY();
//                 wj = _placement.module(k).width();
//                 hj = _placement.module(k).height();
//                 c = _placement.module(k).area();
                
//                 dx = abs(xb-xj);
//                 dy = abs(yb-yj);

//                 ax = 4/((4*wb+wj)*(2*wb+wj));
//                 bx = 4/(wb*(4*wb+wj));
//                 ay = 4/((4*hb+hj)*(2*hb+hj));
//                 by = 4/(hb*(4*hb+hj));

//                 if(0<=dx && dx<=wb+wj/2){
//                     xov = 1-ax*(dx*dx);
//                 }else if(wb+wj/2<=dx && dx<=2*wb+wj/2){
//                     xov = bx*((dx-2*wb-wj/2)*(dx-2*wb-wj/2));
//                 }else if(2*wb+wj/2<=dx){
//                     xov = 0;
//                 }
//                 if(0<=dy && dy<=hb+hj/2){
//                     yov = 1-ay*(dy*dy);
//                 }else if(hb+hj/2<=dy && dy<=2*hb+hj/2){
//                     yov = by*((dy-2*hb-hj/2)*(dy-2*hb-hj/2));
//                 }else if(2*hb+hj/2<=dy){
//                     yov = 0;
//                 }
//                 db += xov*yov*(1/c);     //一個bin density    
//             }
//             if(db>tb){
//                 Db = db-tb; 
//                 //printf("Db=%lf\n",Db);
//             }else{
//                 Db = 0;
//             }
//            // printf("Db=%lf\n",Db);
//             //算gradx
//             d=d1=0;
//             for(int m=0;m<_placement.numModules(); ++m){
//                 if (_placement.module(m).isFixed()){
//                     g[m*2] = 0;
//                     g[m*2+1] = 0;
//                     continue;
//                 }
//                 xi = _placement.module(m).centerX();
//                 yi = _placement.module(m).centerY();
//                 wi = _placement.module(m).width();
//                 hi = _placement.module(m).height();
//                 c = _placement.module(m).area();
//                 //printf("xi=%lf,yi=%lf\n",xi,yi);
//                 dx = abs(xb-xi);
//                 dy = abs(yb-yi);
//                 if(xb<xi) sign = -1;
//                 else sign = 1;

                
//                 if(yb<yi) sign1 = -1;
//                 else sign1 = 1;

            
//                 ax = 4/((4*wb+wi)*(2*wb+wi));
//                 bx = 4/(wb*(4*wb+wi));
//                 ay = 4/((4*hb+hi)*(2*hb+hi));
//                 by = 4/(hb*(4*hb+hi));

//                 if(0<=dx && dx<=wb/2+wi/2){
//                 xov = (-2)*ax*dx;
//                 }else if(wb/2+wi/2<=dx && dx<=wb+wi/2){
//                 xov = 2*bx*(dx-2*wb-wi/2);
//                 }else if(wb+wi/2<=dx){
//                 xov = 0;
//                 }
//                 if(0<=dy && dy<=hb+hi/2){
//                     yov = 1-ay*(dy*dy);
//                 }else if(hb+hi/2<=dy && dy<=2*hb+hi/2){
//                     yov = by*((dy-2*hb-hi/2)*(dy-2*hb-hi/2));
//                 }else if(2*hb+hi/2<=dy){
//                     yov = 0;
//                 }
//                 d = xov*yov*(1);  //gradx

//                 if(0<=dx && dx<=wb+wi/2){
//                     xov1 = 1-ax*(dx*dx);
//                 }else if(wb+wi/2<=dx && dx<=2*wb+wi/2){
//                     xov1 = bx*((dx-2*wb-wi/2)*(dx-2*wb-wi/2));
//                 }else if(2*wb+wi/2<=dx){
//                     xov1 = 0;
//                 }
//                 if(0<=dy && dy<=hb+hi/2){
//                     yov1 = (-2)*ay*dy;
//                 }else if(hb+hi/2<=dy && dy<=2*hb+hi/2){
//                     yov1 = 2*by*(dy-2*hb-hi/2);
//                 }else if(2*hb+hi/2<=dy){
//                     yov1 = 0;
//                 }
//                 //printf("xov=%lf,yov=%lf\n",xov1,yov1);
//                 d1 = xov1*yov1*(1);  
//                 //printf("d=%lf,d1=%lf\n",d,d1);
//                 g[m*2] += 2*Db*d*0.6;
//                 g[m*2+1] += 2*Db*d1*0.6;
//             }
            
//         }
//     }
   
    
// }
// double ExampleFunction::BinDensity(){
//     //printf("e3\n");
//     double xsum,xnsum,ysum,ynsum,d,D,c,ax,bx,ay,by;
//     double xb,yb,wb,hb, xj,yj,dx,dy,xov,yov,wj,hj;
//     xsum=xnsum=ysum=ynsum=d=D=c=ax=bx=ay=by=0.0;
//     xb=yb=wb=hb= xj=yj=dx=dy=xov=yov=wj=hj=0.0;
//     double coreWidth = _placement.boundryRight() - _placement.boundryLeft();
//     double coreHeight = _placement.boundryTop() - _placement.boundryBottom();
//     unsigned binnum = ceil(sqrt(_placement.numModules()));
//     //printf("binnum =%lu\n",binnum);
//     double binwidth = coreWidth/binnum;  
//     //printf("binwidth=%lf\n",binwidth);
//     double binheight = coreHeight/binnum;
//     //printf("binheight=%lf\n",binheight);
//     double tb = getTb();
//     //printf("tb=%lf\n", tb);
//     //c = 2;
//     for(int i=0;i<binnum; ++i){
//         for(int j=0;j<binnum; ++j){
//             wb = binwidth;
//             hb = binheight;
//             xb = _placement.boundryLeft() + i*wb + 0.5*wb;
//             yb = _placement.boundryBottom() + j*hb + 0.5*hb;
//             //printf("wb=%lf,hb=%lf\n", wb,hb);
//             //printf("xb=%lf,yb=%lf\n", xb,yb);
//             d=0;
//             for(int k=0;k<_placement.numModules(); ++k){
//                 //printf("id=%d\n", k);
//                 //printf("wb=%lf,hb=%lf\n", wb,hb);
//                 //printf("xb=%lf,yb=%lf\n", xb,yb);
//                 xj = _placement.module(k).centerX();
//                 yj = _placement.module(k).centerY();
//                 wj = _placement.module(k).width();
//                 hj = _placement.module(k).height();
//                 c = _placement.module(k).area();
//                 //printf("xj=%lf,yj=%lf\n", xj,yj);
//                 //printf("wj=%lf,hj=%lf\n", wj,hj);
//                 dx = abs(xb-xj);
//                 dy = abs(yb-yj);
//                 //printf("dx=%lf,dy=%lf\n", dx,dy);
//                 ax = 4/((wb+wj)*(2*wb+wj));
//                 bx = 4/(wb*(2*wb+wj));
//                 ay = 4/((hb+hj)*(2*hb+hj));
//                 by = 4/(hb*(2*hb+hj));
//                 //printf("ax=%lf,bx=%lf,ay=%lf,by=%lf\n",ax,bx,ay,by);

//                 if(0<=dx && dx<=wb/2+wj/2){
//                     xov = 1-ax*(dx*dx);
//                 }else if(wb/2+wj/2<=dx && dx<=wb+wj/2){
//                     xov = bx*((dx-wb-wj/2)*(dx-wb-wj/2));
//                 }else if(wb+wj/2<=dx){
//                     xov = 0;
//                 }
//                 if(0<=dy && dy<=hb/2+hj/2){
//                     yov = 1-ay*(dy*dy);
//                 }else if(hb/2+hj/2<=dy && dy<=hb+hj/2){
//                     yov = by*((dy-hb-hj/2)*(dy-hb-hj/2));
//                 }else if(hb+hj/2<=dy){
//                     yov = 0;
//                 }
//                 //printf("xov=%lf,yov=%lf\n", xov,yov);
//                 d += xov*yov*(1/c);  //一個bin的density
//                 //printf("d=%lf\n", d);
//             }
//             if(d>tb){
//                 //printf("d=%lf\n", d);
//                 D += (d-tb)*(d-tb);  //加總bins den
//             }
//             //printf("D=%lf\n", D);
    
//         }
//     } 
//     //printf("Dsum=%lf\n", D);
//     return D;
// }
// void ExampleFunction::grad_BinDensity(vector<double> &g){
//     //printf("e4\n");
//     double xsum,xnsum,ysum,ynsum,d,D,c,ax,bx,ay,by,axi,bxi,sign,db,Db;
//     double xb,yb,wb,hb, xj,yj,dx,dy,xov,yov,wj,hj,xi,yi,wi,hi,dxi;
//     double d1,dyi,xov1,yov1,ayi,byi,sign1;
//     d1=dyi=xov1=yov1=ayi=byi=sign1=0.0;
//     xsum=xnsum=ysum=ynsum=d=db=D=Db=c=ax=bx=ay=by=axi=bxi=sign=0.0;
//     xb=yb=wb=hb= xj=yj=dx=dy=xov=yov=wj=hj=xi=yi=wi=hi=dxi=0.0;
//     double coreWidth = _placement.boundryRight() - _placement.boundryLeft();
//     double coreHeight = _placement.boundryTop() - _placement.boundryBottom();
//     unsigned binnum = ceil(sqrt(_placement.numModules()));
//     //printf("binnum =%lu\n",binnum);
//     double binwidth = coreWidth/binnum;  
//     //printf("binwidth=%lf\n",binwidth);
//      double binheight = coreHeight/binnum;
//     //printf("binheight=%lf\n",binheight);
//     double result =0;
//     double tb = getTb();
//     //c = 2;
//     for(int i=0;i<binnum; ++i){
//         for(int j=0;j<binnum; ++j){
//             wb = binwidth;
//             hb = binheight;
//             xb = _placement.boundryLeft() + i*wb + 0.5*wb;
//             yb = _placement.boundryBottom() + j*hb + 0.5*hb;
//             //printf("xb=%lf,yb=%lf\n",xb,yb);
//             db=0;
//             for(int k=0;k<_placement.numModules(); ++k){
//                 xj = _placement.module(k).centerX();
//                 yj = _placement.module(k).centerY();
//                 wj = _placement.module(k).width();
//                 hj = _placement.module(k).height();
                
//                 dx = abs(xb-xj);
//                 dy = abs(yb-yj);

//                 ax = 4/((wb+wj)*(2*wb+wj));
//                 bx = 4/(wb*(2*wb+wj));
//                 ay = 4/((hb+hj)*(2*hb+hj));
//                 by = 4/(hb*(2*hb+hj));

//                 if(0<=dx && dx<=wb/2+wj/2){
//                     xov = 1-ax*(dx*dx);
//                 }else if(wb/2+wj/2<=dx && dx<=wb+wj/2){
//                     xov = bx*((dx-wb-wj)*(dx-wb-wj));
//                 }else if(wb+wj/2<=dx){
//                     xov = 0;
//                 }
//                 if(0<=dy && dy<=hb/2+hj/2){
//                     yov = 1-ay*(dy*dy);
//                 }else if(hb/2+hj/2<=dy && dy<=hb+hj/2){
//                     yov = by*((dy-hb-hj)*(dy-hb-hj));
//                 }else if(hb+hj/2<=dy){
//                     yov = 0;
//                 }
//                 db += xov*yov;     //一個bin density    
//             }
//             if(db>tb){
//                 Db = db-tb; 
//                 //printf("Db=%lf\n",Db);
//             }else{
//                 Db = 0;
//             }
//            // printf("Db=%lf\n",Db);
//             //算gradx
//             d=d1=0;
//             for(int m=0;m<_placement.numModules(); ++m){
//                 xi = _placement.module(m).centerX();
//                 yi = _placement.module(m).centerY();
//                 wi = _placement.module(m).width();
//                 hi = _placement.module(m).height();
//                 //printf("xi=%lf,yi=%lf\n",xi,yi);
//                 dx = abs(xb-xi);
//                 dy = abs(yb-yi);
//                 if(xb<xi) sign = -1;
//                 else sign = 1;

                
//                 if(yb<yi) sign1 = -1;
//                 else sign1 = 1;

            
//                 ax = 4/((wb+wi)*(2*wb+wi));
//                 bx = 4/(wb*(2*wb+wi));
//                 ay = 4/((hb+hj)*(2*hb+hj));
//                 by = 4/(hb*(2*hb+hj));

//                 if(0<=dx && dx<=wb/2+wj/2){
//                 xov = (-2)*ax*(xb-xi)*sign;
//                 }else if(wb/2+wj/2<=dx && dx<=wb+wj/2){
//                 xov = 2*bx*((xb-xi)-wb-wj/2)*sign;
//                 }else if(wb+wj/2<=dx){
//                 xov = 0;
//                 }
//                 if(0<=dy && dy<=hb/2+hj/2){
//                 yov = 1-ay*(dy*dy);
//                 }else if(hb/2+hj/2<=dy && dy<=hb+hj/2){
//                 yov = by*((dy-hb-hj/2)*(dy-hb-hj/2));
//                 }else if(hb+hj/2<=dy){
//                 yov = 0;
//                 }
//                 d = xov*yov;  //gradx

//                 if(0<=dx && dx<=wb/2+wj/2){
//                     xov1 = 1-ax*(dx*dx);
//                 }else if(wb/2+wj/2<=dx && dx<=wb+wj/2){
//                     xov1 = bx*((dx-wb-wj/2)*(dx-wb-wj/2));
//                 }else if(wb+wj/2<=dx){
//                     xov1 = 0;
//                 }
//                 if(0<=dy && dy<=hb/2+hj/2){
//                     yov1 = (-2)*ay*(yb-yi)*sign1;
//                 }else if(hb/2+hj/2<=dy && dy<=hb+hj/2){
//                     yov1 = 2*by*((yb-yi)-hb-hj/2)*sign1;
//                 }else if(hb+hj/2<=dy){
//                     yov1 = 0;
//                 }
//                 //printf("xov=%lf,yov=%lf\n",xov1,yov1);
//                 d1 = xov1*yov1;  
//                 //printf("d=%lf,d1=%lf\n",d,d1);
//                 g[m*2] += 2*Db*d;
//                 g[m*2+1] += 2*Db*d1;
//             }
            
//         }
//     }
    
// }
// double ExampleFunction::gradx_BinDensity(int id){
//     //printf("e4\n");
//     double wl;
//     double a = 0.15;
//     double xsum,xnsum,ysum,ynsum,d,D,c,ax,bx,ay,by,axi,bxi,sign,db,Db;
//     double xb,yb,wb,hb, xj,yj,dx,dy,xov,yov,wj,hj,xi,yi,wi,hi,dxi;
//     xsum=xnsum=ysum=ynsum=d=db=D=Db=c=ax=bx=ay=by=axi=bxi=sign=0.0;
//     xb=yb=wb=hb= xj=yj=dx=dy=xov=yov=wj=hj=xi=yi=wi=hi=dxi=0.0;
//     double coreWidth = _placement.boundryRight() - _placement.boundryLeft();
//     double coreHeight = _placement.boundryTop() - _placement.boundryBottom();
//     unsigned binnum = ceil(sqrt(_placement.numModules()));
//     //printf("binnum =%lu\n",binnum);
//     double binwidth = coreWidth/binnum;  
//     //printf("binwidth=%lf\n",binwidth);
//      double binheight = coreHeight/binnum;
//     //printf("binheight=%lf\n",binheight);
//     double result =0;
//     double tb = getTb();
//     //c = 2;
//     for(int i=0;i<binnum; ++i){
//         for(int j=0;j<binnum; ++j){
//             wb = binwidth;
//             hb = binheight;
//             xb = _placement.boundryLeft() + i*wb + 0.5*wb;
//             yb = _placement.boundryBottom() + j*hb + 0.5*hb;
//             db=0;
//             for(int k=0;k<_placement.numModules(); ++k){
//                 xj = _placement.module(k).centerX();
//                 yj = _placement.module(k).centerY();
//                 wj = _placement.module(k).width();
//                 hj = _placement.module(k).height();
                
//                 dx = abs(xb-xj);
//                 dy = abs(yb-yj);

//                 ax = 4/((wb+wj)*(2*wb+wj));
//                 bx = 4/(wb*(2*wb+wj));
//                 ay = 4/((hb+hj)*(2*hb+hj));
//                 by = 4/(hb*(2*hb+hj));

//                 if(0<=dx && dx<=wb/2+wj/2){
//                     xov = 1-ax*(dx*dx);
//                 }else if(wb/2+wj/2<=dx && dx<=wb+wj/2){
//                     xov = bx*((dx-wb-wj)*(dx-wb-wj));
//                 }else if(wb+wj/2<=dx){
//                     xov = 0;
//                 }
//                 if(0<=dy && dy<=hb/2+hj/2){
//                     yov = 1-ay*(dy*dy);
//                 }else if(hb/2+hj/2<=dy && dy<=hb+hj/2){
//                     yov = by*((dy-hb-hj)*(dy-hb-hj));
//                 }else if(hb+hj/2<=dy){
//                     yov = 0;
//                 }
//                 db += xov*yov;     //一個bin density    
//             }
//             Db = db-tb; 
//             //算gradx

//             xi = _placement.module(id).centerX();
//             yi = _placement.module(id).centerY();
//             wi = _placement.module(id).width();
//             hi = _placement.module(id).height();

//             dxi = abs(xb-xi);
//             dy = abs(yb-yi);
//             if(xb<xi) sign = -1;
//             else sign = 1;

            
//             axi = 4/((wb+wi)*(2*wb+wi));
//             bxi = 4/(wb*(2*wb+wi));
//             ay = 4/((hb+hj)*(2*hb+hj));
//             by = 4/(hb*(2*hb+hj));

//             if(0<=dx && dx<=wb/2+wj/2){
//                 xov = (-2)*axi*(xb-xi)*sign;
//             }else if(wb/2+wj/2<=dx && dx<=wb+wj/2){
//                 xov = 2*bxi*((xb-xi)-wb-wj/2)*sign;
//             }else if(wb+wj/2<=dx){
//                 xov = 0;
//             }
//             if(0<=dy && dy<=hb/2+hj/2){
//                 yov = 1-ay*(dy*dy);
//             }else if(hb/2+hj/2<=dy && dy<=hb+hj/2){
//                 yov = by*((dy-hb-hj/2)*(dy-hb-hj/2));
//             }else if(hb+hj/2<=dy){
//                 yov = 0;
//             }
//             d = xov*yov;  //gradx
           
//             result += 2*Db*d;
//         }
//     }
//     return result;
// }
// double ExampleFunction::grady_BinDensity(int id){
//     //printf("e4\n");
//     double wl;
//     double a = 0.15;
//     double xsum,xnsum,ysum,ynsum,d,D,c,ax,bx,ay,by,ayi,byi,sign,db,Db;
//     double xb,yb,wb,hb, xj,yj,dx,dy,xov,yov,wj,hj,xi,yi,wi,hi,dyi;
//     xsum=xnsum=ysum=ynsum=d=db=D=Db=c=ax=bx=ay=by=ayi=byi=sign=0.0;
//     xb=yb=wb=hb= xj=yj=dx=dy=xov=yov=wj=hj=xi=yi=wi=hi=dyi=0.0;
//     double coreWidth = _placement.boundryRight() - _placement.boundryLeft();
//     double coreHeight = _placement.boundryTop() - _placement.boundryBottom();
//     unsigned binnum = ceil(sqrt(_placement.numModules()));
//     //printf("binnum =%lu\n",binnum);
//     double binwidth = coreWidth/binnum;  
//     //printf("binwidth=%lf\n",binwidth);
//      double binheight = coreHeight/binnum;
//     //printf("binheight=%lf\n",binheight);
//     double result =0;
//     double tb = getTb();
//     //c = 2;
//     for(int i=0;i<binnum; ++i){
//         for(int j=0;j<binnum; ++j){
//             wb = binwidth;
//             hb = binheight;
//             xb = _placement.boundryLeft() + i*wb + 0.5*wb;
//             yb = _placement.boundryBottom() + j*hb + 0.5*hb;
//             db=0;
//             for(int k=0;k<_placement.numModules(); ++k){
//                 xj = _placement.module(k).centerX();
//                 yj = _placement.module(k).centerY();
//                 wj = _placement.module(k).width();
//                 hj = _placement.module(k).height();
                
//                 dx = abs(xb-xj);
//                 dy = abs(yb-yj);

//                 ax = 4/((wb+wj)*(2*wb+wj));
//                 bx = 4/(wb*(2*wb+wj));
//                 ay = 4/((hb+hj)*(2*hb+hj));
//                 by = 4/(hb*(2*hb+hj));

//                 if(0<=dx && dx<=wb/2+wj/2){
//                     xov = 1-ax*(dx*dx);
//                 }else if(wb/2+wj/2<=dx && dx<=wb+wj/2){
//                     xov = bx*((dx-wb-wj)*(dx-wb-wj));
//                 }else if(wb+wj/2<=dx){
//                     xov = 0;
//                 }
//                 if(0<=dy && dy<=hb/2+hj/2){
//                     yov = 1-ay*(dy*dy);
//                 }else if(hb/2+hj/2<=dy && dy<=hb+hj/2){
//                     yov = by*((dy-hb-hj)*(dy-hb-hj));
//                 }else if(hb+hj/2<=dy){
//                     yov = 0;
//                 }
//                 db += xov*yov;     //一個bin density    
//             }
//             Db = db-tb; 
//             //算gradx

//             xi = _placement.module(id).centerX();
//             yi = _placement.module(id).centerY();
//             wi = _placement.module(id).width();
//             hi = _placement.module(id).height();

//             dyi = abs(yb-yi);
//             dx = abs(xb-xi);
//             if(yb<yi) sign = -1;
//             else sign = 1;

            
//             ayi = 4/((hb+hi)*(2*hb+hi));
//             byi = 4/(hb*(2*hb+hi));
//             ax = 4/((wb+wj)*(2*wb+wj));
//             bx = 4/(wb*(2*wb+wj));

//             if(0<=dx && dx<=wb/2+wj/2){
//                 xov = 1-ax*(dx*dx);
//             }else if(wb/2+wj/2<=dx && dx<=wb+wj/2){
//                 xov = bx*((dx-wb-wj/2)*(dx-wb-wj/2));
//             }else if(wb+wj/2<=dx){
//                 xov = 0;
//             }
//             if(0<=dy && dy<=hb/2+hj/2){
//                 yov = (-2)*ayi*(yb-yi)*sign;
//             }else if(hb/2+hj/2<=dy && dy<=hb+hj/2){
//                 yov = 2*byi*((yb-yi)-hb-hj/2)*sign;
//             }else if(hb+hj/2<=dy){
//                 yov = 0;
//             }
//             d = xov*yov; //grady
           
//             result += 2*Db*d;
//         }
//     }
//     return result;
// }
// double ExampleFunction::grady_BinDensity( int id){
//     printf("e5\n");
//     double wl=0.0;
//     double a = 0.15;
//     double xsum,xnsum,ysum,ynsum,d,D,c,ax,bx,ay,by,ayi,byi,sign;
//     xsum=xnsum=ysum=ynsum=d=D=c=ax=bx=ay=by=ayi=byi=sign=0.0;
//     double xb,yb,wb,hb, xj,yj,dx,dy,xov,yov,wj,hj,xi,yi,wi,hi,dyi;
//     xb=yb=wb=hb= xj=yj=dx=dy=xov=yov=wj=hj=xi=yi=wi=hi=dyi=0.0;
//     double coreWidth = _placement.boundryRight() - _placement.boundryLeft();
//     double coreHeight = _placement.boundryTop() - _placement.boundryBottom();
//     unsigned binnum = ceil(sqrt(_placement.numModules()));
//     //printf("binnum =%lu\n",binnum);
//     double binwidth = coreWidth/binnum;  
//     //printf("binwidth=%lf\n",binwidth);
//      double binheight = coreHeight/binnum;
//     //printf("binheight=%lf\n",binheight);
//     double tb = getTb();
//     //c = 2;
//     for(int i=0;i<binnum; ++i){
//         for(int j=0;j<binnum; ++j){
//             wb = binwidth;
//             hb = binheight;
//             xb = _placement.boundryLeft() + i*wb + 0.5*wb;
//             yb = _placement.boundryBottom() + j*hb + 0.5*hb;

//             xi = _placement.module(id).centerX();
//             yi = _placement.module(id).centerY();
//             wi = _placement.module(id).width();
//             hi = _placement.module(id).height();

//             dyi = abs(yb-yi);
//             dx = abs(xb-xi);
//             if(yb<yi) sign = -1;
//             else sign = 1;

            
//             ayi = 4/((hb+hi)*(2*hb+hi));
//             byi = 4/(hb*(2*hb+hi));
//             ax = 4/((wb+wj)*(2*wb+wj));
//             bx = 4/(wb*(2*wb+wj));

//             if(0<=dx && dx<=wb/2+wj/2){
//                 xov = 1-ax*(dx*dx);
//             }else if(wb/2+wj/2<=dx && dx<=wb+wj/2){
//                 xov = bx*((dx-wb-wj/2)*(dx-wb-wj/2));
//             }else if(wb+wj/2<=dx){
//                 xov = 0;
//             }
//             if(0<=dy && dy<=hb/2+hj/2){
//                 yov = (-2)*ayi*(yb-yi)*sign;
//             }else if(hb/2+hj/2<=dy && dy<=hb+hj/2){
//                 yov = 2*byi*((yb-yi)-hb-hj/2)*sign;
//             }else if(hb+hj/2<=dy){
//                 yov = 0;
//             }
//             d = xov*yov-tb;
//             D += d;

//             for(int k=0;k<_placement.numModules(); ++k){
//                 if(k!=id){
//                 xj = _placement.module(k).centerX();
//                 yj = _placement.module(k).centerY();
//                 wj = _placement.module(k).width();
//                 hj = _placement.module(k).height();
                
//                 dx = abs(xb-xj);
//                 dy = abs(yb-yj);

//                 ax = 4/((wb+wj)*(2*wb+wj));
//                 bx = 4/(wb*(2*wb+wj));
//                 ay = 4/((hb+hj)*(2*hb+hj));
//                 by = 4/(hb*(2*hb+hj));

//                 if(0<=dx && dx<=wb/2+wj/2){
//                     xov = 1-ax*(dx*dx);
//                 }else if(wb/2+wj/2<=dx && dx<=wb+wj/2){
//                     xov = bx*((dx-wb-wj)*(dx-wb-wj));
//                 }else if(wb+wj/2<=dx){
//                     xov = 0;
//                 }
//                 if(0<=dy && dy<=hb/2+hj/2){
//                     yov = 1-ay*(dy*dy);
//                 }else if(hb/2+hj/2<=dy && dy<=hb+hj/2){
//                     yov = by*((dy-hb-hj)*(dy-hb-hj));
//                 }else if(hb+hj/2<=dy){
//                     yov = 0;
//                 }
//                 d = xov*yov-tb;
//                 D += d;
//                 }
//             }
//         }
//     }
//     return D;
// }
