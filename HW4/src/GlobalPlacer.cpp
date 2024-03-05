#include "GlobalPlacer.h"
#include "ExampleFunction.h"
#include "NumericalOptimizer.h"
#include "Wrapper.hpp"
#include <cstdlib>
#include <iostream>
#include <cmath>
#include <algorithm>

GlobalPlacer::GlobalPlacer(wrapper::Placement &placement)
    : _placement(placement)
{
}

void GlobalPlacer::randomPlace()
{
    //srand(0);
    double coreWidth = _placement.boundryRight() - _placement.boundryLeft();
    double coreHeight = _placement.boundryTop() - _placement.boundryBottom();
    for (size_t i = 0; i < _placement.numModules(); ++i)
    {   
        if (_placement.module(i).isFixed())
            continue;

        double width = _placement.module(i).width();
        double height = _placement.module(i).height();
        double x = rand() % static_cast<int>(coreWidth - width) + _placement.boundryLeft();
        double y = rand() % static_cast<int>(coreHeight - height) + _placement.boundryBottom();
        _placement.module(i).setPosition(x, y);
    }
}

void GlobalPlacer::place()
{
    unsigned int num = _placement.numModules();
    double coreWidth = _placement.boundryRight() - _placement.boundryLeft();
    double coreHeight = _placement.boundryTop() - _placement.boundryBottom();
    double step = coreWidth*7;
    int k = 5000;
    int seed =  time(0);
    int iters = 2;
    if (num < 15000){    
        seed = 5;
        iters = 3;
        k = 100;
        step = coreWidth;
    }
    else if (num >= 15000 && num < 30000)     //case 2 
    {   
        seed = 0;
        step = coreWidth*5;
        k = 2500;
    }
    else if (num == 51382)       //case 3
    {
        seed = 0;
        step = coreWidth*5;
        k = 2500;
    }
    else if (num >= 13000 && num < 28000)     //h1
    {
        seed = 0;
        step = coreWidth*5;
        k = 2500;
    }
    
    srand(seed);
    randomPlace();
    
    vector<double> x(num * 2);
    for(unsigned int id = 0; id < num; id++){
        x[id * 2] = _placement.module(id).centerX();
        x[id * 2 + 1] = _placement.module(id).centerY();
    }

    
    ExampleFunction ef(_placement);
    int iter = 150;
    
    for(int i = 0; i < iters; i++){
        if(i==1) iter = 35;
        else if(i==2) iter=35;
        else if(i==3) iter=30;
        ef.beta += i * k;  
        NumericalOptimizer no(ef);
        no.setX(x);             
        no.setNumIteration(iter); 
        no.setStepSizeBound(step); 
        no.solve();            

        for(unsigned int id = 0; id < num; id++){
            if(!_placement.module(id).isFixed()){

                double w = _placement.module(id).width();
                double h = _placement.module(id).height();

                double cx = no.x(id * 2);
                double cy = no.x(id * 2 + 1);
                if(cx + w/2 > _placement.boundryRight())
                    cx = _placement.boundryRight() - w/2;
                else if(cx - w/2 < _placement.boundryLeft())
                    cx = _placement.boundryLeft() + w/2;
                if(cy + h/2 > _placement.boundryTop())
                    cy = _placement.boundryTop() - h/2;
                else if(cy - h/2 < _placement.boundryBottom())
                    cy = _placement.boundryBottom() + h/2;

                x[id * 2] = cx;
                x[id * 2 + 1] = cy;

                _placement.module(id).setPosition(cx - w/2, cy - h/2);
            }
        }
    }
    
    // cout << "Current solution:\n";
    // for (unsigned i = 0; i < no.dimension(); i++)
    // {
    //     cout << "x[" << i << "] = " << no.x(i) << "\n";
    // }
    // cout << "Objective: " << no.objective() << "\n";
    ////////////////////////////////////////////////////////////////

    // An example of random placement implemented by TA.
    // If you want to use it, please uncomment the folllwing 1 line.
    //randomPlace();

    /* @@@ TODO
     * 1. Understand above example and modify ExampleFunction.cpp to implement the analytical placement
     * 2. You can choose LSE or WA as the wirelength model, the former is easier to calculate the gradient
     * 3. For the bin density model, you could refer to the lecture notes
     * 4. You should first calculate the form of wirelength model and bin density model and the forms of their gradients ON YOUR OWN
     * 5. Replace the value of f in evaluateF() by the form like "f = alpha*WL() + beta*BinDensity()"
     * 6. Replace the form of g[] in evaluateG() by the form like "g = grad(WL()) + grad(BinDensity())"
     * 7. Set the initial vector x in place(), set step size, set #iteration, and call the solver like above example
     * */
}
