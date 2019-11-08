//2D, 5-equation (Allaire), two material, inert model with Mie Gruneisen EoS

#include <iostream>
#include <cmath>
#include <string>
#include <vector>
#include <fstream>
#include <assert.h>
#include <iomanip>
using namespace std;


//denote floats with decimals or auto interpreted as integer


typedef vector<double> Vector;

template<typename T>
void printVector(vector<T> vec){
    int len = vec.size();
    for(int i; i<len; ++i){
        cout << vec[i] << endl;
    }
}//vector-array print function templated to print for any vector type


#include "Classes.H"

//e_ref and p_ref functions:

double p_ref_calc(MaterialProperties M, string EoS, double rho){

    double p_ref;
    
    if(EoS == "ideal"){
        p_ref = 0;
    }else if(EoS == "stiffened"){
        p_ref = -(M.varGamma+1)*M.p_inf;
    }else if(EoS == "JWL"){
        p_ref = M.A*exp(-M.R1*M.rho_0/rho) + M.B*exp(-M.R2*M.rho_0/rho);
    }else if(EoS == "CC"){
        p_ref = M.A*pow((M.rho_0/rho),-M.xi1) - M.B*pow((M.rho_0/rho),-M.xi2); 
    }else if(EoS == "Hugoniot"){
       if(M.c_0 == 0){
            p_ref = 0.0;
        }else{
            p_ref = M.p_0 + (M.rho_0*rho*pow(M.c_0,2)*(rho-M.rho_0))/pow(rho - M.s*(rho-M.rho_0),2);
        } 
    }else{
        cout << "undefined EoS" << endl;
    }

    return p_ref;
}


double e_ref_calc(MaterialProperties M, string EoS, double rho){
    
    double e_ref;
    
    if(EoS == "ideal" || EoS == "stiffened"){
        e_ref = M.e_0;
    }else if(EoS == "JWL"){
        e_ref = (M.A/(M.rho_0*M.R1))*exp(-M.R1*M.rho_0/rho) + \
                (M.B/(M.rho_0*M.R2))*exp(-M.R2*M.rho_0/rho);
    }else if(EoS == "CC"){
        e_ref = -(M.A/(M.rho_0*(1-M.xi1)))*(pow((rho/M.rho_0),M.xi1-1)-1) + \
                (M.B/(M.rho_0*(1-M.xi2)))*(pow((rho/M.rho_0),M.xi2-1)-1);
    }else if(EoS == "Hugoniot"){
       if(M.c_0 == 0){
            e_ref = 0.0;
        }else{
            e_ref = M.e_0 + (p_ref_calc(M, EoS, rho) + M.p_0)*(rho-M.rho_0)/(2*rho*M.rho_0);
        } 
    }else{
        cout << "undefined EoS" << endl;
    }

    return e_ref;
}
        

#include "TestCases.H"

//~~~~ USER DEFINED SETTINGS HERE ~~~~//
TestCase CASE("WeakShockCavityCollapse");   //test cases defined in TestCases.H file
    
bool MUSCL = 1;                             //second order extension scheme: on/off
string scheme = "UltraBee";                 //slope limiter types - 
                                            //"Superbee" "MinBee" "UltraBee" or "VanLeer"

int view_resolution = 200;                  //number of time steps written to file for visualisation

bool time_stepping = 1;                     //either run full time or just initial construction
bool surface_tension_dt = 0;                //restrict time step based on ST terms
bool curvature = 1;

string material1 = CASE.material1;
string material2 = CASE.material2;
string EoS = CASE.EoS;
double sigma = CASE.sigma;
MaterialProperties M1(material1, EoS);
MaterialProperties M2(material2, EoS);
double varGamma1 = M1.varGamma;
double varGamma2 = M2.varGamma;
double k_avr;

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//


class Properties{
    public:
        Prim W;
        ConsU U;
        ConsF F;
        ConsG G;
        double rho;
        double e;
        double E;
        double eref;
        double pref;
        double eref1, eref2;
        double pref1, pref2;
        double zeta;
        double z2;
        double S; //sound speed
        bool bk;  //binary curvature parameter 
        double k; //curvature parameter
        double z1_avr; //used in curvature calc
        double n_x, n_y; //x and y normal components at surface
        double w;        //weighting parameter for surface tension
        double n_xx, n_yy, n_xy, n_yx; //for ST computation
        double n_x_avr, n_y_avr;
        double k_temp;          //for smoothing/weighting iteration
        Properties(Prim _W):
            W(0,0,0,0,0,0), U(0,0,0,0,0,0), F(0,0,0,0,0,0), G(0,0,0,0,0,0){
            eref1 = e_ref_calc(M1, EoS, _W.rho1); //check this is rho1-rho2 and not rho
            eref2 = e_ref_calc(M2, EoS, _W.rho2);
            pref1 = p_ref_calc(M1, EoS, _W.rho1);
            pref2 = p_ref_calc(M2, EoS, _W.rho2);
            W = _W;
            z2 = 1-W.z1;
            rho = W.z1*W.rho1 + z2*W.rho2;
            zeta = W.z1/(varGamma1) + z2/(varGamma2);
            eref = (1/rho)*(W.z1*W.rho1*eref1 + z2*W.rho2*eref2);
            pref = (1/zeta)*(W.z1*pref1/(varGamma1) + z2*pref2/(varGamma2));
            e = eref + zeta*(W.p - pref)/rho;
            E = 0.5*rho*(pow(W.v,2) + pow(W.u,2)) + rho*e;
            U.z1rho1 = W.z1*W.rho1;
            U.z2rho2 = z2*W.rho2;
            U.rhou = rho*W.u;
            U.rhov = rho*W.v;
            U.E = E;
            U.z1 = W.z1;
            F.mass1 = W.z1*W.rho1*W.u;
            F.mass2 = z2*W.rho2*W.u;
            F.xmom = rho*pow(W.u,2) + W.p;
            F.ymom = rho*W.u*W.v;
            F.en = W.u*(E+W.p);
            F.uz1 = W.u*W.z1;
            G.mass1 = W.z1*W.rho1*W.v;
            G.mass2 = z2*W.rho2*W.v;
            G.xmom = rho*W.u*W.v;
            G.ymom = rho*pow(W.v,2) + W.p;
            G.en = W.v*(E+W.p);
            G.vz1 = W.v*W.z1;
            }
};

//initial vector construction:
vector<Properties> initial_construction(TestCase CASE){ 
                
    vector<Properties> IC_Vector;   //initial state vector
                
    if(CASE.construction == "default2D"){
        //Default construction
        Properties propLHS(CASE.WL);
        Properties propRHS(CASE.WR);
        IC_Vector.resize(CASE.n*CASE.m, propRHS);
        for(int j = 0; j < CASE.m; ++j){
            for(int i = 0; i < CASE.nG + (int)(CASE.x0/CASE.L * (CASE.n-2*CASE.nG)); ++i){
                IC_Vector[j*CASE.n + i] = propLHS;
            }
        }
    }else if(CASE.construction == "3part"){
        //Default construction
        Properties propLHS(CASE.WL);
        Properties propRHS(CASE.WR);
        IC_Vector.resize(CASE.n*CASE.m, propRHS);
        for(int j = 0; j < CASE.m; ++j){
            for(int i = CASE.nG + (int)(CASE.x0/CASE.L * (CASE.n-2*CASE.nG));\
                    i < CASE.nG + (int)(CASE.x1/CASE.L * (CASE.n-2*CASE.nG)); ++i){
                IC_Vector[j*CASE.n + i] = propLHS;
            }
        }
    }else if(CASE.construction == "square"){
        //Default construction
        Properties propLHS(CASE.WL);
        Properties propRHS(CASE.WR);
        IC_Vector.resize(CASE.n*CASE.m, propRHS);
        for(int j = 0; j < CASE.nG + (int)(CASE.y0/CASE.L * (CASE.m-2*CASE.nG)); ++j){
            for(int i = 0; i < CASE.nG + (int)(CASE.x0/CASE.L * (CASE.n-2*CASE.nG)); ++i){
                IC_Vector[j*CASE.n + i] = propLHS;
            }
        }
    }else if(CASE.construction == "y-direction"){
        Properties propLHS(CASE.WL);
        Properties propRHS(CASE.WR);
        IC_Vector.resize(CASE.n*CASE.m, propRHS);
        for(int j = 0; j < CASE.nG + (int)(CASE.y0/CASE.L * (CASE.m-2*CASE.nG)); ++j){
            for(int i = 0; i < CASE.n; ++i){
                IC_Vector[j*CASE.n + i] = propLHS;
            }
        }
    }else if(CASE.construction == "bubble"){
        Properties propAir(CASE.WL);
        Properties propWater(CASE.WR);
        //fill with air:
        IC_Vector.resize(CASE.n*CASE.m, propAir);
        //fill bottom portion with water:
        for(int j = 0; j < CASE.nG + (int)(CASE.y0/CASE.Y * (CASE.m-2*CASE.nG)); ++j){
            for(int i = 0; i < CASE.n; ++i){
                IC_Vector[j*CASE.n + i] = propWater;
            }
        }
        //add the pressurised air bubble:
        propAir.W.p = 1.0e9;
        Properties HP_Air(propAir.W); //high pressure air
        for(int j = 0; j < CASE.m; ++j){
            for(int i = 0; i < CASE.n; ++i){
                if((pow(i*CASE.dx-CASE.o_x, 2)+pow(j*CASE.dy - CASE.o_y,2)) <= pow(CASE.R,2)){
                    IC_Vector[j*CASE.n + i] = HP_Air;
                }
            }
        }
    }else if(CASE.construction == "air-bubble"){
        Properties propAir(CASE.WL);
        Properties propWater(CASE.WR);
        //fill with water:
        IC_Vector.resize(CASE.n*CASE.m, propWater);
        //fill bottom portion with air:
        for(int j = 0; j < CASE.nG + (int)(CASE.y0/CASE.Y * (CASE.m-2*CASE.nG)); ++j){
            for(int i = 0; i < CASE.n; ++i){
                IC_Vector[j*CASE.n + i] = propAir;
            }
        }
        //add the pressurised air bubble:
        propAir.W.p = 1.0e9;
        Properties HP_Air(propAir.W); //high pressure air
        for(int j = 0; j < CASE.m; ++j){
            for(int i = 0; i < CASE.n; ++i){
                if((pow(i*CASE.dx-CASE.o_x, 2)+pow(j*CASE.dy - CASE.o_y,2)) <= pow(CASE.R,2)){
                    IC_Vector[j*CASE.n + i] = HP_Air;
                }
            }
        }
    }else if(CASE.construction == "shocked-cavity"){
        Properties propAir(CASE.WL);
        Properties propWater(CASE.WR);
        //fill with water:
        IC_Vector.resize(CASE.n*CASE.m, propWater);
        //create shocked water state:
        if(CASE.title == "CavityCollapse"){
            propWater.W.rho1 = 1.59;
            propWater.W.rho2 = 1325;
            propWater.W.u = 680.525;
            propWater.W.p = 1.9153e9;
        }else if(CASE.title == "WeakShockCavityCollapse"){
            propWater.W.rho1 = 1.32;
            propWater.W.rho2 = 1117.73;
            propWater.W.u = 150.0;
            propWater.W.p = 0.28918e9;
        }else if(CASE.title == "CavityCollapseST"){
            propWater.W.rho1 = 1.2;
            propWater.W.rho2 = 1033.3;
            propWater.W.u = 5.0;
            propWater.W.p = 8.1536e6;
        }else if(CASE.title == "MicroBubble"){
            propWater.W.rho1 = 1.23;
            propWater.W.rho2 = 1030.0;
            propWater.W.u = 50.0;
            propWater.W.p = 8.483e7;
        }else if(CASE.title == "CavityCollapse1" || CASE.title == "CavityCollapse3" || CASE.title == "CavityCollapse9"){
            propWater.W.rho1 = 1.32;
            propWater.W.rho2 = 1117.0;
            propWater.W.u = 150.0;
            propWater.W.p = 0.28918e9;
        }else if(CASE.title == "CavityCollapse4" || CASE.title == "CavityCollapse7"){
            propWater.W.rho1 = 1.35;
            propWater.W.rho2 = 1162.21;
            propWater.W.u = 250.0;
            propWater.W.p = 0.5435e9;
        }else if(CASE.title == "CavityCollapse2" || CASE.title == "CavityCollapse5" || CASE.title == "CavityCollapse8"){
            propWater.W.rho1 = 1.52;
            propWater.W.rho2 = 1267.0;
            propWater.W.u = 685.0;
            propWater.W.p = 1.9138e9;
        }else if(CASE.title == "CavityCollapse6"){
            propWater.W.rho1 = 1.58;
            propWater.W.rho2 = 1319.69;
            propWater.W.u = 1010.0;
            propWater.W.p = 3.516e9;
        }else{
            throw runtime_error("invalid title for construction type");
        }
        Properties ShockedWater(propWater.W);
        //fill left portion with shocked water:
        for(int j = 0; j < CASE.m; ++j){
            for(int i = 0; i < CASE.nG + (int)(CASE.x0/CASE.L * (CASE.n-2*CASE.nG)); ++i){
                IC_Vector[j*CASE.n + i] = ShockedWater;
            }
        }
        //add the air cavity:
        for(int j = 0; j < CASE.m; ++j){
            for(int i = 0; i < CASE.n; ++i){
                if((pow(i*CASE.dx-CASE.o_x, 2)+pow(j*CASE.dy - CASE.o_y,2)) <= pow(CASE.R,2)){
                    IC_Vector[j*CASE.n + i] = propAir;
                }
            }
        }
    }else if(CASE.construction == "shocked-cavity-nitromethane"){
        Properties propAir(CASE.WL);
        Properties propNitro(CASE.WR);
        //fill with liquid nitromethane:
        IC_Vector.resize(CASE.n*CASE.m, propNitro);
        //create shocked nitromethane state:
        propNitro.W.rho1 = 2.4;
        propNitro.W.rho2 = 1934;
        propNitro.W.u = 2000.0;
        propNitro.W.p = 10.98e9;
        Properties ShockedNitro(propNitro.W);
        //fill left portion with shocked water:
        for(int j = 0; j < CASE.m; ++j){
            for(int i = 0; i < CASE.nG + (int)(CASE.x0/CASE.L * (CASE.n-2*CASE.nG)); ++i){
                IC_Vector[j*CASE.n + i] = ShockedNitro;
            }
        }
        //add the air cavity:
        for(int j = 0; j < CASE.m; ++j){
            for(int i = 0; i < CASE.n; ++i){
                if((pow(i*CASE.dx-CASE.o_x, 2)+pow(j*CASE.dy - CASE.o_y,2)) <= pow(CASE.R,2)){
                    IC_Vector[j*CASE.n + i] = propAir;
                }
            }
        }
    }else if(CASE.construction == "circular"){
        Properties propIn(CASE.WL);
        Properties propOut(CASE.WR);
        //fill with outer state:
        IC_Vector.resize(CASE.n*CASE.m, propOut);
        //add the inner region:
        for(int j = 0; j < CASE.m; ++j){
            for(int i = 0; i < CASE.n; ++i){
                if((pow((i-CASE.nG)*CASE.dx-CASE.o_x, 2)+pow((j-CASE.nG)*CASE.dy - CASE.o_y,2)) <= pow(CASE.R,2)){
                    IC_Vector[j*CASE.n + i] = propIn;
                }
            }
        }
    }else if(CASE.construction == "elliptical"){
        Properties propIn(CASE.WL);
        Properties propOut(CASE.WR);
        //fill with outer state:
        IC_Vector.resize(CASE.n*CASE.m, propOut);
        //add the inner region:
        for(int j = 0; j < CASE.m; ++j){
            for(int i = 0; i < CASE.n; ++i){
                if((pow(i*CASE.dx-CASE.o_x, 2)/pow(0.20, 2)+pow(j*CASE.dy - CASE.o_y,2)/pow(0.12,2)) <= 1.0){
                    IC_Vector[j*CASE.n + i] = propIn;
                }
            }
        }
    }else{
        throw runtime_error("construction type in CASE not valid or not specified");
    }
    
    return IC_Vector;
} 



//Print functions for classes-----------------------------------------------:
void printW(Prim W){
    printf("---------------\n W: \n rho1: %f \n rho2: %f \n u:   %f \n v:   %f \n p:   %f \n z:    %f \n --------------\n",\
            W.rho1, W.rho2, W.u, W.v, W.p, W.z1);
}

void printU(ConsU U){
    printf("---------------\n U: \n z1rho1:  %f \n z2rho2: %f \n rhou: %f \n rhov: %f \nE: %f \n --------------\n"\
            ,U.z1rho1, U.z2rho2, U.rhou, U.rhov, U.E);
}

void printF(ConsF F){
    printf("---------------\n F: \n mass1:%f \n mass2:%f \n xmom: %f \n ymom:%f \n  en:%f \n --------------\n"\
            ,F.mass1, F.mass2, F.xmom, F.ymom, F.en);
}

void printG(ConsG G){
    printf("---------------\n G: \n mass1:%f \n mass2:%f \n xmom: %f \n ymom:%f \n  en:%f \n --------------\n"\
            ,G.mass1, G.mass2, G.xmom, G.ymom, G.en);
}

//--------------------------------------------------------------------------/


//Conversion functions:

void Prop_WtoUF(Properties &Prop){
    //given a properties object, internally update 
    //conservative variables U from primitive
    //variable vector W, by first updating 
    //shared standalone variables:
    double eref1, eref2;
    double pref1, pref2;
    eref1 = e_ref_calc(M1, EoS, Prop.W.rho1); 
    eref2 = e_ref_calc(M2, EoS, Prop.W.rho2);
    pref1 = p_ref_calc(M1, EoS, Prop.W.rho1);
    pref2 = p_ref_calc(M2, EoS, Prop.W.rho2);
    Prop.z2 = 1-Prop.W.z1;
    Prop.rho = Prop.W.z1*Prop.W.rho1 + Prop.z2*Prop.W.rho2;
    Prop.zeta = Prop.W.z1/(varGamma1) + Prop.z2/(varGamma2);
    Prop.eref = (1/Prop.rho)*(Prop.W.z1*Prop.W.rho1*eref1 + Prop.z2*Prop.W.rho2*eref2);
    Prop.pref = (1/Prop.zeta)*(Prop.W.z1*pref1/(varGamma1) + Prop.z2*pref2/(varGamma2));
    Prop.e = Prop.eref + Prop.zeta*(Prop.W.p - Prop.pref)/Prop.rho;
    Prop.E = 0.5*Prop.rho*(pow(Prop.W.u,2)+pow(Prop.W.v,2)) + Prop.rho*Prop.e;
    //Update U:
    Prop.U.z1rho1 = Prop.W.z1*Prop.W.rho1;
    Prop.U.z2rho2 = Prop.z2*Prop.W.rho2;
    Prop.U.rhou = Prop.rho*Prop.W.u;
    Prop.U.rhov = Prop.rho*Prop.W.v;
    Prop.U.E = Prop.E;
    Prop.U.z1 = Prop.W.z1;
    //Update F:
    Prop.F.mass1 = Prop.W.z1*Prop.W.rho1*Prop.W.u;
    Prop.F.mass2 = Prop.z2*Prop.W.rho2*Prop.W.u;
    Prop.F.xmom = Prop.rho*pow(Prop.W.u, 2) + Prop.W.p;
    Prop.F.ymom = Prop.rho*Prop.W.u*Prop.W.v;
    Prop.F.en = Prop.W.u*(Prop.E+Prop.W.p);
    Prop.F.uz1 = Prop.W.u*Prop.W.z1;
    //Update G:
    Prop.G.mass1 = Prop.W.z1*Prop.W.rho1*Prop.W.v;
    Prop.G.mass2 = Prop.z2*Prop.W.rho2*Prop.W.v;
    Prop.G.xmom = Prop.rho*Prop.W.u*Prop.W.v;
    Prop.G.ymom = Prop.rho*pow(Prop.W.v,2) + Prop.W.p; 
    Prop.G.en = Prop.W.v*(Prop.E+Prop.W.p);
    Prop.G.vz1 = Prop.W.v*Prop.W.z1;     
}

void Prop_UtoWF(Properties &Prop){
    //given a properties object, internally update 
    //primitive variable state W from evolved
    //conservative variables U by first updating 
    //shared standalone variables:
    double eref1, eref2;
    double pref1, pref2;
    Prop.z2 = 1-Prop.U.z1; 
    Prop.W.z1 = Prop.U.z1;
    Prop.E = Prop.U.E;
    Prop.W.rho1 = Prop.U.z1rho1/Prop.U.z1;
    Prop.W.rho2 = Prop.U.z2rho2/Prop.z2;
    eref1 = e_ref_calc(M1, EoS, Prop.W.rho1); 
    eref2 = e_ref_calc(M2, EoS, Prop.W.rho2);
    pref1 = p_ref_calc(M1, EoS, Prop.W.rho1);
    pref2 = p_ref_calc(M2, EoS, Prop.W.rho2);
    Prop.rho = Prop.U.z1rho1 + Prop.U.z2rho2;
    Prop.W.u = Prop.U.rhou/Prop.rho;
    Prop.W.v = Prop.U.rhov/Prop.rho;
    Prop.zeta = Prop.W.z1/(varGamma1) + Prop.z2/(varGamma2);
    Prop.eref = (1/Prop.rho)*(Prop.U.z1rho1*eref1 + Prop.U.z2rho2*eref2);
    Prop.pref = (1/Prop.zeta)*(Prop.W.z1*pref1/(varGamma1) + Prop.z2*pref2/(varGamma2));
    Prop.e = (1/Prop.rho)*(Prop.E - 0.5*Prop.rho*(pow(Prop.W.u,2)+pow(Prop.W.v,2)));
    Prop.W.p = Prop.pref + (Prop.rho/Prop.zeta)*(Prop.e - Prop.eref);
    //Subsequently update F:
    Prop.F.mass1 = Prop.W.z1*Prop.W.rho1*Prop.W.u;
    Prop.F.mass2 = Prop.z2*Prop.W.rho2*Prop.W.u;
    Prop.F.xmom = Prop.rho*pow(Prop.W.u, 2) + Prop.W.p;
    Prop.F.ymom = Prop.rho*Prop.W.u*Prop.W.v;
    Prop.F.en = Prop.W.u*(Prop.E+Prop.W.p);
    Prop.F.uz1 = Prop.W.u*Prop.W.z1;
    //Update G:
    Prop.G.mass1 = Prop.W.z1*Prop.W.rho1*Prop.W.v;
    Prop.G.mass2 = Prop.z2*Prop.W.rho2*Prop.W.v;
    Prop.G.xmom = Prop.rho*Prop.W.u*Prop.W.v;
    Prop.G.ymom = Prop.rho*pow(Prop.W.v,2) + Prop.W.p; 
    Prop.G.en = Prop.W.v*(Prop.E+Prop.W.p);
    Prop.G.vz1 = Prop.W.v*Prop.W.z1; 
}

double speed_of_sound(Properties P){
    double c;
    double c1_sqr, c2_sqr;
    if(EoS == "ideal"){
        c1_sqr = (varGamma1+1)*P.W.p/P.W.rho1;
        c2_sqr = (varGamma2+1)*P.W.p/P.W.rho2;
    }
    if(EoS == "stiffened"){
        c1_sqr = (varGamma1+1)*(P.W.p+M1.p_inf)/P.W.rho1;
        c2_sqr = (varGamma2+1)*(P.W.p+M2.p_inf)/P.W.rho2;
    }
    if(EoS == "JWL"){
    c1_sqr = (M1.R1*M1.rho_0/pow(P.W.rho1,2)-(M1.varGamma+1)/P.W.rho1)*M1.A*exp(-M1.R1*M1.rho_0/P.W.rho1)\
            +(M1.R2*M1.rho_0/pow(P.W.rho1,2)-(M1.varGamma+1)/P.W.rho1)*M1.B*exp(-M1.R2*M1.rho_0/P.W.rho1)\
            + (varGamma1+1)*P.W.p/P.W.rho1;
    c2_sqr = (M2.R1*M2.rho_0/pow(P.W.rho2,2)-(M2.varGamma+1)/P.W.rho2)*M2.A*exp(-M2.R1*M2.rho_0/P.W.rho2)\
            +(M2.R2*M2.rho_0/pow(P.W.rho2,2)-(M2.varGamma+1)/P.W.rho2)*M2.B*exp(-M2.R2*M2.rho_0/P.W.rho2)\
            + (varGamma2+1)*P.W.p/P.W.rho2;
    }
    if(EoS == "CC"){
        c1_sqr = (M1.xi1-M1.varGamma-1)*(M1.A/P.W.rho1)*pow(M1.rho_0/P.W.rho1,-M1.xi1) - \
                 (M1.xi2-M1.varGamma-1)*(M1.B/P.W.rho1)*pow(M1.rho_0/P.W.rho1,-M1.xi2) \
                 + (varGamma1+1)*P.W.p/P.W.rho1;
                 
        c2_sqr = (M2.xi1-M2.varGamma-1)*(M2.A/P.W.rho2)*pow(M2.rho_0/P.W.rho2,-M2.xi1) - \
                 (M2.xi2-M2.varGamma-1)*(M2.B/P.W.rho2)*pow(M2.rho_0/P.W.rho2,-M2.xi2) \
                 + (varGamma2+1)*P.W.p/P.W.rho2;
    }
    double dp_ref1, de_ref1;
    double dp_ref2, de_ref2;
    double d_gam1, d_gam2;
    if(EoS == "Hugoniot"){
       if(M1.c_0 == 0){
            dp_ref1 = 0;
            de_ref1 = 0;
        }else{
            dp_ref1 = pow(M1.c_0*M1.rho_0,2)*(M1.s*(P.W.rho1-M1.rho_0)+P.W.rho1)/\
                      pow(P.W.rho1 - M1.s*(P.W.rho1-M1.rho_0),3);
            de_ref1 = 1/(2*pow(P.W.rho1,2))*(p_ref_calc(M1, EoS,P.W.rho1) + M1.p_0) + \
                      (P.W.rho1 - M1.rho_0)/(2*P.W.rho1*M1.rho_0)*dp_ref1;
        }
        if(M2.c_0 == 0){
            dp_ref2 = 0;
            de_ref2 = 0;
        }else{
            dp_ref2 = pow(M2.c_0*M2.rho_0,2)*(M2.s*(P.W.rho2-M2.rho_0)+P.W.rho2)/\
                      pow(P.W.rho2 - M2.s*(P.W.rho2-M2.rho_0),3);
            de_ref2 = 1/(2*pow(P.W.rho2,2))*(p_ref_calc(M2, EoS,P.W.rho2) + M2.p_0) + \
                      (P.W.rho2 - M2.rho_0)/(2*P.W.rho2*M2.rho_0)*dp_ref2;
        }

        c1_sqr = dp_ref1 + ((varGamma1+1)*P.W.p-p_ref_calc(M1, EoS, P.W.rho1))/P.W.rho1 - \
                 varGamma1*P.W.rho1*de_ref1 ; 
        c2_sqr = dp_ref2 + ((varGamma2+1)*P.W.p-p_ref_calc(M2, EoS, P.W.rho2))/P.W.rho2 - \
                 varGamma2*P.W.rho2*de_ref2; 
    }
    
    c = pow(fabs((1/P.zeta)*(P.W.z1*P.W.rho1*c1_sqr/(P.W.rho1*varGamma1) + \
                            (P.z2)*P.W.rho2*c2_sqr/(P.W.rho2*varGamma2))), 0.5);

    //non-physical states are permitted... 
    //...as an intermediate step...
    //... therefore fabs corrects for negative pressure and/or density

    return c;
}

double r_calc(double n_val, double l_val){
    //alternative conventions for dividing by zero RHS slope:
    // -3 encoded to represent -ve infinity
    // +3 encoded to represent +ve infinity 
    double r = l_val;
    if(n_val == 0){
        if(l_val== 0){
            r = 0;           //both left and right slopes are zero
        }else if(l_val < 0){
            r = -3;          //negative slope on left, zero slope on right
        }else{
            r = 3;           //positive slope on left, zero slope on right
        }
    }else{ r = r/n_val;}     //permissible division by delta_i_n

    return r;
}

double delta_bar_calc(string scheme, double r, double delta_i, double sigma_R, double sigma_L){
    double delta_bar;
    if(scheme == "SuperBee"){
        if(r < 0){
            delta_bar = 0;
        }else if(0 <= r && r < 0.5){
            delta_bar = 2*r*delta_i;
        }else if(0.5 <= r && r < 1){
            delta_bar = delta_i;
        }else if(1 <= r && r < 2){
            delta_bar = ((r < sigma_R)? 
                    r*delta_i : sigma_R*delta_i);
        }else{
            delta_bar = ((sigma_R < 2)? 
                    sigma_R*delta_i : 2*delta_i);
        }
    }else if(scheme == "MinBee"){

        if(r <= 0){
            delta_bar = 0;
        }else if(0 <= r && r <= 1.0){
            delta_bar = r*delta_i;
        }else{
            delta_bar = ((1 < sigma_R)? delta_i : sigma_R*delta_i);
        }
    }else if(scheme == "VanLeer"){
        if(r < 0){
            delta_bar = 0;
        }else{
            delta_bar = ((2*r/(1+r) < sigma_R)? 
                    (2*r/(1+r))*delta_i : sigma_R*delta_i);
        }
    }else if(scheme == "UltraBee"){
        if(r <= 0){
            delta_bar = 0;
        }else{
            delta_bar = ((sigma_L < sigma_R)? 
                    sigma_L*delta_i : sigma_R*delta_i);
        }
    }else{
        return delta_i;
    }

    return delta_bar;
}



Prim delta(double omega, Prim Wl, Prim Wi, Prim Wn, string scheme){
    //omega = [-1, 1]
    //Wl = W_{i-1}, Wi = Ui, Wn = W_{i+1}
    Prim delta_i_l(Wi.rho1 - Wl.rho1, Wi.rho2 - Wi.rho2, Wi.u - Wl.u, Wi.v - Wl.v, 
                    Wi.p - Wl.p, Wi.z1 - Wl.z1);
    Prim delta_i_n(Wn.rho1 - Wi.rho1, Wn.rho2 - Wi.rho2, Wn.u - Wi.u, Wn.v - Wi.v,
                    Wn.p - Wi.p, Wn.z1 - Wi.z1);
    
    
    double r_rho1 = r_calc(delta_i_n.rho1, delta_i_l.rho1);
    double r_rho2 = r_calc(delta_i_n.rho2, delta_i_l.rho2);
    double r_u = r_calc(delta_i_n.u, delta_i_l.u);
    double r_v = r_calc(delta_i_n.v, delta_i_l.v);
    double r_p = r_calc(delta_i_n.p, delta_i_l.p);
    double r_z1 = r_calc(delta_i_n.z1, delta_i_l.z1);

    Prim delta_i(0.5*(1+omega)*delta_i_l.rho1+0.5*(1-omega)*delta_i_n.rho1,
                  0.5*(1+omega)*delta_i_l.rho2+0.5*(1-omega)*delta_i_n.rho2,
                  0.5*(1+omega)*delta_i_l.u+0.5*(1-omega)*delta_i_n.u,
                  0.5*(1+omega)*delta_i_l.v+0.5*(1-omega)*delta_i_n.v,
                  0.5*(1+omega)*delta_i_l.p+0.5*(1-omega)*delta_i_n.p,
                  0.5*(1+omega)*delta_i_l.z1+0.5*(1-omega)*delta_i_n.z1);

    Prim sigma_R(2/(1-omega+(1+omega)*r_rho1), 2/(1-omega+(1+omega)*r_rho2), 2/(1-omega+(1+omega)*r_u), 
                    2/(1-omega+(1+omega)*r_v), 2/(1-omega+(1+omega)*r_p), 2/(1-omega+(1+omega)*r_z1));

    Prim sigma_L(2*r_rho1/(1-omega+(1+omega)*r_rho1), 2*r_rho1/(1-omega+(1+omega)*r_rho1), 
                    2*r_u/(1-omega+(1+omega)*r_u), 2*r_v/(1-omega+(1+omega)*r_v),
                     2*r_p/(1-omega+(1+omega)*r_p), 2*r_z1/(1-omega+(1+omega)*r_z1));

    //adjustments for infinity encoded cases:
    if(r_rho1 == -3 || r_rho1 == 3){
        sigma_R.rho1 = 0;
    }
    if(r_rho2 == -3 || r_rho2 == 3){
        sigma_R.rho2 = 0;
    }
    if(r_u == -3 || r_u == 3){
        sigma_R.u = 0;
    }
    if(r_v == -3 || r_v == 3){
        sigma_R.v = 0;
    }
    if(r_p == -3 || r_p == 3){
        sigma_R.p = 0;
    }
    if(r_z1 == -3 || r_z1 == 3){
        sigma_R.z1 = 0;
    }
    
    Prim delta_bar(0,0,0,0,0,0);

    delta_bar.rho1 = delta_bar_calc(scheme, r_rho1, delta_i.rho1, sigma_R.rho1, sigma_L.rho1);
    delta_bar.rho2 = delta_bar_calc(scheme, r_rho2, delta_i.rho2, sigma_R.rho2, sigma_L.rho2);
    delta_bar.u = delta_bar_calc(scheme, r_u, delta_i.u, sigma_R.u, sigma_L.u);
    delta_bar.v = delta_bar_calc(scheme, r_v, delta_i.v, sigma_R.v, sigma_L.v);
    delta_bar.p = delta_bar_calc(scheme, r_p, delta_i.p, sigma_R.p, sigma_L.p);
    delta_bar.z1 = delta_bar_calc(scheme, r_z1, delta_i.z1, sigma_R.z1, sigma_L.z1);

    return delta_bar;

}




class HLLC_x{
    public:
        double aL, aR;                  //left and right sound speeds
        double SL, SR, Sstar, Splus;    //wave speeds
        ConsU UL, UR, ULstar, URstar;   //4 cons. states
        ConsF FL, FR, FLstar, FRstar;   //4 cons. fluxes
        ConsU U;                        //tbd U state
        ConsF Fout;                     //flux out (F_i+1/2)
        double wL, wR;                       //boundary weighting factor
        
        //must declare constructor here:
        HLLC_x(Properties PL, Properties PR):
            //must initialize all class objects here:
            UL(0,0,0,0,0,0), UR(0,0,0,0,0,0), ULstar(0,0,0,0,0,0), URstar(0,0,0,0,0,0),
            FL(0,0,0,0,0,0), FR(0,0,0,0,0,0), FLstar(0,0,0,0,0,0), FRstar(0,0,0,0,0,0),
            U(0,0,0,0,0,0), Fout(0,0,0,0,0,0) {
            //direct wave speed estimates:
            //sound speed calcs:
            aL = speed_of_sound(PL);
            aR = speed_of_sound(PR);
            double k = 0.5*(PL.k+PR.k);
            double sign = 1.0; //((PL.n_x > 0)? 1.0 : -1.0);
            wL = PL.w;
            wR = PR.w;
            //wave speed calcs:    
            SL = ((PL.W.u - aL < PR.W.u - aR)? PL.W.u - aL : PR.W.u - aR);
            SR = ((PL.W.u + aL > PR.W.u + aR)? PL.W.u + aL : PR.W.u + aR);
            Sstar = (PR.W.p - PL.W.p + PL.rho*PL.W.u*(SL-PL.W.u) - PR.rho*PR.W.u*(SR-PR.W.u)\
                    - sigma*k*sign*(PR.w-PL.w))/(PL.rho*(SL-PL.W.u)-PR.rho*(SR-PR.W.u));
            Splus = ((fabs(PL.W.u)+aL > fabs(PR.W.u)+aR)? fabs(PL.W.u)+aL : fabs(PR.W.u)+aR);
            //W-->U and W-->F conversions:
            UL = PL.U;
            UR = PR.U;
            FL = PL.F;
            FR = PR.F;
            //Star region calculations:
            ULstar.z1rho1 = PL.W.z1*PL.W.rho1*((SL-PL.W.u)/(SL-Sstar));
            ULstar.z2rho2 = (1-PL.W.z1)*PL.W.rho2*((SL-PL.W.u)/(SL-Sstar));
            ULstar.rhou = PL.rho*Sstar*((SL-PL.W.u)/(SL-Sstar));
            ULstar.rhov = PL.rho*PL.W.v*((SL-PL.W.u)/(SL-Sstar)); 
            ULstar.E = PL.rho*((SL-PL.W.u)/(SL-Sstar))*(UL.E/PL.rho+(Sstar-PL.W.u)*(Sstar+(PL.W.p - sigma*k*sign*PL.w)\
                        /(PL.rho*(SL-PL.W.u))));
            ULstar.z1 = PL.W.z1;
            //assert(ULstar.E >=0); //check that total energy is non-negative
            URstar.z1rho1 = PR.W.z1*PR.W.rho1*((SR-PR.W.u)/(SR-Sstar));
            URstar.z2rho2 = (1-PR.W.z1)*PR.W.rho2*((SR-PR.W.u)/(SR-Sstar));
            URstar.rhou = PR.rho*Sstar*((SR-PR.W.u)/(SR-Sstar));
            URstar.rhov = PR.rho*PR.W.v*((SR-PR.W.u)/(SR-Sstar)); 
            URstar.E = PR.rho*((SR-PR.W.u)/(SR-Sstar))*(UR.E/PR.rho+(Sstar-PR.W.u)*(Sstar+(PR.W.p - sigma*k*sign*PR.w)\
                        /(PR.rho*(SR-PR.W.u))));
            URstar.z1 = PR.W.z1;
            //assert(URstar.E >=0); //check that total energy is non-negative
            FLstar.mass1 = FL.mass1 + SL*(ULstar.z1rho1 - UL.z1rho1);
            FLstar.mass2 = FL.mass2 + SL*(ULstar.z2rho2 - UL.z2rho2);
            FLstar.xmom = FL.xmom + SL*(ULstar.rhou - UL.rhou);
            FLstar.ymom = FL.ymom + SL*(ULstar.rhov - UL.rhov); //CHECK
            FLstar.en = FL.en + SL*(ULstar.E - UL.E);
            FLstar.uz1 = FLstar.mass1/PL.W.rho1;

            FRstar.mass1 = FR.mass1 + SR*(URstar.z1rho1 - UR.z1rho1);
            FRstar.mass2 = FR.mass2 + SR*(URstar.z2rho2 - UR.z2rho2);
            FRstar.xmom = FR.xmom + SR*(URstar.rhou - UR.rhou);
            FRstar.ymom = FR.ymom + SR*(URstar.rhov - UR.rhov); //CHECK
            FRstar.en = FR.en + SR*(URstar.E - UR.E);
            FRstar.uz1 = FRstar.mass1/PR.W.rho1; 
            
            //state and flux Godunov-type evaulation:
            if(0 < SL){
                U = UL;
                Fout = FL;
            }else if(SL <= 0 && 0 < Sstar){
                U = ULstar;
                Fout = FLstar;
            }else if(Sstar <= 0 && 0 < SR){
                U = URstar;
                Fout = FRstar;
            }else if(SR <= 0){
                U = UR;
                Fout = FR;
            }
        }
};

class HLLC_y{
    public:
        double aL, aR;                  //left and right sound speeds
        double SL, SR, Sstar, Splus;    //wave speeds
        ConsU UL, UR, ULstar, URstar;   //4 cons. states
        ConsG GL, GR, GLstar, GRstar;   //4 cons. fluxes
        ConsU U;                        //tbd U state
        ConsG Gout;                     //flux out (F_i+1/2)
        double wL, wR;                       //boundary weighting factor
        
        //must declare constructor here:
        HLLC_y(Properties PL, Properties PR):
            //must initialize all class objects here:
            UL(0,0,0,0,0,0), UR(0,0,0,0,0,0), ULstar(0,0,0,0,0,0), URstar(0,0,0,0,0,0),
            GL(0,0,0,0,0,0), GR(0,0,0,0,0,0), GLstar(0,0,0,0,0,0), GRstar(0,0,0,0,0,0),
            U(0,0,0,0,0,0), Gout(0,0,0,0,0,0) {
            //direct wave speed estimates:
            //sound speed calcs:
            aL = speed_of_sound(PL);
            aR = speed_of_sound(PR);
            double k = 0.5*(PL.k+PR.k);
            double sign = 1.0; //((PL.n_y > 0)? 1.0 : -1.0);
            wL = PL.w;
            wR = PR.w;
            //wave speed calcs:    
            SL = ((PL.W.v - aL < PR.W.v - aR)? PL.W.v - aL : PR.W.v - aR);
            SR = ((PL.W.v + aL > PR.W.v + aR)? PL.W.v + aL : PR.W.v + aR);
            Sstar = (PR.W.p - PL.W.p + PL.rho*PL.W.v*(SL-PL.W.v) - PR.rho*PR.W.v*(SR-PR.W.v) \
                    - sigma*k*sign*(PR.w-PL.w))/(PL.rho*(SL-PL.W.v)-PR.rho*(SR-PR.W.v));
            Splus = ((fabs(PL.W.v)+aL > fabs(PR.W.v)+aR)? fabs(PL.W.v)+aL : fabs(PR.W.v)+aR);           
            //W-->U and W-->F conversions:
            UL = PL.U;
            UR = PR.U;
            GL = PL.G;
            GR = PR.G;
            //Star region calculations:
            ULstar.z1rho1 = PL.W.z1*PL.W.rho1*((SL-PL.W.v)/(SL-Sstar));
            ULstar.z2rho2 = (1-PL.W.z1)*PL.W.rho2*((SL-PL.W.v)/(SL-Sstar));
            ULstar.rhou = PL.rho*PL.W.u*((SL-PL.W.v)/(SL-Sstar));
            ULstar.rhov = PL.rho*Sstar*((SL-PL.W.v)/(SL-Sstar)); 
            ULstar.E = PL.rho*((SL-PL.W.v)/(SL-Sstar))*(UL.E/PL.rho+(Sstar-PL.W.v)*(Sstar+(PL.W.p-sigma*k*sign*PL.w)\
                        /(PL.rho*(SL-PL.W.v))));
            ULstar.z1 = PL.W.z1;
            //assert(ULstar.E >=0); //check that total energy is non-negative
            URstar.z1rho1 = PR.W.z1*PR.W.rho1*((SR-PR.W.v)/(SR-Sstar));
            URstar.z2rho2 = (1-PR.W.z1)*PR.W.rho2*((SR-PR.W.v)/(SR-Sstar));
            URstar.rhou = PR.rho*PR.W.u*((SR-PR.W.v)/(SR-Sstar));
            URstar.rhov = PR.rho*Sstar*((SR-PR.W.v)/(SR-Sstar)); 
            URstar.E = PR.rho*((SR-PR.W.v)/(SR-Sstar))*(UR.E/PR.rho+(Sstar-PR.W.v)*(Sstar+(PR.W.p - sigma*k*sign*PR.w)\
                        /(PR.rho*(SR-PR.W.v))));
            URstar.z1 = PR.W.z1;
            //assert(URstar.E >=0); //check that total energy is non-negative
            GLstar.mass1 = GL.mass1 + SL*(ULstar.z1rho1 - UL.z1rho1);
            GLstar.mass2 = GL.mass2 + SL*(ULstar.z2rho2 - UL.z2rho2);
            GLstar.xmom = GL.xmom + SL*(ULstar.rhou - UL.rhou);
            GLstar.ymom = GL.ymom + SL*(ULstar.rhov - UL.rhov); //CHECK
            GLstar.en = GL.en + SL*(ULstar.E - UL.E);
            GLstar.vz1 = GLstar.mass1/PL.W.rho1;

            GRstar.mass1 = GR.mass1 + SR*(URstar.z1rho1 - UR.z1rho1);
            GRstar.mass2 = GR.mass2 + SR*(URstar.z2rho2 - UR.z2rho2);
            GRstar.xmom = GR.xmom + SR*(URstar.rhou - UR.rhou);
            GRstar.ymom = GR.ymom + SR*(URstar.rhov - UR.rhov); //CHECK
            GRstar.en = GR.en + SR*(URstar.E - UR.E);
            GRstar.vz1 = GRstar.mass1/PR.W.rho1; 
            
            //state and flux Godunov-type evaulation:
            if(0 < SL){
                U = UL;
                Gout = GL;
            }else if(SL <= 0 && 0 < Sstar){
                U = ULstar;
                Gout = GLstar;
            }else if(Sstar <= 0 && 0 < SR){
                U = URstar;
                Gout = GRstar;
            }else if(SR <= 0){
                U = UR;
                Gout = GR;
            }
        }
};


void updateBoundary(vector<Properties> &P_vec, int n, int m, int nG, string condition){
    assert(n > 2*nG);
    assert(m > 2*nG);
    
    //transmissive component of all BCs:
    for(int g = 0; g < nG; ++g){
        for(int j = 0; j < m; ++j){
            for(int i = 0; i < n; ++i){
                P_vec[g*n+i] = P_vec[nG*n+i];               //bottom
                P_vec[(m-1-g)*n+i] = P_vec[(m-1-nG)*n+i];   //top
                P_vec[j*n+g] = P_vec[j*n+nG];               //left
                P_vec[j*n+(n-1-g)] = P_vec[j*n+(n-1-nG)];   //right
            }
        }
    }

    if(condition == "transmissive"){
        //as above
    }else if(condition == "bottom-reflective"){
        for(int g = 0; g < nG; ++g){
            for(int j = 0; j < m; ++j){
                for(int i = 0; i < n; ++i){
                    P_vec[g*n+i] = P_vec[nG*n+i];               //bottom
                    //reflected y-velocity:
                    //P_vec[(nG-1-g)*n+i] = P_vec[(nG+g)*n+i];
                    P_vec[(nG-1-g)*n+i].W.v = -P_vec[(nG+g)*n+i].W.v;
                    P_vec[(nG-1-g)*n+i].W.u = -P_vec[(nG+g)*n+i].W.u;
                    Prop_WtoUF(P_vec[(nG-1-g)*n+i]);  
                    P_vec[(m-1-g)*n+i] = P_vec[(m-1-nG)*n+i];   //top
                    P_vec[j*n+g] = P_vec[j*n+nG];               //left
                    P_vec[j*n+(n-1-g)] = P_vec[j*n+(n-1-nG)];   //right
                }
            }
        }   
    }else if(condition == "all-reflective"){
        //this BC seems to be causing errors
        for(int g = 0; g < nG; ++g){
            for(int j = 0; j < m; ++j){
                for(int i = 0; i < n; ++i){
                    //bottom surface
                    P_vec[(nG-1-g)*n+i].W.v = -P_vec[(nG+g)*n+i].W.v;
                    P_vec[(nG-1-g)*n+i].W.u = -P_vec[(nG+g)*n+i].W.u;
                    Prop_WtoUF(P_vec[(nG-1-g)*n+i]); 
                    //top surface
                    P_vec[(m-nG+g)*n+i].W.v = -P_vec[(m-nG-1-g)*n+i].W.v;
                    P_vec[(m-nG+g)*n+i].W.u = -P_vec[(m-nG-1-g)*n+i].W.u;
                    Prop_WtoUF(P_vec[(m-nG+g)*n+i]);
                    //left wall:
                    P_vec[j*n + nG-1-g].W.u = -P_vec[j*n + nG+g].W.u;
                    P_vec[j*n + nG-1-g].W.v = -P_vec[j*n + nG+g].W.v;
                    Prop_WtoUF(P_vec[j*n + nG-1-g]);
                    //right wall:
                    P_vec[j*n + m-nG+g].W.u = -P_vec[j*n + m-nG-1-g].W.u;  
                    P_vec[j*n + m-nG+g].W.v = -P_vec[j*n + m-nG-1-g].W.v;
                    Prop_WtoUF(P_vec[j*n + m-nG+g]);  
                }
            }
        }   
    }else{
        cout << "UNDEFINED BOUNDARY CONDITION" << endl;
    } 

}

void reflective_boundaries(vector<Properties> &P_vec, int n, int m, int nG, string condition){
    assert(n > 2*nG);
    assert(m > 2*nG);
    if(condition == "bottom-reflective"){
        for(int g = 0; g < nG; ++g){
            for(int j = 0; j < m; ++j){
                for(int i = 0; i < n; ++i){
                    //reflected y-velocity:
                    P_vec[(nG-1-g)*n+i].W.v = -P_vec[(nG+g)*n+i].W.v;
                    P_vec[(nG-1-g)*n+i].W.u = -P_vec[(nG+g)*n+i].W.u;
                    Prop_WtoUF(P_vec[(nG-1-g)*n+i]);  
                }
            }
        }   
    }
    if(condition == "all-reflective"){
        for(int g = 0; g < nG; ++g){
            for(int j = 0; j < m; ++j){
                for(int i = 0; i < n; ++i){
                    //bottom surface
                    P_vec[(nG-1-g)*n+i].W.v = -P_vec[(nG+g)*n+i].W.v;
                    P_vec[(nG-1-g)*n+i].W.u = -P_vec[(nG+g)*n+i].W.u;
                    Prop_WtoUF(P_vec[(nG-1-g)*n+i]); 
                    //top surface
                    P_vec[(m-nG+g)*n+i].W.v = -P_vec[(m-nG-1-g)*n+i].W.v;
                    P_vec[(m-nG+g)*n+i].W.u = -P_vec[(m-nG-1-g)*n+i].W.u;
                    Prop_WtoUF(P_vec[(m-nG+g)*n+i]);
                    //left wall:
                    P_vec[j*n + nG-1-g].W.u = -P_vec[j*n + nG+g].W.u;
                    P_vec[j*n + nG-1-g].W.v = -P_vec[j*n + nG+g].W.v;
                    Prop_WtoUF(P_vec[j*n + nG-1-g]);
                    //right wall:
                    P_vec[j*n + m-nG+g].W.u = -P_vec[j*n + m-nG-1-g].W.u;  
                    P_vec[j*n + m-nG+g].W.v = -P_vec[j*n + m-nG-1-g].W.v;
                    Prop_WtoUF(P_vec[j*n + m-nG+g]);  
                }
            }
        }   
    }

}
           
//SURFACE TENSION:

void curvature_calc(vector<Properties> &P_vec, TestCase CASE, double &k_avr){
    int nG = CASE.nG;
    int m = CASE.m;
    int n = CASE.n;
    double dx = CASE.dx;
    double dy = CASE.dy;
    double tol = 0.01; //diffusion tolerence
    double z1_avr;
    double z1_sum;

    for(int j = nG; j < m - nG; ++j){
        for(int i = nG; i < n - nG; ++i){

            P_vec[j*n+i].bk = 0;
            
            //if cell lies in diffuse region or borders z1=0.5 level, then 9-point avr should lie between 0-1:
            z1_avr =P_vec[(j+1)*n+i-1].W.z1 + P_vec[(j+1)*n+i].W.z1 + P_vec[(j+1)*n+i+1].W.z1 + \
                    P_vec[j*n+i-1].W.z1     + P_vec[j*n+i].W.z1     + P_vec[j*n+i+1].W.z1 +\
                    P_vec[(j-1)*n+i-1].W.z1 + P_vec[(j-1)*n+i].W.z1 + P_vec[(j-1)*n+i+1].W.z1;
            z1_avr = z1_avr/9.0;
            P_vec[j*n+i].z1_avr = z1_avr; 
            
            if(tol < z1_avr && z1_avr < 1-tol){
                P_vec[j*n+i].bk = 1;
            }
                                                
        }
    }
    
    double nx, ny;

    for(int j = nG; j < m - nG; ++j){
        for(int i = nG; i < n - nG; ++i){
            if(P_vec[j*n+i].bk == 1){
            //calculate normals and weighting factor:
            //fourth order central differences:
               nx = -(1.0/(12*dx))*(-(P_vec[(j+1)*n+i+2].z1_avr + 2*P_vec[(j)*n+i+2].z1_avr + P_vec[(j-1)*n+i+2].z1_avr) +\
                            +  8*(P_vec[(j+1)*n+i+1].z1_avr + 2*P_vec[(j)*n+i+1].z1_avr + P_vec[(j-1)*n+i+1].z1_avr) \
                            -  8*(P_vec[(j+1)*n+i-1].z1_avr + 2*P_vec[(j)*n+i-1].z1_avr + P_vec[(j-1)*n+i-1].z1_avr) \
                            +    (P_vec[(j+1)*n+i-2].z1_avr + 2*P_vec[(j)*n+i-2].z1_avr + P_vec[(j-1)*n+i-2].z1_avr));
                
               ny = -(1.0/(12*dy))*(-(P_vec[(j+2)*n+i+1].z1_avr + 2*P_vec[(j+2)*n+i].z1_avr + P_vec[(j+2)*n+i-1].z1_avr)
                            + 8*(P_vec[(j+1)*n+i+1].z1_avr + 2*P_vec[(j+1)*n+i].z1_avr + P_vec[(j+1)*n+i-1].z1_avr) \
                            - 8*(P_vec[(j-1)*n+i+1].z1_avr + 2*P_vec[(j-1)*n+i].z1_avr + P_vec[(j-1)*n+i-1].z1_avr)
                            +   (P_vec[(j-2)*n+i+1].z1_avr + 2*P_vec[(j-2)*n+i].z1_avr + P_vec[(j-2)*n+i-1].z1_avr));

                P_vec[j*n+i].n_x = nx/pow(pow(nx,2) + pow(ny,2), 0.5);
                P_vec[j*n+i].n_y = ny/pow(pow(nx,2) + pow(ny,2), 0.5); 

                P_vec[j*n+i].w = -1*(6*pow(P_vec[j*n+i].z1_avr,5)-15*pow(P_vec[j*n+i].z1_avr,4)+10*pow(P_vec[j*n+i].z1_avr,3));
            }
        }
    }

    double nR0x = (int)(0.25*(1.0/k_avr)/dx); 
    double nR0y = (int)(0.25*(1.0/k_avr)/dy);
    //cout << "nR0x = " << nR0x << endl;
    double local_nx_sum, local_ny_sum;
    int counter;
    int k_counter = 0;
    double k_sum = 0;
    double sum_w, w, z_hat; 
    double gs = 4.0; //gaussian function steepness
    
    
    for(int j = nR0y; j < m - nR0y; ++j){
        for(int i = nR0x; i < n - nR0x; ++i){
            if(P_vec[j*n+i].bk == 1 ){
                counter = 0;
                sum_w = 0;
                local_nx_sum = 0;
                local_ny_sum = 0;
                for(int Rx = -nR0x; Rx <= nR0x; ++Rx){
                    for(int Ry = -nR0y; Ry <= nR0y; ++Ry){
                        if(P_vec[(j+Ry)*n+i+Rx].bk == 1){
                            ++counter;
                            //counter must at least be 1
                            if(P_vec[(j+Ry)*n+i+Rx].n_x*P_vec[(j)*n+i].n_x >= 0 || \
                                P_vec[(j+Ry)*n+i+Rx].n_y*P_vec[(j)*n+i].n_y >= 0){
                                //check boundary overlap zones are not included
                                z_hat = P_vec[(j+Ry)*n+i+Rx].z1_avr - P_vec[j*n+i].z1_avr;
                                w = exp(-gs*pow(z_hat,2));
                                sum_w += w;
                                local_nx_sum += w*P_vec[(j+Ry)*n+i+Rx].n_x;
                                local_ny_sum += w*P_vec[(j+Ry)*n+i+Rx].n_y;
                            }
                        }
                    }
                }
                P_vec[j*n+i].n_x_avr = local_nx_sum/sum_w;
                P_vec[j*n+i].n_y_avr = local_ny_sum/sum_w;                
            }
        }
    }
    

    for(int j = nG; j < m - nG; ++j){
        for(int i = nG; i < n - nG; ++i){
            if(P_vec[j*n+i].bk == 1){
            //calculate normals and weighting factor:
                if(P_vec[j*n+i-1].bk == 0){
                    P_vec[j*n+i].n_xx = (1.0/(dx))*(P_vec[j*n+i+1].n_x_avr - P_vec[j*n+i].n_x_avr);
                    P_vec[j*n+i].n_yx = (1.0/(dx))*(P_vec[j*n+i+1].n_y_avr - P_vec[j*n+i].n_y_avr);
                }else if(P_vec[j*n+i+1].bk == 0){
                    P_vec[j*n+i].n_xx = (1.0/(dx))*(P_vec[j*n+i].n_x_avr - P_vec[j*n+i-1].n_x_avr);
                    P_vec[j*n+i].n_yx = (1.0/(dx))*(P_vec[j*n+i].n_y_avr - P_vec[j*n+i-1].n_y_avr);
                }else{
                    P_vec[j*n+i].n_xx = (1.0/(2*dx))*(P_vec[j*n+i+1].n_x_avr - P_vec[j*n+i-1].n_x_avr);
                    P_vec[j*n+i].n_yx = (1.0/(2*dx))*(P_vec[j*n+i+1].n_y_avr - P_vec[j*n+i-1].n_y_avr);
                }
                if(P_vec[(j-1)*n+i].bk == 0){
                    P_vec[j*n+i].n_yy =  (1.0/(dy))*(P_vec[(j+1)*n+i].n_y_avr - P_vec[j*n+i].n_y_avr);
                    P_vec[j*n+i].n_xy =  (1.0/(dy))*(P_vec[(j+1)*n+i].n_x_avr - P_vec[j*n+i].n_x_avr);
                }else if(P_vec[(j+1)*n+i].bk == 0){
                    P_vec[j*n+i].n_yy =  (1.0/(dy))*(P_vec[j*n+i].n_y_avr - P_vec[(j-1)*n+i].n_y_avr);
                    P_vec[j*n+i].n_xy =  (1.0/(dy))*(P_vec[j*n+i].n_x_avr - P_vec[(j-1)*n+i].n_x_avr);
                }else{
                    P_vec[j*n+i].n_yy =  (1.0/(2*dy))*(P_vec[(j+1)*n+i].n_y_avr - P_vec[(j-1)*n+i].n_y_avr);
                    P_vec[j*n+i].n_xy =  (1.0/(2*dy))*(P_vec[(j+1)*n+i].n_x_avr - P_vec[(j-1)*n+i].n_x_avr);
                }
                
                P_vec[j*n+i].k = (P_vec[j*n+i].n_y_avr*(P_vec[j*n+i].n_y_avr*P_vec[j*n+i].n_xx - \
                                    P_vec[j*n+i].n_x_avr*P_vec[j*n+i].n_xy) - \
                                    P_vec[j*n+i].n_x_avr*(P_vec[j*n+i].n_y_avr*P_vec[j*n+i].n_yx - \
                                        P_vec[j*n+i].n_x_avr*P_vec[j*n+i].n_yy)) \
                                  /pow(pow(P_vec[j*n+i].n_x_avr,2) + pow(P_vec[j*n+i].n_y_avr,2),1.5);
                 
                ++k_counter;
                k_sum += fabs(P_vec[j*n+i].k);
                
            }else{
                 P_vec[j*n+i].k = 0;
            }
        }
    }
    
    k_avr = k_sum/k_counter;

    
    //correction/weighting iterations:

    double num_its = 1;
    double local_k_sum;
    //nR0x = 3;
    //nR0y = 3;

    for(int iter = 0; iter < num_its; ++iter){
        for(int j = nR0y; j < m - nR0y; ++j){
            for(int i = nR0x; i < n - nR0x; ++i){
                if(P_vec[j*n+i].bk == 1 ){
                    sum_w = 0;
                    local_k_sum = 0;
                    for(int Rx = -nR0x; Rx <= nR0x; ++Rx){
                        for(int Ry = -nR0y; Ry <= nR0y; ++Ry){
                            if(P_vec[(j+Ry)*n+i+Rx].bk == 1){
                                if(P_vec[(j+Ry)*n+i+Rx].n_x*P_vec[(j)*n+i].n_x >= 0 || \
                                P_vec[(j+Ry)*n+i+Rx].n_y*P_vec[(j)*n+i].n_y >= 0){
                                //z_hat = P_vec[(j+Ry)*n+i+Rx].z1_avr - 0.5;
                                w = pow(P_vec[(j+Ry)*n+i+Rx].z1_avr*(1-P_vec[(j+Ry)*n+i+Rx].z1_avr),2); //exp(-gs*pow(z_hat,2));
                                sum_w += w;
                                local_k_sum += w*P_vec[(j+Ry)*n+i+Rx].k;
                                }
                            }
                        }
                    }
                    P_vec[j*n+i].k_temp = local_k_sum/sum_w;               
                }
            }
        }
        for(int j = nR0y; j < m - nR0y; ++j){
            for(int i = nR0x; i < n - nR0x; ++i){
                if(P_vec[j*n+i].bk == 1 ){
                    P_vec[j*n+i].k = P_vec[j*n+i].k_temp;
                }
            }
        }

    }
    
    
} 


vector<Properties> updateP(vector<Properties> &P_vec, TestCase CASE, double &dt, double CFL, 
        double &Smax, bool MUSCL, string scheme){

    //MUSCL-Hancock parameters:
    Prim WL(0,0,0,0,0,0);            
    Prim Wl(0,0,0,0,0,0);
    Prim Wi(0,0,0,0,0,0);
    Prim Wn(0,0,0,0,0,0);
    Prim Wnn(0,0,0,0,0,0);
    Prim WR(0,0,0,0,0,0);
    Prim WbarL(0,0,0,0,0,0);
    Prim WbarR(0,0,0,0,0,0);
    Properties PbarL = P_vec[0];
    Properties PbarR = P_vec[0];   
    double omega = 0;
    Prim delta_i(0,0,0,0,0,0);
    Prim delta_n(0,0,0,0,0,0);
    double cL, cR;
    double rho_i;
    double rho_n;
    vector<Properties> PbarL_vec_x = P_vec;
    vector<Properties> PbarR_vec_x = P_vec;
    vector<Properties> PbarL_vec_y = P_vec;
    vector<Properties> PbarR_vec_y = P_vec;
    int n = CASE.n;
    int m = CASE.m;
    int nG = CASE.nG;
    double dx = CASE.dx;
    double dy = CASE.dy;
    string BC = CASE.BCs;

    //update parameters 
    vector<Properties> P_vec_new = P_vec; //initialised variable for subsequent updating  
    HLLC_x cell0_x(P_vec[0], P_vec[0]);
    HLLC_y cell0_y(P_vec[0], P_vec[0]);
    vector<HLLC_x> cell_vec_x(n*m, cell0_x);
    vector<HLLC_y> cell_vec_y(n*m, cell0_y);

    
    for(int j = nG; j<(m-nG); ++j){    
        for(int i = nG; i<(n-nG); ++i){
            if(MUSCL == 1){
                //MUSCL-Hancock steps 1+2:
                
                //X-sweep:
                Wl = P_vec[j*n+i-1].W; 
                Wi = P_vec[j*n+i].W;
                Wn = P_vec[j*n+i+1].W; 
                Wnn = P_vec[j*n+i+2].W; 
                
                //to copy over other parameters including ST terms
                PbarL = P_vec[j*n+i];
                PbarR = P_vec[j*n+i+1];             
                       
                delta_i = delta(omega, Wl, Wi, Wn, scheme);
                delta_n = delta(omega, Wi, Wn, Wnn, scheme);
                
                
                //combined MUSCL-Hancock steps 1+2:
                //Left state:
                cL = speed_of_sound(P_vec[j*n+i]);
                rho_i = Wi.z1*Wi.rho1 + (1-Wi.z1)*Wi.rho2;
                PbarL.W.rho1 = Wi.rho1 + 0.5*(1-(dt/dx)*Wi.u)*delta_i.rho1 - 0.5*(dt/dx)*Wi.rho1*delta_i.u;            
                PbarL.W.rho2 = Wi.rho2 + 0.5*(1-(dt/dx)*Wi.u)*delta_i.rho2 - 0.5*(dt/dx)*Wi.rho2*delta_i.u;
                PbarL.W.u = Wi.u + 0.5*(1-(dt/dx)*Wi.u)*delta_i.u - 0.5*(dt/dx)*(1/rho_i)*delta_i.p;
                PbarL.W.v = Wi.v + 0.5*(1-(dt/dx)*Wi.u)*delta_i.v;
                PbarL.W.p = Wi.p + 0.5*(1-(dt/dx)*Wi.u)*delta_i.p - 0.5*(dt/dx)*rho_i*pow(cL,2)*delta_i.u;
                PbarL.W.z1 = Wi.z1 + 0.5*(1-(dt/dx)*Wi.u)*delta_i.z1;
                Prop_WtoUF(PbarL);

                //Right state:
                cR = speed_of_sound(P_vec[j*n+i+1]);
                rho_n = Wn.z1*Wn.rho1 + (1-Wn.z1)*Wn.rho2;
                PbarR.W.rho1 = Wn.rho1 - 0.5*(1+(dt/dx)*Wn.u)*delta_n.rho1 - 0.5*(dt/dx)*Wn.rho1*delta_n.u;
                PbarR.W.rho2 = Wn.rho2 - 0.5*(1+(dt/dx)*Wn.u)*delta_n.rho2 - 0.5*(dt/dx)*Wn.rho2*delta_n.u;
                PbarR.W.u = Wn.u - 0.5*(1+(dt/dx)*Wn.u)*delta_n.u - 0.5*(dt/dx)*(1/rho_n)*delta_n.p;
                PbarR.W.v = Wn.v - 0.5*(1+(dt/dx)*Wn.u)*delta_n.v;
                PbarR.W.p = Wn.p - 0.5*(1+(dt/dx)*Wn.u)*delta_n.p - 0.5*(dt/dx)*rho_n*pow(cR,2)*delta_n.u;
                PbarR.W.z1 = Wn.z1 - 0.5*(1+(dt/dx)*Wn.u)*delta_n.z1;
                Prop_WtoUF(PbarR);
           
                PbarL_vec_x[j*n+i] = PbarL;
                PbarR_vec_x[j*n+i] = PbarR;
                
                
                //Y-sweep:
                Wl = P_vec[(j-1)*n+i].W; 
                Wi = P_vec[j*n+i].W;
                Wn = P_vec[(j+1)*n+i].W; 
                Wnn = P_vec[(j+2)*n+i].W;

                //to copy over other parameters including ST terms
                PbarL = P_vec[j*n+i];
                PbarR = P_vec[(j+1)*n+i]; 
                
                       
                delta_i = delta(omega, Wl, Wi, Wn, scheme);
                delta_n = delta(omega, Wi, Wn, Wnn, scheme);
                
                
                //combined MUSCL-Hancock steps 1+2:
                //Left state:
                cL = speed_of_sound(P_vec[j*n+i]);
                rho_i = Wi.z1*Wi.rho1 + (1-Wi.z1)*Wi.rho2;
                PbarL.W.rho1 = Wi.rho1 + 0.5*(1-(dt/dy)*Wi.v)*delta_i.rho1 - 0.5*(dt/dy)*Wi.rho1*delta_i.v;        
                PbarL.W.rho2 = Wi.rho2 + 0.5*(1-(dt/dy)*Wi.v)*delta_i.rho2 - 0.5*(dt/dy)*Wi.rho2*delta_i.v;
                PbarL.W.u = Wi.u + 0.5*(1-(dt/dy)*Wi.v)*delta_i.u;
                PbarL.W.v = Wi.v + 0.5*(1-(dt/dy)*Wi.v)*delta_i.v - 0.5*(dt/dy)*(1/rho_i)*delta_i.p;
                PbarL.W.p = Wi.p + 0.5*(1-(dt/dy)*Wi.v)*delta_i.p - 0.5*(dt/dy)*rho_i*pow(cL,2)*delta_i.v;
                PbarL.W.z1 = Wi.z1 + 0.5*(1-(dt/dy)*Wi.v)*delta_i.z1;
                Prop_WtoUF(PbarL);

                //Right state:
                cR = speed_of_sound(P_vec[(j+1)*n+i]);
                rho_n = Wn.z1*Wn.rho1 + (1-Wn.z1)*Wn.rho2;
                PbarR.W.rho1 = Wn.rho1 - 0.5*(1+(dt/dy)*Wn.v)*delta_n.rho1 - 0.5*(dt/dy)*Wn.rho1*delta_n.v;
                PbarR.W.rho2 = Wn.rho2 - 0.5*(1+(dt/dy)*Wn.v)*delta_n.rho2 - 0.5*(dt/dy)*Wn.rho2*delta_n.v;
                PbarR.W.u = Wn.u - 0.5*(1+(dt/dy)*Wn.v)*delta_n.u;
                PbarR.W.v = Wn.v - 0.5*(1+(dt/dy)*Wn.v)*delta_n.v - 0.5*(dt/dy)*(1/rho_n)*delta_n.p;
                PbarR.W.p = Wn.p - 0.5*(1+(dt/dy)*Wn.v)*delta_n.p - 0.5*(dt/dy)*rho_n*pow(cR,2)*delta_n.v;
                PbarR.W.z1 = Wn.z1 - 0.5*(1+(dt/dy)*Wn.v)*delta_n.z1;
                Prop_WtoUF(PbarR);
           
                PbarL_vec_y[j*n+i] = PbarL;
                PbarR_vec_y[j*n+i] = PbarR;
                

            }else{
                PbarL_vec_x[j*n+i] = P_vec[j*n+i];
                PbarR_vec_x[j*n+i] = P_vec[j*n+i+1];
                PbarL_vec_y[j*n+i] = P_vec[j*n+i];
                PbarR_vec_y[j*n+i] = P_vec[(j+1)*n+i];
            }
        }
    }
    
    reflective_boundaries(PbarL_vec_x, n, m, nG, BC);
    reflective_boundaries(PbarR_vec_x, n, m, nG, BC);
    reflective_boundaries(PbarL_vec_y, n, m, nG, BC);
    reflective_boundaries(PbarR_vec_y, n, m, nG, BC);
    
    double rho_min = P_vec[0].rho;
    double ST_dt;
    
    for(int j = 0; j<m; ++j){
        for(int i = 0; i<n; ++i){
            //HLLC evaluation (within HLLC_x and HLLC_y classes):
            cell_vec_x[j*n+i] = HLLC_x(PbarL_vec_x[j*n+i], PbarR_vec_x[j*n+i]);
            cell_vec_y[j*n+i] = HLLC_y(PbarL_vec_y[j*n+i], PbarR_vec_y[j*n+i]);             

            if(cell_vec_x[j*n+i].Splus > Smax){
                Smax = cell_vec_x[j*n+i].Splus;
            }
            if(cell_vec_y[j*n+i].Splus > Smax){
                Smax = cell_vec_y[j*n+i].Splus;
            }
            if(P_vec[j*n+i].rho < rho_min){
                rho_min = P_vec[j*n+i].rho;
            }
            //max wave speed in x or y direction limits max time step

        }
    }

    assert(Smax > 0);
    //assert(dx==dy);
    dt = ((dx <= dy)? CFL*dx/Smax : CFL*dy/Smax);
    assert(dt > 0);     //check this is true rather than enforcing fabs(dt)
    
    if(surface_tension_dt){
        ST_dt = ((dx <= dy)? 0.5*pow(rho_min/(8*pow(M_PI,3)*CASE.sigma),0.5)*pow(dx,1.5) :\
                0.5*pow(rho_min/(8*pow(M_PI,3)*CASE.sigma),0.5)*pow(dy,1.5));

        assert(ST_dt > 0); 
        dt = ((dt <= ST_dt)? dt : ST_dt);
    }
    HLLC_x cell_x = cell0_x;  //any initialisation
    HLLC_x prev_cell_x = cell0_x;
    HLLC_y cell_y = cell0_y;  //any initialisation
    HLLC_y prev_cell_y = cell0_y;
    //place-holders for update variables
    ConsU Unew(0,0,0,0,0,0);
    ConsU Ui(0,0,0,0,0,0);
    ConsF zeroF(0,0,0,0,0,0);
    ConsG zeroG(0,0,0,0,0,0);
    Properties Pnew = P_vec[0];
    double Fw_x, Fw_y;

    double plusE;
    
    for(int j = 0; j<m; ++j){
        for(int i = 0; i<n; ++i){

            cell_x = cell_vec_x[j*n+i];
            cell_y = cell_vec_y[j*n+i];
            Ui = P_vec[j*n+i].U;

            if(i == 0){
                prev_cell_x = cell_x;
            }else{
                prev_cell_x = cell_vec_x[j*n+i-1];
            }

            if(j == 0){
                prev_cell_y = cell_y;
            }else{
                prev_cell_y = cell_vec_y[(j-1)*n+i];
            }

            if(P_vec[j*n+i].n_x >= 0){
                Fw_x = cell_x.wR - cell_x.wL;
            }else{
                Fw_x = prev_cell_x.wR - prev_cell_x.wL;
            }

            if(P_vec[j*n+i].n_y >= 0){
                Fw_y = cell_y.wR - cell_y.wL;
            }else{
                Fw_y = prev_cell_y.wR - prev_cell_y.wL;
            }

            //conservative update formula:
            Unew.z1rho1 = Ui.z1rho1 + dt/dx * (prev_cell_x.Fout.mass1 - cell_x.Fout.mass1) \
                                    + dt/dy * (prev_cell_y.Gout.mass1 - cell_y.Gout.mass1);
            Unew.z2rho2 = Ui.z2rho2 + dt/dx * (prev_cell_x.Fout.mass2 - cell_x.Fout.mass2) \
                                    + dt/dy * (prev_cell_y.Gout.mass2 - cell_y.Gout.mass2);
            Unew.rhou = Ui.rhou + dt/dx * (prev_cell_x.Fout.xmom - cell_x.Fout.xmom) \
                                + dt/dy * (prev_cell_y.Gout.xmom - cell_y.Gout.xmom) \
                                + (dt)*(-sigma*P_vec[j*n+i].k*(1.0/dx)*(Fw_x));
            Unew.rhov = Ui.rhov + dt/dx * (prev_cell_x.Fout.ymom - cell_x.Fout.ymom) \
                                + dt/dy * (prev_cell_y.Gout.ymom - cell_y.Gout.ymom) \
                                + (dt)*(-sigma*P_vec[j*n+i].k*(1.0/dy)*(Fw_y));
            Unew.E = Ui.E + dt/dx * (prev_cell_x.Fout.en - cell_x.Fout.en) \
                          + dt/dy * (prev_cell_y.Gout.en - cell_y.Gout.en)
                          + dt*(-sigma*P_vec[j*n+i].k*\
                                  ((cell_x.Fout.uz1/cell_x.U.z1)*(1.0/dx)*(Fw_x)\
                                  + (cell_y.Gout.vz1/cell_y.U.z1)*(1.0/dy)*(Fw_y)));
            //pseudo-conservative z1 update
            Unew.z1 = Ui.z1 - dt/dx*(cell_x.Fout.uz1 - prev_cell_x.Fout.uz1 \
                     - Ui.z1*(cell_x.Fout.uz1/cell_x.U.z1 - prev_cell_x.Fout.uz1/prev_cell_x.U.z1))\
                     - dt/dy*(cell_y.Gout.vz1 - prev_cell_y.Gout.vz1 \
                     - Ui.z1*(cell_y.Gout.vz1/cell_y.U.z1 - prev_cell_y.Gout.vz1/prev_cell_y.U.z1));
            

            //variable conversion
            Pnew = P_vec[j*n+i]; //carry across other property info incl. curvature
            Pnew.U = Unew;
            Prop_UtoWF(Pnew);   //reference variable &Pnew - equates properties internally
            Pnew.S = cell_x.aL; //for plotting only
            P_vec_new[j*n+i] = Pnew;
                        
        }
    }
    

    //curvature_calc(P_vec_new, CASE, k_avr); 
    updateBoundary(P_vec_new, n, m, nG, BC);

    return P_vec_new;
}


vector< vector<Properties> > generate_solution(TestCase CASE, Vector &time_vec, 
        bool MUSCL, string scheme){

    int n = CASE.n;
    int m = CASE.m;
    int nG = CASE.nG;
    double CFL = CASE.CFL;
    double T0 = CASE.T0;
    double Tf = CASE.Tf;
    double dx = CASE.dx;
    double dy = CASE.dy;
    double T = T0;
    double L = CASE.L;
    double x0 = CASE.x0;
    string BC = CASE.BCs;
    vector<Properties> P_vec = initial_construction(CASE);    
    vector<Properties> P_vec_new = P_vec; //initialised variable for subsequent updating
    

    HLLC_x cell0_x(P_vec[0], P_vec[0]);
    vector<HLLC_x> cell_vec_x(n*m, cell0_x);
    HLLC_y cell0_y(P_vec[0], P_vec[0]);
    vector<HLLC_y> cell_vec_y(n*m, cell0_y);
    double Smax=0;
    double dt = 0;
    double view_dt = Tf/view_resolution;
    int view_num = 0;
    
    vector< vector<Properties> > MATRIX;  //solution is generated as a space*time matrix 
                                    //of the primitive variables 
    
    k_avr = (1.0/CASE.R);
    if(curvature){
        curvature_calc(P_vec, CASE, k_avr);
    }
    updateBoundary(P_vec, n, m, nG, BC);
    
    MATRIX.push_back(P_vec);
    ++view_num; 

    double rho_min = P_vec[0].rho;
    double ST_dt; //surface tension restricted dt


    for(int j = 1; j < m; ++j){
        for(int i = 1; i < n; ++i){
            cell0_x = HLLC_x(P_vec[j*n+i-1], P_vec[j*n+i]);
            cell0_y = HLLC_y(P_vec[(j-1)*n+i], P_vec[j*n+i]);
            if(cell0_x.Splus > Smax){
                Smax = cell0_x.Splus;
            }
            if(cell0_y.Splus > Smax){
                Smax = cell0_y.Splus;
            }
            if(P_vec[j*n+i].rho < rho_min){
                rho_min = P_vec[j*n+i].rho;
            }
        }
    }
    cout << "Smax0 = " << Smax << endl;
    cout << "rho_min = " << rho_min << endl;
    assert(Smax > 0);
    dt = ((dx <= dy)? CFL*dx/Smax : CFL*dy/Smax);
    
    if(surface_tension_dt){
        ST_dt = ((dx <= dy)? 0.5*pow(rho_min/(8*pow(M_PI,3)*CASE.sigma),0.5)*pow(dx,1.5) :\
                0.5*pow(rho_min/(8*pow(M_PI,3)*CASE.sigma),0.5)*pow(dy,1.5));
        
        cout << "dt = " << dt << endl;
        cout << "ST required dt <= " << ST_dt << endl;

        dt = ((dt <= ST_dt)? dt : ST_dt);
        cout << "actual dt = " << dt << endl;
    }
    time_vec.push_back(0);
    int nT = 1; //zero time is the first time step
    
    double percentage;

    if(time_stepping){
    while(T < Tf){
        Smax = 0;
        if(curvature){
            curvature_calc(P_vec, CASE, k_avr);
        }
        P_vec_new = updateP(P_vec, CASE, dt, CFL, Smax, MUSCL, scheme);
        P_vec = P_vec_new; //save time by not copying?
        T += dt;
        percentage = T/CASE.Tf * 100;
        cout << "T = " << T << setw(10) << " (" << percentage << "%) " << endl;
        //current simulation time printed to terminal
        time_vec.push_back(T);        
        nT += 1; //time step counter
        if(T >= view_num*view_dt){
            MATRIX.push_back(P_vec);
            ++view_num;
        }
        //if(nT >= 2){
        //    break;
        //}
    }
    }
    
    cout << "(actual final time = " << T << "s)" << endl; 
    
    return MATRIX;
}


//(beginning) DATA AND VISUALISATION -----------------------------------------------------//

Vector writetofile(TestCase CASE, vector< vector<Properties> > MATRIX, int view_resolution, double x0, double x1, double y0, double y1){

    Vector output(14); //returning: {u_min, u_max, rho_min, rho_max,
                                // p_min, p_max, e_min, e_max, c_min, c_max, rho_grad_min/max}
                                // for use in automatic plot scaling later

    string name = CASE.datfile;
    int n = CASE.n;     //number of spatial nodes
    int nT = CASE.nT;   //number of time steps
    double dx = CASE.dx;
    double dy = CASE.dy;
    assert(MATRIX.size() == nT);
    assert(MATRIX[0].size() == n*CASE.m);
    ofstream ufile, rhofile, pfile, efile, z1file, kfile, finalT, slice1D, MSrho;
    ofstream k_dist;
    ufile.open("./data/u_" + name);
    rhofile.open("./data/rho_" + name);
    pfile.open("./data/p_" + name);
    efile.open("./data/e_" + name);
    z1file.open("./data/z1_" + name);
    kfile.open("./data/k_" + name);
    finalT.open("./data/finalT" + name);
    slice1D.open("./data/slice1D" + name);
    MSrho.open("./data/MSrho_" + name);
    k_dist.open("./data/k_distribution.txt");
    //just initialising these to possible max/min values
    double u_min=CASE.WL.u, u_max=u_min, rho_min=CASE.WL.rho1, rho_max=rho_min; 
    double p_min=CASE.WL.p, p_max=p_min, e_min=MATRIX[0][0].e, e_max=e_min;
    double z1_min=1, z1_max=0, MSrho_min=1, MSrho_max=0, k_min=1, k_max=0;
    double V;
    double rho_grad;
    double x,y;
    //1D slice variables:
    assert(x0 == x1 || y0 == y1);

    for(int j = 0; j < CASE.m; ++j){
        //first column of each file is the space vector
        for(int i = 0; i < n; ++i){
            x = (i-CASE.nG)*dx;
            y = (j-CASE.nG)*dy;
            ufile << x << ' ' << y << ' ';
            rhofile << x << ' ' << y << ' ';
            pfile << x << ' ' << y << ' ';
            efile << x << ' ' << y << ' ';
            z1file << x << ' ' << y << ' ';
            kfile << x << ' ' << y << ' ';
            finalT << x << ' ' << y << ' ';
            if(i < n-1 && j < CASE.m-1){
                MSrho << x << ' ' << y << ' ';
            }
            for(int k = 0; k < nT; ++k){

                V = pow(pow(MATRIX[k][j*n+i].W.u, 2) + pow(MATRIX[k][j*n+i].W.v, 2),0.5);
                ufile << V << ' ';
                if(k== nT-1){
                    finalT << V << ' ';
                    if(x >= x0 && x < x1+dx && y >= y0 && y< y1+dy){
                    slice1D << x << ' ' << y << ' ' << V << ' ';
                    }
                }
                if(V < u_min){
                    u_min = V;
                }
                if(V > u_max){
                    u_max = V;
                }

                rhofile << MATRIX[k][j*n+i].rho << ' ';
                if(k== nT-1){
                    finalT << MATRIX[k][j*n+i].rho << ' ';
                    if(x >= x0 && x < x1+dx && y >= y0 && y< y1+dy){
                    slice1D << MATRIX[k][j*n+i].rho << ' ';
                    }
                }
                if(MATRIX[k][j*n+i].rho < rho_min){
                    rho_min = MATRIX[k][j*n+i].rho;
                }
                if(MATRIX[k][j*n+i].rho > rho_max){
                    rho_max = MATRIX[k][j*n+i].rho;
                }

                if(i < n-1 && j < CASE.m-1){
                    rho_grad = pow(pow((MATRIX[k][j*n+i+1].rho - MATRIX[k][j*n+i].rho)/dx,2) + \
                                    pow((MATRIX[k][(j+1)*n+i].rho - MATRIX[k][j*n+i].rho)/dy,2), 0.5);
                    MSrho << rho_grad << ' ';
                    if(rho_grad < MSrho_min){
                        MSrho_min = rho_grad;
                    }
                    if(rho_grad > MSrho_max){
                        MSrho_max = rho_grad;
                    }
                }


                pfile << MATRIX[k][j*n+i].W.p << ' ';
                if(k== nT-1){
                    finalT << MATRIX[k][j*n+i].W.p << ' ';
                    if(x >= x0 && x < x1+dx && y >= y0 && y< y1+dy){
                    slice1D << MATRIX[k][j*n+i].W.p << ' ';
                    }
                }
                if(MATRIX[k][j*n+i].W.p < p_min){
                    p_min = MATRIX[k][j*n+i].W.p;
                }
                if(MATRIX[k][j*n+i].W.p > p_max){
                    p_max = MATRIX[k][j*n+i].W.p;
                }
            

                efile << MATRIX[k][j*n+i].e << ' ';
                if(k== nT-1){
                    finalT << MATRIX[k][j*n+i].e << ' ';
                    if(x >= x0 && x < x1+dx && y >= y0 && y< y1+dy){
                    slice1D << MATRIX[k][j*n+i].e << ' ';
                    }
                }
                if(MATRIX[k][j*n+i].e < e_min){
                    e_min = MATRIX[k][j*n+i].e;
                }
                if(MATRIX[k][j*n+i].e > e_max){
                    e_max = MATRIX[k][j*n+i].e;
                }
                

                z1file << MATRIX[k][j*n+i].U.z1 << ' ';
                if(k== nT-1){
                    finalT << MATRIX[k][j*n+i].U.z1 << ' ';
                    finalT << MATRIX[k][j*n+i].S << ' ';
                    finalT << MATRIX[k][j*n+i].k << ' ';
                    if(x >= x0 && x < x1+dx && y >= y0 && y< y1+dy){
                        slice1D << MATRIX[k][j*n+i].k << ' ';
                        slice1D << MATRIX[k][j*n+i].U.z1 << ' ';
                    }
                    if(MATRIX[k][j*n+i].k != 0){
                        k_dist << MATRIX[k][j*n+i].k;
                        k_dist << ' ';
                    }
                    if(x >= x0 && x < x1+dx && y >= y0 && y< y1+dy){
                    slice1D << MATRIX[k][j*n+i].W.z1 << ' ';
                    }
                }
                if(MATRIX[k][j*n+i].W.z1 < z1_min){
                    z1_min = MATRIX[k][j*n+i].W.z1;
                }
                if(MATRIX[k][j*n+i].W.z1 > z1_max){
                    z1_max = MATRIX[k][j*n+i].W.z1;
                }

                kfile << MATRIX[k][j*n+i].k << ' ';
                if(MATRIX[k][j*n+i].k < k_min){
                    k_min = MATRIX[k][j*n+i].k;
                }
                if(MATRIX[k][j*n+i].k > k_max){
                    k_max = MATRIX[k][j*n+i].k;
                }
                //essentially transposing matrix to file write s.t.
                //vector downwards in space, with each column stepping in time
            }
            ufile << '\n';
            rhofile << '\n';
            pfile << '\n';
            efile << '\n';
            z1file << '\n';
            kfile << '\n';
            finalT << '\n';
            slice1D << '\n';
            if(i < n-1 && j < CASE.m-1){
                MSrho << '\n';
            }
        }
    }
    ufile.close();
    rhofile.close();
    pfile.close();
    efile.close();
    z1file.close();
    kfile.close();
    finalT.close();
    slice1D.close();
    MSrho.close();
    k_dist.close();
    double tol = 1e-3;
    if(fabs(u_max - u_min) < tol){ u_min += -0.1; u_max += 0.1;}
    if(fabs(rho_max - rho_min) < tol){ rho_min += -0.1; rho_max += 0.1;}
    if(fabs(p_max - p_min) < tol){ p_min += -0.1; p_max += 0.1;}
    if(fabs(e_max - e_min) < tol){ e_min += -0.1; e_max += 0.1;}
    if(fabs(z1_max - z1_min) < tol){ z1_min += -0.1; z1_max += 0.1;}
    if(fabs(k_max - k_min) < tol){ k_min += -0.1; k_max += 0.1;}
    if(fabs(MSrho_max - MSrho_min) < tol){ MSrho_min += -0.1; MSrho_max += 0.1;}    

    double arr[14] = {u_min, u_max, rho_min, rho_max, p_min, p_max, e_min, e_max, 
        z1_min, z1_max, MSrho_min, MSrho_max, k_min, k_max};
    output.assign(arr, &arr[14]);
    return output;
}

void GNUplot_gif(TestCase CASE, Vector domain){
    string gpname = CASE.gpname;
    string datfile = CASE.datfile;
    string title = CASE.title;
    string giffile = CASE.giffile;
    int nT = CASE.nT;

    string u_datfile = "./data/u_" + datfile;
    string rho_datfile = "./data/rho_" + datfile;
    string p_datfile = "./data/p_" + datfile;
    string e_datfile = "./data/e_" + datfile;
    string z1_datfile = "./data/z1_" + datfile;
    string k_datfile = "./data/k_" + datfile;

    ofstream file;
    file.open(gpname);
    file << "#!/usr/local/bin/gnuplot -persist \n";
    file << "set terminal gif animate delay 8 \n";
    file << "n = " << nT << '\n';
    file << "set out \'" << giffile << "\' \n";
    file << "set key off \n";
    file << "set size " << CASE.L << ", " << CASE.Y << "\n";
    file << "set view map \n";
    file << "unset surface \n";
    file << "set contour \n";
    file << "set palette rgbformulae 33,13,10 \n";
    
    file << "do for [i=0:" << nT-1 << "] { \n    ";
    file << "j = i+3 \n";
    file << "set multiplot layout 2,2 \n";
    //rho plot:
    file << "set title \"" << "Density - rho" << "\" \n";
    file << "set cbrange[" << domain[2] - 0.3*fabs(domain[3]-domain[2]) << \
        ":" << domain[3] + 0.3*fabs(domain[3]-domain[2]) << "] \n";
    file << "plot \"" << rho_datfile << "\" using 1:2:j with image\n";
    //p plot:
    file << "set title \"" << "Pressure - p" << "\" \n";
    file << "set cbrange[" << domain[4] - 0.3*fabs(domain[5]-domain[4]) << \
        ":" << domain[5] + 0.3*fabs(domain[5]-domain[4]) << "] \n";
    file << "plot \"" << p_datfile << "\" using 1:2:j with image\n";
    /*    
    //u plot:
    file << "set title \"" << "Absolute velocity - V" << "\" \n";
    file << "set cbrange[" << domain[0] - 0.3*fabs(domain[1]-domain[0]) << \
        ":" << domain[1] + 0.3*fabs(domain[1]-domain[0]) << "] \n";
    file << "plot \"" << u_datfile << "\" using 1:2:j with image\n";
    */
    
    //z1 plot:
    file << "set title \"" << "Component mixture level - z1" << "\" \n";
    file << "set cbrange[" << domain[8] - 0.3*fabs(domain[9]-domain[8]) << \
        ":" << domain[9] + 0.3*fabs(domain[9]-domain[8]) << "] \n";
    file << "plot \"" << z1_datfile << "\" using 1:2:j with image\n";
    
    //k plot:
    file << "set title \"" << "Curvature - k" << "\" \n";
    file << "set cbrange[" << domain[12] - 0.3*fabs(domain[13]-domain[12]) << \
        ":" << domain[13] + 0.3*fabs(domain[13]-domain[12]) << "] \n";
    file << "plot \"" << k_datfile << "\" using 1:2:j with image\n";
    
    /*
    //e plot:
    file << "set title \"" << "Internal energy - e" << "\" \n";
    file << "set cbrange[" << domain[6] - 0.3*fabs(domain[7]-domain[6]) << \
        ":" << domain[7] + 0.3*fabs(domain[7]-domain[6]) << "] \n";
    file << "plot \"" << e_datfile << "\" using 1:2:j with image\n";
    */
    
    file << "unset multiplot \n } \n";
    file << "exit";
    file.close();

}

void GNUplot_density_plots(TestCase CASE, Vector domain, string type){
    //type must be: MS = mock-schlieren or rho = normal
    string gpname = type + CASE.gpname;
    string datfile = "rho_" + CASE.datfile;
    string title = type + CASE.title;
    string giffile = CASE.giffile;
    int nT = CASE.nT;

    string rho_datfile;
    if(type == "MS"){
        rho_datfile = "./data/MS" + datfile;
    }else if(type == "rho"){
        rho_datfile = "./data/" + datfile;
    }

    ofstream file;
    file.open(gpname);
    file << "#!/usr/local/bin/gnuplot -persist \n";
    file << "set terminal gif animate delay 8 \n";
    file << "n = " << nT << '\n';
    file << "set out \'./visualisation/" << title << ".gif\' \n";
    file << "set key off \n";
    file << "set title \"" << "Density ";
    if(type == "MS"){
       file << "gradient - Mock-Schlieren";
    } 
    file << "\" \n";
    file << "set size " << 1 << ", " << CASE.Y/CASE.L << "\n";
    file << "set view map \n";
    file << "unset surface \n";
    file << "set contour \n";
    file << "set palette rgbformulae 33,13,10 \n";
    
    file << "do for [i=0:" << nT-1 << "] { \n    ";
    file << "j = i+3 \n";
    //rho plot:
    file << "plot \"" << rho_datfile << "\" using 1:2:j with image\n";

    file << "} \n";
    
    file << "exit";
    file.close();

}

void GNUplot_curvature_plot(TestCase CASE){
    //type must be: MS = mock-schlieren or rho = normal
    string gpname = "k_" + CASE.gpname;
    string datfile = "./data/k_" + CASE.datfile;
    string title = "k_" + CASE.title;
    string giffile = CASE.giffile;
    int nT = CASE.nT;

    ofstream file;
    file.open(gpname);
    file << "#!/usr/local/bin/gnuplot -persist \n";
    file << "set terminal gif animate delay 8 \n";
    file << "n = " << nT << '\n';
    file << "set out \'./visualisation/" << title << ".gif\' \n";
    file << "set key off \n";
    file << "set title \"Curvature - k\" \n";
    file << "set size " << 1 << ", " << CASE.Y/CASE.L << "\n";
    file << "set view map \n";
    file << "unset surface \n";
    file << "set contour \n";
    file << "set palette rgbformulae 33,13,10 \n";
    
    file << "do for [i=0:" << nT-1 << "] { \n    ";
    file << "j = i+3 \n";
    //rho plot:
    file << "plot \"" << datfile << "\" using 1:2:j with image\n";

    file << "} \n";
    
    file << "exit";
    file.close();

}




void GNUplot_finalT(TestCase CASE, string filetype){
    ofstream file;
    file.open(CASE.title + " finalT.gp");
    string datfile = "./data/finalT" + CASE.datfile;
    file << "#!/usr/local/bin/gnuplot -persist \n";
    file << "set terminal " << filetype << " \n";
    file << "set out './visualisation/"<< CASE.title << " finalT." << filetype <<"\' \n";
    file << "set view map \n";
    file << "unset surface \n";
    file << "set contour \n";
    file << "set dgrid3d \n";
    file << "set key off \n";
    file << "set size " << CASE.L << ", " << CASE.Y << "\n";
    file << "set palette rgbformulae 33,13,10 \n";
    file << "set multiplot layout 2,2 \n";
    file << "set title \"Density - rho\" \n";
    //file << "set title \"x-component normal - Fx\" \n";
    file << "plot '"<< datfile << "\' using 1:2:4 with image";
    if(CASE.exact_sol_file == "none"){
        file << " \n";
    }else{
        file << ", '" << CASE.exact_sol_file << "rho.dat' using 1:2 with lines lc rgb 'black' \n";
    }

    file << "set title \"Pressure - p\" \n";
    //file << "set title \"y-component normal - Fy\" \n";
    file << "plot '"<< datfile << "\' using 1:2:5 with image";
    if(CASE.exact_sol_file == "none"){
        file << " \n";
    }else{
        file << ", '" << CASE.exact_sol_file << "p.dat' using 1:2 with lines lc rgb 'black' \n";
    }
    
    file << "set title \"Velocity - V\" \n";
    file << "plot '"<< datfile << "\' using 1:2:3 with image";
    if(CASE.exact_sol_file == "none"){
        file << " \n";
    }else{
        file << ", '" << CASE.exact_sol_file << "u.dat' using 1:2 with lines lc rgb 'black' \n";
    }
    
    file << "set title \"Curvature - k\" \n";
    file << "plot '"<< datfile << "\' using 1:2:9 with image";
    if(CASE.exact_sol_file == "none"){
        file << " \n";
    }else{
        file << ", '" << CASE.exact_sol_file << "k.dat' using 1:2 with lines lc rgb 'black' \n";
    }
    /*
    file << "set title \"z1\" \n";
    file << "plot '"<< datfile << "\' using 1:2:7 with image";
    if(CASE.exact_sol_file == "none"){
        file << " \n";
    }else{
        file << ", '" << CASE.exact_sol_file << "z1.dat' using 1:2 with lines lc rgb 'black' \n";
    }
    */
    /*
    file << "set title \"Internal energy - e\" \n";
    file << "plot '"<< datfile << "\' using 1:2:6 with image";
    if(CASE.exact_sol_file == "none"){
        file << " \n";
    }else{
        file << ", '" << CASE.exact_sol_file << "e.dat' using 1:2 with lines lc rgb 'black' \n";
    }
    */
    /*
    file << "set title \" Sound Speed\" \n";
    file << "splot '"<< datfile << "\' using 1:2:8";
    if(CASE.exact_sol_file == "none"){
        file << " \n";
    }else{
        file << ", '" << CASE.exact_sol_file << "S.dat' using 1:2 with lines lc rgb 'black' \n";
    }
    */
    file << "exit";
}

void GNUplot_slice1D(TestCase CASE, string filetype, string slice){
    ofstream file;
    file.open(CASE.title + " slice1D.gp");
    string datfile = "./data/slice1D" + CASE.datfile;
    file << "#!/usr/local/bin/gnuplot -persist \n";
    file << "set terminal " << filetype << " \n";
    file << "set out './visualisation/"<< CASE.title << " slice1D." << filetype <<"\' \n";
    file << "set key off \n";
    file << "set multiplot layout 3,2 \n";
    file << "set title \"Density - rho\" \n";
    if(slice == "x-slice"){
        file << "plot '"<< datfile << "\' using 1:4 with linespoints pt 6 ps 0.3 lc rgb \'blue\'";
    }else if(slice == "y-slice"){
        file << "plot '"<< datfile << "\' using 2:4 with linespoints pt 6 ps 0.3 lc rgb \'blue\'";
    }
    if(CASE.exact_sol_file == "none"){
        file << " \n";
    }else{
        file << ", '" << CASE.exact_sol_file << "rho.dat' using 1:2 with lines lc rgb 'black' \n";
    }
    file << "set title \"Absolute velocity - V\" \n";
    if(slice == "x-slice"){
        file << "plot '"<< datfile << "\' using 1:3 with linespoints pt 6 ps 0.3 lc rgb \'blue\'";
    }else if(slice == "y-slice"){
        file << "plot '"<< datfile << "\' using 2:3 with linespoints pt 6 ps 0.3 lc rgb \'blue\'";
    }
    if(CASE.exact_sol_file == "none"){
        file << " \n";
    }else{
        file << ", '" << CASE.exact_sol_file << "u.dat' using 1:2 with lines lc rgb 'black' \n";
    }
    file << "set title \"Pressure - p\" \n";
    if(slice == "x-slice"){
        file << "plot '"<< datfile << "\' using 1:5 with linespoints pt 6 ps 0.3 lc rgb \'blue\'";
    }else if(slice == "y-slice"){
        file << "plot '"<< datfile << "\' using 2:5 with linespoints pt 6 ps 0.3 lc rgb \'blue\'";
    }
    if(CASE.exact_sol_file == "none"){
        file << " \n";
    }else{
        file << ", '" << CASE.exact_sol_file << "p.dat' using 1:2 with lines lc rgb 'black' \n";
    }
    
    file << "set title \"Internal energy - e\" \n";
    if(slice == "x-slice"){
        file << "plot '"<< datfile << "\' using 1:6 with linespoints pt 6 ps 0.3 lc rgb \'blue\'";
    }else if(slice == "y-slice"){
        file << "plot '"<< datfile << "\' using 2:6 with linespoints pt 6 ps 0.3 lc rgb \'blue\'";
    }
    if(CASE.exact_sol_file == "none"){
        file << " \n";
    }else{
        file << ", '" << CASE.exact_sol_file << "e.dat' using 1:2 with lines lc rgb 'black' \n";
    }
    file << "set title \"Curvature - k\" \n";
    if(slice == "x-slice"){
        file << "plot '"<< datfile << "\' using 1:7 with linespoints pt 6 ps 0.3 lc rgb \'blue\' \n";
    }else if(slice == "y-slice"){
        file << "plot '"<< datfile << "\' using 2:7 with linespoints pt 6 ps 0.3 lc rgb \'blue\' \n";
    }
    file << "set title \"mixture level - z1\" \n";
    if(slice == "x-slice"){
        file << "plot '"<< datfile << "\' using 1:8 with linespoints pt 6 ps 0.3 lc rgb \'blue\' \n";
    }else if(slice == "y-slice"){
        file << "plot '"<< datfile << "\' using 2:8 with linespoints pt 6 ps 0.3 lc rgb \'blue\' \n";
    }



    /*
    file << "set title \" Sound Speed\" \n";
    if(slice == "x-slice"){
        file << "plot '"<< datfile << "\' using 1:8 with linespoints pt 6 ps 0.3 lc rgb \'blue\'";
    }else if(slice == "y-slice"){
        file << "plot '"<< datfile << "\' using 2:8 with linespoints pt 6 ps 0.3 lc rgb \'blue\'";
    }
    if(CASE.exact_sol_file == "none"){
        file << " \n";
    }else{
        file << ", '" << CASE.exact_sol_file << "S.dat' using 1:2 with lines lc rgb 'black' \n";
    }
    */
    file << "exit";
}

//--------------------------------------------------------------(end) DATA AND VISUALISATION //

void tellmethings(TestCase CASE, bool MUSCL, string scheme){
    cout << "|| " << CASE.title << " || " << endl;
    cout << "Test case description: " << CASE.description << endl;
    cout << "Test case models materials: \"" << material1 << "\" and \"" << material2;
    cout << "\" under the \"" << EoS << "\" equation of state. " << endl;
    if(surface_tension_dt){
        cout << "with surface tension time-step restriction" << endl;
    }
    cout << "spatial discretisation x * y = " << CASE.n - 2*CASE.nG << " * " << CASE.m - 2*CASE.nG;
    cout << " = " << (CASE.n - 2*CASE.nG)*(CASE.m - 2*CASE.nG) << " cells (Cartesian grid)" << endl;
    cout << "simulated time = " << CASE.Tf << "s" << endl;
    cout << "HLLC solver- run as ";
    if(MUSCL){
        cout << "second order scheme: MUSCL, slope limiter = " << scheme << endl;
    }else{ cout << "first order scheme " << endl;
    }
    cout << "with " << CASE.BCs << " boundary conditions" << endl;
    cout << "number of time steps stored in simulated time = " << CASE.nT << endl;
    cout << "max. number of time steps written for visualisation = " << view_resolution << endl;
    cout << "5 data files generated in ./data/: \n  *" << CASE.datfile << "(x5 properties)\n";
    cout << "6 .gp files generated which can be compiled for visualisation: \n";
    cout << "for gif animation of all properties: \n  *" << CASE.gpname << endl;
    cout << "for final time solution image: \n  *" << CASE.title << " finalT.gp" << endl;
    cout << "for 1D cross-sectional profile: \n  *" << CASE.title << " slice1D.gp" << endl;
    cout << "for mock-schlieren gif animation: \n  *MS" << CASE.gpname << endl;
    cout << "for density plot gif animation: \n  *rho" << CASE.gpname << endl;
    cout << "for curvature plot gif animation: \n  *k_" << CASE.gpname << endl;
    cout << "compile with: \"chmod u+x <name.gp> \" " << \
        "then run \"./<name>.gp\", to generate file located in ./visualisation/" << endl;
    cout << "(compiling a gif with more than 300 time steps will take a super long time and is not recommended)\n";
    cout << "done." << endl;
}


int main(){
        
    Vector time_vec;
    
    vector< vector<Properties> > MATRIX = generate_solution(CASE, time_vec, MUSCL, scheme);
    CASE.nT = MATRIX.size();
    cout << "dimensions: " << MATRIX[0].size() << endl;

    tellmethings(CASE, MUSCL, scheme); 

    cout << "generating visualisation files ..." << endl;

    //slice1D variables:
    //x or y must be constant
    //check values lie in valid spatial range
    double x0 = 0.00;
    double x1 = 1.0;
    double y0 = 0.5;
    double y1 = 0.5;
    string slice = "x-slice"; //must be "x-slice" or "y-slice"

    Vector domain = writetofile(CASE, MATRIX, view_resolution, x0, x1, y0, y1);
    cout << "writetofile executed" << endl;
    
    GNUplot_gif(CASE, domain);
    cout << "GNUplot_gif file written" << endl;
    GNUplot_finalT(CASE, "svg");
    cout << "GNUplot_finalT file written" << endl;
    GNUplot_slice1D(CASE, "svg", slice);
    cout << "GNUplot_slice1D file written" << endl;
    GNUplot_density_plots(CASE, domain, "MS");
    cout << "GNUplot mock-schlieren gif file written" << endl;
    GNUplot_density_plots(CASE, domain, "rho");
    cout << "GNUplot density gif file written" << endl;
    if(curvature){
        GNUplot_curvature_plot(CASE);
        cout << "GNUplot curvature gif file written" << endl;
    }

    //cout << "post-shock pressure = " << MATRIX[CASE.nT-1][(int)(0.1*CASE.m)*CASE.n+(int)(0.5*CASE.n)].W.p << endl;
    //cout << "post-shock density = " << MATRIX[CASE.nT-1][(int)(0.1*CASE.m)*CASE.n+(int)(0.5*CASE.n)].rho << endl;
    //cout << "post-shock velocity = " << MATRIX[CASE.nT-1][(int)(0.1*CASE.m)*CASE.n+(int)(0.5*CASE.n)].W.v << endl;



    return 0;
}







