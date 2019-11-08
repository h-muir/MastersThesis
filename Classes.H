//three class definitions for: 
//   - Primitive variables (W) 
//   - Conservative variables (U)
//   - Conservative fluxes (F) 

class Prim{
    public:
        Prim(double, double, double, double, double, double);
        double rho1;
        double rho2;
        double u;
        double v;
        double p;
        double z1;
};
Prim::Prim(double _rho1, double _rho2, double _u, double _v, double _p, double _z1){
    rho1 = _rho1; rho2 = _rho2, u = _u; v = _v; p=_p; z1=_z1;
}

class ConsU{
    public:
        ConsU(double, double, double, double, double, double);
        double z1rho1;
        double z2rho2;
        double rhou;
        double rhov;
        double E;
        double z1;
};
ConsU::ConsU(double _z1rho1, double _z2rho2, double _rhou, double _rhov, double _E, double _z1){
    z1rho1 = _z1rho1; z2rho2 = _z2rho2; rhou = _rhou; rhov = _rhov; E = _E; z1 = _z1;
}

class ConsF{
    public:
        ConsF(double, double, double, double, double, double);
        double mass1;
        double mass2;
        double xmom;
        double ymom;
        double en;
        double uz1;
};
ConsF::ConsF(double _mass1, double _mass2, double _xmom, double _ymom, double _en, double _uz1){
    mass1 = _mass1; mass2 = _mass2; xmom = _xmom; ymom = _ymom; en = _en; uz1 = _uz1;
}

class ConsG{
    public:
        ConsG(double, double, double, double, double, double);
        double mass1;
        double mass2;
        double xmom;
        double ymom;
        double en;
        double vz1;
};
ConsG::ConsG(double _mass1, double _mass2, double _xmom, double _ymom, double _en, double _vz1){
    mass1 = _mass1; mass2 = _mass2; xmom = _xmom; ymom = _ymom; en = _en; vz1 = _vz1;
}


//MATERIAL PROPERTIES class:

class MaterialProperties{
    public:
        string material, EoS;
        double varGamma;                //all EoS
        double A, B, R1, R2, rho_0;     //additional JWL
        double xi1, xi2, Cv;            //additional cochran-chan
        double c_0, s, e_0, p_0;     //additional Hugoniot
        //error with use of a and s in equation
        double p_inf;                   //stiffened gas
        
        MaterialProperties(string _material, string _EoS){
            material = _material;
            EoS = _EoS;
            if(_material == "air" && _EoS == "ideal"){
                varGamma = 1.4 - 1;
            }
            else if(_material == "helium" && _EoS == "ideal"){
                varGamma = 1.67 - 1;
            }
            else if(_material == "TNT" && _EoS == "JWL"){
                varGamma = 0.25;
                A = 854.5e9;
                B = 20.5e9;
                R1 = 4.6;
                R2 = 1.35;
                rho_0 = 1840;
            }
            else if(_material == "TNT" && _EoS == "CC"){
                varGamma = 0.93;
                A = 12.87e9;
                B = 13.42e9;
                xi1 = 4.1;
                xi2 = 3.1;
                rho_0 = 1840;
                Cv = 1087;
            }
            else if(_material == "copper" && _EoS == "CC"){
                varGamma = 2.0;
                A = 145.67e9;
                B = 147.75e9;
                xi1 = 2.99;
                xi2 = 1.99;
                rho_0 = 8900;
                Cv = 393;
            }else if(_material == "nitromethane" && _EoS == "CC"){
                varGamma = 1.19;
                A = 0.819181e9;
                B = 1.50835e9;
                xi1 = 4.52969;
                xi2 = 1.42144;
                rho_0 = 1134;
                Cv = 1714;
            }else if(_material == "air" && _EoS == "CC"){
                varGamma = 0.4;
                A = 0.0;
                B = 0.0;
                xi1 = 0.0;
                xi2 = 0.0;
                rho_0 = 1134;
                Cv = 718;
            }else if(_material == "HMX" && _EoS == "Hugoniot"){
                varGamma = 0.7;
                c_0 = 3070;
                s = 1.79;
                rho_0 = 1891;
                e_0 = 0;        //CAN'T FIND ACTUAL VALUE
                p_0 = 0;        //CAN'T FIND ACTUAL VALUE
            }else if(_material == "air" && _EoS == "stiffened"){
                varGamma = 0.4;
                p_inf = 0.0;
            }else if(_material == "water" && _EoS == "stiffened"){
                varGamma = 3.4;
                p_inf = 6.0e8;
            }else if(_material == "gelatin10" && _EoS == "stiffened"){
                //10% gelatin in water mixture
                varGamma = 5.72;
                p_inf = 3.7e8;
            }else if(_material == "LS_water" && _EoS == "stiffened"){
                varGamma = 1.35;
                p_inf = 1.0e9;
                e_0 = -1167.0e3;
            }else if(_material == "air" && _EoS == "Hugoniot"){
                varGamma = 0.40;
                c_0 = 0.0;
                s = 0.0;
                rho_0 = 1.2;
                e_0 = 0;        
                p_0 = 0;        
            }else if(_material == "gel" && _EoS == "Hugoniot"){
                varGamma = 0.17;
                c_0 = 1520.0;
                s = 1.87;
                rho_0 = 1060.0;
                e_0 = 0;        
                p_0 = 0;        
            }else if(_material == "fake-liquid" && _EoS == "ideal"){
                varGamma = 2.1 - 1;
            }
            else{
                throw runtime_error("non-specified material + EoS pair: " + _material + " + " + _EoS);
            }
        }
};
