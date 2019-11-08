
class TestCase{
    public:
        int n;          //x-dimension spatial discretisation
        int m;          //y-dimension spatial discretisation
        double L, Y;    //1D spatial Length
        double Tf;      //final time
        double T0;      //initial time
        double CFL;     //Courant number
        double dx;      //x spatial step size
        double dy;      //y spatial step size
        double x0, x1, y0;      //property discontinuity position
        double o_x, o_y, R, k0; //some additional "bubble" construction paramaters
        double gamma1, gamma2;      
        string material1, material2;    //handles 2 materials
        string EoS;                     //single EoS 
        Prim WL, WR;
        string gpname, datfile, title, giffile, description;
        string exact_sol_file;
        string construction, BCs;
        int nG;         //number of ghost cells at the boundary: on each side
                        //2 required for reflected shock BCs
        int nT;         //number of time steps,
                        //not assigned until after solution and number of time steps are known
        double sigma;   //surface tension coefficient in N/m, with default 0 = no surface tension effects
        TestCase(string CASE):
            WL(0,0,0,0,0,0), WR(0,0,0,0,0,0){
                exact_sol_file = "none"; //auto declaration if no solution file exists
                if(CASE == "TestA"){
                    title = "TestA";
                    description = "constant properties test";
                    datfile = "euler_eqs_testA.dat";
                    gpname = "gifplot_testA.gp";
                    giffile = "./visualisation/testA.gif";
                    BCs = "transmissive";
                    WL.rho1 = 1.0; WL.rho2 = 1.0; WL.u = 1.0; WL.p = 1.0; WL.z1 = 1 - 1e-6;
                    WR.rho1 = 1.0; WR.rho2 = 1.0; WR.u = 1.0; WR.p = 1.0; WR.z1 = 1 - 1e-6;
                    nG = 2;
                    n = 100 + 2*nG;
                    m = n;
                    L = 1.0;    //space domain 0 < x < 1
                    Tf = 0.5;
                    T0 = 0;
                    CFL = 0.4;
                    dx = L/(n-1-2*nG);
                    dy = L/(m-1-2*nG);
                    x0 = 0.25;
                    material1 = "air";
                    material2 = "air";
                    EoS = "ideal";
                    construction = "default";
                }
                if(CASE == "TestA2"){
                    title = "TestA2";
                    description = "z1 discontinuity";
                    datfile = "euler_eqs_testA2.dat";
                    gpname = "gifplot_testA2.gp";
                    giffile = "./visualisation/testA2.gif";
                    BCs = "transmissive";
                    WL.rho1 = 1.0; WL.rho2 = 1.0; WL.u = 1.0; WL.p = 1.0; WL.z1 = 1 - 1e-6;
                    WR.rho1 = 1.0; WR.rho2 = 1.0; WR.u = 1.0; WR.p = 1.0; WR.z1 = 1e-6;
                    nG = 2;
                    n = 100 + 2*nG;
                    m = n;
                    L = 1.0;    //space domain 0 < x < 1
                    Tf = 0.5;
                    T0 = 0;
                    CFL = 0.4;
                    dx = L/(n-1-2*nG);
                    dy = L/(m-1-2*nG);
                    x0 = 0.25;
                    material1 = "air";
                    material2 = "air";
                    EoS = "ideal";
                    construction = "default";
                }
                if(CASE == "TestA3"){
                    title = "TestA3";
                    description = "z1 discontinuity with different gammas"; 
                    datfile = "euler_eqs_testA3.dat";
                    gpname = "gifplot_testA3.gp";
                    giffile = "./visualisation/testA3.gif";
                    BCs = "transmissive";
                    WL.rho1 = 1.0; WL.rho2 = 1.0; WL.u = 1.0; WL.p = 1.0; WL.z1 = 1 - 1e-6;
                    WR.rho1 = 1.0; WR.rho2 = 1.0; WR.u = 1.0; WR.p = 1.0; WR.z1 = 1e-6;
                    nG = 2;
                    n = 200 + 2*nG;
                    m = n;
                    L = 1.0;    //space domain 0 < x < 1
                    Tf = 0.15;
                    T0 = 0;
                    CFL = 0.4;
                    dx = L/(n-1-2*nG);
                    dy = L/(m-1-2*nG);
                    x0 = 0.25;
                    material1 = "air";
                    material2 = "helium";
                    EoS = "ideal";
                    construction = "default";
                }
                if(CASE == "Test1"){
                    title = "Test1";
                    description = "Toro test 1 equivalent";
                    datfile = "test1.dat";
                    gpname = "gifplot_test1.gp";
                    giffile = "./visualisation/test1.gif";
                    BCs = "transmissive";
                    //exact_sol_file = "./ExactSols/Toro411";
                    WL.rho1 = 1.0;   WL.rho2 = 1.0;   WL.u = 0.0; WL.v = 0; WL.p = 1.0; WL.z1 = 1 - 1e-6;
                    WR.rho1 = 0.125; WR.rho2 = 0.125; WR.u = 0.0; WR.v = 0; WR.p = 0.1; WR.z1 = 1 - 1e-6;
                    nG = 2;
                    n = 50 + 2*nG;
                    m = n;
                    L = 1.0;    //space domain 0 < x < 1
                    Tf = 0.25;
                    T0 = 0;
                    CFL = 0.4;
                    dx = L/(n-1-2*nG);
                    dy = L/(m-1-2*nG);
                    x0 = 0.5;
                    material1 = "air";
                    material2 = "air";
                    EoS = "ideal";
                    construction = "default2D";
                }
                if(CASE == "Test1_y"){
                    title = "Test1_y";
                    description = "Toro test 1 y-direction equivalent";
                    datfile = "test1_y.dat";
                    gpname = "gifplot_test1_y.gp";
                    giffile = "./visualisation/test1_y.gif";
                    BCs = "bottom-reflective";
                    //exact_sol_file = "./ExactSols/Toro411";
                    WL.rho1 = 0.125;   WL.rho2 = 0.125;   WL.u = 0.0; WL.v = 0.0; WL.p = 0.1; WL.z1 = 1 - 1e-6;
                    WR.rho1 = 1.0; WR.rho2 = 1.0; WR.u = 0.0; WR.v = 0.0; WR.p = 1.0; WR.z1 = 1 - 1e-6;
                    nG = 2;
                    n = 50 + 2*nG;
                    m = n;
                    L = 1.0;    //space domain 0 < x < 1
                    Tf = 0.50;
                    T0 = 0;
                    CFL = 0.4;
                    dx = L/(n-1-2*nG);
                    dy = L/(m-1-2*nG);
                    x0 = 0.5;
                    y0 = 0.5;
                    material1 = "air";
                    material2 = "air";
                    EoS = "ideal";
                    construction = "y-direction";
                }
                if(CASE == "Test1_V2"){
                    title = "Test1_V2";
                    description = "Toro test 1 equivalent";
                    datfile = "test1_V2.dat";
                    gpname = "gifplot_test1_V2.gp";
                    giffile = "./visualisation/test1_V2.gif";
                    BCs = "transmissive";
                    //exact_sol_file = "./ExactSols/Toro411";
                    WL.rho1 = 1.0;   WL.rho2 = 0.125;   WL.u = 0.0; WL.p = 1.0; WL.z1 = 1 - 1e-6;
                    WR.rho1 = 1.0; WR.rho2 = 0.125; WR.u = 0.0;  WR.p = 0.1; WR.z1 = 1e-6;
                    nG = 2;
                    n = 50 + 2*nG;
                    m = n;
                    L = 1.0;    //space domain 0 < x < 1
                    Tf = 0.25;
                    T0 = 0;
                    CFL = 0.4;
                    dx = L/(n-1-2*nG);
                    dy = L/(m-1-2*nG);
                    x0 = 0.5;
                    material1 = "air";
                    material2 = "air";
                    EoS = "ideal";
                    construction = "default2D";
                }
                if(CASE == "Test2"){
                    title = "Test2";
                    description = "Toro test 2 equivalent";
                    datfile = "test2.dat";
                    gpname = "gifplot_test2.gp";
                    giffile = "./visualisation/test2.gif";
                    BCs = "transmissive";
                    //exact_sol_file = "./ExactSols/Toro412";
                    WL.rho1 = 1.0;   WL.rho2 = 1.0;   WL.u = -2.0; WL.p = 0.4; WL.z1 = 1 - 1e-6;
                    WR.rho1 = 1.0; WR.rho2 = 1.0; WR.u = 2.0;  WR.p = 0.4; WR.z1 = 1 - 1e-6;
                    nG = 2;
                    n = 50 + 2*nG;
                    m = n;
                    L = 1.0;    //space domain 0 < x < 1
                    Tf = 0.50;
                    T0 = 0;
                    CFL = 0.4;
                    dx = L/(n-1-2*nG);
                    dy = L/(m-1-2*nG);
                    x0 = 0.5;
                    material1 = "air";
                    material2 = "air";
                    EoS = "ideal";
                    construction = "default2D";
                }
                if(CASE == "Test2_y"){
                    title = "Test2_y";
                    description = "Toro test 2 y-dimension equivalent";
                    datfile = "test2_y.dat";
                    gpname = "gifplot_test2_y.gp";
                    giffile = "./visualisation/test2_y.gif";
                    BCs = "bottom-reflective";
                    //exact_sol_file = "./ExactSols/Toro412";
                    WL.rho1 = 1.0;   WL.rho2 = 1.0;   WL.v = -2.0; WL.p = 0.4; WL.z1 = 1 - 1e-6;
                    WR.rho1 = 1.0; WR.rho2 = 1.0; WR.v = 2.0;  WR.p = 0.4; WR.z1 = 1 - 1e-6;
                    nG = 2;
                    n = 50 + 2*nG;
                    m = n;
                    L = 1.0;    //space domain 0 < x < 1
                    Tf = 0.50;
                    T0 = 0;
                    CFL = 0.4;
                    dx = L/(n-1-2*nG);
                    dy = L/(m-1-2*nG);
                    x0 = 0.5;
                    y0 = 0.5;
                    material1 = "air";
                    material2 = "air";
                    EoS = "ideal";
                    construction = "y-direction";
                }
                if(CASE == "Test3"){
                    title = "Test3";
                    description = "Toro test 3 equivalent";
                    datfile = "test3.dat";
                    gpname = "gifplot_test3.gp";
                    giffile = "./visualisation/test3.gif";
                    BCs = "transmissive";
                    //exact_sol_file = "./ExactSols/Toro413";
                    WL.rho1 = 1.0;   WL.rho2 = 1.0;   WL.u = 0.0; WL.p = 1000.0; WL.z1 = 1 - 1e-6;
                    WR.rho1 = 1.0; WR.rho2 = 1.0; WR.u = 0.0;  WR.p = 0.01; WR.z1 = 1 - 1e-6;
                    nG = 2;
                    n = 50 + 2*nG;
                    m = n;
                    L = 1.0;    //space domain 0 < x < 1
                    Tf = 0.012;
                    T0 = 0;
                    CFL = 0.4;
                    dx = L/(n-1-2*nG);
                    dy = L/(m-1-2*nG);
                    x0 = 0.5;
                    material1 = "air";
                    material2 = "air";
                    EoS = "ideal";
                    construction = "default2D";
                }
                if(CASE == "Test4"){
                    title = "Test4";
                    description = "Toro test 4 equivalent";
                    datfile = "test4.dat";
                    gpname = "gifplot_test4.gp";
                    giffile = "./visualisation/test4.gif";
                    BCs = "transmissive";
                    //exact_sol_file = "./ExactSols/Toro414";
                    WL.rho1 = 1.0;   WL.rho2 = 1.0;   WL.u = 0.0; WL.p = 0.01; WL.z1 = 1 - 1e-6;
                    WR.rho1 = 1.0; WR.rho2 = 1.0; WR.u = 0.0;  WR.p = 100.0; WR.z1 = 1 - 1e-6;
                    nG = 2;
                    n = 50 + 2*nG;
                    m = n;
                    L = 1.0;    //space domain 0 < x < 1
                    Tf = 0.035;
                    T0 = 0;
                    CFL = 0.4;
                    dx = L/(n-1-2*nG);
                    dy = L/(m-1-2*nG);
                    x0 = 0.5;
                    material1 = "air";
                    material2 = "air";
                    EoS = "ideal";
                    construction = "default2D";
                }
                if(CASE == "Test5"){
                    title = "Test5";
                    description = "Toro test 5 equivalent";
                    datfile = "euler_eqs_test5.dat";
                    gpname = "test5.gp";
                    giffile = "./visualisation/test5.gif";
                    BCs = "transmissive";
                    //exact_sol_file = "./ExactSols/Toro415";
                    WL.rho1 = 5.99924;   WL.rho2 = 5.99924;   WL.u = 19.5975; WL.p = 460.894; WL.z1 = 1 - 1e-6;
                    WR.rho1 = 5.99924; WR.rho2 = 5.99924; WR.u = -6.19633;  WR.p = 46.0950; WR.z1 = 1 - 1e-6;
                    nG = 2;
                    n = 50 + 2*nG;
                    m = n;
                    L = 1.0;    //space domain 0 < x < 1
                    Tf = 0.035;
                    T0 = 0;
                    CFL = 0.4;
                    dx = L/(n-1-2*nG);
                    dy = L/(m-1-2*nG);
                    x0 = 0.5;
                    material1 = "air";
                    material2 = "air";
                    EoS = "ideal";
                    construction = "default2D";
                }
                if(CASE == "Test1_2D"){
                    title = "Test1_2D";
                    description = "Toro test 1 2D equivalent";
                    datfile = "test1_2D.dat";
                    gpname = "gifplot_test1_2D.gp";
                    giffile = "./visualisation/test1_2D.gif";
                    BCs = "transmissive";
                    //exact_sol_file = "./ExactSols/Toro411";
                    WL.rho1 = 1.0;   WL.rho2 = 1.0;   WL.u = 0.0; WL.v = 0; WL.p = 1.0; WL.z1 = 1 - 1e-6;
                    WR.rho1 = 0.125; WR.rho2 = 0.125; WR.u = 0.0; WR.v = 0; WR.p = 0.1; WR.z1 = 1 - 1e-6;
                    nG = 2;
                    n = 50 + 2*nG;
                    m = 50 + 2*nG;
                    L = 1.0;    //space domain 0 < x < 1
                    Tf = 0.25;
                    T0 = 0;
                    CFL = 0.4;
                    dx = L/(n-1-2*nG);
                    dy = L/(m-1-2*nG);
                    x0 = 0.5;
                    y0 = 0.5;
                    material1 = "air";
                    material2 = "air";
                    EoS = "ideal";
                    construction = "square";
                }
                if(CASE == "Test2_2D"){
                    title = "Test2_2D";
                    description = "Toro test 2 2D equivalent";
                    datfile = "test2_2D.dat";
                    gpname = "gifplot_test2_2D.gp";
                    giffile = "./visualisation/test2_2D.gif";
                    BCs = "transmissive";
                    //exact_sol_file = "./ExactSols/Toro412";
                    WL.rho1 = 1.0;   WL.rho2 = 1.0;   WL.u = -2.0; WL.v = 2.0; WL.p = 0.4; WL.z1 = 1 - 1e-6;
                    WR.rho1 = 1.0; WR.rho2 = 1.0; WR.u = 2.0;  WR.v = -2.0; WR.p = 0.4; WR.z1 = 1 - 1e-6;
                    nG = 2;
                    n = 50 + 2*nG;
                    m = n;
                    L = 1.0;    //space domain 0 < x < 1
                    Tf = 0.15;
                    T0 = 0;
                    CFL = 0.4;
                    dx = L/(n-1-2*nG);
                    dy = L/(m-1-2*nG);
                    x0 = 0.5;
                    y0 = 0.5;
                    material1 = "air";
                    material2 = "air";
                    EoS = "ideal";
                    construction = "square";
                }
                if(CASE == "UnderwaterExplosion"){
                    title = "Underwater Explosion";
                    description = "2D inert underwater explosion test";
                    datfile = "underwater_explosion.dat";
                    gpname = "gifplot_underwater_explosion.gp";
                    giffile = "./visualisation/underwater_explosion.gif";
                    BCs = "bottom-reflective";
                    //exact_sol_file = "./ExactSols/Toro412";
                    //WL -> air, WR -> water
                    WL.rho1 = 1.225; WL.rho2 = 0.001; WL.u = 0.0; WL.v = 0.0; WL.p = 1.0e5; WL.z1 = 1 - 1e-6;
                    WR.rho1 = 0.001;  WR.rho2 = 1250; WR.u = 0.0; WR.v = 0.0; WR.p = 1.0e5; WR.z1 = 1e-6;
                    nG = 2;
                    n = 400 + 2*nG;
                    m = 250 + 2*nG;
                    L = 4.0;    //space domain 0 < x < 1
                    Y = 2.5;
                    Tf = 1.0e-3;
                    T0 = 0;
                    CFL = 0.4;
                    dx = L/(n-1-2*nG);
                    dy = Y/(m-1-2*nG);
                    x0 = 0.5;
                    y0 = 1.5;
                    material1 = "air";
                    material2 = "water";
                    EoS = "stiffened";
                    construction = "bubble";
                    o_x = 2;
                    o_y = 1.2;
                    R = 0.12;
                }
                if(CASE == "CircularExplosion"){
                    title = "Circular Explosion";
                    description = "2D inert circular gas explosion test";
                    datfile = "circular_explosion.dat";
                    gpname = "gifplot_circular_explosion.gp";
                    giffile = "./visualisation/circular_explosion.gif";
                    BCs = "transmissive";
                    //exact_sol_file = "./ExactSols/Toro412";
                    //WL -> air, WR -> water
                    WL.rho1 = 1.0; WL.rho2 = 1.0; WL.u = 0.0; WL.v = 0.0; WL.p = 1.0; WL.z1 = 1 - 1e-6;
                    WR.rho1 = 0.125; WR.rho2 = 0.125; WR.u = 0.0; WR.v = 0.0; WR.p = 0.1; WR.z1 = 1 - 1e-6;
                    nG = 2;
                    n = 200 + 2*nG;
                    m = 200 + 2*nG;
                    L = 2.0;    //space domain 0 < x < 1
                    Y = 2.0;
                    Tf = 0.25;
                    T0 = 0;
                    CFL = 0.4;
                    dx = L/(n-1-2*nG);
                    dy = Y/(m-1-2*nG);
                    x0 = 1.0;
                    y0 = 1.0;
                    material1 = "air";
                    material2 = "air";
                    EoS = "ideal";
                    construction = "circular";
                    o_x = 1.0;
                    o_y = 1.0;
                    R = 0.4;
                }
                if(CASE == "reflected_shock"){
                    title = "reflected_shock";
                    description = "reflected shock test for bottom wall BC";
                    datfile = "reflected_shock.dat";
                    gpname = "gifplot_reflected_shock.gp";
                    giffile = "./visualisation/reflected_shock.gif";
                    BCs = "bottom-reflective";
                    //exact_sol_file = "./ExactSols/Toro412";
                    WL.rho1 = 1.0;   WL.rho2 = 1.0;   WL.u = 0.0; WL.v = -3.0; WL.p = 0.4; WL.z1 = 1 - 1e-6;
                    WR.rho1 = 1.0; WR.rho2 = 1.0; WR.u = 0.0;  WR.v = -3.0; WR.p = 0.4; WR.z1 = 1 - 1e-6;
                    nG = 2;
                    n = 50 + 2*nG;
                    m = n;
                    L = 1.0;    //space domain 0 < x < 1
                    Tf = 0.15;
                    T0 = 0;
                    CFL = 0.4;
                    dx = L/(n-1-2*nG);
                    dy = L/(m-1-2*nG);
                    x0 = 0.5;
                    y0 = 0.5;
                    material1 = "air";
                    material2 = "air";
                    EoS = "ideal";
                    construction = "default2D";
                }
                if(CASE == "AirExplosion"){
                    title = "Air Explosion";
                    description = "2D inert air + water explosion test";
                    datfile = "air_explosion.dat";
                    gpname = "gifplot_air_explosion.gp";
                    giffile = "./visualisation/air_explosion.gif";
                    BCs = "all-reflective";
                    //exact_sol_file = "./ExactSols/Toro412";
                    //WL -> air, WR -> water
                    WL.rho1 = 1.225; WL.rho2 = 0.001; WL.u = 0.0; WL.v = 0.0; WL.p = 1.0e5; WL.z1 = 1 - 1e-6;
                    WR.rho1 = 0.001;  WR.rho2 = 1250; WR.u = 0.0; WR.v = 0.0; WR.p = 1.0e5; WR.z1 = 1e-6;
                    nG = 2;
                    n = 200 + 2*nG;
                    m = 125 + 2*nG;
                    L = 4.0;    //space domain 0 < x < 1
                    Y = 2.5;
                    Tf = 0.5e-3;
                    T0 = 0;
                    CFL = 0.4;
                    dx = L/(n-1-2*nG);
                    dy = Y/(m-1-2*nG);
                    x0 = 0.5;
                    y0 = 1.5;
                    material1 = "air";
                    material2 = "water";
                    EoS = "stiffened";
                    construction = "air-bubble";
                    o_x = 2;
                    o_y = 1.2;
                    R = 0.20;
                }
                if(CASE == "CavityCollapse"){
                    title = "CavityCollapse";
                    description = "2D inert single cavity collapse";
                    datfile = "cavity_collapse.dat";
                    gpname = "gifplot_cavity_collapse.gp";
                    giffile = "./visualisation/cavity_collapse.gif";
                    BCs = "transmissive";
                    //exact_sol_file = "./ExactSols/Toro412";
                    //WL -> air bubble, WR -> ambient water
                    WL.rho1 = 1.2; WL.rho2 = 1000; WL.u = 0.0; WL.v = 0.0; WL.p = 1.0e5; WL.z1 = 1 - 1e-6;
                    WR.rho1 = 1.2;  WR.rho2 = 1000; WR.u = 0.0; WR.v = 0.0; WR.p = 1.0e5; WR.z1 = 1e-6;
                    nG = 2;
                    n = 700 + 2*nG;
                    m = 500 + 2*nG;
                    L = 14.0;    //space domain 0 < x < 1
                    Y = 12.0;
                    Tf = 5.0e-3;
                    T0 = 0;
                    CFL = 0.4;
                    dx = L/(n-1-2*nG);
                    dy = Y/(m-1-2*nG);
                    x0 = 0.6;
                    y0 = 0.0;
                    material1 = "air";
                    material2 = "water";
                    EoS = "stiffened";
                    construction = "shocked-cavity";
                    o_x = 6.0;
                    o_y = 6.0;
                    R = 3.0;
                }
                if(CASE == "CavityCollapse_Nitromethane"){
                    title = "Cavity Collapse - Nitromethane";
                    description = "2D inert single air cavity collapse in Nitromethane";
                    datfile = "cavity_collapse_nitromethane.dat";
                    gpname = "gifplot_cavity_collapse_nitromethane.gp";
                    giffile = "./visualisation/cavity_collapse_nitromethane.gif";
                    BCs = "transmissive";
                    //exact_sol_file = "./ExactSols/Toro412";
                    //WL -> air bubble, WR -> ambient nitromethane
                    WL.rho1 = 1.2; WL.rho2 = 1134; WL.u = 0.0; WL.v = 0.0; WL.p = 1.0e5; WL.z1 = 1 - 1e-6;
                    WR.rho1 = 1.2;  WR.rho2 = 1134; WR.u = 0.0; WR.v = 0.0; WR.p = 1.0e5; WR.z1 = 1e-6;
                    nG = 2;
                    n = 600 + 2*nG;
                    m = 450 + 2*nG;
                    L = 14.0;    //space domain 0 < x < 1
                    Y = 12.0;
                    Tf = 2.7e-3;
                    T0 = 0;
                    CFL = 0.4;
                    dx = L/(n-1-2*nG);
                    dy = Y/(m-1-2*nG);
                    x0 = 0.6;
                    y0 = 0.0;
                    material1 = "air";
                    material2 = "nitromethane";
                    EoS = "CC";
                    construction = "shocked-cavity-nitromethane";
                    o_x = 6.0;
                    o_y = 6.0;
                    R = 3.0;
                }
                if(CASE == "ReflectedGelShock"){
                    title = "ReflectedGelShock";
                    description = "testing for post-shock conditions";
                    datfile = "ReflectedGelShock.dat";
                    gpname = "gifplot_ReflectedGelShock.gp";
                    giffile = "./visualisation/ReflectedGelShock.gif";
                    BCs = "bottom-reflective";
                    //exact_sol_file = "./ExactSols/Toro412";
                    //WL -> air bubble, WR -> ambient water
                    WL.rho1 = 1000; WL.rho2 = 1.2; WL.u = 0.0; WL.v = -50; WL.p = 1.0e5; WL.z1 = 1 - 1e-6;
                    WR.rho1 = 1000;  WR.rho2 = 1.2; WR.u = 0.0; WR.v = -50; WR.p = 1.0e5; WR.z1 = 1 - 1e-6;
                    nG = 2;
                    n = 100 + 2*nG;
                    m = 100 + 2*nG;
                    L = 1.0;    //space domain 0 < x < 1
                    Y = 1.0;
                    Tf = 3.0e-4;
                    T0 = 0;
                    CFL = 0.4;
                    dx = L/(n-1-2*nG);
                    dy = Y/(m-1-2*nG);
                    x0 = 0.5;
                    y0 = 0.5;
                    material1 = "water";
                    material2 = "water";
                    EoS = "stiffened";
                    construction = "y-direction";
                }
                if(CASE == "WeakShockCavityCollapse"){
                    title = "WeakShockCavityCollapse";
                    description = "2D inert single cavity collapse for \
                                   weak shock case in Bourne and Field's paper";
                    datfile = "weak_shock_cavity_collapse.dat";
                    gpname = "gifplot_weak_shock_cavity_collapse.gp";
                    giffile = "./visualisation/weak_shock_cavity_collapse.gif";
                    BCs = "transmissive";
                    //exact_sol_file = "./ExactSols/Toro412";
                    //WL -> air bubble, WR -> ambient water
                    WL.rho1 = 1.2; WL.rho2 = 1030; WL.u = 0.0; WL.v = 0.0; WL.p = 1.0e5; WL.z1 = 1 - 1e-6;
                    WR.rho1 = 1.2;  WR.rho2 = 1030; WR.u = 0.0; WR.v = 0.0; WR.p = 1.0e5; WR.z1 = 1e-6;
                    nG = 2;
                    n = 200 + 2*nG;
                    m = 200 + 2*nG;
                    L = 0.03;    //space domain 0 < x < 1
                    Y = 0.03;
                    Tf = 17.0e-6; //usually 0.0275e-3
                    T0 = 0;
                    CFL = 0.4;
                    dx = L/(n-1-2*nG);
                    dy = Y/(m-1-2*nG);
                    x0 = 0.004;
                    y0 = 0.0;
                    material1 = "air";
                    material2 = "gelatin10";
                    EoS = "stiffened";
                    construction = "shocked-cavity";
                    o_x = 0.012;
                    o_y = 0.015;
                    R = 0.006;
                    k0 = 1.0/R;
                    sigma = 0.0728; //0.0728;                 
                }
                if(CASE == "CircularSTtest"){
                    title = "Circular ST Test";
                    description = "2D circular ST test - expanding circular bubble";
                    datfile = "circular_ST.dat";
                    gpname = "gifplot_circular_ST.gp";
                    giffile = "./visualisation/circular_ST.gif";
                    BCs = "transmissive";
                    //exact_sol_file = "./ExactSols/Toro412";
                    //WL -> air, WR -> water
                    WL.rho1 = 1.0; WL.rho2 = 0.125; WL.u = 0.0; WL.v = 0.0; WL.p = 1.0; WL.z1 = 1-1e-6;
                    WR.rho1 = 1.0; WR.rho2 = 0.125; WR.u = 0.0; WR.v = 0.0; WR.p = 0.1; WR.z1 = 1e-6;
                    nG = 2;
                    n = 300 + 2*nG;
                    m = 300 + 2*nG;
                    L = 1.0;    //space domain 0 < x < 1
                    Y = 1.0;
                    Tf = 0.10;
                    T0 = 0;
                    CFL = 0.4;
                    dx = L/(n-1-2*nG);
                    dy = Y/(m-1-2*nG);
                    x0 = 1.0;
                    y0 = 1.0;
                    material1 = "air";
                    material2 = "air";
                    EoS = "ideal";
                    construction = "circular";
                    o_x = 0.5;
                    o_y = 0.5;
                    R = 0.2;
                    k0 = 1.0/R;
                    sigma = 0.06;
                }
                if(CASE == "DropletTest"){
                    title = "DropletTest";
                    description = "2D circular ST test - expanding circular bubble";
                    datfile = "DropletTest.dat";
                    gpname = "DropletTest.gp";
                    giffile = "./visualisation/DropletTest.gif";
                    BCs = "transmissive";
                    //exact_sol_file = "./ExactSols/Toro412";
                    //WR -> air bubble, WL -> ambient water
                    WL.rho1 = 1000; WL.rho2 = 1.2; WL.u = 0.0; WL.v = 0.0; WL.p = 1.36e5; WL.z1 = 1 - 1e-6;
                    WR.rho1 = 1000;  WR.rho2 = 1.2; WR.u = 0.0; WR.v = 0.0; WR.p = 1.0e5; WR.z1 = 1e-6;              
                    nG = 2;
                    n = 150 + 2*nG;
                    m = 150 + 2*nG;
                    L = 1.0;    //space domain 0 < x < 1
                    Y = 1.0;
                    Tf = 0.1;
                    T0 = 0;
                    CFL = 0.4;
                    dx = L/(n-1-2*nG);
                    dy = Y/(m-1-2*nG);
                    x0 = 1.0;
                    y0 = 1.0;
                    material1 = "fake-liquid";
                    material2 = "air";
                    EoS = "ideal";
                    construction = "elliptical";
                    o_x = 0.5;
                    o_y = 0.5;
                    R = 0.2;
                    k0 = 1.0/R;
                    sigma = 0.072;
                }

                if(CASE == "ST_Test1"){
                    title = "ST_Test1";
                    description = "Sod test";
                    datfile = "ST_test1.dat";
                    gpname = "gifplot_ST_test1.gp";
                    giffile = "./visualisation/ST_test1.gif";
                    BCs = "transmissive";
                    //exact_sol_file = "./ExactSols/Toro411";
                    WL.rho1 = 1.0; WL.rho2 = 0.125; WL.u = 0.0; WL.v = 0.0; WL.p = 1.0; WL.z1 = 1 - 1e-6;
                    WR.rho1 = 1.0; WR.rho2 = 0.125; WR.u = 0.0; WL.v = 0.0; WR.p = 0.1; WR.z1 = 1e-6;
                    nG = 2;
                    n = 200 + 2*nG;
                    m = 200 + 2*nG;
                    L = 1.0;    //space domain 0 < x < 1
                    Y = 1.0;
                    Tf = 0.20;
                    T0 = 0;
                    CFL = 0.4;
                    dx = L/(n-1-2*nG);
                    dy = Y/(m-1-2*nG);
                    x0 = 0.50;
                    material1 = "air";
                    material2 = "air";
                    EoS = "ideal";
                    construction = "default2D";
                    R = 1.0;
                    k0 = 1.0/R;
                    sigma = 0.1;
                }

                if(CASE == "ST_Test2"){
                    title = "ST_Test2";
                    description = "double direction Sod test";
                    datfile = "ST_test2.dat";
                    gpname = "gifplot_ST_test2.gp";
                    giffile = "./visualisation/ST_test2.gif";
                    BCs = "transmissive";
                    WL.rho1 = 1.0; WL.rho2 = 0.125; WL.u = 0.0; WL.v = 0.0; WL.p = 1.0; WL.z1 = 1 - 1e-6;
                    WR.rho1 = 1.0; WR.rho2 = 0.125; WR.u = 0.0; WL.v = 0.0; WR.p = 0.1; WR.z1 = 1e-6;
                    nG = 2;
                    n = 150 + 2*nG;
                    m = 150 + 2*nG;
                    L = 2.0;    //space domain 0 < x < 1
                    Y = 2.0;
                    Tf = 0.20; //0.20
                    T0 = 0;
                    CFL = 0.4;
                    dx = L/(n-1-2*nG);
                    dy = Y/(m-1-2*nG);
                    x0 = 0.50;
                    x1 = 1.50;
                    material1 = "air";
                    material2 = "air";
                    EoS = "ideal";
                    construction = "3part";
                    R = 1.0;
                    k0 = 1.0/R;
                    sigma = 0.1;
                }

                if(CASE == "CircularSTgel"){
                    title = "CircularSTgel";
                    description = "2D circular ST test - collapsing circular air bubble in gelatin";
                    datfile = "CircularSTgel.dat";
                    gpname = "CircularSTgel.gp";
                    giffile = "./visualisation/CircularSTgel.gif";
                    BCs = "transmissive";
                    //exact_sol_file = "./ExactSols/Toro412";
                    //WL -> air, WR -> gel
                    WL.rho1 = 1.2; WL.rho2 = 1031; WL.u = 0.0; WL.v = 0.0; WL.p = 1.0e5; WL.z1 = 1-1e-6;//air bubble
                    WR.rho1 = 1.2; WR.rho2 = 1031; WR.u = 0.0; WR.v = 0.0; WR.p = 14.6e5; WR.z1 = 1e-6; //gel
                    nG = 2;
                    n = 200 + 2*nG;
                    m = 200 + 2*nG;
                    L = 3.0e-3;    //space domain 0 < x < 1
                    Y = 3.0e-3;
                    Tf = 30.0e-6;
                    T0 = 0;
                    CFL = 0.4;
                    dx = L/(n-1-2*nG);
                    dy = Y/(m-1-2*nG);
                    material1 = "air";
                    material2 = "gelatin10";
                    EoS = "stiffened";
                    construction = "circular";
                    o_x = 1.5e-3;
                    o_y = 1.5e-3;
                    R = 1.0e-3;
                    k0 = 1.0/R;
                    sigma = 72.8;
                }

                if(CASE == "ExpandingAirBubble"){
                    title = "ExpandingAirBubble";
                    description = "2D circular ST test - expanding circular air bubble in water";
                    datfile = "ExpandingAirBubble.dat";
                    gpname = "ExpandingAirBubble.gp";
                    giffile = "./visualisation/ExpandingAirBubble.gif";
                    BCs = "transmissive";
                    //exact_sol_file = "./ExactSols/Toro412";
                    //WL -> air, WR -> water
                    WL.rho1 = 1.2; WL.rho2 = 1000; WL.u = 0.0; WL.v = 0.0; WL.p = 1.0e7; WL.z1 = 1-1e-6;//air bubble
                    WR.rho1 = 1.2; WR.rho2 = 1000; WR.u = 0.0; WR.v = 0.0; WR.p = 1.0e5; WR.z1 = 1e-6; //water
                    nG = 2;
                    n = 250 + 2*nG;
                    m = 250 + 2*nG;
                    L = 10.0e-3;    //space domain 0 < x < 1
                    Y = 10.0e-3;
                    Tf = 12.0e-6;
                    T0 = 0;
                    CFL = 0.4;
                    dx = L/(n-1-2*nG);
                    dy = Y/(m-1-2*nG);
                    material1 = "air";
                    material2 = "water";
                    EoS = "stiffened";
                    construction = "circular";
                    o_x = 5.0e-3;
                    o_y = 5.0e-3;
                    R = 1.0e-3;
                    k0 = 1.0/R;
                    sigma = 0.0;
                }

                if(CASE == "OscillatingBubble"){
                    title = "OscillatingBubble";
                    description = "bubble oscillating under ST affects";
                    datfile = "OscillatingBubble.dat";
                    gpname = "OscillatingBubble.gp";
                    giffile = "./visualisation/OscillatingBubble.gif";
                    BCs = "transmissive";
                    //exact_sol_file = "./ExactSols/Toro412";
                    //WR -> air bubble, WL -> ambient water
                    WL.rho1 = 100; WL.rho2 = 1.0; WL.u = 0.0; WL.v = 0.0; WL.p = 1.0e5; WL.z1 = 1 - 1e-4;
                    WR.rho1 = 100;  WR.rho2 = 1.0; WR.u = 0.0; WR.v = 0.0; WR.p = 1.0e5; WR.z1 = 1e-4;              
                    nG = 2;
                    n = 250 + 2*nG;
                    m = 250 + 2*nG;
                    L = 0.8;    //space domain 0 < x < 1
                    Y = 0.8;
                    Tf = 0.14;
                    T0 = 0;
                    CFL = 0.25;
                    dx = L/(n-1-2*nG);
                    dy = Y/(m-1-2*nG);
                    x0 = 1.0;
                    y0 = 1.0;
                    material1 = "fake-liquid";
                    material2 = "air";
                    EoS = "ideal";
                    construction = "elliptical";
                    o_x = 0.4;
                    o_y = 0.4;
                    R = 0.15825;
                    k0 = 1.0/R;
                    sigma = 341.64; 
 
                }

                //EXPERIMENTAL TEST CASES 1-9
                
                if(CASE == "CavityCollapse1"){
                    title = "CavityCollapse1";
                    description = "2D inert single cavity collapse";
                    datfile = "cavity_collapse1.dat";
                    gpname = "gifplot_cavity_collapse1.gp";
                    giffile = "./visualisation/cavity_collapse1.gif";
                    BCs = "transmissive";
                    //exact_sol_file = "./ExactSols/Toro412";
                    //WL -> air bubble, WR -> ambient water
                    WL.rho1 = 1.2; WL.rho2 = 1000; WL.u = 0.0; WL.v = 0.0; WL.p = 1.0e5; WL.z1 = 1 - 1e-6;
                    WR.rho1 = 1.2;  WR.rho2 = 1000; WR.u = 0.0; WR.v = 0.0; WR.p = 1.0e5; WR.z1 = 1e-6;
                    nG = 2;
                    n = 230 + 2*nG;
                    m = 200 + 2*nG;
                    L = 0.03;    //space domain 0 < x < 1
                    Y = 0.03;
                    Tf = 0.020e-3; //usually 0.0275e-3
                    T0 = 0;
                    CFL = 0.2;
                    dx = L/(n-1-2*nG);
                    dy = Y/(m-1-2*nG);
                    x0 = 0.004;
                    y0 = 0.0;
                    material1 = "air";
                    material2 = "water";
                    EoS = "stiffened";
                    construction = "shocked-cavity";
                    o_x = 0.012;
                    o_y = 0.015;
                    R = 0.006; 
                }
                if(CASE == "CavityCollapse2"){
                    title = "CavityCollapse2";
                    description = "2D inert single cavity collapse, with post shock pressure: 1.9 GPa";
                    datfile = "cavity_collapse2.dat";
                    gpname = "gifplot_cavity_collapse2.gp";
                    giffile = "./visualisation/cavity_collapse2.gif";
                    BCs = "transmissive";
                    //exact_sol_file = "./ExactSols/Toro412";
                    //WL -> air bubble, WR -> ambient water
                    WL.rho1 = 1.2; WL.rho2 = 1000; WL.u = 0.0; WL.v = 0.0; WL.p = 1.0e5; WL.z1 = 1 - 1e-6;
                    WR.rho1 = 1.2;  WR.rho2 = 1000; WR.u = 0.0; WR.v = 0.0; WR.p = 1.0e5; WR.z1 = 1e-6;
                    nG = 2;
                    n = 230 + 2*nG;
                    m = 200 + 2*nG;
                    L = 0.03;    //space domain 0 < x < 1
                    Y = 0.03;
                    Tf = 0.020e-3; //usually 0.0275e-3
                    T0 = 0;
                    CFL = 0.2;
                    dx = L/(n-1-2*nG);
                    dy = Y/(m-1-2*nG);
                    x0 = 0.004;
                    y0 = 0.0;
                    material1 = "air";
                    material2 = "water";
                    EoS = "stiffened";
                    construction = "shocked-cavity";
                    o_x = 0.012;
                    o_y = 0.015;
                    R = 0.006;
                }
                if(CASE == "CavityCollapse3"){
                    title = "CavityCollapse3";
                    description = "2D single cavity collapse, with post shock pressure: 0.26 GPa and D=6mm";
                    datfile = "cavity_collapse3.dat";
                    gpname = "gifplot_cavity_collapse3.gp";
                    giffile = "./visualisation/cavity_collapse3.gif";
                    BCs = "transmissive";
                    //exact_sol_file = "./ExactSols/Toro412";
                    //WL -> air bubble, WR -> ambient water
                    WL.rho1 = 1.2; WL.rho2 = 1030; WL.u = 0.0; WL.v = 0.0; WL.p = 1.0e5; WL.z1 = 1 - 1e-6;
                    WR.rho1 = 1.2;  WR.rho2 = 1030; WR.u = 0.0; WR.v = 0.0; WR.p = 1.0e5; WR.z1 = 1e-6;
                    nG = 2;
                    n = 230 + 2*nG;
                    m = 200 + 2*nG;
                    L = 0.015;    //space domain 0 < x < 1
                    Y = 0.015;
                    Tf = 12.0e-6; //usually 0.0275e-3
                    T0 = 0;
                    CFL = 0.2;
                    dx = L/(n-1-2*nG);
                    dy = Y/(m-1-2*nG);
                    x0 = 0.002;
                    y0 = 0.0;
                    material1 = "air";
                    material2 = "gelatin10";
                    EoS = "stiffened";
                    construction = "shocked-cavity";
                    o_x = 0.006;
                    o_y = 0.0075;
                    R = 0.003;
                }
                if(CASE == "CavityCollapse4"){
                    title = "CavityCollapse4";
                    description = "2D single cavity collapse, with post shock pressure: 0.5 GPa and D=6mm";
                    datfile = "cavity_collapse4.dat";
                    gpname = "gifplot_cavity_collapse4.gp";
                    giffile = "./visualisation/cavity_collapse4.gif";
                    BCs = "transmissive";
                    //exact_sol_file = "./ExactSols/Toro412";
                    //WL -> air bubble, WR -> ambient water
                    WL.rho1 = 1.2; WL.rho2 = 1000; WL.u = 0.0; WL.v = 0.0; WL.p = 1.0e5; WL.z1 = 1 - 1e-6;
                    WR.rho1 = 1.2;  WR.rho2 = 1000; WR.u = 0.0; WR.v = 0.0; WR.p = 1.0e5; WR.z1 = 1e-6;
                    nG = 2;
                    n = 200 + 2*nG;
                    m = 200 + 2*nG;
                    L = 0.015;    //space domain 0 < x < 1
                    Y = 0.015;
                    Tf = 8.0e-6; //usually 0.0275e-3
                    T0 = 0;
                    CFL = 0.2;
                    dx = L/(n-1-2*nG);
                    dy = Y/(m-1-2*nG);
                    x0 = 0.002;
                    y0 = 0.0;
                    material1 = "air";
                    material2 = "water";
                    EoS = "stiffened";
                    construction = "shocked-cavity";
                    o_x = 0.006;
                    o_y = 0.0075;
                    R = 0.003;
                }
                if(CASE == "CavityCollapse5"){
                    title = "CavityCollapse5";
                    description = "2D single cavity collapse, with post shock pressure: 1.9 GPa and D=6mm";
                    datfile = "cavity_collapse5.dat";
                    gpname = "gifplot_cavity_collapse5.gp";
                    giffile = "./visualisation/cavity_collapse5.gif";
                    BCs = "transmissive";
                    //exact_sol_file = "./ExactSols/Toro412";
                    //WL -> air bubble, WR -> ambient water
                    WL.rho1 = 1.2; WL.rho2 = 1000; WL.u = 0.0; WL.v = 0.0; WL.p = 1.0e5; WL.z1 = 1 - 1e-6;
                    WR.rho1 = 1.2;  WR.rho2 = 1000; WR.u = 0.0; WR.v = 0.0; WR.p = 1.0e5; WR.z1 = 1e-6;
                    nG = 2;
                    n = 230 + 2*nG;
                    m = 200 + 2*nG;
                    L = 0.015;    //space domain 0 < x < 1
                    Y = 0.015;
                    Tf = 5.0e-6; //usually 0.0275e-3
                    T0 = 0;
                    CFL = 0.2;
                    dx = L/(n-1-2*nG);
                    dy = Y/(m-1-2*nG);
                    x0 = 0.002;
                    y0 = 0.0;
                    material1 = "air";
                    material2 = "water";
                    EoS = "stiffened";
                    construction = "shocked-cavity";
                    o_x = 0.006;
                    o_y = 0.0075;
                    R = 0.003;
                }
                if(CASE == "CavityCollapse6"){
                    title = "CavityCollapse6";
                    description = "2D single cavity collapse, with post shock pressure: 3.5 GPa and D=6mm";
                    datfile = "cavity_collapse6.dat";
                    gpname = "CavityCollapse6.gp";
                    giffile = "./visualisation/CavityCollapse6.gif";
                    BCs = "transmissive";
                    //exact_sol_file = "./ExactSols/Toro412";
                    //WL -> air bubble, WR -> ambient water
                    WL.rho1 = 1.2; WL.rho2 = 1000; WL.u = 0.0; WL.v = 0.0; WL.p = 1.0e5; WL.z1 = 1 - 1e-6;
                    WR.rho1 = 1.2;  WR.rho2 = 1000; WR.u = 0.0; WR.v = 0.0; WR.p = 1.0e5; WR.z1 = 1e-6;
                    nG = 2;
                    n = 230 + 2*nG;
                    m = 200 + 2*nG;
                    L = 0.015;    //space domain 0 < x < 1
                    Y = 0.015;
                    Tf = 3.0e-6; //usually 0.0275e-3
                    T0 = 0;
                    CFL = 0.2;
                    dx = L/(n-1-2*nG);
                    dy = Y/(m-1-2*nG);
                    x0 = 0.002;
                    y0 = 0.0;
                    material1 = "air";
                    material2 = "water";
                    EoS = "stiffened";
                    construction = "shocked-cavity";
                    o_x = 0.006;
                    o_y = 0.0075;
                    R = 0.003;
                }
                if(CASE == "CavityCollapse7"){
                    title = "CavityCollapse7";
                    description = "2D single cavity collapse, with post shock pressure: 0.5 GPa and D=3mm";
                    datfile = "cavity_collapse7.dat";
                    gpname = "CavityCollapse7.gp";
                    giffile = "./visualisation/CavityCollapse7.gif";
                    BCs = "transmissive";
                    //exact_sol_file = "./ExactSols/Toro412";
                    //WL -> air bubble, WR -> ambient water
                    WL.rho1 = 1.2; WL.rho2 = 1000; WL.u = 0.0; WL.v = 0.0; WL.p = 1.0e5; WL.z1 = 1 - 1e-6;
                    WR.rho1 = 1.2;  WR.rho2 = 1000; WR.u = 0.0; WR.v = 0.0; WR.p = 1.0e5; WR.z1 = 1e-6;
                    nG = 2;
                    n = 230 + 2*nG;
                    m = 200 + 2*nG;
                    L = 0.006;    //space domain 0 < x < 1
                    Y = 0.006;
                    Tf = 5.0e-6; //usually 0.0275e-3
                    T0 = 0;
                    CFL = 0.2;
                    dx = L/(n-1-2*nG);
                    dy = Y/(m-1-2*nG);
                    x0 = 0.0006;
                    y0 = 0.0;
                    material1 = "air";
                    material2 = "water";
                    EoS = "stiffened";
                    construction = "shocked-cavity";
                    o_x = 0.0025;
                    o_y = 0.003;
                    R = 0.0015;
                }
                if(CASE == "CavityCollapse8"){
                    title = "CavityCollapse8";
                    description = "2D single cavity collapse, with post shock pressure: 1.9 GPa and D=3mm";
                    datfile = "cavity_collapse8.dat";
                    gpname = "CavityCollapse8.gp";
                    giffile = "./visualisation/CavityCollapse8.gif";
                    BCs = "transmissive";
                    //exact_sol_file = "./ExactSols/Toro412";
                    //WL -> air bubble, WR -> ambient water
                    WL.rho1 = 1.2; WL.rho2 = 1000; WL.u = 0.0; WL.v = 0.0; WL.p = 1.0e5; WL.z1 = 1 - 1e-6;
                    WR.rho1 = 1.2;  WR.rho2 = 1000; WR.u = 0.0; WR.v = 0.0; WR.p = 1.0e5; WR.z1 = 1e-6;
                    nG = 2;
                    n = 230 + 2*nG;
                    m = 200 + 2*nG;
                    L = 0.006;    //space domain 0 < x < 1
                    Y = 0.006;
                    Tf = 2.0e-6; //usually 0.0275e-3
                    T0 = 0;
                    CFL = 0.2;
                    dx = L/(n-1-2*nG);
                    dy = Y/(m-1-2*nG);
                    x0 = 0.0006;
                    y0 = 0.0;
                    material1 = "air";
                    material2 = "water";
                    EoS = "stiffened";
                    construction = "shocked-cavity";
                    o_x = 0.0025;
                    o_y = 0.003;
                    R = 0.0015;
                }
                if(CASE == "CavityCollapse9"){
                    title = "CavityCollapse9";
                    description = "2D single cavity collapse, with post shock pressure: 0.26 GPa and D=3mm";
                    datfile = "cavity_collapse9.dat";
                    gpname = "CavityCollapse9.gp";
                    giffile = "./visualisation/CavityCollapse9.gif";
                    BCs = "transmissive";
                    //exact_sol_file = "./ExactSols/Toro412";
                    //WL -> air bubble, WR -> ambient water
                    WL.rho1 = 1.2; WL.rho2 = 1000; WL.u = 0.0; WL.v = 0.0; WL.p = 1.0e5; WL.z1 = 1 - 1e-6;
                    WR.rho1 = 1.2;  WR.rho2 = 1000; WR.u = 0.0; WR.v = 0.0; WR.p = 1.0e5; WR.z1 = 1e-6;
                    nG = 2;
                    n = 200 + 2*nG;
                    m = 200 + 2*nG;
                    L = 0.006;    //space domain 0 < x < 1
                    Y = 0.006;
                    Tf = 12.0e-6; //usually 0.0275e-3
                    T0 = 0;
                    CFL = 0.4;
                    dx = L/(n-1-2*nG);
                    dy = Y/(m-1-2*nG);
                    x0 = 0.0006;
                    y0 = 0.0;
                    material1 = "air";
                    material2 = "water";
                    EoS = "stiffened";
                    construction = "shocked-cavity";
                    o_x = 0.0025;
                    o_y = 0.003;
                    R = 0.0015;
                }

                if(CASE == "CavityCollapseST"){
                    title = "CavityCollapseST";
                    description = "2D single cavity collapse, with post shock pressure: 3.3 MPa and D=6mm";
                    datfile = "cavity_collapseST.dat";
                    gpname = "CavityCollapseST.gp";
                    giffile = "./visualisation/CavityCollapseST.gif";
                    BCs = "transmissive";
                    //exact_sol_file = "./ExactSols/Toro412";
                    //WL -> air bubble, WR -> ambient water
                    WL.rho1 = 1.2; WL.rho2 = 1030; WL.u = 0.0; WL.v = 0.0; WL.p = 1.0e5; WL.z1 = 1 - 1e-6;
                    WR.rho1 = 1.2;  WR.rho2 = 1030; WR.u = 0.0; WR.v = 0.0; WR.p = 1.0e5; WR.z1 = 1e-6;
                    nG = 2;
                    n = 200 + 2*nG;
                    m = 200 + 2*nG;
                    L = 0.015;    //space domain 0 < x < 1
                    Y = 0.015;
                    Tf = 100.0e-6; //usually 0.0275e-3
                    T0 = 0;
                    CFL = 0.4;
                    dx = L/(n-1-2*nG);
                    dy = Y/(m-1-2*nG);
                    x0 = 0.002;
                    y0 = 0.0;
                    material1 = "air";
                    material2 = "gelatin10";
                    EoS = "stiffened";
                    construction = "shocked-cavity";
                    o_x = 0.006;
                    o_y = 0.0075;
                    R = 0.003;
                    sigma = 0;
                }

                if(CASE == "MicroBubble"){
                    title = "MicroBubble";
                    description = "1 micrometer bubble collapse, under very weak loading";
                    datfile = "MicroBubble.dat";
                    gpname = "MicroBubble.gp";
                    giffile = "./visualisation/MicroBubble.gif";
                    BCs = "transmissive";
                    //exact_sol_file = "./ExactSols/Toro412";
                    //WL -> air bubble, WR -> ambient water
                    WL.rho1 = 1.2; WL.rho2 = 1000; WL.u = 0.0; WL.v = 0.0; WL.p = 1.0e5; WL.z1 = 1 - 1e-6;
                    WR.rho1 = 1.2;  WR.rho2 = 1000; WR.u = 0.0; WR.v = 0.0; WR.p = 1.0e5; WR.z1 = 1e-6;
                    nG = 2;
                    n = 300 + 2*nG;
                    m = 300 + 2*nG;
                    L = 2.0e-6;    //space domain 0 < x < 1
                    Y = 2.0e-6;
                    Tf = 0.006e-6; //usually 0.0275e-3
                    T0 = 0;
                    CFL = 0.4;
                    x0 = 0.2e-6;
                    dx = L/(n-1-2*nG);
                    dy = Y/(m-1-2*nG);
                    material1 = "air";
                    material2 = "water";
                    EoS = "stiffened";
                    construction = "shocked-cavity";
                    o_x = 0.8e-6;
                    o_y = 1.0e-6;
                    R = 0.5e-6;
                    sigma = 0.00;
                }




            }
};



