
#include "calc.h"
#include "MersenneTwister.h"

using std::cin;

int main(int argc, char **argv) {


  if (argc != 3 && argc != 5) {
    std::cerr 
      << "calc s <m> < <cfg> # plotting shape\n"
      << "     p <m> < <cfg> # plotting path\n"
      << "     S <phi:l> <px:py> <t> # plotting the trajectory\n"
      << "     P <phi:sign> <N:seed> <n> # Poincare portraits\n"

      << "note:\n"
      << "  m - number of points per a plotted curve\n"
      << "  phi:l - initial point on a cross-section\n"
      << "  px:py - initial moment of the particle\n"
      << "  t - length of trajectory in time\n"
      << "  N:seed - number of samples points, seed for the random generator\n"
      << "  n - length of trajectory in number of iterations\n"
      << "  phi:sign - crossection and direction of travel\n";
    
    exit(EXIT_FAILURE);
  }

  
  double d;

  cin >> d;

  Tmode a;
  std::vector <Tmode> c;

  while (cin >> a.s >> a.n >> a.val) c.push_back(a);

  Tbill bill(d, c);

  if (argv[1][0] == 's'){
  
    int m = atoi(argv[2]);

    std::cout << "#m="<< m << '\n'
              << "#coef.size=" <<  c.size() << '\n';

    #if 0
    bill.plot_shape(m, std::cout);
    #else
    
    int i,j;
    Tpoint **curve = matrix<Tpoint>(2,m+1);

    bill.plot_shape(m, curve);

    for (i = 0; i < 2; ++i) {
      for (j = 0; j <= m; ++j)
        std::cout << curve[i][j].x << '\t' << curve[i][j].y << '\n';
      std::cout << '\n';
    }
    
    free_matrix <Tpoint> (curve);
    
    #endif

  } else if (argv[1][0] == 'p'){

    int m = atoi(argv[2]);

    double dphi = M_2PI/m,
           phi = 0;

    std::cout << "#m="<< m << '\n'
              << "#coef.size=" <<  c.size() << '\n';

    for (int i = 0; i <= m; ++i, phi += dphi)
      std::cout << phi << '\t' << bill.s(phi) << '\n';


  } else if (argv[1][0] == 'S') { 

    double phi, l, p[2], t = atof(argv[4]);

    sscanf(argv[2], "%lf:%lf", &phi, &l);
    sscanf(argv[3], "%lf:%lf", p,p+1);

    bill.plot_traj(phi, l, p, t, 0, &std::cout);

    switch (bill.status){
      case 1: std::cerr << "#There was an error in INSIDE-TO-WALL\n"; break;
      case 2: std::cerr << "#There was an error in WALL-TO-WALL\n"; break;
    }


  } else if (argv[1][0] == 'P'){

    int i, j, sign, N, n;
    
    unsigned long seed;

    double phi, t, J0[2], J1[2];
    
    sscanf(argv[2], "%lf:%d", &phi, &sign);
    sscanf(argv[3], "%d:%lu", &N, &seed);

    n = atoi(argv[4]);

    MTRand rnd(seed);

    #define BORDER 1e-4  
    for (i = 0; i < N; ++i) {
      J0[0] = (d-2*BORDER)*rnd.rand53() + BORDER;
      J0[1] = (2-2*BORDER)*rnd.rand53() - 1+BORDER;

      std::cout << J0[0] << '\t' << J0[1] << '\t' << 0 << '\n';
      
      for (j = 0; j < n; ++j){
      //  t = bill.poincare_map(sign, phi, J0, J1, &std::cerr);
        t = bill.poincare_map(sign, phi, J0, J1);

        if (bill.status != 0){ 
          switch (bill.status){
            case 1: 
              std::cerr << "#There was an error in INSIDE-TO-WALL\n"; 
            break;
            case 2: 
              std::cerr << "#There was an error in WALL-TO-WALL\n"; 
            break;
            case 3: 
              std::cerr << "#Length along Out-of-bounds " << J1[0] << '\n'; 
            break;
            case 4: 
              std::cerr << "#Time is negative \n";
            break;
          }
          std::cerr << "#Point=" << J0[0]<< ':' << J0[1] << '\n'; 
          break;    
        } else {
          J0[0] = J1[0]; J0[1] = J1[1];
          std::cout << J0[0] << '\t' << J0[1] << '\t' << t << '\n'; 
        }
      }
     }
    #undef BORDER
  }


  return EXIT_SUCCESS;
}
