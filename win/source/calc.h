#if !defined(__calc_h)
#define __calc_h

#include <cmath>
#include <vector>
#include <iostream>

#define M_2PI 6.2831853071795865
#define M_2PIL 6.2831853071795864769
#define GRID 128
#define INFTY 1e300

#if defined(WIN32)
void sincos(const double &x, double *s, double *c){
  *s = sin(x);
  *c = cos(x);
}
#endif

/* allocate a double matrix with subscript range m[nrl..nrh][ncl..nch] */
template <class T> T** matrix(int nrow, int ncol) {

        T **m = new T* [nrow];
        m[0] = (T *) new char [nrow*ncol*sizeof(T)];
        for( int i = 1; i < nrow; i++) m[i] = m[i-1]+ ncol;

        return m;       
}

/* free a double matrix allocated by dmatrix() */
template <class T> void free_matrix(T **m) {
        delete [] (char*) m[0];
        delete [] m;
}

/*
   Structure for storing information about a Fourier mode

   val*( (1-s)*cos(n*x)+ s*sin(n*x) )
*/

struct Tmode {

  int s,            // sc in {0,1}    
      n;            // n in {0,1, 2,...} 

  double val;
};


/*
  Two dimensional point    
*/

struct Tpoint {
  double x,y;
  Tpoint(const double &x, const double &y):x(x),y(y){}
  Tpoint(double v[2]):x(v[0]),y(v[1]) {}
};

/*

  Considering billards with inner boundary given with

    vec{r}_1(phi) = r(phi) (cos(phi), sin(phi)) phi in [0, 2 pi]

  and outer boundary

    vec{r}_1(phi) = vec{r}_1(phi) + d n(phi)

    n(phi) -- normal vector on the curve vec{r}_1(phi) of unit size
*/

#define OPT(x,y)(std::abs(1 - ((x)[0]*(y)[0] + (x)[1]*(y)[1])))

class Tbill{
  
  int M;

  double d;


  int cof_size;
  Tmode *cof, *cof_end;
  
  int path_size;
  double *path;

  struct Tresult {
    int side;

    double phi, len, x[2];
    
    Tresult()
    :side(0), phi(0), len(INFTY) {
      x[0]=x[1]=0;
    }
    void store(int side, const double & phi, const double &len, double x[2]){
      this->side = side;
      this->phi = phi;
      this->len = len;
      this->x[0] = x[0]; this->x[1] = x[1];
    }
  };

  /*
    Distance of the torus [0,2PI], shortest arc between angles x,y 
    devided by 2 Pi
  */

  double dist(const double & x, const double & y)  { 
    double t = (x - y)/M_2PI;

    if (t > 1) t -= int(t);
    if (t < 0.0) t -= int(t)-1;

    return std::min(t, 1 - t);
    //return M_2PI*std::min(t,1-t);
  }

  public:

  int status; // status = 0  -- OK
              // status = 1  -- INSIDE-TO-WALL
              // status = 2  -- WALL-TO-WALL

  Tbill(double d, std::vector<Tmode> & coef)
    : d(d)
  {
    M = 0; 
    
    cof = new Tmode [cof_size = coef.size()];
    
    for (int i = 0; i < cof_size; ++i) {
      cof[i].s = coef[i].s;
      cof[i].n = coef[i].n;
      cof[i].val = coef[i].val;
      M = std::max(M, cof[i].n);
    }
    
    cof_end = cof + cof_size;
    path_size = 0;
  }
  
  ~Tbill(){
    delete [] cof;
    if (path_size) delete [] path;
  }
  
  double get_width() const {
    return d;
  }
  
  void init_path(){

    if (path_size != 0) return;

    int N = GRID*(M + 1);

     // number of points in path  = (maximal mode number M)  x factor + 1
    path = new double [path_size = N+1];

    path[0] = 0.0;

    double phi = 0, dphi = M_2PI/N, 
           r, dr, sum;
    
    for (int i = 0; i < N; ++i, phi += dphi) {

      // Boole's rule integration of 5 order, error of the 6 order
      rad(phi, r, dr);            sum = 7*sqrt(r*r + dr*dr);
      rad(phi + 0.2*dphi, r, dr); sum += 32*sqrt(r*r + dr*dr);
      rad(phi + 0.4*dphi, r, dr); sum += 12*sqrt(r*r + dr*dr);
      rad(phi + 0.8*dphi, r, dr); sum += 32*sqrt(r*r + dr*dr);
      rad(phi + dphi, r, dr);     sum += 7*sqrt(r*r + dr*dr);

      path[i+1] = dphi*sum/90 + path[i];
    }
  }
  /*
    The length arc in the angle interval [0,phi]
  */
  double s(const double &phi) {

    if (path_size == 0 ) init_path();

    double theta = phi*(path_size - 1)/M_2PI;

    int n = int(theta);

    if (n >= int(path_size - 1)) return path[n];

    theta -= n;  
    
    return (1 - theta)*path[n] + theta*path[n+1];
  }
  
  /*
    Gives the point on the curve
  */  
  void curve(const int& side, const double &phi, double x[2]){

    double r, nx, ny;

    r = rad(phi);
    sincos(phi, x+1, x);
        
    x[0] *= r;
    x[1] *= r;
        
    if (side){
      normal(phi, nx, ny);
      x[0] += d*nx;
      x[1] += d*ny;
    } 
  }


  /*
    Gives the point on the curve x and pointer from x0 -> x stored in p
    input: side, phi, x0
    output: x, p, len
  */  
  void pointer(const int& side, const double &phi, double x0[2], 
              double x[2], double p[2], double &len){

    double r;

    sincos(phi, x+1, x);
        
    x[0] *= (r = rad(phi));
    x[1] *= r;
        
    if (side){
      double nx, ny;

      normal(phi, nx, ny);

      x[0] += d*nx;
      x[1] += d*ny;
    }

    p[0] = x[0] - x0[0];
    p[1] = x[1] - x0[1];

    
    p[0] /= (len = sqrt(p[0]*p[0]+p[1]*p[1]));
    p[1] /= len; 
  }


  /*
    Radius as a function of angle
  */
  double rad(const double &phi) {

    double sum = 0.0;
    
    for (Tmode *it = cof; it != cof_end; ++it)
      if ((*it).s) 
        sum += (*it).val*sin((*it).n*phi);
      else 
        sum += (*it).val*cos((*it).n*phi);

    return sum;
  }

  /*
    Radius as a function of angle
  */
  void rad(const double &phi, double &r, double &dr) {

    double s, c;

    r = dr = 0.0;
    
    for (Tmode *it = cof; it != cof_end; ++it) {
    
      sincos((*it).n*phi, &s, &c);

      if ((*it).s) {
        r += (*it).val*s;
        dr += (*it).n*(*it).val*c;
     } else {
        r += (*it).val*c;
        dr -= (*it).n*(*it).val*s;     
     }
    }
     
  }

  /*
    Normal vector to the inner boundary
  */
  void normal(const double &phi, double &nx, double &ny) {

    double r, dr, s, c;

    rad(phi, r, dr);

    sincos(phi, &s, &c);
    
    nx = dr*s + r*c;
    ny = -dr*c + r*s;
 
    nx /= (r = sqrt(nx*nx + ny*ny));
    ny /= r;
  }

/*
    Normal vector to side-th boundary, side  in {0,1}
  */
  void normal(const int &side, const double &phi, double &nx, double &ny) {

    double r, dr, s, c;

    rad(phi, r, dr);

    sincos(phi, &s, &c);
    
    nx = dr*s + r*c;
    ny = -dr*c + r*s;
 
    nx /= (r = sqrt(nx*nx + ny*ny));
    ny /= r;
    
    if (side) {
      nx = -nx;
      ny = -ny;
    }
  }

  /* Plotting the shape of the billiard with m equally distributed points in
     gnuplot format: First outputing the inner boundary and then outer 
      x_1(phi_i) y_1(phi_i)     i = 0
      ....    
      x_1(phi_i) y_1(phi_i)     i = m-1
      (empty line)
      x_2(phi_i) y_2(phi_i)     i = 0
      ....    
      x_2(phi_i) y_2(phi_i)     i = m-1
  */
  void plot_shape(const int m, std::ostream & os) {

    int i;

    double dphi = M_2PI/m,  phi = 0, 
           r, s, c, nx,ny;

    for (i = 0; i <= m; ++i, phi += dphi) {
      r = rad(phi);
      sincos(phi, &s, &c);
      os  <<  r*c << '\t' << r*s << '\n';
    } 

    os << '\n';
  
    phi = 0;
    for (i = 0; i <= m; ++i, phi += dphi) {
      r = rad(phi);
      sincos(phi, &s, &c);
      normal(phi, nx, ny);
      os  <<  r*c + d*nx << '\t' << r*s + d*ny << '\n';
    } 

  } 
  
  void plot_shape(const int m, Tpoint **curve) {

    double dphi = M_2PI/m, phi = 0, 
           r, s, c, nx, ny;

    for (int i = 0; i <= m; ++i, phi += dphi) {
      r = rad(phi);
      sincos(phi, &s, &c);
      normal(phi, nx, ny);
      curve[1][i].x = (curve[0][i].x = r*c) + d*nx;
      curve[1][i].y = (curve[0][i].y = r*s) + d*ny;
    }
  } 

  void plot_shape(const int m, std::vector<Tpoint> curve[2]) {

    double dphi = M_2PI/m, phi = 0, 
           r, s, c, nx, ny;

    for (int i = 0; i <= m; ++i, phi += dphi) {
      r = rad(phi);
      sincos(phi, &s, &c);
      normal(phi, nx, ny);
      c *= r; s *= r;
      nx *=d; ny *=d;

      curve[0].push_back(Tpoint(c, s));
      curve[1].push_back(Tpoint(c + nx, s + ny));
    }
  } 

  /*
    Trajectories starting on the cross section
  */
  void cross_section(const double  &phi, const double & l, double x[2]) {
    double r = rad(phi), nx, ny;
    
    sincos(phi, x+1, x);
    normal(phi, nx, ny);
    x[0] *=r;
    x[1] *=r;
    x[0] += l*nx;
    x[1] += l*ny;
  }

  /*
    Calculate the intersection of a the line
      x0 + v0 t, t > 0
    starting in on one bondary of the billiard (side0) and 
    ending on the boundary (side0). 
    
    Input: side0 -- side of the boundary
           phi0 -- angle from the  origin 
           x0[2] -- starting point on the boundary
           v0[2] -- starting direction of travel
    
    Output: side -- side of the boundary
            phi -- angle from the  origin
            x1[2] -- ending point on the boundary
            v1[2] -- new direction of travel
    Return: time        
  */  
  #define EPS_PHI 1e-12
  #define EPS_OPT 1e-11
  #define EPS_ERR 1e-12
  //#define DEBUG
  double walltowall(int side0, const double & phi0, double x0[2], double v0[2],
                    int & side1, double & phi1, double x1[2], double v1[2]) {
        
    #if defined(DEBUG)
    std::cerr << "walltowall:begin\n";
    #endif

    bool ok;

    int dir, i, j, N = GRID*(M + 1);
  
    double nx, ny, len1, len0,
           ph, dph = M_2PI/N, dph_ = 1e-2*dph, 
           xr[2], y[2], e[2];   

    Tresult result;
    
    e[1] = len1 = 1;

    normal(phi0, nx, ny);

    // determine direction of travel w.r.t. to the angle axis
    if (-v0[0]*ny + v0[1]*nx > 0) { 
      dir = 1; 
    } else {
      dir =-1;
      dph = -dph;
    }
    
    for (i = 0; i < 2; ++i) {

      // move initial angle a bit away from starting point
      ok = true;
      ph = phi0 + dir*EPS_PHI;
      
      for (j = 0; j < N && ok; ++j, ph += dph) {
        // calculate the error
        len0 = len1;
        pointer(i, ph, x0, xr, y, len1);

        e[0] = e[1];        
        e[1] = v0[1]*y[0] - v0[0]*y[1]; 

        // detect if the interval contains the solution
        // the first solution away from phi0 is already correct one
        if (j > 0 && e[0]*e[1] <= 0.0 ) {

          double ph_,
                 ph0 = ph - dph, 
                 ph1 = ph,
                 err0 = e[0], 
                 err1 = e[1],
                 err, len; 

          #if defined(DEBUG)
          std::cerr << "Hit: i = " << i << " j=" << j 
                    << " phi0 =" << ph0 << " len0 =" << len0
                    << " phi1 =" << ph1 << " len1 =" << len1  << '\n'
                    << "e0 ="    << e[0] << " e1 =" << e[0] << '\n';
          #endif

      
          if (std::abs(err0) < EPS_ERR && j > 1) {

            #if defined(DEBUG)
            std::cerr << "1: err0=" << err0 << '\n';
            #endif

            normal(i, ph0, nx, ny);
            pointer(i, ph0, x0, xr, y, len0);
                
            if (OPT(v0, y) < EPS_OPT && v0[0]*nx + v0[1]*ny < 0.0)  {
              if (len0 < result.len) result.store(i, ph0, len0, xr);
              ok = false;
              break;
            }   
          } else if (std::abs(err1) < EPS_ERR) {
            #if defined(DEBUG)
            std::cerr << "2: err1=" << err1 << '\n';
            #endif

            normal(i, ph1, nx, ny);
            pointer(i, ph1, x0, xr, y, len1);

            if (OPT(v0,y) < EPS_OPT && v0[0]*nx + v0[1]*ny < 0.0) {
               if (len1 < result.len) result.store(i, ph1, len1, xr);
               ok = false;
               break;
            }

          } else {
            #if defined(DEBUG)
            std::cerr << "3:\n";
            #endif
            do {
              #if defined(DEBUG)
              std::cerr << "(" << ph0<< ':'<< err0 << "), (" 
                               << ph1<< ':'<< err1 << ")\n";
              #endif                                              
              ph_ = 0.5*(ph0 + ph1);  
              pointer(i, ph_, x0, xr, y, len);

              err = v0[1]*y[0] - v0[0]*y[1];

              if (std::abs(err) < EPS_ERR) {  
                normal(i, ph_, nx, ny);

                if (OPT(v0,y) < EPS_OPT && v0[0]*nx + v0[1]*ny < 0.0) {
                  if (len < result.len) result.store(i, ph_, len, xr);
                  ok = false;
                }  
                break;
              }
            
              if (err*err0 < 0) {
                err1 = err;
                ph1 = ph_;
              } else if (err*err1 < 0) {
                err0 = err;
                ph0 = ph_;
              } else {
                std::cerr << "This should not happen\n";
                exit(EXIT_FAILURE);
              }

              if (std::abs(ph0 - ph1) < dph_ && v0[0]*y[0] + v0[1]*y[1] < 0)
                break; 

              if (std::abs(ph0 - ph1) < EPS_PHI) {
                ph_ = 0.5*(ph1 + ph0);
                normal(i, ph_, nx, ny);
                pointer(i, ph_, x0, xr, y, len);
                
                if (OPT(v0,y) < EPS_OPT && v0[0]*nx + v0[1]*ny < 0.0) {
                  if (len < result.len) result.store(i, ph_, len, xr);
                  ok = false;
                }
                break;                
              }
            } while (1);
          }   
        }
      }
    }

    if (result.len == INFTY){
      status = 2;
      return 0;
    }
    
    #if defined(DEBUG)
    std::cerr << "Processing result\n";
    #endif
   
    side1 = result.side;
    phi1 = result.phi;    
    
    // making modulus w.r.t 2*PI of outputing angles
    long double phi_ = phi1;
    phi_ /= M_2PIL;
    
    if (phi_ > 1) 
      phi_ -= int(phi_);
    else if ( phi_ < 0) 
      phi_ -= int(phi_)-1;
      
    phi1 = (phi_ *= M_2PIL);  
    
    x1[0] = result.x[0];
    x1[1] = result.x[1];
     
    normal(side1, phi1, nx, ny);
        
    double fac = 2*(nx*v0[0] + ny*v0[1]);
    v1[0] = v0[0] - fac*nx;
    v1[1] = v0[1] - fac*ny;

    #if defined(DEBUG)
    std::cerr << "walltowall:end\n";
    #endif

    return result.len;
  }
  #undef DEBUG
  
  /*
    Calculate the intersection between a the line
      x0 + v0 t, t> 0
    starting in inside the billiard and the boundary (side). 
    
    Input: x0[2] 
           v0[2]
    
    Output: side -- side of the boundary
            phi -- angle from the  origin
            x1[2] -- point on the boundary
            v1[2] -- new direction of travel
    Return: time
  */

  #define EPS_PHI 1e-12
  #define EPS_ERR 1e-12
  //#define DEBUG
  double insidetowall(double x0[2], double v0[2], 
                      int &side, double &phi, double x1[2], double v1[2]) {
    #if defined(DEBUG)
    std::cerr << "insidetowall:begin\n";
    #endif
    
    int i, j, N = GRID*(M + 1);
  
    double nx, ny, len0, len1,
           ph, dph = M_2PI/N, dph_ = 1e-2*dph, 
           y[2], e[2], xr[2];

    Tresult result;
    
    e[1] = len1 = 1;

    for (i = 0; i < 2; ++i) {
      
      ph = 0;      
      for (j = 0; j < N; ++j, ph += dph) {
        // calculate the error
        len0 = len1;
        pointer(i, ph, x0, xr, y, len1);

        e[0] = e[1];
        e[1] = v0[1]*y[0] - v0[0]*y[1]; 

        // detect if the interval contains the solution
        if (j > 0 && e[0]*e[1] <= 0.0 ) {
      
               // research first in more detail 
          double ph_,
                 ph0 = ph - dph, 
                 ph1 = ph,
                 err0 = e[0], 
                 err1 = e[1], 
                 err, len; 
        
          #if defined(DEBUG)    
          std::cerr << "Hit: i = " << i   << " j="     << j 
                    << " phi0 ="  << ph0  << " len0 =" << len0
                    << " phi1 ="  << ph1  << " len1 =" << len1  << '\n'
                    << " e0 ="    << e[0] << " e1 ="   << e[1] << '\n';
          #endif          
        
          if (std::abs(err0) < EPS_ERR) {
            #if defined(DEBUG)
            std::cerr << "1: err0=" << err0 << '\n';
            #endif

            normal(i, ph0, nx, ny);
            pointer(i, ph0, x0, xr, y, len0);
                
            if (len0 < result.len &&
                OPT(v0,y) < EPS_OPT && v0[0]*nx + v0[1]*ny < 0.0) 
              result.store(i, ph0, len0, xr);   
            
            continue;
          } else if (std::abs(err1) < EPS_ERR) {
            #if defined(DEBUG)
            std::cerr << "2: err1=" << err1 << '\n';
            #endif

            normal(i, ph1, nx, ny);
            pointer(i, ph1, x0, xr, y, len1);
                
            if (len1 < result.len &&
                OPT(v0,y) < EPS_OPT && v0[0]*nx + v0[1]*ny < 0.0)
              result.store(i, ph1, len1, xr); 

            continue;
          } else {
            #if defined(DEBUG)
            std::cerr << "3 fine polish\n";
            #endif
            do {
              #if defined(DEBUG)
              std::cerr << "(" << ph0<< ':'<< err0 << "), (" 
                               << ph1<< ':'<< err1 << ")\n";
              #endif                                              

              ph_ = 0.5*(ph0 + ph1);
              pointer(i, ph_, x0, xr, y, len);

              err = v0[1]*y[0] - v0[0]*y[1];

              if (std::abs(err) < EPS_ERR) {
                normal(i, ph_, nx, ny);

                if (len < result.len && 
                    OPT(v0,y) < EPS_OPT && v0[0]*nx + v0[1]*ny < 0.0) 
                  result.store(i, ph_, len, xr); 
                
                break;
              }
            
              if (err*err0 < 0) {
                err1 = err;
                ph1 = ph_;
              } else if (err*err1 < 0) {
                err0 = err;
                ph0 = ph_;
              } else {
                std::cerr << "This should not happen\n";
                exit(EXIT_FAILURE);
              }

              if (std::abs(ph0 - ph1) < dph_ && v0[0]*y[0] + v0[1]*y[1] < 0)
                break;

              if (std::abs(ph0 - ph1) < EPS_PHI) {
                ph_ = 0.5*(ph1 + ph0);

                normal(i, ph_, nx, ny);
                pointer(i, ph_, x0, xr, y, len);
              
                if (len < result.len &&
                    OPT(v0,y) < EPS_OPT && v0[0]*nx + v0[1]*ny < 0.0 )
                  result.store(i, ph_, len, xr);
                break;
              }
            } while (1);
          }
        }  
      }
    }
    
    if (result.len == INFTY){
      status = 1;
      return 0;
    }

    #if defined(DEBUG)
    std::cerr << "Processing result: size=" << results.size() << '\n';
    #endif
   
    // calculating the solution    

    side = result.side;
    phi = result.phi;
    x1[0] = result.x[0];
    x1[1] = result.x[1];
     
    normal(side, phi, nx, ny);
        
    double fac = 2*(nx*v0[0] + ny*v0[1]);
    v1[0] = v0[0] - fac*nx;
    v1[1] = v0[1] - fac*ny;

    #if defined(DEBUG)
    std::cerr << "insidetowall:end\n";
    #endif
    
    return result.len;
  }
  #undef DEBUG
  #undef EPS_PHI
  #undef EPS_ERR
  

  /*
    Calculate the intersection between a the line
      x0 + v0 t, t> 0
    starting on a cross-section at phi0, more precisely at length l, 
    and the boundary (side). 
    
    Input: x0[2] 
           v0[2]
    
    Output: phi0 -- angle of the cross-section
            x0[2] -- initial point in Cartesian coordinates
            v0[2] -- initial direction of travel
            x1[2] -- final point on the boundary in Cartesian coordinates
            v1[2] -- new direction of travel

    Return: time
  */
  #define EPS_PHI 1e-12
  #define EPS_ERR 1e-12
  //#define DEBUG

  double cross_sectiontowall(double phi0, double x0[2], double v0[2], 
                             int &side1, double &phi1, double x1[2], double v1[2]){
    
    #if defined(DEBUG)
    std::cerr << "walltowall:begin\n";
    #endif

    bool ok;

    int dir, i, j, N = GRID*(M + 1);
  
    double nx, ny, len1, len0,
           ph, dph = M_2PI/N, dph_ = 1e-2*dph, 
           y[2], e[2], xr[2];   
    
    Tresult result;
    
    e[1] = len1 = 1;

    normal(phi0, nx, ny);

    // determine direction of travel w.r.t. to the angle axis
    if (-v0[0]*ny + v0[1]*nx > 0) { 
      dir = 1; 
    } else {
      dir =-1;
      dph = -dph;
    }
    
    for (i = 0; i < 2; ++i) {

      // move initial angle a bit away from starting point
      ok = true;
      ph = phi0 + dir*EPS_PHI;
      for (j = 0; j < N && ok; ++j, ph += dph) {
        // calculate the error
        len0 = len1;
        pointer(i, ph, x0, xr, y, len1);

        e[0] = e[1];
        e[1] = v0[1]*y[0] - v0[0]*y[1]; 

        // detect if the interval contains the solution
        // the first solution away from phi0 is already correct one
        if (j > 0 && e[0]*e[1] <= 0.0 ) {

          double ph_,
                 ph0 = ph - dph, 
                 ph1 = ph,
                 err0 = e[0], 
                 err1 = e[1],
                 err, len; 

          #if defined(DEBUG)
          std::cerr << "Hit: i = " << i << " j=" << j 
                    << " phi0 =" << ph0 << " len0 =" << len0
                    << " phi1 =" << ph1 << " len1 =" << len1  << '\n'
                    << " e0 ="   << e[0] << " e1 =" << e[1] << '\n';
          #endif

      
          if (std::abs(err0) < EPS_ERR && j > 1) {
            #if defined(DEBUG)
            std::cerr << "1: err0=" << err0 << '\n';
            #endif

            normal(i, ph0, nx, ny);
            pointer(i, ph0, x0, xr, y, len0);
                
            if (OPT(v0,y) < EPS_OPT  && v0[0]*nx + v0[1]*ny < 0.0) {
              if (len0 < result.len) result.store(i, ph0, len0, xr);
              ok = false;
              break;
            }
          } else if (std::abs(err1) < EPS_ERR) {
            #if defined(DEBUG)
            std::cerr << "2: err1=" << err1 << '\n';
            #endif

            normal(i, ph1, nx, ny);
            pointer(i, ph1, x0, xr, y, len1);
                
            if (OPT(v0,y) < EPS_OPT  &&  v0[0]*nx + v0[1]*ny < 0.0) {
              if (len1 < result.len) result.store(i, ph1, len1, xr);
              ok = false;
              break;
            }
          } else {
            #if defined(DEBUG)
            std::cerr << "3:\n";
            #endif
            do {
              #if defined(DEBUG)
              std::cerr << "(" << ph0<< ':'<< err0 << "), (" 
                               << ph1<< ':'<< err1 << ")\n";
              #endif                                              
              ph_ = 0.5*(ph0 + ph1);
              pointer(i, ph_, x0, xr, y, len);
              
              err = v0[1]*y[0] - v0[0]*y[1];

              if (std::abs(err) < EPS_ERR) {
                normal(i, ph_, nx, ny);        

                if (OPT(v0,y) < EPS_OPT &&  v0[0]*nx + v0[1]*ny < 0.0 ) {
                  if (len < result.len) result.store(i, ph_, len, xr);
                  ok = false;
                }  
                break;
              }
            
              if (err*err0 < 0) {
                err1 = err;
                ph1 = ph_;
              } else if (err*err1 < 0) {
                err0 = err;
                ph0 = ph_;
              } else {
                std::cerr << "This should not happen\n";
                exit(EXIT_FAILURE);
              }

              if (std::abs(ph0 - ph1) < dph_ && v0[0]*y[0] + v0[1]*y[1] < 0)
                break; 

              if (std::abs(ph0 - ph1) < EPS_PHI) {
                ph_ = 0.5*(ph1 + ph0);
                normal(i, ph_, nx, ny);           
                pointer(i, ph_, x0, xr, y, len);

                if (OPT(v0,y) < EPS_OPT && v0[0]*nx + v0[1]*ny < 0.0) {
                  if (len < result.len) result.store(i, ph_, len, xr);
                  ok = false;
                }
                break;
              }
            } while (1);
          }   
        }
      }
    }

    if (result.len == INFTY){
      status = 2;
      return 0;
    }
    
    #if defined(DEBUG)
    std::cerr << "Processing result: size=" << results.size() << '\n';
    #endif
   
    // calculating the solution    
    side1 = result.side;
    phi1 = result.phi;
    
    // making modulus w.r.t 2*PI of outputing angles
    long double phi_ = phi1;
    phi_ /= M_2PIL;
    
    if (phi_ > 1) 
      phi_ -= int(phi_);
    else if ( phi_ < 0) 
      phi_ -= int(phi_)-1;
      
    phi1 = (phi_ *= M_2PIL);  
    
    x1[0] = result.x[0];
    x1[1] = result.x[1];
     
    normal(side1, phi1, nx, ny);
        
    double fac = 2*(nx*v0[0] + ny*v0[1]);
    v1[0] = v0[0] - fac*nx;
    v1[1] = v0[1] - fac*ny;

    #if defined(DEBUG)
    std::cerr << "walltowall:end\n";
    #endif

    return result.len;


  }
  #undef DEBUG
  #undef EPS_PHI
  #undef EPS_ERR


  #define EPS_WALL 1e-10
  //#define DEBUG
  
  void plot_traj(double phi, double l, double p[2], 
                 double t_end, std::vector<Tpoint> *traj = 0, std::ostream *os = 0){

    #if defined(DEBUG)
    std::cerr << "#plot_traj:start\n"; 
    #endif

    // reset the status
    status = 0;

    int side0, side1;

    double t, fac,
           V = sqrt(p[0]*p[0] + p[1]*p[1]),
           v0[2] = {p[0]/V, p[1]/V}, v1[2],
           x0[2], x1[2], 
           phi0 = phi, phi1;

    cross_section(phi0, l, x0);
    
    t_end *= V;     
    if (traj) traj->push_back(Tpoint(x0));
     
    if (os)
      (*os) << x0[0] << '\t' << x0[1] << '\t' 
           << v0[0] << '\t' << v0[1] << '\n';

    if (std::abs(l)< EPS_WALL) {
      side0 = 0;
      t = walltowall(0, phi0, x0, v0, side1, phi1, x1, v1);
    } else if (std::abs(l-d) < EPS_WALL) {
      side0 = 1;
      t = walltowall(1, phi0, x0, v0, side1, phi1, x1, v1);
    } else  {     
      side0 = -1;
      t = cross_sectiontowall(phi0, x0, v0, side1, phi1, x1, v1);
    }
             
    if (t > t_end) {
      fac = t_end/t; 
      x1[0] = x0[0] + (x1[0] - x0[0])*fac;
      x1[1] = x0[1] + (x1[1] - x0[1])*fac;
      
      if (traj) traj->push_back(Tpoint(x1));
      
      if (os)
        (*os) << x1[0] << '\t' << x1[1] << '\t' 
              << v1[0] << '\t' << v1[1] << '\t'
              << side1 << '\t' << phi1  << '\t' 
              << t << '\n';
         
      return;
    }
    
    
    double dt;
    
     while (status == 0) {

      side0 = side1; phi0 = phi1;
      x0[0] = x1[0]; x0[1] = x1[1]; 
      v0[0] = v1[0]; v0[1] = v1[1];

      if (traj) traj->push_back(Tpoint(x0));
     
      if (os)
        (*os) << x0[0] << '\t' << x0[1] << '\t' 
              << v0[0] << '\t' << v0[1] << '\t'
              << side0 << '\t' << phi0  << '\t' 
              << t << '\n';

      t += (dt = walltowall(side0, phi0, x0, v0, side1, phi1, x1, v1));
      
      if (t > t_end) {
        fac = (t - t_end)/dt; 
        x1[0] = x0[0] + (x1[0] - x0[0])*fac;
        x1[1] = x0[1] + (x1[1] - x0[1])*fac;
        break;
      } else if (t == t_end) break;
      
    }
    
    if (traj) traj->push_back(Tpoint(x1));
     
    if (os)
      (*os) << x1[0] << '\t' << x1[1] << '\t' 
            << v1[0] << '\t' << v1[1] << '\t'
            << side1 << '\t' << phi1  << '\t' 
            << t << '\n';

    #if defined(DEBUG)
    std::cerr << "#plot_traj:end\n"; 
    #endif
  }
  #undef EPS_WALL
  #undef DEBUG

  /*
    Poincare mapping on the crossection at phi:
    
    J=(l,v_l)  -> (l,v_l)' = J'     v_l = cos(angle)
    
    where l in [0,1] is relative coordinate along the crossection and 
          v_l in [-1,1]  is projection of the speed on the crossection. 
    The total speed is of unit size.

    Input:  sign in {-1, 1}  -- sign of the travel  .
            phi in [0, 2pi]  -- angle of the cross-section
            J0[2] = (l, v_l) -- the initial point
      
    Output: J1[2] = (l, v_l) -- the final point

    Returns: length

    NOTE: The particle on the initial cross-section can not be 
          shot over the initial cross-section.
  */
  #define EPS_D 1e-10
  #define EPS_WALL 1e-10
  //#define DEBUG
  double poincare_map(int sign, double phi, double J0[2], double J1[2],
                      std::ostream *os = 0) {
    
    #if defined(DEBUG)
    std::cerr << "#poincare_map:start\n";    
    #endif

    // reset the status
    status = 0;

    int side0, side1;

    double v0[2], v1[2], 
           x0[2], x1[2], 
           phi0 = phi, phi1,
           t, fac, dt,  
           nx, ny;
           
    cross_section(phi0, J0[0], x0);

    normal(phi0, nx, ny);

    fac = sign*sqrt(1.0 - J0[1]*J0[1]);

    v0[0] = nx*J0[1] - fac*ny;
    v0[1] = ny*J0[1] + fac*nx;
   
    if (os) 
      (*os) << x0[0] << '\t' << x0[1] << '\t' 
            << v0[0] << '\t' << v0[1] << '\n';

    if (std::abs(J0[0])< EPS_WALL) {
      side0 = 0;
      t = walltowall(0, phi0, x0, v0, side1, phi1, x1, v1);
    } else if (std::abs(J0[0] - d)< EPS_WALL) {
      side0 = 1;
      t = walltowall(1, phi0, x0, v0, side1, phi1, x1, v1);
    } else  {     
      side0 = -1;
      t = cross_sectiontowall(phi0, x0, v0, side1, phi1, x1, v1);
    }
    
    while (status == 0){
      side0 = side1; phi0 = phi1;
      x0[0] = x1[0]; x0[1] = x1[1]; 
      v0[0] = v1[0]; v0[1] = v1[1];
    
      if (os)  
        (*os) << x0[0] << '\t' << x0[1] << '\t' 
              << v0[0] << '\t' << v0[1] << '\t'
              << side0 << '\t' << phi0  << '\n';
     
      t += (dt = walltowall(side0, phi0, x0, v0, side1, phi1, x1, v1));

      // is crossing the initial cross-section?
      if (status == 0 &&
          std::abs(dist(phi0, phi)+dist(phi,phi1)-dist(phi0, phi1)) < EPS_D) {

        t -= dt;
        #if defined(DEBUG)
        std::cerr << "#phi0="  << phi0
                  << " phi1=" << phi1
                  << " phi="  << phi 
                  << '\n'; 
        #endif
        double r = rad(phi), s, c, dx, dy, det;

        sincos(phi, &s, &c);
        normal(phi, nx, ny);

        dx = x0[0] - r*c;
        dy = x0[1] - r*s;
        det = ny*v0[0] - nx*v0[1];

        J1[0] = (dy*v0[0] - dx*v0[1])/det;

        if (J1[0] < 0 || J1[0] > d) {
          status = 3;
          return 0.0;
        }

        J1[1] = v0[0]*nx + v0[1]*ny;
        
        t += (dt = (dy*nx - dx*ny)/det);

        if (dt < 0) {
          status = 4;
          return 0.0;
        }

        break;
      }
      
    }
    
    if (os)  
      (*os) << x1[0] << '\t' << x1[1] << '\t' 
            << v1[0] << '\t' << v1[1] << '\t'
            << side1 << '\t' << phi1  << '\n';

    #if defined(DEBUG)
    std::cerr << "#poincare_map:end\n";    
    #endif

    return t; 
  }
  #undef DEBUG
  #undef EPS_WALL
  #undef EPS_D
  
};
#endif // if !defined(__calc_h)
