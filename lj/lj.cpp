#include <iostream>
#include <vector>
#include <array>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <random>
#define PI 3.1415926

using namespace std;
typedef std::array<double, 3> cart_coord;
template<typename T> inline T POW3(T x){
    return x*x*x;
}
template<typename T> inline T SQR(T x){
    return (x*x);
}

class md_cell{
    std::vector<int> link;
    std::vector<std::vector<std::vector<int> > > cell;
    std::vector<cart_coord> force;
    std::vector<double> pot;
    static inline double potential_w(double r){
        double r6 = 1.0/POW3(SQR(r));
        return (24*r6*(1 - 2*r6));
    }
    inline double get_rho() const{
        return (np/POW3(box));
    }
public:
    double critical_r, box, rc, delg;
    double KE, TE, ecor, ecut, Tb;
    double dt, dt2, dt_2, dt_4, dt_8;
    int nl, np;
    double *xi, *vxi, *Q;
    cart_coord psum;
    std::vector<cart_coord> sites_cart, xm, trace, vel;
    std::vector<int> g_r;
    std::vector<double> pcf;
    std::size_t nhis;
    std::size_t ngr;
    int nc;
    md_cell (double, double, double, double, double, std::size_t, std::vector<cart_coord>, std::vector<cart_coord>);
    void move(bool);
    void new_cell();
    void force_calculation(bool);
    void chain(double);
    void bootstrap(double);
    void print_force();
    void sample();
    void ComputePCF();
    double ComputePressure(double, double);
    void nvemove(bool);
    void verlet();
    void f_ver(bool);
};

md_cell::md_cell (double T, double KE_, double c_r, double b, double r_c, std::size_t nh, std::vector<cart_coord> carts, std::vector<cart_coord> velocities){
    critical_r = c_r, box = b, rc = r_c, nhis = nh;
    KE = KE_;
    sites_cart = carts;
    xm = carts;
    vel = velocities;
    delg = box/(2.*nhis);
    ngr = 0;
    nc = ceil(box/rc);
    cout << "the number of cells is " << nc << endl;
    nl = 2;
    np = sites_cart.size();
    cout << "the number of particles is "<< np << endl;
    xi = new double[nl];
    vxi = new double[nl];
    Q = new double[nl];
    for(std::size_t i = 0; i < nl; i++){
        xi[i] = vxi[i] = 0.0;
    }
    Q[1] = 2e-6; 
    Q[0] = 3.*Q[1]*np;
    Tb = T;
    dt = 0.001;
    dt2 = dt*dt;
    dt_2 = 0.5*dt;
    dt_4 = 0.5*dt_2;
    dt_8 = 0.5*dt_4;
    for(std::size_t i = 0; i < np; i++){
        cart_coord f_i = {0.0, 0.0, 0.0};
        cart_coord t_i = {0.0, 0.0, 0.0};
        cart_coord xm_i = {0.0, 0.0, 0.0};
        force.push_back(f_i);
        trace.push_back(t_i);
        xm.push_back(xm_i);
    }
    for(std::size_t i = 0; i < nc; i++){
        std::vector<std::vector<int> > cx;
        for(std::size_t j = 0; j < nc; j++){
            std::vector<int> cy;
            for(std::size_t k = 0; k < nc; k++){
                cy.push_back(-1);
            }
            cx.push_back(cy);
        }
        cell.push_back(cx);
    }
    for(std::size_t i = 0; i < np; i++){
        link.push_back(-1);
        pot.push_back(0.0);
    }
    f_ver(false);
//    new_cell();
//    force_calculation(false);
    for(std::size_t i = 0; i < nhis; i++){
        g_r.push_back(0);
        pcf.push_back(0.0);
    }
}

void md_cell::new_cell(){
    for(std::size_t i = 0; i < nc; i++){
        for(std::size_t j = 0; j < nc; j++){
            for(std::size_t k = 0; k < nc; k++){
                cell[i][j][k] = -1;
            }
        }
    }
    for(std::size_t i = 0; i < sites_cart.size(); i++){
        link[i] = -1;
    }
    for(std::size_t i = 0; i < sites_cart.size(); i++){
        int cellx = sites_cart[i][0]/rc;
        int celly = sites_cart[i][1]/rc;
        int cellz = sites_cart[i][2]/rc;
        link[i] = cell[cellx][celly][cellz];
        cell[cellx][celly][cellz] = i;
    }
}
void md_cell::print_force(){
    for(std::size_t i = 0; i < sites_cart.size(); i++){
        cout<< force[i][0] << "," << force[i][1] << "," << force[i][2] << endl;
    }
}
void md_cell::f_ver(bool Switch){
    if(Switch) ngr++;
    for(int i = 0; i < np; i++){
        for(int m = 0; m < 3; m++){
            force[i][m] = 0.0;
        }
    }
    for(int i = 0; i < np - 1; i++){
        for(int j = i+1; j < sites_cart.size(); j++){
            cart_coord xr = {0., 0., 0.};
            double r2 = 0.;
            for(int m = 0; m < 3; m++){
                xr[m] = sites_cart[i][m] - sites_cart[j][m];
                if(xr[m] > box/2.) xr[m] -= box;
                else if(xr[m] < -box/2.) xr[m] += box;
                r2 += xr[m]*xr[m];
            }
            if( r2 < critical_r*critical_r){
                double r2i = 1/r2;
                double r6i = POW3(r2i);
                double ff = 48*r2i*r6i*(r6i - 0.5);
                for(int m = 0; m < 3; m++){
                    force[i][m] += xr[m]*ff;
                    force[j][m] -= xr[m]*ff;
                }
            }
            double r = sqrt(r2);
            if(Switch and r < box/2.){
                int ig = r/delg;
                g_r[ig] += 2;
            }
        }
    }
}

void md_cell::force_calculation(bool Switch){
    if(Switch) ngr++;
//    cart_coord pressure = {.0, .0, .0};
    for(std::size_t i = 0; i < sites_cart.size(); i++){
        int cellx = sites_cart[i][0]/rc;
        int celly = sites_cart[i][1]/rc;
        int cellz = sites_cart[i][2]/rc;
        pot[i] = 0.0;
        for(std::size_t m = 0; m < 3; m++){
            force[i][m] = 0.0;
        }
        for(int ic = cellx -1; ic <= cellx + 1; ic++){
            for(int jc = celly- 1; jc <= celly + 1; jc++){
                for(int kc = cellz - 1; kc <= cellz + 1; kc++){
                    cart_coord pos_tmp = {0., 0., 0.};
                    cart_coord distance = {0., 0., 0.};
                    std::array<int, 3> period = {0, 0, 0};
                    std::array<int, 3> index = {ic, jc, kc};
                    for(std::size_t m = 0; m < 3; m++){
                        if(index[m] < 0){
                            index[m] += nc;
                            period[m] = 1;
                        }
                        else if(index[m] > nc -1){
                            index[m] -= nc;
                            period[m] = 2;
                        }
                    }
                    int id = cell[index[0]][index[1]][index[2]];
                    while(id != -1){
                        double radius2 = 0.0;
                        for(std::size_t m = 0; m < 3; m++){
                            if(period[m] == 0){
                                pos_tmp[m] = sites_cart[id][m];
                            }
                            else if(period[m] == 1){
                                pos_tmp[m] = sites_cart[id][m] - box;
                            }
                            else if(period[m] == 2){
                                pos_tmp[m] = sites_cart[id][m] + box;
                            }
//                            cout << period[m] << floor((period[m] - 1)/double(nc - 2)) << endl;
                            distance[m] = sites_cart[i][m] - pos_tmp[m];
                            radius2 += distance[m]*distance[m];
                        }
                        if(radius2 < critical_r*critical_r and id != i){
                            double radius6 = POW3(radius2);
                            double radius6i = 1./radius6;
                            double uij = 4*radius6i*(radius6i - 1);
                            pot[i] += uij;
                            for(std::size_t m = 0; m < 3; m++){
                               force[i][m] += distance[m]*6*(uij+4*radius6i*radius6i)/radius2;
                            }
                        }
                        double r = sqrt(radius2);
                        if(Switch and r < box/2. and id != i){
                            int ig = r/delg;
                            g_r[ig] += 1;
                        }
                        id = link[id];
                    }
                }
            }
        }
    }
}
void md_cell::ComputePCF(){
    for(std::size_t i = 0; i < g_r.size(); i++){
        double vb = (POW3(i+1) - POW3(i))*POW3(delg);
        double nid = 4./3*PI*vb*get_rho();
        pcf[i] = g_r[i]/(ngr*np*nid);
    }
}

double md_cell::ComputePressure(double rmax, double T_curr){
    double sum = 0.0;
    int n = g_r.size();
    int np = sites_cart.size();
    for(int i = 0; i < n; i++){
        double ri = delg*(i+0.5);
        sum += potential_w(ri)*pcf[i]*SQR(ri);
    }
    sum *= (-np*PI/1.5*get_rho()*delg);

    sum += np*T_curr;

    double rm_3=1.0/POW3(rmax);
	sum += 16.0/3.0*PI*np*get_rho()*rm_3*(SQR(rm_3)/1.5-1);

    return (sum/POW3(box));
}
void md_cell::chain(double T){
    int N;
    double G1, G2, s;
    N = sites_cart.size();
    G2 = Q[0]*vxi[0]*vxi[0] - T;
    vxi[1] += G2*dt_4;
    vxi[0] += exp(-vxi[1]*dt_8);
    G1 = (2*KE - 3*N*T)/Q[0];
    vxi[0] += G1*dt_4;
    vxi[0] *= exp(-vxi[1]*dt_8);
    xi[0] += vxi[0]*dt_2;
    xi[1] += vxi[1]*dt_2;
    s = exp(-vxi[0]*dt_2);
    for(std::size_t i = 0; i < N; i++){
        for(std::size_t j = 0; j < 3; j++){
            vel[i][j] *=s;
        }
    }
    KE *=(s*s);
    vxi[0] *= exp(-vxi[1]*dt_8);
    G1 = (2*KE - 3*N*T)/Q[0];
    vxi[0] += G1*dt_4;
    vxi[0] *= exp(-vxi[1]*dt_8);
    G2 = (Q[0]*vxi[0]*vxi[0] - T)/Q[1];
    vxi[1] += G2*dt_4;
}

void md_cell::verlet(){
    psum[0] = 0.0, psum[1] = .0, psum[2] =  0.0;
    double sumv2 = 0.0;
    for(std::size_t i = 0; i < np; i++){
        for(std::size_t j = 0; j < 3; j++){
            double xx = 2*sites_cart[i][j] - xm[i][j] + SQR(dt)*force[i][j];
            xx -= box*floor(xx/box);
            double delta_x = xx - xm[i][j];
            if(delta_x > box/2.) delta_x -= box;
            if(delta_x < -box/2.) delta_x += box;
            vel[i][j] = delta_x/(2*dt);
            trace[i][j] += (sites_cart[i][j] - xm[i][j]);
            psum[j] += vel[i][j];
            sumv2 += SQR(vel[i][j]);
            xm[i][j] = sites_cart[i][j];
            sites_cart[i][j] = xx;
        }
    }
    KE = sumv2*0.5;
}
void md_cell::nvemove(bool Switch){
    new_cell();
    force_calculation(Switch);
    verlet();
}
   
void md_cell::move(bool Switch){
    chain(Tb);
    KE = 0.0;
    int N = sites_cart.size();
    for(std::size_t i = 0; i < N; i++){
        for(std::size_t j = 0; j < 3; j++){
            sites_cart[i][j] += vel[i][j]*dt_2;
            trace[i][j] += vel[i][j]*dt_2;
            sites_cart[i][j] -= box*floor(sites_cart[i][j]/box);
        }
    }
    f_ver(Switch);
//    new_cell();
//    force_calculation(Switch);
    for(std::size_t i = 0; i < N; i++){
        for(std::size_t j = 0; j < 3; j++){
            vel[i][j] += dt*force[i][j];
            sites_cart[i][j] += vel[i][j]*dt_2;
            trace[i][j] += vel[i][j]*dt_2;
            KE += vel[i][j]*vel[i][j];
        }
    }
    KE *= 0.5;
//    cout << KE << endl;
    chain(Tb);
}

void md_cell::sample(){
    cart_coord pressure = {0.0, 0.0, 0.0};
    std::size_t N = sites_cart.size();
    for(std::size_t i = 0; i < N - 1; i++){
        for(std::size_t j = i+1; j < N; j++){
            double radius = 0.0;
            for(std::size_t m = 0; m < 3; m++){
                radius += (sites_cart[i][m] - sites_cart[j][m])*(sites_cart[i][m] - sites_cart[j][m]);
            }
            double f = 4.0/(radius*radius*radius*radius)*(12./(radius*radius*radius) - 6);
            for(std::size_t m = 0; m < 3; m++){
                pressure[m] += (sites_cart[i][m] - sites_cart[j][m])*(sites_cart[i][m] - sites_cart[j][m])*f;
            }
        }
    }
    double v = box*box*box;
    cout << (pressure[0]+pressure[1]+pressure[2])/(3.0*v) + N/v*Tb << endl;
}

void md_cell::bootstrap(double T){
    cart_coord sumv = {0.0, 0.0, 0.0};
    double sumv2 = 0.0;
    int N = sites_cart.size();
    for(std::size_t i = 0; i < sites_cart.size(); i++){
        for(std::size_t j = 0; j < 3; j++){
//            double xx = 2*sites_cart[i][j] - xm[i][j] + force[i][j]*dt*dt;
            double delta_pos = vel[i][j]*dt + force[i][j]/2.0*dt*dt;
//            vel[i][j] = (xx - xm[i][j])/(2*dt);
//            xm[i][j] = sites_cart[i][j];
            xm[i][j] = sites_cart[i][j];
            sites_cart[i][j] += delta_pos;
            trace[i][j] += delta_pos;
            vel[i][j] += force[i][j]*dt;
            sumv[j] += vel[i][j];
            sumv2 += vel[i][j]*vel[i][j];
            sites_cart[i][j] -= floor(sites_cart[i][j]/box)*box;
        }
    }
    for(std::size_t i = 0; i < 3; i++){
        sumv[i] /= N;
    }
    sumv2 /= N;
    double sf = sqrt(3.*T/sumv2);
    sumv2 = 0.;
    for(std::size_t i = 0; i < N; i++){
        for(std::size_t j = 0; j < 3; j++){
            vel[i][j] = (vel[i][j] - sumv[j])*sf;
            sumv2 += vel[i][j]*vel[i][j];
//            xm[i][j] -= vel[i][j]*dt;
        }
    }
}

void init_md(std::vector<cart_coord> & sites_cart, std::vector<cart_coord> & velocities, std::size_t total_atom, double box, double T, double &KE){
    double sigma = sqrt(T);
    std::default_random_engine generator;
    std::normal_distribution<double> distribution(0.0, sigma);
    int n3 = 2;
    std::size_t counter = 0;
    while(n3*n3*n3 < total_atom) n3++;
    int iix = 0, iiy = 0, iiz = 0;
    for(std::size_t i = 0; i < total_atom; i++){
        cart_coord r = {0., 0., 0.};
        r[0] = (double(iix) + 0.5)*box/n3;
        r[1] = (double(iiy) + 0.5)*box/n3;
        r[2] = (double(iiz) + 0.5)*box/n3;
        sites_cart.push_back(r);
        iix++;
        if(iix == n3){
            iix = 0;
            iiy++;
            if(iiy == n3){
                iiy = 0;
                iiz++;
            }
        }
    }
    for(std::size_t i = 0; i < total_atom; i++){
        cart_coord vel = {0.0, 0.0, 0.0};
        for(std::size_t j = 0; j < 3; j++){
            vel[j] = distribution(generator);
        }
        velocities.push_back(vel);
    }
    double cmx = 0., cmy = 0., cmz = 0., v2_sum = 0.;
    for(std::size_t i = 0; i < total_atom; i++){
        cmx += velocities[i][0];
        cmy += velocities[i][1];
        cmz += velocities[i][2];
        for(std::size_t j = 0; j < 3; j++){
            v2_sum += velocities[i][j]*velocities[i][j];
        }
    }
    v2_sum /= double(total_atom);
    double sf = sqrt(3.0*T/v2_sum);
    KE = 0.;
    for(std::size_t i = 0; i < total_atom; i++){
        velocities[i][0] -= cmx/double(total_atom);
        velocities[i][1] -= cmy/double(total_atom);
        velocities[i][2] -= cmz/double(total_atom);
        velocities[i][0] *= sf;
        velocities[i][1] *= sf;
        velocities[i][2] *= sf;
        KE += velocities[i][0]*velocities[i][0] + velocities[i][1]*velocities[i][1] + velocities[i][2]*velocities[i][2];
    }
    KE *= 0.5;
}


//void md_cell::sample(){
//    ngr++;
//    for(std::size_t i = 0; i < sites_cart.size(); i++){
//        for(std::size_t k = 0; k < vlist[i].size(); k++){
//            std::size_t j = vlist[i][k];
//            radius = sqrt(radius);
//            if(radius <= rv){
//                int ig = int(radius/delg);
//                g_r[ig] += 1;
//            }
//        }
//    }
//}

int main(int argc, char* argv[]){
    double total_atom, box;
    std::size_t ncycl;
    for(std::size_t i = 1; i < argc; i++){
        if(strcmp(argv[i], "atom") == 0){
            i++;
            total_atom = atof(argv[i]);
        }
        else if(strcmp(argv[i], "box") == 0){
            i++;
            box = atof(argv[i]);
        }
        else if(strcmp(argv[i], "ncycl") == 0){
            i++;
            ncycl = atof(argv[i]);
        }
    }
//    box = 10;
//    total_atom = 800;
//    ncycl = 10000;
    double critical_r = 2.5;
    double delg = 0.02;
    double rc = 5.0;
    double KE = 0.0;
    double Tb = 2.0;
    double sum_pressure = 0.0;
    bool Switch = false;
    std::size_t nsampl = 10;
    std::size_t nhis = box/(2*delg);
    std::vector<cart_coord> sites_cart;
    std::vector<cart_coord> velocities;
    std::vector<double> gr;
    init_md(sites_cart, velocities, total_atom, box, Tb, KE);
    cout << "initial KE is "<< KE << endl;
    md_cell md(Tb, KE, critical_r, box, rc, nhis, sites_cart, velocities);
    for(std::size_t i = 0; i < 1; i++){
        md.bootstrap(Tb);
        md.f_ver(false);
//        md.new_cell();
//        md.force_calculation(false);
    }
    int counter = 0;
    for(std::size_t i = 0; i < ncycl; i++){
        if( i%nsampl == 0 && i > 200) Switch = true;
        else Switch = false;
//        md.nvemove(Switch);
        md.move(Switch);
        double trace2 = 0.;
        for(std::size_t j = 0; j < md.np; j++){
            trace2 += md.trace[j][0]*md.trace[j][0] + md.trace[j][1]*md.trace[j][1] + md.trace[j][2]*md.trace[j][2];
        }
        trace2 /= md.np;
        double Tc = 2.*md.KE/(3.*md.np);
        cout << i << " " << trace2 << endl;
        if( i%nsampl == 0 && i > 200000){
            md.ComputePCF();
            double Tc = 2.*md.KE/(3.*md.np);
//            double pres = md.ComputePressure(md.box/2., Tb);
            double pres = md.ComputePressure(md.box/2., Tc);
//            cout << Tc << " " << pres << endl;
            sum_pressure += pres;
            counter++;
        }
//        if( i < 300 ){
//            for(std::size_t m = 0; m < 3; m++){
//                md.psum[m] /= md.np;
//            }
//            double sf = sqrt(3.*Tb/(2.*md.KE));
//            for(std::size_t j = 0; j < md.np; j++){
//                for(std::size_t k = 0; k < 3; k++){
//                    md.vel[j][k] = (md.vel[j][k] - md.psum[k])*sf;
//            xm[i][j] -= vel[i][j]*dt;
//                }
//            }
//        }
    }
//    for(std::size_t i = 0; i < md.g_r.size(); i++){
//        double vb = (POW3(i+1) - POW3(i))*POW3(md.delg);
//        double vb = (i+0.5)*(i+0.5)*md.delg*md.delg*md.delg;
//        double nid = 4./3*PI*vb*sites_cart.size()/POW3(box);
//        gr.push_back(md.g_r[i]/(md.ngr*total_atom*nid));
//    }
//    for(std::size_t i = 0; i < nhis; i++){
//       cout << i*md.delg << " " << md.pcf[i] << endl;
//    }
//    cout << sum_pressure/counter << endl;
}
