#include <iostream>
#include <vector>
#include <array>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "mt19937.h"

#define PI 3.1415926

using namespace std;
typedef std::array<double, 3> cart_coord;

bool ener(std::vector<cart_coord> & sites_cart, double critical_r, double box){
    std::size_t total_atom = sites_cart.size();
    int k = 1;
//    for(int x = -k; x <= k; x++){
//        for(int y = -k; y <= k; y++){
//            for(int z = -k; z <= k; z++){
                for(std::size_t i = 0; i < total_atom; i++){
                    for(std::size_t j = 0; j < i; j++){
                        double radius = 0.0;
                        cart_coord new_cart = {0,0,0};
                        cart_coord new_cart_1 = sites_cart[i];
                        cart_coord new_cart_2 = sites_cart[j];
                        for(std::size_t k = 0; k < 3; k++){
                            new_cart_1[k] -= box*int(floor(new_cart_1[k]/box));
                            new_cart_2[k] -= box*int(floor(new_cart_2[k]/box));
                        }
                        new_cart[0] = new_cart_1[0] - new_cart_2[0];
                        new_cart[1] = new_cart_1[1] - new_cart_2[1];
                        new_cart[2] = new_cart_1[2] - new_cart_2[2];
                        radius = new_cart[0]*new_cart[0] + new_cart[1]*new_cart[1] + new_cart[2]*new_cart[2];
                        if(radius < critical_r*critical_r){
                            return false;
                        }
                    }
                }
    return true;
}
class mc_verlet{
    std::vector<std::vector<int> > vlist;
    std::vector<cart_coord> xv;
public:
    double critical_r, box, del_x, rv, delg;
    std::vector<cart_coord> sites_cart;
    std::vector<int> g_r;
    std::size_t nhis;
    std::size_t ngr;
    mc_verlet (double, double, double, double, std::size_t, std::vector<cart_coord>);
    bool move();
    void new_vlist();
    bool ener(std::size_t );
    double radius_square(cart_coord, cart_coord);
    void optimize_delx();
    void sample();
};

mc_verlet::mc_verlet (double c_r, double b, double d_x, double r_v, std::size_t nh, std::vector<cart_coord> carts){
    critical_r = c_r, box = b, del_x = d_x, rv = r_v, nhis = nh;
    sites_cart = carts;
    delg = rv/nhis;
    ngr = 0;
    new_vlist();
    for(std::size_t i = 0; i < nhis; i++){
        g_r.push_back(0);
    }
}

void mc_verlet::new_vlist(){
    xv = sites_cart;
    if(vlist.size() != 0){
        vlist.erase(vlist.begin(), vlist.end());
    }
    for(std::size_t i = 0; i < sites_cart.size(); i++){
        std::vector<int> vlist_i;
        vlist.push_back(vlist_i);
    }
    for(std::size_t i = 0; i < sites_cart.size(); i++){
//        std::vector<int> vlist_i, vlist_j;
        for(std::size_t j = 0; j < i; j++){
            double radius = 0.0;
            cart_coord xr = {0, 0, 0};
            for(std::size_t k = 0; k < 3; k++){
                xr[k] = sites_cart[i][k] - sites_cart[j][k];
                if( xr[k] > box/2 ){
                    xr[k] -= box;
                }
                else if( xr[k] < -box/2 ){
                    xr[k] += box;
                }
            }
            radius = xr[0]*xr[0] + xr[1]*xr[1] + xr[2]*xr[2];
            if(radius < rv*rv){
                vlist[i].push_back(j);
                vlist[j].push_back(i);
            }
        }
    }
}

double mc_verlet::radius_square(cart_coord site1, cart_coord site2){
    double radius = 0.0;
    cart_coord xr = {0, 0, 0};
    for(std::size_t k = 0; k < 3; k++){
        xr[k] = site1[k] - site2[k];
        if( xr[k] > box*0.5){
            xr[k] -= box;
        }
        else if( xr[k] < -box*0.5 ){
            xr[k] += box;
        }
    }
    radius = xr[0]*xr[0] + xr[1]*xr[1] + xr[2]*xr[2];
    return radius;
}

bool mc_verlet::move(){
    std::size_t total_atom = sites_cart.size();
    std::size_t o = genrand_real2()*total_atom;
    double radius = radius_square(sites_cart[o], xv[o]);
//    for(std::size_t i = 0; i < 3; i++){
//        radius += (sites_cart[o][i] - xv[o][i])*(sites_cart[o][i] - xv[o][i]);
//    }
    if(radius > (rv - critical_r)*(rv - critical_r)/4){
        new_vlist();
    }
    cart_coord old_site = sites_cart[o];
    cart_coord new_site = sites_cart[o];
    new_site[0] += (genrand_real1() - 0.5)*del_x;
    new_site[1] += (genrand_real1() - 0.5)*del_x;
    new_site[2] += (genrand_real1() - 0.5)*del_x;
    for(std::size_t i = 0; i < 3; i++){
        new_site[i] -= box*int(floor(new_site[i]/box));
    }
    sites_cart[o] = new_site;
    radius = radius_square(sites_cart[o], xv[o]);
//    for(std::size_t i = 0; i < 3; i++){
//        radius += (new_site[i] - xv[o][i])*(new_site[i] - xv[o][i]);
//    }
    if(radius > (rv - critical_r)*(rv - critical_r)/4){
        new_vlist();
    }
    if(ener(o)){
        return true;
    }
    else{
        sites_cart[o] = old_site;
        return false;
    }
}
   
bool mc_verlet::ener(std::size_t o){
    for(std::size_t j = 0; j < vlist[o].size(); j++){
        double radius = radius_square(sites_cart[o], sites_cart[vlist[o][j]]);
        if(radius < critical_r*critical_r){
            return false;
        }
    }
    return true;
}

void mc_verlet::optimize_delx(){
    std::vector<double> prob;
//    for(std::size_t i = 0; i < 10*sites_cart.size(); i++){
//        move();
//    }
    for(std::size_t i = 1; i <= 20; i++){
        double del_x = critical_r*i/20.;
        std::size_t counter = 0;
        for(std::size_t j = 0; j < 10000; j++){
            if(move()){
                counter++;
            }
        }
        prob.push_back(counter/10000.);
    }
    double tmp = fabs(prob[0] - 0.2);
    del_x = critical_r/10.;
    for(std::size_t i = 0; i < prob.size(); i++){
        cout<< i <<", " <<  prob[i] << endl;
        double tmp1 = fabs(0.2 - prob[i]);
        if(tmp1 < tmp){
            tmp = tmp1;
            del_x = critical_r*(i+1)/10.;
        }
    }
    cout <<"del_x is " << del_x <<endl;
}

bool mcmove(std::vector<cart_coord> & sites_cart, double critical_r, double box, double del_x){
    std::size_t total_atom = sites_cart.size();
    std::size_t o = genrand_real2()*total_atom;
    std::vector<cart_coord> new_sites = sites_cart;
    new_sites[o][0] += (genrand_real1() - 0.5)*del_x;
    new_sites[o][1] += (genrand_real1() - 0.5)*del_x;
    new_sites[o][2] += (genrand_real1() - 0.5)*del_x;
    if(ener(new_sites, critical_r, box)){
        sites_cart = new_sites;
        return true;
    }
    return false;
}

void init_mc(std::vector<cart_coord> & sites_cart, std::size_t total_atom, double box, double critical_r){
    std::vector<std::array<int, 3> > index;
    int max_index = box/critical_r;
    std::size_t counter = 0;
//    for(std::size_t i = 1; i < max_index; i++){
//        for(std::size_t j = 1; j < max_index; j++){
//            for(std::size_t k = 1; k < max_index; k++){
//                counter++;
//                std::array<int, 3> index_i = {i, j, k};
//                index.push_back(index_i);
//                if(counter == total_atom){
//                    break;
//                }
//            }
//            if(counter == total_atom)
//                break;
//        }
//        if(counter == total_atom)
//            break;
//    }
    for(std::size_t i = 0; i < total_atom; i++){
        if(i == 0){
            std::array<int,3> index_i = {1,1,1};
            index.push_back(index_i);
        }
        else{
            std::array<int, 3> index_i = {0,0,0};
            while(1){
                index_i[0] = genrand_real3()*max_index;
                index_i[1] = genrand_real3()*max_index;
                index_i[2] = genrand_real3()*max_index;
                std::size_t counter = 0;
                for(std::size_t j = 0; j < index.size(); j++){
                    if(index_i[0] == index[j][0] and index_i[1] == index[j][1] and index_i[2] == index[j][2]){
                        break;
                    }
                    else{
                        counter++;
                    }
                }
                if(counter == index.size()){
                    index.push_back(index_i);
//                    cout << i << " successfully place a hard ball" << endl;
                    break;
                }
            }
        }
    }
    for(std::size_t i = 0; i < total_atom; i++){
        cart_coord site = {index[i][0]*critical_r, index[i][1]*critical_r, index[i][2]*critical_r};
        sites_cart.push_back(site);
    }
    cout << "the number of atoms in simulation is " << sites_cart.size() <<endl;
}

void mc_verlet::sample(){
    ngr++;
    for(std::size_t i = 0; i < sites_cart.size(); i++){
        for(std::size_t k = 0; k < vlist[i].size(); k++){
            std::size_t j = vlist[i][k];
//            double radius = radius_square(sites_cart[i], sites_cart[j]);
            cart_coord new_cart = {0,0,0};
            new_cart[0] = sites_cart[i][0] - sites_cart[j][0];
            new_cart[1] = sites_cart[i][1] - sites_cart[j][1];
            new_cart[2] = sites_cart[i][2] - sites_cart[j][2];
            double radius = new_cart[0]*new_cart[0] + new_cart[1]*new_cart[1] + new_cart[2]*new_cart[2];
            radius = sqrt(radius);
            if(radius < rv){
                int ig = int(radius/delg);
                g_r[ig] += 1;
            }
        }
    }
}

int main(int argc, char* argv[]){
    unsigned long init[4]={0x123, 0x234, 0x345, 0x456}, length=4;
    init_by_array(init, length);
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
//    double box = 20;
//    std::size_t total_atom = 4000;
    double critical_r = 1;
    double del_x = critical_r;
    double delg = 0.1*critical_r;
//    double rv = 0.25*box;
    double rv = 5*critical_r;
    std::size_t nsampl = 10;
    std::size_t nhis = rv/delg;
    std::vector<cart_coord> sites_cart;
    std::vector<double> gr;
    init_mc(sites_cart, total_atom, box, critical_r);
    mc_verlet mc(critical_r, box, del_x, rv, nhis, sites_cart);
    mc.optimize_delx();
    for(std::size_t i = 0; i < ncycl; i++){
        bool acc = mc.move();
        if( i%nsampl == 0){
            mc.sample();
        }
    }
    for(std::size_t i = 0; i < mc.g_r.size(); i++){
        double vb = ((i+1)*(i+1)*(i+1) - i*i*i)*mc.delg*mc.delg*mc.delg;
//        double vb = (i+0.5)*(i+0.5)*mc.delg*mc.delg*mc.delg;
        double nid = 4./3*PI*vb*sites_cart.size()/(box*box*box);
        gr.push_back(mc.g_r[i]/(mc.ngr*total_atom*nid));
    }
    for(std::size_t i = 0; i < gr.size(); i++){
       cout << i*mc.delg << " " << gr[i] << endl;
    }
}
