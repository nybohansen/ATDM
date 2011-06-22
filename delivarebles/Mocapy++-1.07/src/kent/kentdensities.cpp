/*
 *  kentdensities.cpp
 *
 *  Copyright (C) 2008, Kasper Stovgaard, The Bioinformatics Centre, University of Copenhagen.
 *
 *  This file is part of Mocapy++.
 *
 *  Mocapy++ is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Mocapy++ is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with Mocapy++.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <sstream>
#include "kentdensities.h"

using namespace std;

namespace mocapy {


KentDensities::KentDensities( vector<double> new_user_kappas,
                              vector<double> new_user_betas,
                              vector<MDArray<double> > new_user_axes,
                              double new_kappa_max,
                              double new_beta_max,
                              double new_beta_factor,
                              double new_kappa_init,
                              double new_beta_init) {

     user_kappas = new_user_kappas;
     user_betas = new_user_betas;
     user_es = new_user_axes;

     kappa_max = new_kappa_max;
     beta_max = new_beta_max;
     beta_factor = new_beta_factor;
     kappa_init = new_kappa_init;
     beta_init = new_beta_init;

     initialize();
};

void KentDensities::initialize() {
     type = KENT;
     output_size = 3;
}

// Initialize your parameters
void KentDensities::construct(vector<uint> & parent_sizes) {
	node_size = parent_sizes[0];
     // Check/construct parameters
     if (!user_kappas.empty() && !user_betas.empty()) {

          assert(user_kappas.size() == node_size);
          kappas = user_kappas;
          assert(user_betas.size() == node_size);
          betas = user_betas;
     }
     else if (user_kappas.empty() && user_betas.empty()) {
          double kappa, beta;

          for(uint i = 0; i < node_size; i++) {
               kappa=0; beta=0;
               // Create random kappa, beta pair
               while (!(kappa>beta_factor*beta)) {
            	   double d1 = randomGen->get_rand();
            	   double d2 = randomGen->get_rand();
            	   kappa = kappa_init*d1;
            	   beta = beta_init*d2;
               }
               kappas.push_back(kappa);
               betas.push_back(beta);
          }

     }
     else{
         assert(0); //Kappa/beta list mismatch
     }

     if (!user_es.empty()) {
          assert(user_es.size() == node_size);
          es = user_es;
     }
     else {
          MDArray<double> tmp;
          tmp.set_shape(3,3);

          for(uint i = 0; i < node_size; i++) {
               create_rand_axes(tmp);
               es.push_back(tmp);
          }
     }

     MDArray<double> rho_t;
     rho_t.set_shape(3,3);
     double ikappa,ibeta;

     for(unsigned int i=0; i<node_size; i++) {

          ikappa=kappas[i];
          ibeta=betas[i];

          // Normalization constant for the density function
          logc.push_back(log(2*M_PI)+ikappa-log(sqrt((ikappa-2*ibeta)*(ikappa+2*ibeta))));

          // Variables for the sampler
          a.push_back(4*(ikappa-2*ibeta));
          b.push_back(4*(ikappa+2*ibeta));
          gamma.push_back((b[i]-a[i])/2);

          if(b[i]>0) {
               c2.push_back(0.5*b[i]/(b[i]-gamma[i])); }
          else {
               c2.push_back(0.0); }

          lam1.push_back(sqrt(a[i]+2*sqrt(gamma[i])));
          lam2.push_back(sqrt(b[i]));

          //Current e1, e2, e3
          rho_t.clear();
          rho_t=es[i];
          rho_t.transpose();
          rho.push_back(rho_t);


     }
}

// Estimate the parameters of the node based on the ESS
// Calculate the parameters (kappa, beta, axes) given the expected
// sufficient statistics ybar (average vector) and S (dispersion).
//
void KentDensities::estimate(vector<MDArray<double> > & ess) {
     assert(!ess.empty());

     MDArray<double> ybar, S, v, H, B, K, lowerB,rho_t;
     rho_t.set_shape(3,3);
     K.set_shape(3,3);
     double r1, r2, ca, cb, psi, cpsi, spsi, curr_kappa, curr_beta;

     for(unsigned int p = 0; p < node_size; p++) {

          // Too little data for estimation
          if(ess[K_N][p] < 10) {
               continue; }

          ybar=ess[K_YBAR].get_view(p)/ess[K_N][p];
          S=ess[K_S].get_view(p)/ess[K_N][p];

           // normalize axis
          r1=sqrt(ybar[0]*ybar[0]+ybar[1]*ybar[1]+ybar[2]*ybar[2]);

          ybar.div_inplace(r1);
          v.set_shape(1,3);
          v.set(0,0,1);v.set(0,1,0);v.set(0,2,0);
          MDArray<double>& vref=v.get_view(0);

          H.set_shape(3,3);
          H.makeRotationMatrix(vref,ybar);

          MDArray<double> & tmp = S;
          tmp.mult3(H);

          B=H;
          B.transpose();
          B.mult3(tmp);

          psi=0.5*atan2(2*B.get(1,2),(B.get(1,1)-B.get(2,2)));

          //Matrix K
          cpsi=cos(psi);
          spsi=sin(psi);
          K.set(0,0,1.0); K.set(0,1,0.0);  K.set(0,2,0.0);   // [1     0         0   ]
          K.set(1,0,0.0); K.set(1,1,cpsi); K.set(1,2,-spsi); // [0 cos(psi) -sin(psi)]
          K.set(2,0,0.0); K.set(2,1,spsi); K.set(2,2,cpsi);  // [0 sin(psi)  cos(psi)]

          //Matrix G
          MDArray<double> G=H;
          G.mult3(K);
          G.transpose();
          es[p]=G;


          lowerB.set_shape(2,2);
          lowerB.set(0,0,B.get(1,1));lowerB.set(0,1,B.get(1,2));
          lowerB.set(1,0,B.get(2,1));lowerB.set(1,1,B.get(2,2));

          // Get (eigenVal, eigenVector)
          pair<vector<double>, vector<MDArray<double> > > eig = lowerB.eigen();

          r2=abs(eig.first[1]-eig.first[0]);
          ca=1/(2-2*r1-r2); // large kappa limit approximation
          cb=1/(2-2*r1+r2);
          curr_kappa=ca+cb;
          curr_beta=0.5*(ca-cb);

          // Check overflow
          if(curr_kappa > kappa_max) {
               double f=kappa_max/curr_kappa;
               curr_kappa=kappa_max;
               curr_beta=f*curr_beta; }

          if(curr_kappa < beta_factor*curr_beta) {
               curr_beta=curr_kappa/beta_factor; }

          kappas[p]=curr_kappa;
          betas[p]=curr_beta;

          // Update normalization const for density function
          logc[p]=log(2*M_PI)+curr_kappa-log(sqrt((curr_kappa-2*curr_beta)*(curr_kappa+2*curr_beta)));

          // Update variables for the sampler
          a[p]=4*(curr_kappa-2*curr_beta);
          b[p]=4*(curr_kappa+2*curr_beta);
          gamma[p]=(b[p]-a[p])/2;

          if(b[p]>0) {
               c2[p]=0.5*b[p]/(b[p]-gamma[p]); }
          else {
               c2[p]=0; }

          lam1[p]=sqrt(a[p]+2*sqrt(gamma[p]));
          lam2[p]=sqrt(b[p]);

          rho_t.clear();
          rho_t=es[p];
          rho_t.transpose();
          rho[p]=rho_t;

     }

}


// Sample from the Kent distribution desribed in
// Kent, JT (1982). The Fisher-Bingham distribution on the sphere. J Royal Statist. Soc. B 44, 71\ufffd80.
vector<double> KentDensities::sample(vector<double> & pv) {

     uint par = (uint) pv.front();

     vector<double> sample;
     double x1, x2, v1, v2, u1, u2,ratio1, ratio2, sign1, sign2, theta, phi, ct;
     double st, cp, sp, x, y ,z;

     while(1) {

          v1=randomGen->get_rand();
          v2=randomGen->get_rand();
          u1=randomGen->get_rand();
          u2=randomGen->get_rand();

          try {
               x1=-log(1-v1*(1-exp(-lam1[par])))/lam1[par]; }
          catch(...) {
               x1=v1; }

          try {
               x2=-log(1-v2*(1-exp(-lam2[par])))/lam2[par]; }
          catch(...) {
               x2=v2; }

          if ((x1*x1+x2*x2)>1) {
               continue; }

          ratio1=exp(-0.5*(a[par]*x1*x1+gamma[par]*x1*x1*x1*x1)-1+lam1[par]*x1);

          if (u1>ratio1) {
               continue; }

          ratio2=exp(-0.5*(b[par]*x2*x2-gamma[par]*x2*x2*x2*x2)-c2[par]+lam2[par]*x2);

          if (u2>ratio2) {
               continue; }

          sign1= (randomGen->get_rand()-0.5) < 0 ? -1.0 : 1.0;
          sign2= (randomGen->get_rand()-0.5) < 0 ? -1.0 : 1.0;

          x1=x1*sign1;
          x2=x2*sign2;

          theta=acos(1-2*(x1*x1+x2*x2));
          phi=atan2(x2, x1);
          ct=cos(theta);
          st=sin(theta);
          cp=cos(phi);
          sp=sin(phi);
          x=ct;
          y=st*cp;
          z=st*sp;

          MDArray<double> p;
          p.set_shape(1,3);p.set(0,0,x);p.set(0,1,y);p.set(0,2,z);
          p.mult3(rho[par]); // rotate sample back to mean direction

          vector<double> tt;
          tt = cart_to_polar(p);

          if(tt[0]<0) {
               continue; }

          sample.push_back(p.get(0,0));sample.push_back(p.get(0,1));sample.push_back(p.get(0,2));

          break;
     }

     return sample;
}

// Return likelihood, P(child|parent)
double KentDensities::get_lik(vector<double> & ptv, bool log_flag) {

     uint p = (uint) ptv.front();

     MDArray<double> v;
     v.set_shape(1,3);
     v.set(0,0,ptv.at(1));v.set(0,1,ptv.at(2));v.set(0,2,ptv.at(3));

     MDArray<double>& v_ref = v.get_view(0);
     MDArray<double>& e1 = es[p].get_view(0);
     MDArray<double>& e2 = es[p].get_view(1);
     MDArray<double>& e3 = es[p].get_view(2);

     double da=kappas[p]*(v_ref.dot(e1));
     double db=v_ref.dot(e2);
     double dc=v_ref.dot(e3);
     double dd=betas[p]*(db*db-dc*dc);

     double log_density = -logc[p]+da+dd;

     if (log_flag) {
          return log_density; }
     else {
          return exp(log_density);
     }
}

// Return the distribtion's parameters
vector<MDArray<double> > KentDensities::get_parameters() {
     vector<MDArray<double> > parms;

     //parms.push_back(kappas);
     //parms.push_back(betas);
     //parms.push_back(es);

     return parms;
}

ostream& operator<<(ostream& out, const KentDensities& a)  {
    for(uint i=0; i<a.node_size; i++){
         out<<a.kappas[i];
         out<< " ";
         out<<a.betas[i];
         out<< "\n";
         out<<a.es[i];
         out<< "\n";
    }

    return out;
}


// Unit vector to array of theta/tau
vector<double> KentDensities::cart_to_polar(MDArray<double> &x) {

    vector<double> tt;
    tt.push_back(acos(x.get(0,0)));
    tt.push_back(atan2(x.get(0,2), x.get(0,1)));

    return tt;
}

// Unit theta/tau angles to x,y,z array
void KentDensities::polar_to_cart(double theta, double tau, MDArray<double> &v, uint dim) {

    v.set(dim,0,cos(theta));
    v.set(dim,1,sin(theta)*cos(tau));
    v.set(dim,2,sin(theta)*sin(tau));
}

MDArray<double> KentDensities::sphere_rand(uint dim){
     assert(dim>0);

     MDArray<double> res;
     res.set_shape(1,dim);
     double l=0;
     for (uint i=0; i<dim; i++) {
          double r = randomGen->get_rand_normal();
          res.set(0,i,r);
          l+=r*r;
     }

     res.div_inplace(sqrt(l));
     return res;
}

// Create three orthonormal vectors e1,e2,e3 with e1 along sphere_rand
void KentDensities::create_rand_axes(MDArray<double> &a) {

     // random vector on the sphere
     MDArray<double> srand=sphere_rand();

     MDArray<double> e;
     e.eye(3);

     // Make random rotation matrix
     MDArray<double> v1=sphere_rand(2); //random x,y

     MDArray<double> rz; //rot around z
     rz.set_shape(3,3);
     rz.set(0,0,v1.get(0,0)); rz.set(0,1,v1.get(0,1)); rz.set(0,2,0);
     rz.set(1,0,-v1.get(0,1)); rz.set(1,1,v1.get(0,0)); rz.set(1,2,0);
     rz.set(2,0,0); rz.set(2,1,0); rz.set(2,2,1);

     MDArray<double> v2=sphere_rand(2);
     MDArray<double> ry; //rot around y
     ry.set_shape(3,3);
     ry.set(0,0,v2.get(0,0)); ry.set(0,1,0); ry.set(0,2,v2.get(0,1));
     ry.set(1,0,0); ry.set(1,1,1); ry.set(1,2,0);
     ry.set(2,0,-v2.get(0,1)); ry.set(2,1,0); ry.set(2,2,v2.get(0,0));

     MDArray<double> v3=sphere_rand(2);
     MDArray<double> rx; //rot around x
     rx.set_shape(3,3);
     rx.set(0,0,1); rx.set(0,1,0); rx.set(0,2,0);
     rx.set(1,0,0); rx.set(1,1,v3.get(0,0)); rx.set(1,2,v3.get(0,1));
     rx.set(2,0,0); rx.set(2,1,-v3.get(0,1)); rx.set(2,2,v3.get(0,0));

     ry.mult3(rz);
     rx.mult3(ry); //Combine the x, y,z rotations
     rx.mult3(e);

     //Find matrix that rotates e1 to srand
     a.set_shape(3,3);
     MDArray<double>& e1=rx.get_view(0);
     MDArray<double>& sr=srand.get_view(0);
     a.makeRotationMatrix(e1,sr);

     //Rotate e1,e2,e3 so that e1 is along p
     a.mult3(rx);
     a.transpose();
}


}
