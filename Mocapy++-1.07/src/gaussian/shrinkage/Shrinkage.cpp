/*
 *  Shrinkage.cpp
 *
 *  Copyright (C) 2008, Pedro, The Bioinformatics Centre,
 *  University of Copenhagen.
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

#include "Shrinkage.h"

using namespace std;

namespace mocapy
{
Shrinkage::Shrinkage(MDArray<double> data)
{
        N=data.get_shape()[0];
        P=data.get_shape()[1];
        this->data=data;
        wk=wk_StDM=vec(N,P,P);
        S_MLE=vec(P,P);
        StDM=vec(N,P);
	generate_MLE();
}
Shrinkage::Shrinkage(MDArray<double> data, MDArray<double> MLE)
{	
        N=data.get_shape()[0];
        P=data.get_shape()[1];
	assert(MLE.get_shape()[0]==P && MLE.get_shape()[1]==P);
        this->data=data;
        wk=wk_StDM=vec(N,P,P);
        S_MLE=vec(P,P);
        StDM=vec(N,P);
	set_S_MLE(MLE);
}
void Shrinkage::generate_wk(MDArray<double> mat,int SDM_type)
{
        int n=mat.get_shape()[0];
        int p=mat.get_shape()[1];
        vector<uint> index;
        double sum,term1,term2;
        MDArray<double> sub_means(get_sub_means(mat));
        for(int i=0;i<p;i++)
                for(int j=0;j<p;j++){
                        sum=0;
                        for(int k=0;k<n;k++){
                                term1=mat.get(k,i)-sub_means[i];
                                term2=mat.get(k,j)-sub_means[j];
                                index=vec((uint)k, (uint)i, (uint)j);
                                if(SDM_type)
                                        wk_StDM[index]=term1*term2;
                                else
                                        wk[index]=term1*term2;
                        }
        }
}
MDArray<double> Shrinkage::combine(MDArray<double> target,MDArray<double> MLE,double lambda)
{
	assert(target.get_shape() == MLE.get_shape());

	//        assert(lambda>=0 && lambda<=1);
	lambda = std::max(lambda, 0.0);
	lambda = std::min(lambda, 1.0);

	return add_matrices(mult_num_to_matrix(target,lambda),mult_num_to_matrix(MLE,1-lambda));
}
MDArray<double> Shrinkage::get_data()
{
        return data;
}
MDArray<double> Shrinkage::get_S_MLE()
{
        return S_MLE;
}
MDArray<double> Shrinkage::get_StDM()
{
        return StDM;
}
MDArray<double> Shrinkage::get_wk()
{
        return wk;
}
MDArray<double> Shrinkage::get_wk_StDM()
{
        return wk_StDM;
}
void Shrinkage::set_S_MLE(MDArray<double> mat)
{
	S_MLE=mat;
}
MDArray<double> Shrinkage::compute_f()
{
        double term1,term2;
        MDArray<double> f(vec(P,P));
        for(int i=0;i<P;i++)
                for(int j=0;j<P;j++){
                        term1=sqrt(S_MLE.get(j,j)*S_MLE.get(i,i));
                        term1*=w_covar(wk,i,i,i,j);
                        term2=sqrt(S_MLE.get(i,i)*S_MLE.get(j,j));
                        term2*=w_covar(wk,j,j,i,j);
                        f.set(i,j,0.5*(term1+term2));

                }
        return f;
}
void Shrinkage::generate_MLE()
{
        double sum;
        vector<uint> index;
        generate_wk(data,0);
        for(int i=0;i<P;i++)
                for(int j=0;j<P;j++){
                        sum=0;
                        for(int k=0;k<N;k++){
                                index=vec((uint)k, (uint)i, (uint)j);
                                sum+=wk[index];
                        }
                        S_MLE.set(i,j,sum/(N-1));
        }
}

double Shrinkage::get_lambda_cor()
{
        double sum_up=0,sum_d=0;
        StDM=get_SDM(data);
        generate_wk(StDM,1);

        MDArray<double> part_cor(get_partial_correlation(S_MLE));
        for(int i=0;i<P;i++)
                for(int j=0;j<P;j++)
                        if (i!=j){
                                sum_up+=w_var(wk_StDM,i,j);
                                sum_d+=part_cor.get(i,j)*part_cor.get(i,j);
                        }
        return round_lambda(sum_up/sum_d);
}
double Shrinkage::get_lambda_var()
{
        generate_wk(data,0);
        double sum_up=0,sum_d=0;
        double median=get_median(get_diag(S_MLE));
        for (int i=0;i<P;i++){
                sum_up+=w_var(wk,i,i);
                sum_d+=(S_MLE.get(i,i)-median)*(S_MLE.get(i,i)-median);
        }
        return round_lambda(sum_up/sum_d);
}
double Shrinkage::get_lambda_D()
{
        double sum_up=0,sum_d=0;
        generate_wk(data,0);

        for(int i=0;i<P;i++)
                for(int j=0;j<P;j++)
                        if (i!=j){
                                sum_up+=w_var(wk,i,j);
                                sum_d+=S_MLE.get(i,j)*S_MLE.get(i,j);
                        }
        return round_lambda(sum_up/sum_d);

}
double Shrinkage::get_lambda_A()
{
        double sum_up=0,sum_d=0;
        generate_wk(data,0);

        for(int i=0;i<P;i++)
                for(int j=0;j<P;j++)
                        if (i!=j){
                                sum_up+=w_var(wk,i,j);
                                sum_d+=S_MLE.get(i,j)*S_MLE.get(i,j);
                        }
                        else{
                                sum_up+=w_var(wk,i,i);
                                sum_d+=(S_MLE.get(i,i)-1)*(S_MLE.get(i,i)-1);

                        }
        return round_lambda(sum_up/sum_d);

}
double Shrinkage::get_lambda_B()
{
        double sum_up=0,sum_d=0,avg=get_mean(get_diag(S_MLE));
        generate_wk(data,0);

        for(int i=0;i<P;i++)
                for(int j=0;j<P;j++)
                        if (i!=j){
                                sum_up+=w_var(wk,i,j);
                                sum_d+=S_MLE.get(i,j)*S_MLE.get(i,j);
                        }
                        else{
                                sum_up+=w_var(wk,i,i);
                                sum_d+=(S_MLE.get(i,i)-avg)*(S_MLE.get(i,i)-avg);

                        }
        return round_lambda(sum_up/sum_d);
}
double Shrinkage::get_lambda_C()
{
        double sum_up=0,sum_d=0,avg_diag=get_mean(get_diag(S_MLE)),avg_off_diag(get_mean(get_off_diag(S_MLE)));
        generate_wk(data,0);
        for(int i=0;i<P;i++)
                for(int j=0;j<P;j++)
                        if (i!=j){
                                sum_up+=w_var(wk,i,j);
                                sum_d+=(S_MLE.get(i,j)-avg_off_diag)*(S_MLE.get(i,j)-avg_off_diag);
                        }
                        else{
                                sum_up+=w_var(wk,i,i);
                                sum_d+=(S_MLE.get(i,i)-avg_diag)*(S_MLE.get(i,i)-avg_diag);

                        }
        return round_lambda(sum_up/sum_d);
}
MDArray<double> Shrinkage::get_target_A()
{
        MDArray<double> target(vec(P,P));
        for(int i=0;i<P;i++)
                for(int j=0;j<P;j++)
                        if (i==j)
                                target.set(i,i,1);
        return target;
}
MDArray<double> Shrinkage::get_target_B()
{
        MDArray<double> target(vec(P,P));
        double avg=get_mean(get_diag(S_MLE));
        for(int i=0;i<P;i++)
                for(int j=0;j<P;j++)
                        if (i==j)
                                target.set(i,i,avg);
        return target;
}
MDArray<double> Shrinkage::get_target_C()
{
        MDArray<double> target(vec(P,P));
        double avg_diag=get_mean(get_diag(S_MLE)),avg_off_diag(get_mean(get_off_diag(S_MLE)));
        for(int i=0;i<P;i++)
                for(int j=0;j<P;j++)
                        if (i==j)
                                target.set(i,i,avg_diag);
                        else
                                target.set(i,j,avg_off_diag);

        return target;
}
MDArray<double> Shrinkage::get_target_D()
{
        MDArray<double> target(vec(P,P));
        for(int i=0;i<P;i++)
                for(int j=0;j<P;j++)
                        if (i==j)
                                target.set(i,i,S_MLE.get(i,i));
        return target;
}
}
