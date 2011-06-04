/*
 *  ShrinkageHelpers.cpp
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

#include "ShrinkageHelpers.h"

using namespace std;

namespace mocapy
{
double get_mean(MDArray<double> vec)
{
        double mean=0;
        int len=vec.get_shape()[0];
        for (int i=0;i<len;i++){
                mean+=vec[i];
        }
        return mean/len;

}
double get_var(MDArray<double> vec)
{
        double var=0;
        int len=vec.get_shape()[0];
        double mean=get_mean(vec);

        for (int i=0;i<len;i++){
                var+=(vec[i]-mean)*(vec[i]-mean);
        }
        return var/(len-1);

}
double get_median(MDArray<double> vec)
{
        MDArray<double> theValues;
        theValues = Sort(vec);
        int len=theValues.get_shape()[0];

        if ((len % 2) == 1)

                return theValues[(len+1)/2-1];
        else{

                double lower = theValues[len/2-1];
                double upper = theValues[len/2];
                return (lower + upper)/2;
        }
}
MDArray<double> Sort(MDArray<double> vec)
{
        MDArray<double> sorted(vec);
        int index_of_next_smallest;
        int len=vec.get_shape()[0];

        for (int index = 0; index < len - 1; index++)
        {
                index_of_next_smallest = index_of_smallest(sorted, index,len);
                Swap(sorted[index], sorted[index_of_next_smallest]);
        }
        return sorted;
}

void Swap(double& v1, double& v2)
{
        double temp;
        temp = v1;
        v1 = v2;
        v2 = temp;
}

int index_of_smallest(MDArray<double> vec, int start_index, int number_used)
{
        double min = vec[start_index],
        index_of_min = start_index;
        for (int index = start_index + 1; index < number_used; index++)
                if(vec[index] < min)
                {
                        min = vec[index];
                        index_of_min = index;
                }
        return index_of_min;
}
double get_tot_mean(MDArray<double> vec)
{
        int len=vec.get_shape()[0];
        double mean=0;
        for (int i=0;i<len;i++)
                for (int j=0;j<len;j++)
                        mean+=vec.get(i,j);
        return mean/len;
}
double unb_emp_var(MDArray<double> mat,int i)
{
        double wii_mean=0;
        int n=mat.get_shape()[0];
        vector<uint> index;

        for (int k=0;k<n;k++){
                index=vec((uint)k, (uint)i, (uint)i);
                wii_mean+=mat[index];
        }
        wii_mean/=n;
        
	return wii_mean*n/(n-1);
}
double unb_emp_covar(MDArray<double> mat,int i,int j)
{
        double wij_mean=0;
        int n=mat.get_shape()[0];
        vector<uint> index;

        for (int k=0;k<n;k++){
                index=vec((uint)k, (uint)i, (uint)j);
                wij_mean+=mat[index];
        }
        wij_mean/=n;
        
        return wij_mean*n/(n-1);
}
double w_var(MDArray<double> mat,int i,int j)
{
        double wij_mean=0,var_w=0;
        int n=mat.get_shape()[0];
        vector<uint> index;

        for (int k=0;k<n;k++){
                index=vec((uint)k, (uint)i, (uint)j);
                wij_mean+=mat[index];
        }
        wij_mean/=n;
        for(int k=0;k<n;k++){
                index=vec((uint)k, (uint)i, (uint)j);
                var_w+=(mat[index]-wij_mean)*(mat[index]-wij_mean);
        }
        return var_w*n/((n-1)*(n-1)*(n-1));
}
double w_covar(MDArray<double> mat,int i,int j,int l,int m)
{
        double wij_mean,wlm_mean,sum1=0,sum2=0;
        int n=mat.get_shape()[0];
        vector<uint> index1,index2;

        for (int k=0;k<n;k++){
                index1=vec((uint)k, (uint)i, (uint)j);
                index2=vec((uint)k, (uint)l, (uint)m);
                sum1+=mat[index1];
                sum2+=mat[index2];
        }
        wij_mean=sum1/n;
        wlm_mean=sum2/n;

        sum1=0;
        for(int k=0;k<n;k++){
                index1=vec((uint)k, (uint)i, (uint)j);
                index2=vec((uint)k, (uint)l, (uint)m);
                sum1+=(mat[index1]-wij_mean)*(mat[index2]-wlm_mean);
        }
        return sum1*n/((n-1)*(n-1)*(n-1));
}

MDArray<double> get_diag(MDArray<double> mat)
{
        int n=mat.get_shape()[0];
        int p=mat.get_shape()[1];

	assert(n==p);

        MDArray<double> diag(vec(n));
        for (int i=0;i<p;i++)
                diag.set(i,mat.get(i,i));

        return diag;
}
MDArray<double> get_off_diag(MDArray<double> mat)
{
        int n=mat.get_shape()[0];
        int p=mat.get_shape()[1];

	assert(n==p);
	
        int index=0;
        MDArray<double> diag(vec(n*p-n));
        for (int i=0;i<n;i++){
                for (int j=0;j<p;j++)
                        if (i!=j){
                                diag.set(index,mat.get(i,j));
                                index++;
                        }
        }
        return diag;
}
MDArray<double> get_column(MDArray<double> mat, int col)
{
        int rows=mat.get_shape()[0];
        int columns=mat.get_shape()[1];
        
	assert(col<=columns-1);
        
	MDArray<double> col_vec(vec(rows));
        for (int i=0;i<rows;i++){
                col_vec[i]=mat.get(i,col);
        }
        return col_vec;
}

MDArray<double> get_SDM(MDArray<double> samples)
{
        int n=samples.get_shape()[0];
        int p=samples.get_shape()[1];
        MDArray<double> SDM(vec(n,p));
        MDArray<double> means_j(vec(p));
        MDArray<double> sd_j(vec(p));

        for (int j=0;j<p;j++){
                means_j[j]=get_mean(get_column(samples,j));
                sd_j[j]=sqrt(get_var(get_column(samples,j)));
        }
        for (int i=0;i<n;i++)
                for(int j=0;j<p;j++)
                        SDM.set(i,j,((samples.get(i,j)-means_j[j])/sd_j[j]));

        return SDM;
}
MDArray<double> get_partial_correlation(MDArray<double> MLE)
{
        int n=MLE.get_shape()[0];
        int p=MLE.get_shape()[1];
        MDArray<double> part_cor(vec(n,p));

        for (int i=0;i<p;i++)
                for (int j=0;j<n;j++)
                        part_cor.set(i,j,MLE.get(i,j)/sqrt(MLE.get(i,i)*MLE.get(j,j)));
	
        return part_cor;

}

double round_lambda(double lambda)
{
        if (lambda<0)
                return 0;
        else if (lambda>1)
                return 1;
        return lambda;
}

MDArray<double> get_sub_means(MDArray<double> mat)
{
        int n=mat.get_shape()[0];
        int p=mat.get_shape()[1];
        MDArray<double> sub_means(vec(p));
        double sum;
        for (int i=0;i<p;i++){
                sum=0;
                for (int j=0;j<n;j++)
                        sum+=mat.get(j,i);
                sub_means[i]=sum/n;
        }
        return sub_means;
}
MDArray<double> add_matrices(MDArray<double> m1,MDArray<double> m2)
{
        assert(m1.get_shape() == m2.get_shape());
	
	int n=m1.get_shape()[0];
        int p=m1.get_shape()[1];
	MDArray<double> sum(vec(n,p));
	for(int i=0;i<n;i++)
		for(int j=0;j<p;j++)
			sum.set(i,j,m1.get(i,j)+m2.get(i,j));
	return sum;
}
MDArray<double> mult_num_to_matrix(MDArray<double> mat,double num)
{
        int n=mat.get_shape()[0];
        int p=mat.get_shape()[1];
        MDArray<double> prod(vec(n,p));
        for(int i=0;i<n;i++)
                for(int j=0;j<p;j++)
                        prod.set(i,j,num*mat.get(i,j));
        return prod;
}


}
