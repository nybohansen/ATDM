/* pythoninterface.i */

%module wrapper
%include "std_vector.i"


namespace std {
   %template(IntVector) vector<int>;
   %template(DoubleVector) vector<double>;
};


%{
extern void makeDBN();
extern void loadDBN(char* filename);
extern std::vector<double> sample(int i); 
%}




extern void makeDBN();
extern void loadDBN(char* filename);
extern std::vector<double> sample(int i);



