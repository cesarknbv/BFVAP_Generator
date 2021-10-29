#include <ilcplex/ilocplex.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <tuple>
#include <iomanip>
#include <string>
#include <cstdlib>
#include <new>
#include <vector>
#include <map>
#include <tuple>
#include <iterator>
#include <utility> 
#include <cstddef>
#include <type_traits>
#include <chrono>
#include <cfloat>
#include <time.h>



using namespace std;

#ifndef vectors_typedef
#define vectors_typedef
typedef std::vector<std::vector<std::vector<std::vector<double> > > >        vec4Dou;
typedef std::vector<std::vector<std::vector<double> > >                      vec3Dou;
typedef std::vector<std::vector<double> >                                    vec2Dou;
typedef std::vector<double>                                                  vec1Dou;
typedef std::vector<std::vector<std::vector<std::vector<int> > > >           vec4Int;
typedef std::vector<std::vector<std::vector<int> > >                         vec3Int;
typedef std::vector<std::vector<int> >                                       vec2Int;
typedef std::vector<int>                                                     vec1Int;
typedef vector<bool>                                      	                 vec1Bool;

typedef tuple<int,int,int>                  	                 			 tup3Int;
typedef tuple<int,int,int,int>                  	             			 tup4Int;
typedef vector<tup4Int>                	                     				 vec_tup4Int;
typedef vector<tup3Int>                	                     				 vec_tup3Int;
typedef vector< pair<int,int> >                	                             vec_pair;

typedef map<tup3Int,double,bool(*)(const tup3Int&,const tup3Int&)>           map_tup3Int_Dou;


typedef vec1Int::iterator                                          		 iteInt;
typedef vec1Int::const_iterator                                    		 citeInt;
typedef vec1Int::reverse_iterator                                  		 riteInt;
typedef vec1Int::const_reverse_iterator                            		 criteInt;
typedef vec1Dou::iterator                                          		 iteDou;
typedef vec1Dou::const_iterator                                    		 citeDou;
typedef vec1Dou::reverse_iterator                                  		 riteDou;
typedef vec1Dou::const_reverse_iterator                            		 criteDou;
                                                                      
typedef vec_pair::const_iterator                        vec_pair_cit;
typedef vec_pair::iterator                              vec_pair_it;
typedef vec_tup3Int::const_iterator						vec_tup3Int_cit;
typedef vec_tup3Int::iterator							vec_tup3Int_it;
typedef vec_tup4Int::const_iterator						vec_tup4Int_cit;
typedef vec_tup4Int::iterator							vec_tup4Int_it;
typedef map_tup3Int_Dou::const_iterator						map_tup3Int_Dou_cit;
typedef map_tup3Int_Dou::iterator							map_tup3Int_Dou_it;

                                                                 
 #endif


#ifndef DVAP
#define DVAP





struct Instance_data
{
    int I;              
	int J;              
    int T;              
    int V;              
	vec2Int  tau;   
	vec3Dou  cos;   
	vec3Dou  cosf;   
	vec3Dou  pro;   
	vec3Int  m;	  
	vec3Int  dem;	
	vec3Int  A;	    
	vec3Dou  In;
	vec2Int  Kap;
	
	
};



#endif

bool tup_comp (const tup3Int& lhs, const tup3Int& rhs) {

	return ( get<2>(lhs) < get<2>(rhs) )                                                                            || 
		   ((get<2>(lhs) == get<2>(rhs)) && (( get<0>(lhs) < get<0>(rhs) )))                                        ||
		   ((get<2>(lhs) == get<2>(rhs)) && (( get<0>(lhs) == get<0>(rhs) )) && (( get<1>(lhs) < get<1>(rhs) )))     ;
}

bool between (const int& lhs, const int& val, const int& rhs) {

	return ( (val >= lhs) && ( val <= rhs) );
}

 #ifndef graph_def
#define graph_def
class Graph 
{
	

	public:
	    /* The following members are static since its values and dimensions are the same for each objecto of graph*/
	    static int I,T,V;
		static int          n;
		static int          m;
		static int          oriented;
		static vec1Int   orr;
		static vec1Int   ex;
		
		/* The following are adjacency list of arcs, adjacency list of nodes, and tracker of values for each node*/
		vec1Int          la;
		vec1Int          lp;
		vec1Int          ls;
		vec1Dou          cost;
		
		Graph( ): la(),
				  lp(),
				  ls(),
				  cost()
				  {
					//cout << "\nEntering contructor of Graph is called\n" << endl << flush;
					//cout << "\nExiting contructor of Graph is called\n" << endl << flush;
				  }
				  
				  
		friend ostream &operator<<( ostream &output, const Graph &objeto )
	    {
			
			try{
			
				if ( m == 0 || orr.empty() || ex.empty() || I == 0 ){
					output << "\nSome of the paramenters are not initialized yet\n " << flush;
				}else{
				
					output << "\n" << setw(8) << "arc" << setw(10) << "i" << setw(4) << "t" << setw(10) << "j" << setw(4) << "t" << flush;
					
					for ( int k = 0; k < m/*(m-I)*/; k++ ){
						
						output << "\n" << setw(8) << k << setw(10) << orr.at(k)%I 
							   << setw(4) << orr.at(k)/I << setw(10) << ex.at(k)%I 
							   << setw(4) << ex.at(k)/I << flush;
					}
				}
			}
			catch (const exception &exc)
			{
				// catch anything thrown within try block that derives from std::exception
				cerr << "\n\nprint graph - " << exc.what() << flush;
				exit(EXIT_FAILURE);
			}
			
			
		}
				  
		
		~Graph(){
			//cout << "\n  Entering destructor of Graph\n" << flush;
			if ( !la.empty() ) la.clear() ;
			if ( !lp.empty() ) lp.clear() ;
			if ( !ls.empty() ) ls.clear() ;
			if ( !cost.empty() ) cost.clear();
			//cout << "\n  Exiting destructor of Graph\n" << flush;
		}
    
    
    
};

/* Initialize static members of G*/
int         Graph::I         =  0;
int         Graph::V         =  0;
int         Graph::T         =  0;
int         Graph::n         =  0;
int         Graph::m         =  0;
int         Graph::oriented  =  1;
vec1Int  Graph::orr(0);       
vec1Int  Graph::ex(0);  

typedef int Vertex;
typedef int Veh;
typedef int Arc;
typedef double Flow;

unsigned semilla;
string n_Glo; 

static vec1Dou d_tableau;
static vec1Int prec_tableau;

//const double TOL_ZERO = 1.E-6;

#endif 

ostream_iterator< int >   output_arci( cout, " " );
ostream_iterator< double > output_arcd( cout, " " );

const int    zeroInt = 0; 
const int    oneInt =  1; 
const double zeroDou = 0.0;
const double oneDou =  1.0;

const double EPS1 = 1E-1;
const double EPS2 = 1E-2;
const double EPS3 = 1E-3;
const double EPS4 = 1E-4;
const double EPS5 = 1E-5;
const double EPS6 = 1E-6;
const double EPS7 = 1E-7;
const double EPS8 = 1E-8;

/* ************************************************** */


void adjacency_list( 
	Instance_data &instance, 
	Graph &G );
	
void succesor_list( 
	Instance_data &instance, 
	Graph &G );
	
void set_network( 
    Instance_data &instance,
	Graph &G, 
	int oriented);


void DAG_solveshortest( 
	Instance_data &instance, 
	Graph &G, 
	Veh v,
	Vertex orr );
	
void update_Costs_Veh ( 
	Instance_data &instance, 
	Graph &G,
    Veh v );
    
    
inline int flat_ixt(Instance_data *instance, const int t, const int i ){
		
	return (  (instance->I * t)+i < (instance->I*(instance->T+1)) ? (instance->I * t)+i:(instance->I * instance->T)+i  );
	
}

inline int flat_ixt(Instance_data &instance, const int t, const int i ){
		
	return (  (instance.I * t)+i < (instance.I*(instance.T+1)) ? (instance.I * t)+i:(instance.I * instance.T)+i  );
	
}






