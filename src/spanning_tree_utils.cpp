// [[Rcpp::depends(BH)]]
// [[Rcpp::depends(RcppEigen)]]
#define BOOST_DISABLE_ASSERTS
#define EIGEN_NO_DEBUG

// Header files
#include <RcppArmadillo.h>
#include <RcppEigen.h>
#include <math.h>
#include <boost/config.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/random_spanning_tree.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/graph/graph_utility.hpp>
#include <boost/graph/subgraph.hpp>
#include <boost/graph/undirected_graph.hpp>
#include <boost/graph/exterior_property.hpp>
#include <boost/graph/johnson_all_pairs_shortest.hpp>
#include <boost/graph/floyd_warshall_shortest.hpp>
#include <boost/array.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include <random>

using namespace Rcpp;
using namespace boost;

using Eigen::Map;                 // 'maps' rather than copies 
using Eigen::MatrixXd;                  // variable size matrix, double  precision
using Eigen::VectorXd;                  // variable size vector, double precision
using Eigen::SelfAdjointEigenSolver;    // one of the eigenvalue solvers
using Eigen::MatrixXi;
using Eigen::MatrixBase;

// Function to convert a List to a Boost graph object 
adjacency_list<> conv_List_boost(List aList){
  // Convert a List to a Boost graph object 
  adjacency_list<> mygraph; 
  // Loop through vectors in aList
  for(int i = 0; i < aList.size(); i++){
    // Get i'th entry in list
    NumericVector list1;
    list1 = aList(i);    
    // Loop through elements in list1
    for(int j = 0; j < list1.size(); j++){
      add_edge(i,list1(j),mygraph);
    }
  }    
  return(mygraph);
}

// Function to output a random spanning tree from a Boost graph object 
adjacency_list<> uniform_spanning_tree(adjacency_list<> mygraph){
  int ms = num_vertices(mygraph);
  random::uniform_int_distribution<> dist(0, ms-1);
  //std::array<int> predeces;
  array<int, 10000> predeces;
  //std::mt19937 rg( (unsigned)time(NULL) );
  static std::random_device rg;
  int root = dist(rg);
  random_spanning_tree(mygraph, rg,
		       predecessor_map(predeces.begin()).
		       root_vertex(root));
  //std::cout << "root: " << root << '\n'; //Rcpp::Rcout << 
  //std::cout << "predeces.size(): " << predeces.max_size() << '\n';

  adjacency_list<> treegraph;
  for (int i = 0; i < ms; i++) {
    if(predeces[i]!=-1) {  
      int j = predeces[i];
      //std::cout << "node: " << i ;
      //std::cout << " predecessors[i]: " << j << '\n';
      add_edge(i, j, treegraph);
      add_edge(j, i, treegraph);
    }
  }
  return(treegraph);
}

// Function to convert a Boost graph object to a List
List conv_boost_List(adjacency_list<> treegraph){
  int ms = num_vertices(treegraph);
  List tList;  
  // Loop through 
  for(int i = 0; i < ms; i++){
    NumericVector list1;
    //int k = 0; 
    for(int j = 0; j < ms; j++){  
      if(edge(i, j, treegraph).second != 0){
	list1.push_back(j);
      }       
    }
    tList.push_back(list1);
  }  
  return(tList); 
}

/* Sample randon integer in a given range */
IntegerVector sample_int(IntegerVector pool, int n) {
  //IntegerVector pool = seq(min, max);
  std::random_shuffle(pool.begin(), pool.end());
  return pool[Range(0, n - 1)];
}

/* Integrated function generate random spanning tree from
   List object */
List UST(List aList0){
  // Convert aList0 to boost object
  adjacency_list<> boost0=conv_List_boost(aList0);
  adjacency_list<> tboost = uniform_spanning_tree(boost0);
  List aList2 = conv_boost_List(tboost);
  return(aList2); 
}

/* Uniform partition sampling gviven a tree */
IntegerVector UPS(List aList, int cutnum){
  //IntegerMatrix v((aList.size()-1),2);
  IntegerVector vv = sample_int(Range(0,aList.size()-2),cutnum);
  //Rcpp::Rcout << "vv:  " << vv << '\n';
  int k = 0;
  adjacency_list<> mygraph; 
  // Loop through vectors in aList
  //Rcpp::Rcout << "aList.size():  " << aList.size() << '\n';
  for(int i = 0; i < aList.size(); i++){
    add_vertex(mygraph);
    // Get i'th entry in list
    NumericVector list1;
    list1 = aList(i);    
    // Loop through elements in list1
    for(int j = 0; j < list1.size(); j++){      
      if (i > list1(j)){
        int kk = 0;
        for(int l=0; l < vv.size(); l++){
          if (vv[l]==k){
            kk=1;
          }
        }
        if (kk == 0){
	  //Rcpp::Rcout << "add_edge i:" << i << "  list1(j): " <<  list1(j) << '\n';
	  add_edge(i,list1(j),mygraph);
	  add_edge(list1(j),i,mygraph);
        }else{
	  //Rcpp::Rcout << "cut_edge i:" << i << "  list1(j): " <<  list1(j) << '\n';
	}
	//v[k,1]=i;v[k,2]=list1(j);
	//Rcpp::Rcout << k << " v:  " << v[k,1] << "; v:  " << v[k,2] <<'\n';
	k++;
      }
    }
  }    
  // int nm = num_vertices (mygraph);
  // //Rcpp::Rcout << "num_vertices (mygraph):  " << nm << '\n';
  // //Rcpp::Rcout << "aList.size():  " << aList.size() << '\n';
  // if (nm!=aList.size()){
  //   add_vertex(mygraph);
  // }
  IntegerVector component(num_vertices (mygraph));
  //Rcpp::Rcout << "component:  " << component << '\n';
  int num_components = connected_components(mygraph, &component[0]);
  //Rcpp::Rcout << "num_components:  " << num_components << '\n';
  //Rcpp::Rcout << "component:  " << component << '\n';
  // Rcpp::Rcout << "num_components:  " << num_components << '\n';
  // for (int i = 0; i < component.size(); i++){
  //   Rcpp::Rcout << "Vertex " << i <<" is in component " << component[i] << '\n';
  // }
  return(component); 
}

/* Graph Laplacian determinant calculator 
   GL matrix is divided by a constant to prevent overflow*/
double LapDet(NumericMatrix A){

  int N = A.nrow(); 
  NumericMatrix D(N-1,N-1);
  std::fill( D.begin(), D.end(), 0) ;

  for(int i = 0; i < N; i++){
    for(int j = 0; j < N; j++){      
      if(A(i,j)>0){if(i!=(N-1)){
	  D(i,i)=D(i,i)+A(i,j);
	  if(j!=(N-1)){
	    D(i,j)=-A(i,j);
	  }  
	}}
    }
  } 
  D=D/pow(10,1);
  const Eigen::Map<Eigen::MatrixXd> MM(as<Eigen::Map<Eigen::MatrixXd> >(D));
  double d = MM.determinant();
  return(d);
}  

/* Spanning tree number calculator */
double NT(List aList, IntegerVector component){
  //,IntegerVector popvec, 
  int u = max(component)+1; 
  NumericVector NC(u);
  NumericMatrix M1(u,u); //multigraph 
  std::fill( M1.begin(), M1.end(), 0) ;
  //det vector
  NumericVector dv(u+1);
  //accumulate by district 
  for(int k = 0; k < u; k++){    
    NC(k)=sum(component==k);
  }
  for(int k = 0; k < u; k++){ 
    NumericMatrix Mt(aList.size(),aList.size());  
    std::fill( Mt.begin(), Mt.end(), 0); 
    for(int i = 0; i < aList.size(); i++){
      // Get i'th entry in list
      NumericVector list1;
      list1 = aList(i);    
      // Loop through elements in list1
      for(int j = 0; j < list1.size(); j++){      
	if(component(i)==k && component(list1(j))==k){
	  Mt(i,list1(j))=1;//assign by num 
	}
	if(k==0){
	  if(component(i)!=component(list1(j))){
	    M1(component(i),component(list1(j)))=M1(component(i),component(list1(j)))+1;}
	}}}
    if(NC(k)>1){ 
      int ii = -1;int jj = -1;
      NumericMatrix Mt2(NC(k),NC(k));
      std::fill( Mt2.begin(), Mt2.end(), 0); 
      for(int i = 0; i < aList.size(); i++){
	if(component(i)==k){ii=ii+1;
	  for(int j = 0; j < aList.size(); j++){
	    if(component(j)==k){jj=jj+1;
	      if(component(j)==k){
		Mt2(ii,jj)=Mt(i,j);
	      }}}          
	}
	jj = -1;}   
      dv(k)=LapDet(Mt2); 
    } else{
      dv(k)=1;
    }
    //Rcpp::Rcout << "NC(k):  " << NC(k) << '\n';
    //Rcpp::Rcout << "dv(k):  " << dv(k) << '\n';
  }

  dv(u) = LapDet(M1);
  //Rcpp::Rcout << "dv(u):  " << dv(u) << '\n';
  //Rcpp::Rcout << "dv(0):  " << dv(0) << '\n';
  //double dd = std::accumulate(dv.begin(),dv.end(), 1, std::multiplies<int>());
  double dd =1;
  for(int k = 0; k < (u+1); k++){
    dd = dd*dv(k);
    //   Rcpp::Rcout << "dv(k):  " << dv(k) << '\n';
  }
  //dd = std::abs(dd);
  return(dd);
}

// Sample partition function
// [[Rcpp::export]]
List sample_partition(List aList, int num_partitions){

  // UST step
  List tree_out = UST(aList);

  // UPS step
  IntegerVector ups_out = UPS(tree_out, num_partitions - 1);

  // Get NT for inverse probability reweighting
  double nt_out = NT(aList, ups_out);

  // Create output
  List out;
  out["partition"] = ups_out;
  out["prob_sample_partition"] = nt_out;

  return out;
  
}

