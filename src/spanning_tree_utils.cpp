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

/* Graph reduction given a population constraint*/
List TR(List aList, NumericVector popvec, int numdist, double delta){//
  // Convert aList TREE to boost object
  adjacency_list<> mygraph=conv_List_boost(aList);
  int ms = num_vertices(mygraph);
  IntegerVector RV(ms);
  std::fill( RV.begin(), RV.end(), 0) ;
  // arma::vec RV;
  // for(int i = 0; i < ms; i++){
  //   RV(i)=0;
  // }  

  //IntegerVector ind(ms);
  //std::fill( ind.begin(), ind.end(), 1) ;
  //std::cout << "ind(0): " << ind(0) ;
  IntegerVector deg(ms);
  IntegerVector K1(ms);

  //IntegerVector EL(ms);
  int popall;popall=0;
  for(int j = 0; j < ms; j++){
    popall += popvec(j);
  }
  double popmean;popmean=popall/numdist;  
  double popmin;
  double popmax;
  //temporary for unconstrained delta=1
  if(delta==1){
    popmin=0;
    popmax=popall;    
  }else{  
    popmin=popmean*(1-delta);
    popmax=popmean*(1+delta);}

  random::uniform_int_distribution<> dist(0, ms-1);
  //std::array<int> predeces;
  array<int, 30> predeces;//array<int, 2000> predeces;
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
      //EL(i)=j;
      add_edge(i, j, treegraph);
      add_edge(j, i, treegraph);
    }
  }

  List aList2 = conv_boost_List(treegraph);
  
  for(int i = 0; i < aList2.size(); i++){
    // Get i'th entry in list
    NumericVector list1;
    list1 = aList2(i);    
    // Loop through elements in list1
    deg(i) = list1.size();
  }     

  IntegerVector K2(ms);
  IntegerVector K3(ms);

  int k = 0; NumericVector list1;
  for(int i = 0; i < aList2.size(); i++){
    if(deg(i) == 1){      
      k=k+1;  //label
      RV(i) = k;//ind(i)=0;
      int popk = popvec(i);
      
      list1 = aList2(i);   
      int j = list1(0); //i's neighbors 
      if(deg(j)==2){
	RV(j) = k;//ind(j)=0;
	int K=1;
	int dj=deg(j); 
        while(dj == 2 && K>0){//K
          list1 = aList2(j); int KK=1;
          for(int jj = 0; jj < list1.size(); jj++){ 
            //std::cout << "list1(jj): " << list1(jj) << '\n';
            if(RV(list1(jj))==0){ 
	      //if(deg(list1(jj))==2){RV(list1(jj))=k;j=list1(jj);dj=2;KK=0;}
	      int popk2=popk+popvec(list1(jj));
	      if(deg(list1(jj))==2 && popk2<popmin){RV(list1(jj))=k;popk=popk+popvec(list1(jj));j=list1(jj);dj=2;KK=0;}
	      if(deg(list1(jj))>2){K3(k)=list1(jj);K=0;}
	    }
            if(jj+1 == list1.size() && KK==1){K=0;}
          }//for   
        }//while    
	K2(k)=j;//K2.push_back(j);
      }else{K2(k)=i;K3(k)=j;}//else{K2.push_back(i);K3.push_back(j);}
      //if(deg(j)>2){K2(k)=j;}    
    }//if 
  }//END first-stage
  
  //IntegerVector K02;
  //IntegerVector K03;
  //for(int i = 0; i < 10; i++){//max(RV)+1
  //  K02(i) = K2(i);
  //  K03(i) = K3(i);
  //}  
  //K2 = K02;
  //K3 = K03;

  IntegerVector RV2(ms);
  for(int i = 0; i < ms; i++){
    RV2(i)=RV(i);
  }

  //arma::vec RV2;
  //for(int i = 0; i < ms; i++){
  //  RV2(i)=RV(i);
  //} 

  arma::vec RV3 = as<arma::vec>(RV) ;
  //arma::uvec rv0 = find(RV3 == 0);
  //int numrv = max(RV)+rv0.n_elem;// need to fix 

  int mrv = max(RV);
  arma::uvec cid;int pop;int jj;
  for(int k=1; k < mrv+1; k++){ 
    //std::cout << "k: " << k << "; K2(k): " << K2(k) << '\n';  
    for(int k2=1; k2 < mrv+1; k2++){   
      if(K3(k)==K3(k2) && k<k2){
	cid = arma::find(RV3 == k || RV3 == k2);   
	pop = popvec(K3(k));
	for(int j = 0; j < cid.n_elem; j++){
	  pop += popvec(cid(j));
	}   
	cid = arma::find(RV3 == k); 

	if(pop<=popmin){
	  RV2(K3(k)) = k2;
	  for(int j = 0; j < cid.n_elem; j++){
	    RV2(cid(j)) = k2;
	    //numrv=numrv-1;        
	  } 
	}
      }
    }}

  RV3 = as<arma::vec>(RV2) ;
  IntegerVector RV4(ms);

  IntegerVector RVV = sort_unique(RV2);arma::uvec ind;
  for(int k=1; k < RVV.size(); k++){
    ind = find(RV3==RVV(k));
    for(int i=0; i < ind.n_elem; i++){
      RV4(ind(i)) = k;
    }
  } 
  //Rcpp::NumericVector(RV.begin(), RV.end());
  //Rcpp::NumericVector(RV2.begin(), RV2.end());
  //int Rvm = max(RV); 

  List out;
  out["RV"] = RV;
  out["RV2"] = RV4;
  out["K2"] = K2;
  out["K3"] = K3;
  out["aList2"] = aList2; //List
  out["popmean"] = popmean; 
  out["popmin"] = popmin; 
  out["popmax"] = popmax; 
  out["numdist"] = numdist; 
  out["delta"] = delta; 
  out["popvec"] = popvec; 
  //out["numrv"] = numrv; 
  return(out);
  //  return(RV); 
} 

/* Sample partitions on a reduced tree */
List RS(List aList, List out, int iter){//
  IntegerVector IND(iter);
  
  int numdist = out["numdist"]; 
  int cutnum = numdist-1;
  NumericVector popvec = out["popvec"];
  double popmin = out["popmin"];
  double popmax = out["popmax"];
  IntegerVector RV = out["RV2"]; 
  int ms = aList.size();
  //Rcpp::Rcout << "ms:  " << ms << '\n';
  IntegerMatrix component2(ms,iter);
  IntegerVector component3(ms,1);

  int numrv=0;
  NumericVector list1;
  for(int i = 0; i < aList.size(); i++){    
    list1 = aList(i); 
    for(int j = 0; j < list1.size(); j++){ 
      if(i > list1(j)){
	if((RV(i)!=RV(list1(j))) || (RV(i)==0 && RV(list1(j))==0))
	  {numrv = numrv + 1;}
      }}}
  
  arma::uvec cid;int pop;
  //////////////////////////////////////////////////////////
  arma::vec component(ms);//RV3 = as<arma::vec>(RV) ;
  int ii = 0;
  while(ii < iter){//IND < numdist && 
    ii = ii + 1;   
    IntegerVector vv = sample_int(Range(0,numrv-1),cutnum); // 
    //int ntr = numrv-1; 
    //Rcpp::Rcout << "numrv-1:  " << ntr << '\n';
    //for(int l=0; l < vv.size(); l++){
    //Rcpp::Rcout << "vv(l): " << vv(l) << '\n';
    //}
    int k = -1;
    int ll = 0;
    adjacency_list<> mygraph; 
    // Loop through vectors in aList
    for(int i = 0; i < aList.size(); i++){    
      add_vertex(mygraph);
      // Get i'th entry in list
      //NumericVector list1;
      list1 = aList(i);        
      // Loop through elements in list1
      for(int j = 0; j < list1.size(); j++){      
	if((RV(i)!=RV(list1(j))) || (RV(i)==0 && RV(list1(j))==0)){ll=1;if (i > list1(j)){k++;}}else{ll=0;}
        if (i > list1(j)){
	  int kk = 0;
	  for(int l=0; l < vv.size(); l++){
	    if (vv[l]==k && ll==1){
	      kk=1;
	    }
	  }
	  if (kk == 0){
	    //Rcpp::Rcout << "k: " << k << " ll: " << ll << " add_edge i:" << i << "  list1(j): " <<  list1(j) << '\n';
	    add_edge(i,list1(j),mygraph);
	    add_edge(list1(j),i,mygraph);
	  }else{
	    //Rcpp::Rcout << "k: " << k << " ll: " << ll << " cut_edge i:" << i << "  list1(j): " <<  list1(j) << '\n';        
	  }
	}
      }
    }    
    //IntegerVector component(ms);
    //Rcpp::Rcout << "num_vertices (mygraph):  " << num_vertices (mygraph) << '\n';
    //Rcpp::Rcout << "component:  " << component << '\n';
    int num_components = connected_components(mygraph, &component[0]);

  
    IND(ii-1)=0;
    for(int i = 0; i < component.max()+1; i++){
      cid = arma::find(component == i);   
      pop = 0;
      for(int j = 0; j < cid.n_elem; j++){
        pop += popvec(cid(j));
      } 
      if(pop>=popmin && pop<=popmax){IND(ii-1) = IND(ii-1)+1;}
    }
    //Rcpp::Rcout << "IND:  " << IND << '\n';
    component3 = wrap(component);
    component2(_,ii-1) = component3;
  }//WHILE  
  //////////////////////////////////////////////////////////
  
  //  IntegerVector component2(ms);
  //  Rcpp::Rcout << "num_vertices (mygraph):  " << num_vertices (mygraph) << '\n';
  //  int num_components = connected_components(mygraph, &component2[0]);

  List out2;
  out2["RV"] = component2;
  out2["IND"] = IND;
  return(out2);
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

