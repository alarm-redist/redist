// [[Rcpp::depends(BH)]]
// [[Rcpp::depends(RcppEigen)]]
#define BOOST_DISABLE_ASSERTS
#define EIGEN_NO_DEBUG
// #define ARMA_64BIT_WORD 1

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

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;
using namespace boost;

using Eigen::Map;                 // 'maps' rather than copies
using Eigen::MatrixXd;                  // variable size matrix, double  precision
using Eigen::VectorXd;                  // variable size vector, double precision
using Eigen::SelfAdjointEigenSolver;    // one of the eigenvalue solvers
using Eigen::MatrixXi;
using Eigen::MatrixBase;

typedef Eigen::SparseMatrix<double> SpMat;

// Function to convert a List to a Boost graph object
adjacency_list<> conv_List_boost(std::vector<arma::vec> aList){
  // Convert a List to a Boost graph object
  adjacency_list<> mygraph;
  // Loop through vectors in aList
  for(int i = 0; i < aList.size(); i++){
    // Get i'th entry in list
    arma::vec list1 = aList[i];
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
std::vector<arma::vec> conv_boost_List(adjacency_list<> treegraph){
  int ms = num_vertices(treegraph);
  std::vector<arma::vec> tList(ms);
  // Loop through
  for(int i = 0; i < ms; i++){
    arma::vec list1(ms);
    int counter = 0;
    //int k = 0;
    for(int j = 0; j < ms; j++){
      if(edge(i, j, treegraph).second != 0){
	list1(counter) = j;
	counter++;
      }
    }
    list1.resize(counter);
    tList[i] = list1;
  }
  return(tList);
}

/* Sample randon integer in a given range */
arma::vec sample_int(arma::vec pool, int n) {
  //IntegerVector pool = seq(min, max);
  arma::vec pool_shuffle = shuffle(pool);
  return(pool.subvec(0, n - 1));
}

/* Integrated function generate random spanning tree from
   List object */
std::vector<arma::vec> UST(adjacency_list<> aList0){
  adjacency_list<> tboost = uniform_spanning_tree(aList0);
  std::vector<arma::vec> aList1 = conv_boost_List(tboost);
  return(aList1);
}

/* Uniform partition sampling gviven a tree */
arma::vec UPS(std::vector<arma::vec> aList, int cutnum){
  //IntegerMatrix v((aList.size()-1),2);
  arma::vec vv = sample_int(arma::linspace<arma::vec>(0, aList.size() - 2, aList.size() - 1), cutnum);
  //Rcpp::Rcout << "vv:  " << vv << '\n';
  int k = 0;
  adjacency_list<> mygraph;
  // Loop through vectors in aList
  //Rcpp::Rcout << "aList.size():  " << aList.size() << '\n';
  for(int i = 0; i < aList.size(); i++){
    add_vertex(mygraph);
    // Get i'th entry in list
    arma::vec list1;
    list1 = aList[i];
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
  arma::vec component(num_vertices (mygraph));
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
double LapDet(arma::mat AAA){//, double p
  // NumericMatrix AA = wrap(AAA) ;
  double d = 0.0;
  double n = AAA.n_rows;
  //Rcpp::Rcout << "n:  " << n << '\n';
  if(n>1){
    // Construct D = diag(rowSums(A)) - A
    // SpMat A(Rcpp::as<Eigen::MatrixXd>(AA).sparseView());
    SpMat A(Eigen::Map<Eigen::MatrixXd>(AAA.memptr(), AAA.n_rows, AAA.n_cols).sparseView());
    Eigen::VectorXd D_diag = A * Eigen::VectorXd::Ones(A.cols());
    Eigen::MatrixXd D_diag_mat = D_diag.asDiagonal();
    SpMat D = D_diag_mat.sparseView();
    D -= A;
    // Calculate d = det(D[c(1:(dim(A)[1]-1)),c(1:(dim(A)[2]-1))]);
    Eigen::SparseLU<SpMat> solver;
    solver.compute(D.topLeftCorner(A.cols() - 1, A.cols() - 1));
    d = solver.logAbsDeterminant();
  }
  return(d);
}

/* Spanning tree number calculator */
double NT(arma::mat A, arma::vec y){
  int u = max(y) + 1;
  arma::mat L(u,u); L.zeros();
  arma::uvec cid;
  arma::uvec cid2;
  double d = 0.0;
  double d2;
  for(int k = 0; k < u; k++){
    cid = find(y == k);
    arma::mat As=A.submat(cid,cid);
    // Calculate district specific #trees
    d2 = LapDet(As);
    d = d + d2;
    for(int j = 0; j < u; j++){
      if(j > k){
    	cid2 = find(y == j);
    	L(k,j) = accu(A.submat(cid,cid2));
    	L(j,k) = L(k,j);
      }
    }
  }
  // Calculate #trees for the multigraph
  double d1 = LapDet(L);
  d = d + d1;
  return(d);
}

// Sample partition function
// [[Rcpp::export]]
List sample_partition(const std::vector<arma::vec> aList, const arma::mat aMat, const int num_partitions, const int num_samples, const int threads = 1){

  // Create storage - shared
  arma::mat store_parts(aList.size(), num_samples);
  arma::vec store_probs(num_samples);

  // Create objects intermediate - private
  std::vector<arma::vec> tree_out; arma::vec ups_out; double nt_out;

  // Convert to boost - firstprivate
  adjacency_list<> aList_boost = conv_List_boost(aList);

#ifdef _OPENMP
  omp_set_num_threads(threads);
  Rcpp::Rcout << "    Parallelizing partition sampling using OpenMP. "
	<< threads << " threads out of "
	<< omp_get_num_procs() << " are used."
	<< std::endl;
#endif
  for(int i = 0; i < num_samples; i++){
    tree_out = UST(aList_boost);
    ups_out = UPS(tree_out, num_partitions - 1);
    nt_out = NT(aMat, ups_out);
    store_parts.col(i) = ups_out;
    store_probs(i) = nt_out;
  }

  // Create output
  return(List::create(
		      _["partitions"] = store_parts,
		      _["prob_partitions"] = store_probs
		      ));

}

// /* Graph reduction given a population constraint*/
// List TR(List aList, NumericVector popvec, int numdist, double delta){//
//   // Convert aList TREE to boost object
//   adjacency_list<> mygraph=conv_List_boost(aList);
//   int ms = num_vertices(mygraph);
//   IntegerVector RV(ms);
//   std::fill( RV.begin(), RV.end(), 0) ;
//   // arma::vec RV;
//   // for(int i = 0; i < ms; i++){
//   //   RV(i)=0;
//   // }

//   //IntegerVector ind(ms);
//   //std::fill( ind.begin(), ind.end(), 1) ;
//   //std::cout << "ind(0): " << ind(0) ;
//   IntegerVector deg(ms);
//   IntegerVector K1(ms);

//   //IntegerVector EL(ms);
//   int popall;popall=0;
//   for(int j = 0; j < ms; j++){
//     popall += popvec(j);
//   }
//   double popmean;popmean=popall/numdist;
//   double popmin;
//   double popmax;
//   //temporary for unconstrained delta=1
//   if(delta==1){
//     popmin=0;
//     popmax=popall;
//   }else{
//     popmin=popmean*(1-delta);
//     popmax=popmean*(1+delta);}

//   random::uniform_int_distribution<> dist(0, ms-1);
//   //std::array<int> predeces;
//   array<int, 30> predeces;//array<int, 2000> predeces;
//   //std::mt19937 rg( (unsigned)time(NULL) );
//   static std::random_device rg;
//   int root = dist(rg);
//   random_spanning_tree(mygraph, rg,
// 		       predecessor_map(predeces.begin()).
// 		       root_vertex(root));
//   //std::cout << "root: " << root << '\n'; //Rcpp::Rcout <<
//   //std::cout << "predeces.size(): " << predeces.max_size() << '\n';

//   adjacency_list<> treegraph;
//   for (int i = 0; i < ms; i++) {
//     if(predeces[i]!=-1) {
//       int j = predeces[i];
//       //std::cout << "node: " << i ;
//       //std::cout << " predecessors[i]: " << j << '\n';
//       //EL(i)=j;
//       add_edge(i, j, treegraph);
//       add_edge(j, i, treegraph);
//     }
//   }

//   List aList2 = conv_boost_List(treegraph);

//   for(int i = 0; i < aList2.size(); i++){
//     // Get i'th entry in list
//     NumericVector list1;
//     list1 = aList2(i);
//     // Loop through elements in list1
//     deg(i) = list1.size();
//   }

//   IntegerVector K2(ms);
//   IntegerVector K3(ms);

//   int k = 0; NumericVector list1;
//   for(int i = 0; i < aList2.size(); i++){
//     if(deg(i) == 1){
//       k=k+1;  //label
//       RV(i) = k;//ind(i)=0;
//       int popk = popvec(i);

//       list1 = aList2(i);
//       int j = list1(0); //i's neighbors
//       if(deg(j)==2){
// 	RV(j) = k;//ind(j)=0;
// 	int K=1;
// 	int dj=deg(j);
//         while(dj == 2 && K>0){//K
//           list1 = aList2(j); int KK=1;
//           for(int jj = 0; jj < list1.size(); jj++){
//             //std::cout << "list1(jj): " << list1(jj) << '\n';
//             if(RV(list1(jj))==0){
// 	      //if(deg(list1(jj))==2){RV(list1(jj))=k;j=list1(jj);dj=2;KK=0;}
// 	      int popk2=popk+popvec(list1(jj));
// 	      if(deg(list1(jj))==2 && popk2<popmin){RV(list1(jj))=k;popk=popk+popvec(list1(jj));j=list1(jj);dj=2;KK=0;}
// 	      if(deg(list1(jj))>2){K3(k)=list1(jj);K=0;}
// 	    }
//             if(jj+1 == list1.size() && KK==1){K=0;}
//           }//for
//         }//while
// 	K2(k)=j;//K2.push_back(j);
//       }else{K2(k)=i;K3(k)=j;}//else{K2.push_back(i);K3.push_back(j);}
//       //if(deg(j)>2){K2(k)=j;}
//     }//if
//   }//END first-stage

//   //IntegerVector K02;
//   //IntegerVector K03;
//   //for(int i = 0; i < 10; i++){//max(RV)+1
//   //  K02(i) = K2(i);
//   //  K03(i) = K3(i);
//   //}
//   //K2 = K02;
//   //K3 = K03;

//   IntegerVector RV2(ms);
//   for(int i = 0; i < ms; i++){
//     RV2(i)=RV(i);
//   }

//   //arma::vec RV2;
//   //for(int i = 0; i < ms; i++){
//   //  RV2(i)=RV(i);
//   //}

//   arma::vec RV3 = as<arma::vec>(RV) ;
//   //arma::uvec rv0 = find(RV3 == 0);
//   //int numrv = max(RV)+rv0.n_elem;// need to fix

//   int mrv = max(RV);
//   arma::uvec cid;int pop;int jj;
//   for(int k=1; k < mrv+1; k++){
//     //std::cout << "k: " << k << "; K2(k): " << K2(k) << '\n';
//     for(int k2=1; k2 < mrv+1; k2++){
//       if(K3(k)==K3(k2) && k<k2){
// 	cid = arma::find(RV3 == k || RV3 == k2);
// 	pop = popvec(K3(k));
// 	for(int j = 0; j < cid.n_elem; j++){
// 	  pop += popvec(cid(j));
// 	}
// 	cid = arma::find(RV3 == k);

// 	if(pop<=popmin){
// 	  RV2(K3(k)) = k2;
// 	  for(int j = 0; j < cid.n_elem; j++){
// 	    RV2(cid(j)) = k2;
// 	    //numrv=numrv-1;
// 	  }
// 	}
//       }
//     }}

//   RV3 = as<arma::vec>(RV2) ;
//   IntegerVector RV4(ms);

//   IntegerVector RVV = sort_unique(RV2);arma::uvec ind;
//   for(int k=1; k < RVV.size(); k++){
//     ind = find(RV3==RVV(k));
//     for(int i=0; i < ind.n_elem; i++){
//       RV4(ind(i)) = k;
//     }
//   }
//   //Rcpp::NumericVector(RV.begin(), RV.end());
//   //Rcpp::NumericVector(RV2.begin(), RV2.end());
//   //int Rvm = max(RV);

//   List out;
//   out["RV"] = RV;
//   out["RV2"] = RV4;
//   out["K2"] = K2;
//   out["K3"] = K3;
//   out["aList2"] = aList2; //List
//   out["popmean"] = popmean;
//   out["popmin"] = popmin;
//   out["popmax"] = popmax;
//   out["numdist"] = numdist;
//   out["delta"] = delta;
//   out["popvec"] = popvec;
//   //out["numrv"] = numrv;
//   return(out);
//   //  return(RV);
// }

// /* Sample partitions on a reduced tree */
// List RS(List aList, List out, int iter){//
//   IntegerVector IND(iter);

//   int numdist = out["numdist"];
//   int cutnum = numdist-1;
//   NumericVector popvec = out["popvec"];
//   double popmin = out["popmin"];
//   double popmax = out["popmax"];
//   IntegerVector RV = out["RV2"];
//   int ms = aList.size();
//   //Rcpp::Rcout << "ms:  " << ms << '\n';
//   IntegerMatrix component2(ms,iter);
//   IntegerVector component3(ms,1);

//   int numrv=0;
//   NumericVector list1;
//   for(int i = 0; i < aList.size(); i++){
//     list1 = aList(i);
//     for(int j = 0; j < list1.size(); j++){
//       if(i > list1(j)){
// 	if((RV(i)!=RV(list1(j))) || (RV(i)==0 && RV(list1(j))==0))
// 	  {numrv = numrv + 1;}
//       }}}

//   arma::uvec cid;int pop;
//   //////////////////////////////////////////////////////////
//   arma::vec component(ms);//RV3 = as<arma::vec>(RV) ;
//   int ii = 0;
//   while(ii < iter){//IND < numdist &&
//     ii = ii + 1;
//     IntegerVector vv = sample_int(Range(0,numrv-1),cutnum); //
//     //int ntr = numrv-1;
//     //Rcpp::Rcout << "numrv-1:  " << ntr << '\n';
//     //for(int l=0; l < vv.size(); l++){
//     //Rcpp::Rcout << "vv(l): " << vv(l) << '\n';
//     //}
//     int k = -1;
//     int ll = 0;
//     adjacency_list<> mygraph;
//     // Loop through vectors in aList
//     for(int i = 0; i < aList.size(); i++){
//       add_vertex(mygraph);
//       // Get i'th entry in list
//       //NumericVector list1;
//       list1 = aList(i);
//       // Loop through elements in list1
//       for(int j = 0; j < list1.size(); j++){
// 	if((RV(i)!=RV(list1(j))) || (RV(i)==0 && RV(list1(j))==0)){ll=1;if (i > list1(j)){k++;}}else{ll=0;}
//         if (i > list1(j)){
// 	  int kk = 0;
// 	  for(int l=0; l < vv.size(); l++){
// 	    if (vv[l]==k && ll==1){
// 	      kk=1;
// 	    }
// 	  }
// 	  if (kk == 0){
// 	    //Rcpp::Rcout << "k: " << k << " ll: " << ll << " add_edge i:" << i << "  list1(j): " <<  list1(j) << '\n';
// 	    add_edge(i,list1(j),mygraph);
// 	    add_edge(list1(j),i,mygraph);
// 	  }else{
// 	    //Rcpp::Rcout << "k: " << k << " ll: " << ll << " cut_edge i:" << i << "  list1(j): " <<  list1(j) << '\n';
// 	  }
// 	}
//       }
//     }
//     //IntegerVector component(ms);
//     //Rcpp::Rcout << "num_vertices (mygraph):  " << num_vertices (mygraph) << '\n';
//     //Rcpp::Rcout << "component:  " << component << '\n';
//     int num_components = connected_components(mygraph, &component[0]);


//     IND(ii-1)=0;
//     for(int i = 0; i < component.max()+1; i++){
//       cid = arma::find(component == i);
//       pop = 0;
//       for(int j = 0; j < cid.n_elem; j++){
//         pop += popvec(cid(j));
//       }
//       if(pop>=popmin && pop<=popmax){IND(ii-1) = IND(ii-1)+1;}
//     }
//     //Rcpp::Rcout << "IND:  " << IND << '\n';
//     component3 = wrap(component);
//     component2(_,ii-1) = component3;
//   }//WHILE
//   //////////////////////////////////////////////////////////

//   //  IntegerVector component2(ms);
//   //  Rcpp::Rcout << "num_vertices (mygraph):  " << num_vertices (mygraph) << '\n';
//   //  int num_components = connected_components(mygraph, &component2[0]);

//   List out2;
//   out2["RV"] = component2;
//   out2["IND"] = IND;
//   return(out2);
// }
