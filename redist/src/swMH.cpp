/* Swendsen-Wang Algorithm
   Metropolis-Hastings Implementation */

#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
#include <igraph/igraph.h>

using namespace arma;

// Function to get RPI
// [[Rcpp::export]]
Rcpp::NumericVector getRPI(Rcpp::NumericMatrix plans,
			   Rcpp::NumericVector origplan,
			   Rcpp::NumericVector pop,
			   Rcpp::NumericMatrix ssdmat){

  // Get the unique cd labels
  Rcpp::NumericVector labs = unique(plans(Rcpp::_,0));

  // Get original SSD
  // Outer loop over districts
  arma::vec getCds = Rcpp::as<arma::vec> (origplan);
  double ssdorig = 0.0;
  for(int i = 0; i < labs.size(); i++){

    // Initialize sum, find precincts
    double origssd = 0.0;
    arma::uvec findprecs = find(getCds == i);

    // Loop over elements in the district
    for(int j = 0; j < findprecs.size(); j++){

      // Create vector of pwd's relative to precinct findprecs(j)
      for(int k = j + 1; k < findprecs.size(); k++){
	origssd += (double)ssdmat(findprecs(j), findprecs(k)) * 
	  pop(findprecs(j)) * pop(findprecs(k));
      }

    }

    // Add to ssdorig
    ssdorig += origssd;

  }

  // Loop over plans, get RPI
  Rcpp::NumericVector rpivec(plans.ncol());
  for(int i = 0; i < plans.ncol(); i++){

    // Get plan
    arma::vec cdPlan = plans(Rcpp::_,i);

    // Initialize rpi
    double rpi = 0.0;

    // Loop over cd 
    for(int j = 0; j < labs.size(); j++){

      // Initialize rpi
      double ssd = 0.0;

      // Get precincts in cd j
      uvec findcds = find(cdPlan == j);

      // SSD for this precinct
      for(int k = 0; k < findcds.size(); k++){
	for(int l = k + 1; l < findcds.size(); l++){
	  ssd += (double)ssdmat(findcds(k), findcds(l)) * 
	    pop(findcds(k)) * pop(findcds(l));
	}
      }

      // Add ssd to rpi
      rpi += (double)ssd / ssdorig;

    }
  
    // Add to vector
    rpivec(i) = rpi;

  }

  return rpivec;

}

// Function to calculate % deviation from original plan
// [[Rcpp::export]]
Rcpp::NumericVector devOrig(Rcpp::NumericMatrix mat,
			    Rcpp::NumericVector origplan){

  // Get size of original plan
  int nprec = origplan.size();

  // Object to keep track of similarity
  int count;

  // Object to store similarity
  Rcpp::NumericVector devStore(mat.ncol());

  // Loop over plans
  for(int i = 0; i < mat.ncol(); i++){

    // Get a plan
    Rcpp::NumericVector partition = mat(Rcpp::_,i);

    // Initialize count
    count = 0;

    // Loop over precincts to compare
    for(int j = 0; j < nprec; j++){

      if(partition(j) == origplan(j)){
	count++;
      }
      
    }

    // Convert to %
    double devplan = (double)count / origplan.size();

    // Store deviation for the plan
    devStore(i) = devplan;

  }

  return devStore;

}

// Function to calculate deviation from parity from cd matrix
// [[Rcpp::export]]
Rcpp::NumericVector distParity(Rcpp::NumericMatrix mat,
			       Rcpp::NumericVector popvec,
			       double parity){

  // Get the unique cd labels
  Rcpp::NumericVector labs = unique(mat(Rcpp::_,0));

  // Container vector of plan deviations
  Rcpp::NumericVector plandevs(mat.ncol());

  // Loop through plans
  for(int i = 0; i < mat.ncol(); i++){

    // Get plan
    arma::vec plan = mat(Rcpp::_,i);

    // Loop through assignments
    double maxdev = 0.0;
    for(int j = 0; j < labs.size(); j++){

      uvec assignments = find(plan == labs(j));

      // Loop over precincts in plan
      int distpop = 0;
      for(int k = 0; k < assignments.size(); k++){

	distpop += popvec(assignments(k));

      }

      // Get deviation from parity
      double plandev = std::abs((double)((double)distpop - parity) / parity);
      if(plandev > maxdev){
	maxdev = plandev;
      }

    }

    // Store maxdev in plandevs
    plandevs(i) = maxdev;

  }

  return plandevs;

}

// Function to generate adjacency graph and count clusters
int genGraph(Rcpp::List aList){

  // Declare graph, igraph vector
  igraph_t g;
  igraph_vector_t v;

  // Resize vector to accomodate all edges
  int vecsize = 0;
  for(int i = 0; i < aList.size(); i++){
    Rcpp::NumericVector adj = aList(i);
    vecsize += adj.size();
  }
  igraph_vector_init(&v, 2 * vecsize);

  // Create vector of edges
  int ind = 0;
  for(int i = 0; i < aList.size(); i++){
    Rcpp::NumericVector adj = aList(i);
    for(int j = 0; j < adj.size(); j++){
      VECTOR(v)[ind] = i;
      ind++;
      VECTOR(v)[ind] = adj(j);
      ind++;
    }
  }

  // Create full graph
  igraph_create(&g, &v, aList.size(), 0);

  // Objects to calculate clusters
  igraph_integer_t n;

  igraph_clusters(&g, NULL, NULL, &n, IGRAPH_WEAK);

  // Destroy objects to free memory
  igraph_vector_destroy(&v);
  igraph_destroy(&g);

  return(n);

}

// Function to calculate summary statistics
// [[Rcpp::export]]
Rcpp::NumericMatrix sumstat(Rcpp::NumericMatrix distmat,
			    Rcpp::NumericVector hisp,
			    Rcpp::NumericVector afam,
			    Rcpp::NumericVector dem,
			    Rcpp::NumericVector rep,
			    Rcpp::NumericVector pop){

  // Matrix to hold dissimilarity indices
  Rcpp::NumericMatrix diMat(distmat.ncol(), 4);

  // Population parameters
  int T = sum(pop);
  double P_h = (double)sum(hisp) / T;
  double P_a = (double)sum(afam) / T;
  double P_d = (double)sum(dem) / T;
  double P_r = (double)sum(rep) / T;

  // Denominators
  double d_h = (double)1 / (2 * T * P_h * (1 - P_h));
  double d_a = (double)1 / (2 * T * P_a * (1 - P_a));
  double d_d = (double)1 / (2 * T * P_d * (1 - P_d));
  double d_r = (double)1 / (2 * T * P_r * (1 - P_r));

  // Get number of unique plans
  Rcpp::NumericVector cd1 = distmat(Rcpp::_,1);
  arma::vec cdVec1 = Rcpp::as<arma::vec> (cd1);
  arma::vec cdLabs = arma::unique(cdVec1);

  // Range to look over for cds
  int start;
  int end = max(cd1) + 1;
  if(min(cd1) == 1){
    start = 1;
  } else{
    start = 0;
  }

  // Loop over possible plans
  for(int i = 0; i < distmat.ncol(); i++){

    // Create dissimilarity objects
    double dissim_h = 0;
    double dissim_a = 0;
    double dissim_d = 0;
    double dissim_r = 0;

    // Get a plan
    Rcpp::NumericVector cdvec = distmat(Rcpp::_,i);
    arma::vec cds = Rcpp::as<arma::vec> (cdvec);

    // Loop over congressional districts
    for(int j = start; j < end; j++){

      // Initialize counts of groups
      int tpop = 0;
      int hpop = 0;
      int apop = 0;
      int dpop = 0;
      int rpop = 0;

      // Which precincts in the plan are in this cd?
      uvec findCds = find(cds == j);

      // Loop over precincts
      for(int k = 0; k < findCds.size(); k++){

	// Add population counts
	tpop += pop(findCds(k));
	hpop += hisp(findCds(k));
	apop += afam(findCds(k));
	dpop += dem(findCds(k));
	rpop += rep(findCds(k));

      }

      // Get district proportions
      double p_h = (double)hpop / tpop;
      double p_a = (double)apop / tpop;
      double p_d = (double)dpop / tpop;
      double p_r = (double)rpop / tpop;

      // Add to dissimilarity index
      dissim_h += (double)d_h * tpop * std::abs(p_h - P_h);
      dissim_a += (double)d_a * tpop * std::abs(p_a - P_a);
      dissim_d += (double)d_d * tpop * std::abs(p_d - P_d);
      dissim_r += (double)d_r * tpop * std::abs(p_r - P_r);

    }
    
    // Add to matrix
    diMat(i,0) = dissim_h;
    diMat(i,1) = dissim_a;
    diMat(i,2) = dissim_d;
    diMat(i,3) = dissim_r;

  }

  // Return matrix object
  return diMat;

}

// Function that applies the Geyer Thompson algorithm to anneal beta
// [[Rcpp::export]]
double changeBeta(arma::vec betavec,
		  double beta,
		  double constraint,
		  Rcpp::NumericVector weights){

  // Find beta in betavec
  arma::uvec findBetaVec = find(betavec == beta);
  int findBeta = findBetaVec(0);

  // Object to test whether beta is at RHS of vector
  int betaLoc = betavec.size() - 1;

  // Get transition probabilities and propose a new beta
  double qij;
  double qji;
  double wi;
  double wj;
  double propBeta;
  if(findBeta == 0){ // At first element in betavec
    qij = 1;
    qji = .5;
    wi = weights(0);
    wj = weights(1);
    propBeta = betavec(1);
  } else if(findBeta == betaLoc){ // At last element in betavec
    qij = 1;
    qji = .5;
    wi = weights(betaLoc);
    wj = weights(betaLoc - 1);
    propBeta = betavec(betaLoc - 1);
  } else{ // Anywhere in the middle of betavec
    qij = .5;
    qji = .5;
    wi = weights(findBeta);
    vec betaswitch = Rcpp::runif(1);
    if(betaswitch(0) < .5){
      propBeta = betavec(findBeta - 1);
      wj = weights(findBeta - 1);
    }
    if(betaswitch(0) >= .5){
      propBeta = betavec(findBeta + 1);
      wj = weights(findBeta + 1);
    }
  }

  // Accept or reject the proposal
  double mhprobGH = (double)((double)exp((double)propBeta * constraint) / 
			     exp((double)beta * constraint)) * ((double)wj / wi) * ((double)qji / qij);
  if(mhprobGH > 1){
    mhprobGH = 1;
  }
  vec testkeepGH = Rcpp::runif(1);
  if(testkeepGH(0) <= mhprobGH){
    beta = propBeta;
  }

  return(beta);

}

// [[Rcpp::export]]
Rcpp::NumericVector pBias(Rcpp::NumericVector ndem,
			  Rcpp::NumericVector nrep,
			  Rcpp::NumericMatrix cdMat,
			  double swing) {

  // Convert swing object
  double swingpct = std::abs(swing);
  
  // Columns of votes with swing
  Rcpp::NumericVector dems(ndem.size());
  Rcpp::NumericVector reps(nrep.size());
  
  // Output vector
  Rcpp::NumericVector rSeats(cdMat.ncol());

  // Loop over columns of cd assignments
  for(int c = 0; c < cdMat.ncol(); c++){
    
    // Get c'th col of cdmat
    Rcpp::NumericVector cds = cdMat(Rcpp::_,c);

    // Convert to arma objects
    arma::vec cdVec = Rcpp::as<arma::vec> (cds);
    
    // Get vector of unique cd labels
    arma::vec cdLabs = arma::unique(cdVec);
    
    // Induce partisan swing
    // Positive is pro-Democrat
    if(swing > 0){
      for(int i = 0; i < nrep.size(); i++){
	int ssize = round(nrep(i) * swingpct);
	reps(i) = nrep(i) - ssize;
	dems(i) = ndem(i) + ssize;
      }
    }
    if(swing < 0){
      for(int i = 0; i < ndem.size(); i++){
	int ssize = round(ndem(i) * swingpct);
	dems(i) = ndem(i) - ssize;
	reps(i) = nrep(i) + ssize;
      }
    }
    if(swing == 0){
      dems = ndem;
      reps = nrep;
    }
    
    // Calculate district vote totals  
    int nSeats = 0;
    for(int i = 0; i < cdLabs.size(); i++){
      
      // Find districts
      arma::uvec precs = find(cdVec == cdLabs(i));
      int dvotes = 0;
      int rvotes = 0;
      
      // Sum dem, rep votes for that district
      for(int j = 0; j < precs.size(); j++){
	dvotes += dems(precs(j));
	rvotes += reps(precs(j));
      }
      
      // If rvotes > dvotes...
      if(rvotes - dvotes > 0){
	nSeats++;
      }
    }

    rSeats(c) = nSeats;

  }  

  return rSeats;
  
}
  
// Check to see if all ties are reciprocated
// [[Rcpp::export]]
int checkRecip(Rcpp::List aList){
  
  int notie = 0;
  for(int i = 0; i < aList.size(); i++){
    Rcpp::NumericVector adj1 = aList(i);

    for(int j = 0; j < adj1.size(); j++){
      Rcpp::NumericVector adj2 = aList(adj1(j));
      if(is_true(any(adj2 == i)) == FALSE){
	notie++;
      }
    }

  }

  return notie;

}

// Add ties that were left out
// [[Rcpp::export]]
Rcpp::List addTies(Rcpp::List aList){
  
  for(int i = 0; i < aList.size(); i++){
    Rcpp::NumericVector list1 = aList(i);
    
    for(int j = 0; j < list1.size(); j++){
      Rcpp::NumericVector list2 = aList(list1(j));
      
      if(is_true(any(list2 == i)) == FALSE){
	list2.push_back(i);
	aList(list1(j)) = list2;
      }
    
    }
 
  }

  return aList;

}

/* This function can be used to create 
   the connected components from each cutedge list
   when generating workable data */
// [[Rcpp::export]]
Rcpp::List genCC(Rcpp::List aList,
		 Rcpp::List cutList) {
  
  int units = aList.size();
  Rcpp::List alConnected(units);
  
  // Creat adjacency list
  for(int i = 0; i < units; i++){
    
    Rcpp::NumericVector avec = aList(i);
    Rcpp::NumericVector cut = cutList(i);

    Rcpp::NumericVector adj;
    for(int j = 0; j < avec.size(); j++){
      if(is_true(any(cut == avec(j))) == FALSE){
	adj.push_back(avec(j));
      }
    }

    alConnected(i) = adj;

  }

  return alConnected;

}

// Another cc function with different inputs
// [[Rcpp::export]]
Rcpp::List genCC2(Rcpp::List aList,
		  Rcpp::NumericVector cds,
		  int zeroInd = 1) {

  int units = aList.size();
  Rcpp::List alConnected(units);

  // Create adjacency list
  for(int i = 0; i < units; i++){
    
    Rcpp::NumericVector avec = aList(i);
    int cd = cds(i);
    
    Rcpp::NumericVector adj;
    for(int j = 0; j < avec.size(); j++){
      if(cds(avec(j)) == cd){
	if(zeroInd == 1){
	  adj.push_back(avec(j));
	}
	if(zeroInd == 0){
	  adj.push_back(avec(j) + 1);
	}
      }
    }

    alConnected(i) = adj;

  }

  return alConnected;

}

// [[Rcpp::export]]
Rcpp::List bSearchCppBound(Rcpp::List adjList,
			  Rcpp::NumericVector tractnum,
			  Rcpp::NumericVector boundaryindex) {
  
  // Convert vec to arma object
  arma::vec tNum = Rcpp::as<arma::vec> (tractnum);
  arma::vec bInd = Rcpp::as<arma::vec> (boundaryindex);
  arma::uvec findBound1 = find(bInd == 1);
  Rcpp::NumericVector findBound(findBound1.n_elem);
  for(int i = 0; i < findBound1.n_elem; i++){
    findBound(i) = findBound1(i);
  }

  // Set up container objects
  Rcpp::List bSearch;
  int j = 0;
  int lengthT = tNum.n_elem;
  //int lengthB = findBound.size();
  Rcpp::NumericVector mark(lengthT);
  Rcpp::NumericVector partition;

  // Initialize w first element in findBound
  int indexS = findBound(0);
  mark(indexS) = tNum(indexS);
  partition.push_back(tNum(indexS));
  Rcpp::NumericVector q;
  Rcpp::NumericVector initQ = adjList(indexS);
  int sizeQ = initQ.size();
  for(int i = 0; i < sizeQ; i++){
    int push = initQ(i);
    q.push_back(push);
  }

  do{

    while(q.size() > 0){

      int u = q(0);
      arma::uvec vecIndU = find(tNum == u);
      int indexU = vecIndU(0);
      mark(indexU) = u;
      bool inPart = is_true(any(partition == u));
      if(inPart == FALSE){
	partition.push_back(u);
      }  
      
      // Ensure adj vector has only unique elements
      Rcpp::NumericVector adjF = adjList(indexU);
      Rcpp::NumericVector adj;
      int adjFLength = adjF.size();
      for(int i = 0; i < adjFLength; i++){
	int findDup = adjF(i);
	bool dup = is_true(any(adj == findDup));
	if(dup == FALSE){
	  adj.push_back(findDup);
	}
      }
      
      int adjLength = adj.size();
      if(adjLength > 0){
	
	for(int i = 0; i < adjLength; i++){
	  int v = adj(i);
	  arma::uvec vecIndV = find(tNum == v);
	  int indexV = vecIndV(0);
          bool inMark = is_true(any(mark == v));
          if(inMark == FALSE){
	    mark(indexV) = v;
	    partition.push_back(v);
	    q.push_back(v);
          }
	}
	
      }
      
      q.erase(q.begin());

    }

    if(q.size() == 0){

      /* Find elements in partition that are in findBound
	 If they are in partition, remove from findBound */
      for(int i = findBound.size() - 1; i >= 0; i--){
	bool testPart = is_true(any(partition == findBound(i)));
	if(testPart == TRUE){
	  findBound.erase(i);
	}	
      }

      // Store partition
      bSearch.push_back(partition);
      j++;
      partition.erase(partition.begin(), partition.end());

      if(findBound.size() > 0){

	indexS = findBound(0);

	// From new seed, generate new queue
        Rcpp::NumericVector newQ = adjList(indexS);
	int newQL = newQ.size();
	for(int i = 0; i < newQL; i++){
	  int push = newQ(i);
	  bool inMark = is_true(any(mark == push));
	  if(inMark == FALSE){
	    q.push_back(push);
	  }
	}
	
	mark(indexS) = tNum(indexS);
	
	partition.push_back(tNum(indexS));

      }

    }

  }while(findBound.size() > 0);

  return bSearch;

}

// [[Rcpp::export]]
int bSearchPerimCpp(Rcpp::List adjList, 
		    Rcpp::NumericVector perimVec,
		    Rcpp::NumericVector tractnum) {
  
  // Convert vec to arma object
  arma::vec tNum = Rcpp::as<arma::vec> (tractnum);
  
  // Output - 1 if all are reached in one go, 0 if not
  int allreached = 0;
  
  /* Set up container objects
     search.list contains each connected component
     mark keeps track of units that have been reached
     partition stores the current connected component
     until stored in search.list */
  Rcpp::List bSearch;
  int j = 0;
  int length = tNum.n_elem;
  Rcpp::NumericVector mark(length);
  Rcpp::NumericVector partition;

  /* Initializing - start with first element in tNum
     q is a queue of points that have been reached */
  int indexS = 0;
  mark(perimVec(indexS)) = perimVec(indexS);
  partition.push_back(perimVec(indexS));
  Rcpp::NumericVector q;
  Rcpp::NumericVector initQ = adjList(perimVec(indexS));
  for(int i = 0; i < initQ.size(); i++){
    if(is_true(any(perimVec == initQ(i))) == TRUE){
      q.push_back(initQ(i));
    }
  }
    
  // Begin finding connected components
  while(q.size() > 0){
    
    int u = q(0);
    arma::uvec vecIndU = find(tNum == u);
    int indexU = vecIndU(0);
    mark(indexU) = u;
    bool inPart = is_true(any(partition == u));
    if(inPart == FALSE){
      partition.push_back(u);
    }  
    
    // Ensure adj vector has only unique elements of perimeter units
    Rcpp::NumericVector adjF = adjList(indexU);
    Rcpp::NumericVector adj;
    for(int i = 0; i < adjF.size(); i++){
      
      if(is_true(any(perimVec == adjF(i))) == TRUE && 
	 is_true(any(mark == adjF(i))) == FALSE) {
	adj.push_back(adjF(i));
      }
      
    }
    
    // Adding adjacent perimeter units to queue
    if(adj.size() > 0){
      
      for(int i = 0; i < adj.size(); i++){
	int v = adj(i);
	arma::uvec vecIndV = find(tNum == v);
	int indexV = vecIndV(0);
	bool inMark = is_true(any(mark == v));
	if(inMark == FALSE){
	  mark(indexV) = v;
	  partition.push_back(v);
	  q.push_back(v);
	}
      }
      
    }
    
    q.erase(q.begin());
    
  }
  
  // Once queue is empty, count number of adjvec in mark
  int checkAdj = 0;
  for(int i = 0; i < perimVec.size(); i++){
    if(is_true(any(mark == perimVec(i))) == TRUE){
	checkAdj++;
    }
  }
  
  if(checkAdj == perimVec.size()){
    allreached++;
  }
  
  return allreached;

}

// Generates boundary indicator vector from cutedge list
// [[Rcpp::export]]
Rcpp::NumericVector getBoundary(Rcpp::List list) {
  int units = list.size();
  Rcpp::NumericVector boundary(units);
  for(int i = 0; i < units; i++){
    Rcpp::NumericVector testCut = list(i);

    if(testCut.size() > 0){
      boundary(i) = 1;
    }
  }
  
  return boundary;
  
}

// Generates boundary indicator vector by comparing full and connected list
// [[Rcpp::export]]
Rcpp::NumericVector getBoundary2(Rcpp::List fullList,
				 Rcpp::List conList){

  int units = fullList.size();
  Rcpp::NumericVector boundary(units);
  for(int i = 0; i < units; i++){

    Rcpp::NumericVector full = fullList(i);
    Rcpp::NumericVector conn = conList(i);

    if(full.size() > conn.size()){
	boundary(i) = 1;
    }
      
  }
  
  return boundary;
  
}

// Function to get adjacent congressional districts for connected component
// [[Rcpp::export]]
Rcpp::NumericVector getCds(Rcpp::List fullList,
			   Rcpp::NumericVector boundInds,
			   Rcpp::NumericVector cds) {

  int thisCd = cds(boundInds(0));
  Rcpp::NumericVector candCds;
  for(int i = 0; i < boundInds.size(); i++){
    
    Rcpp::NumericVector adjVec = fullList(boundInds(i));
    
    for(int j = 0; j < adjVec.size(); j++){
      
      int cd = cds(adjVec(j));
      if(cd != thisCd && is_true(any(candCds == cd)) == FALSE){
	candCds.push_back(cd);
      }

    }

  }

  return candCds;

}

// [[Rcpp::export]]
Rcpp::List genCutC(Rcpp::List aList, 
		     Rcpp::NumericVector cds) {
  
  int units = cds.size();
  
  Rcpp::List alCutedge(units);

  for(int i = 0; i < units; i++){
    
    // Adj list for tract i
    Rcpp::NumericVector listSub = aList(i);
    int cd = cds(i);
    
    // Vector of edges cut
    Rcpp::NumericVector cutedgeUnit;
    for(int j = 0; j < listSub.size(); j++){
      int testCd = listSub(j);
      if(cds(testCd) != cd){
	cutedgeUnit.push_back(testCd);
      }

    }

    alCutedge(i) = cutedgeUnit;

  }

  return alCutedge;

}

// Feed vector of precincts in, get back out those adjacent to specific congressional district
// [[Rcpp::export]]
Rcpp::NumericVector getAdj(Rcpp::List aList,
			   Rcpp::NumericVector unitVec,
			   Rcpp::NumericVector cdVec,
			   int cd) {
  
  Rcpp::NumericVector adjUnit;
  for(int i = 0; i < unitVec.size(); i++){
    
    Rcpp::NumericVector adjVec = aList(unitVec(i));
    
    for(int j = 0; j < adjVec.size(); j++){
      if(cdVec(adjVec(j)) == cd){
	adjUnit.push_back(unitVec(i));
	break;
      }
    }

  }

  return adjUnit;

}

// [[Rcpp::export]]
Rcpp::NumericVector popCalc(Rcpp::NumericVector switchPartition,
			    Rcpp::NumericVector pops,
			    Rcpp::NumericVector distPops,
			    int oldCd,
			    int newCd) {
  
  Rcpp::NumericVector newDistPops = Rcpp::clone(distPops);
  for(int i = 0; i < switchPartition.size(); i++){
    newDistPops(oldCd) -= pops(switchPartition(i));
    newDistPops(newCd) += pops(switchPartition(i));
  }

  return newDistPops;

}

// Function to get adjacency list of cd's
// [[Rcpp::export]]
Rcpp::List getAlCd(Rcpp::List aList,
		   Rcpp::NumericVector cds){

  // Get num of cds
  Rcpp::NumericVector ncds = unique(cds);

  // Detemrine boundary units to initialize
  Rcpp::List alConnected = genCC2(aList, cds);
  Rcpp::NumericVector boundary = getBoundary2(aList, alConnected);
  arma::vec boundary1 = Rcpp::as<arma::vec> (boundary);
  uvec bIndex = find(boundary1 == 1);
  
  // Create adjacency list for congressional districts
  Rcpp::List alCd(ncds.size());
  for(int i = 0; i < alCd.size(); i++){
    alCd(i) = 1000000;
  }
  for(int i = 0; i < bIndex.size(); i++){

    // Get congressional district and adjacent districts of boundary unit i
    int cdInd = cds(bIndex(i));
    Rcpp::NumericVector adj = aList(bIndex(i));

    // Get current adjacent cds for cdInd
    Rcpp::NumericVector cdAdj = alCd(cdInd);
    
    // Gather the congressional districts of adjacent precincts
    Rcpp::NumericVector cdVec(adj.size());
    for(int j = 0; j < adj.size(); j++){
      cdVec(j) = cds(adj(j));
    }

    // Add unique units not equal to cd of precinct i
    for(int j = 0; j < cdVec.size(); j++){
      if(cdVec(j) != cdInd & is_false(any(cdAdj == cdVec(j))) == TRUE){
	cdAdj.push_back(cdVec(j));
      }
    }

    // Add back to alCd
    alCd(cdInd) = cdAdj;
    
  }

  // Remove entries of 1000000
  for(int i = 0; i < alCd.size(); i++){
    Rcpp::NumericVector adj = alCd(i);
    Rcpp::NumericVector addVec;
    for(int j = 0; j < adj.size(); j++){
      if(adj(j) != 1000000){
	addVec.push_back(adj(j));
      }
    }
    alCd(i) = addVec;
  }

  return(alCd);

}

// [[Rcpp::export]]
Rcpp::List swMH(Rcpp::List aList,
		Rcpp::NumericVector cds, 
		Rcpp::NumericVector origplan,
		int nrep,  
		double eprob, 
		Rcpp::NumericVector pops, 
		Rcpp::NumericVector grouppop,
		int parity, 
		int margin, 
		int dists,
		int lambda, 
		Rcpp::NumericMatrix ssd, 
		double beta, 
		double betadiss, 
		double betapop,
		double betaswitch,
		Rcpp::NumericVector betavec,
		Rcpp::NumericVector betadissvec,
		Rcpp::NumericVector betapopvec,
		Rcpp::NumericVector betaswitchvec,
		Rcpp::NumericVector betaweights,
		int origpartcompact = 0,
		int annealbeta = 0,
		int annealbetadiss = 0,
		int annealbetapop = 0,
		int annealbetaswitch = 0,
		int elimCheck = 1
		) {

  int units = cds.size();

  if(min(cds) == 1){
    for(int i = 0; i < units; i++){
      cds(i)--;
    }
  }

  // Get vector of ssd's for the original district    
  // Outer loop over districts
  arma::vec getCds = Rcpp::as<arma::vec> (origplan);
  Rcpp::NumericVector ssdDists(dists);
  for(int i = 0; i < dists; i++){

    // Initialize sum, find precincts
    double origssd = 0;
    arma::uvec findprecs = find(getCds == i);

    // Loop over elements in the district
    for(int j = 0; j < findprecs.size(); j++){

      // Create vector of pwd's relative to precinct findprecs(j)
      for(int k = j + 1; k < findprecs.size(); k++){
	origssd += ssd(findprecs(j), findprecs(k)) * pops(findprecs(j)) * 
	  pops(findprecs(k));
      }

    }

    // Add to ssdorig
    ssdDists(i) = origssd;

  }

  // Generate al.cutedge
  Rcpp::List alCutedge = genCutC(aList, cds);
  
  // Index for breadth search
  Rcpp::NumericVector index(units);
  for(int i = 0; i < units; i++){
    index(i) = i;
  }

  // Create initial vector of district deviation from parity
  arma::vec getInitPlan = Rcpp::as<arma::vec> (cds);
  Rcpp::NumericVector distParity(dists);
  for(int i = 0; i < dists; i++){
    uvec indices = find(getInitPlan == i);
    int pop = 0;
    for(int j = 0; j < indices.size(); j++){
      pop += pops(indices(j));
    }
    pop -= parity;
    distParity(i) = pop;
  }

  // Create ecuts object, store output of alg here
  Rcpp::NumericMatrix cdDist(units, nrep);
  Rcpp::NumericVector numAcc(nrep);
  Rcpp::NumericVector numCC(nrep);
  Rcpp::NumericVector rejSplit(nrep);
  Rcpp::NumericVector rejPar(nrep);
  Rcpp::NumericVector rejAdj(nrep);
  Rcpp::NumericVector drawLam(nrep);
  Rcpp::NumericVector betastore(nrep);
  Rcpp::NumericVector betaDisstore(nrep);
  Rcpp::NumericVector betaPopstore(nrep);
  Rcpp::NumericVector betaSwitchstore(nrep);
  Rcpp::NumericVector diffBetastore(nrep);
  Rcpp::NumericVector diffDisstore(nrep);
  Rcpp::NumericVector diffPopstore(nrep);
  Rcpp::NumericVector diffSwitchstore(nrep);
  Rcpp::NumericVector decisionCount(nrep);
  int ncc;
  int rsplit;
  int rpar;
  int radj;
  int rlam;
  int k = 0;
  int z = 0;
  int store = 0;

  // Threshold for edge cut
  double cutBound = 1 - eprob;

  // Create denominator of dissimilarity index
  int tAll = 0;
  int tGroup = 0;  
  for(int i = 0; i < pops.size(); i++){
    tAll += pops(i);
    tGroup += grouppop(i);
  }
  double pAll = (double)tGroup / tAll;
  double denom = (double)1 / (2 * tAll * pAll * (1 - pAll));

  // Create objects for constraints
  double diffDissim;
  double distsumssd;
  double newrpi;
  double oldrpi;
  // double diffPop;
  double diffSwitch;
  
  double psiNewPop; // These two are to test whether we've been accepting incorrectly
  double psiOldPop;

  // This do loop iterates through the number of simulations nrep
  while(k < nrep){
    
    //////////////////////////////////////////
    // First - determine the boundary cases //
    //////////////////////////////////////////

    Rcpp::List alConnected = genCC2(aList, cds);

    Rcpp::NumericVector boundary = getBoundary(alCutedge);

    /////////////////////////////////////////////////////////////
    // Second - turn off edges with some probability within CD //
    /////////////////////////////////////////////////////////////

    // Declare objects needed for do loop
    Rcpp::List alConnectedPost(units);
    Rcpp::List alCutedgePost(units);
    Rcpp::NumericVector candCd;
    Rcpp::NumericVector switchPartition;
    Rcpp::NumericVector cdProp;
    Rcpp::NumericVector newCds;
    Rcpp::NumericVector distParityProp;
    int breakOuter;
    int propSwitch;
    int nccCheck;
    int thisPart;
    int actccs;
    double mhprob;	

    // Insert outer repeat loop here (1) - break when breakOuter == 1
    do{
      
      actccs = 0;
      mhprob = 1;

      // This loop cuts edges, generates alConnectedPost and alCutedgePost
      for(int i = 0; i < units; i++){
	
	Rcpp::NumericVector connected1 = alConnected(i);
	Rcpp::NumericVector connected;
	
	// First, subset connected to elements greater than i
	for(int j = 0; j < connected1.size(); j++){
	  if(connected1(j) > i){
	    connected.push_back(connected1(j));
	  }
	}
	
	// Take random draws from [0,1] uniform. If
	// draw is less than 1 - eprob, cut the edge.
	// If greater than 1 - eprob, keep the tie
	vec testkeep = Rcpp::runif(connected.size());
	Rcpp::NumericVector keep;
	Rcpp::NumericVector cut;
	for(int j = 0; j < connected.size(); j++){
	  if(testkeep(j) < cutBound){
	    cut.push_back(connected(j));
	  } else{
	    keep.push_back(connected(j));
	  }
	}
	
	// Here, ensure both directions of edge are cut
	for(int j = 0; j < cut.size(); j++){	    
	  Rcpp::NumericVector alConnCut = alConnected(cut(j));
	  arma::vec findCut = Rcpp::as<arma::vec> (alConnCut);
	  arma::uvec toCut = find(findCut == i);
	  
	  if(toCut.size() > 0){
	    alConnCut.erase(toCut(0));
	    alConnected(cut(j)) = alConnCut;
	  }
	  
	}
	
	// Define keep list
	alConnectedPost(i) = keep;
	  
	// Define cut list
	Rcpp::NumericVector cutedgeI = alCutedge(i);
	for(int j = 0; j < cut.size(); j++){
	  bool testCut = is_true(any(cutedgeI == cut(j)));
	  if(testCut == FALSE){
	    cutedgeI.push_back(cut(j));
	  }
	}
	
	alCutedgePost(i) = cutedgeI;
	
      }
      
      ////////////////////////////////////////////////////////////////
      // Third - search to find connected components within each cd //
      ////////////////////////////////////////////////////////////////

      Rcpp::List searchListBound = bSearchCppBound(alConnectedPost, 
						   index, 
						   boundary);
      ncc = searchListBound.size();

      //////////////////////////////////////////////////////////////////////////
      // Fourth - select a connected component using the uniform distribution //
      //////////////////////////////////////////////////////////////////////////

      // Objects for multiple swaps
      rsplit = 0;
      rpar = 0;
      radj = 0;
      newCds = Rcpp::clone(cds);
      distParityProp = Rcpp::clone(distParity);
      Rcpp::NumericVector accunits;
      int p;
      if(lambda > 0){
	p = R::rpois(lambda);
	p++;
	rlam = p;
      } else{
	p = 1;
	rlam = p;
      }
      // Begin loop over p to collect connected components
      for(int n = 0; n < p; n++){
      
	// Declare objects needed
	nccCheck = searchListBound.size();
	int testPart;
	breakOuter = 0;
	int reps = nccCheck;
	int breakInner = 0;

	// while there are still >0 partitions to select and no swap eliminated cd
	for(int c = 0; c < reps; c++){ 
	  
	  // Sample a test partition to propose switching
	  Rcpp::NumericVector numCCVec(nccCheck);
	  for(int i = 0; i < nccCheck; i++){
	    numCCVec(i) = i;
	  }
	  if(numCCVec.size() > 1){
	    Rcpp::NumericVector testPart1 = Rcpp::RcppArmadillo::sample(numCCVec, 1, TRUE);
	    testPart = testPart1(0);
	  } else{
	    testPart = numCCVec(0);
	  }
	  switchPartition = searchListBound(testPart);
	  searchListBound.erase(testPart);
	  nccCheck--;

	  ///////////////////////////////////////////////
	  // Compare against accunits - is it adjacent?
	  int adjcheck = 0;
	  for(int i = 0; i < switchPartition.size(); i++){

	    Rcpp::NumericVector avec = aList(i);
	    for(int j = 0; j < avec.size(); j++){
	      bool testadj = is_true(any(accunits == avec(j)));
	      
	      if(testadj == TRUE){
		radj++;
		adjcheck = 1;
		break;
	      }

	    }

	    if(adjcheck == 1){
	      break;
	    }
	    
	  }
	  if(adjcheck == 1){
	    continue;
	  }
	  ///////////////////////////////////////////////
	  
	  // Get candidate congressional districts for that partition
	  candCd.erase(candCd.begin(), candCd.end());
	  for(int i = 0; i < switchPartition.size(); i++){
	    
	    Rcpp::NumericVector adjUnitTemp = aList(switchPartition(i));
	    Rcpp::NumericVector candCdI;
	    for(int j = 0; j < adjUnitTemp.size(); j++){
	      
	      int thisCd = newCds(adjUnitTemp(j));
	      bool testCd = is_true(any(candCdI == thisCd));
	      if(testCd == FALSE){
		candCdI.push_back(newCds(adjUnitTemp(j)));
	      }
	      
	    }
	    for(int j = 0; j < candCdI.size(); j++){
	      candCd.push_back(candCdI(j));
	    }
	  }

	  thisPart = cds(switchPartition(0));
	  
	  //////////////////////////////////////////////////////////////////////// 
	  // Break if chosen partition won't eliminate a congressional district //
	  if(elimCheck == 1){
	    int nThisCd = 0;
	    for(int i = 0; i < units; i++){
	      if(newCds(i) == thisPart){
		nThisCd++;
	      }
	    }
	  
	    if(switchPartition.size() == nThisCd){
	      continue;
	    }
	  }
	  /////////////////////////////////////////////////////////////////////////
	  
	  // Define objects
	  Rcpp::NumericVector chooseCd = Rcpp::unique(candCd);
	  arma::vec rmThisCd = Rcpp::as<arma::vec> (chooseCd);
      
	  arma::uvec toRemove = find(rmThisCd == thisPart);
      
	  // Remove old congressional district
	  if(toRemove.n_elem > 0){
	    chooseCd.erase(toRemove(0));
	  }
      
	  int parityCond;
	  int searchListProp;

	  int ccLoop = chooseCd.size();

	  // Inner loop tests possible congressional districts
	  for(int d = 0; d < ccLoop; d++){ 
	  
	    int indProp;
	    Rcpp::NumericVector cdSize;
	    for(int i = 0; i < chooseCd.size(); i++){
	      cdSize.push_back(i);
	    }
	    if(cdSize.size() > 1){
	      Rcpp::NumericVector indProp1 = Rcpp::RcppArmadillo::sample(cdSize, 1, TRUE);
	      indProp = indProp1(0);
	      propSwitch = chooseCd(indProp);
	    } else{
	      indProp = 0;
	      propSwitch = chooseCd(indProp);
	    }

	    chooseCd.erase(indProp);
	  
	    // Here, we test to see if the partition will generate a new congressional district
	    cdProp = Rcpp::clone(newCds);
	    for(int i = 0; i < switchPartition.size(); i++){
	      cdProp(switchPartition(i)) = propSwitch;
	    }
	  
	    arma::vec cdPropTest = Rcpp::as<arma::vec> (cdProp);
	    arma::vec cdCurrTest = Rcpp::as<arma::vec> (newCds);

	    Rcpp::List alConnectedProp = genCC2(aList, cdProp);
	    Rcpp::List alConnectedOld = genCC2(aList, newCds);

	    ////////////////////////////////////
	    // Break if a district is shattered
	    int distCheck = genGraph(alConnectedProp);

	    if(distCheck != dists){
	      rsplit++;
	      continue;
	    }
	    ////////////////////////////////////
	    // Check to see if equal population constraint is violated
	    Rcpp::NumericVector parThisSwap = popCalc(switchPartition, pops, distParityProp, 
						      thisPart, propSwitch);
	    
	    parityCond = 0;
	    for(int i = 0; i < dists; i++){
	      if(margin > std::abs(parThisSwap(i))){
		parityCond++;
	      }
	    }
	    //////////////////////////////////////////////
	    // Break if parity condition not satisfied
	    if(parityCond != dists){
	      rpar++;
	      continue;
	    }
	    //////////////////////////////////////////////

	    actccs++;

	    // Recalculate metropolis-hastings probability
	    int c1 = 0;
	    int c2 = 0;
	    arma::vec cdOrig = Rcpp::as<arma::vec> (newCds); // CD's from last swap
	    arma::uvec vLPIndices = find(cdPropTest == propSwitch);
	    for(int i = 0; i < switchPartition.size(); i++){
	      
	      Rcpp::NumericVector adjVec = aList(switchPartition(i));
	      for(int j = 0; j < adjVec.size(); j++){
		
		// Calculate C(V_0, V_l' \ V_0)
		// Add 1 if you are adjacent to switched partition
		// AND your old cd assignment is the proposed switch
		if(cdOrig(adjVec(j)) == propSwitch){
		  c1++;
		}
		
		// Calculate C(V_0, V_l \ V_0)
		// Add 1 if you are adjacent to the proposed switch,
		// if your cd assignment is the cc's old cong district,
		// AND you are not in the switch partition
		if(cdOrig(adjVec(j)) == cdOrig(switchPartition(0)) && 
		   is_true(any(switchPartition == adjVec(j))) == FALSE){
		  c2++;
		}
		
	      }
	      
	    }

	    // Pair of congressional districts
	    Rcpp::NumericVector cdswitch(2);
	    cdswitch(0) = propSwitch;
	    cdswitch(1) = thisPart;

	    // Calculate sum of squared distances for compactness
	    // Numerator is ssd within proposed district for switched
	    // Denominator is ssd within old district for switched
	    double compRat;
	    distsumssd = 0;
	    newrpi = 0;
	    oldrpi = 0;

	    double ssdDenom = ssdDists(propSwitch) + ssdDists(thisPart);
	      
	    // Loop over congressional district pair
	    for(int i = 0; i < cdswitch.size(); i++){
		
	      // Initialize objects
	      double oldssd = 0;
	      double newssd = 0;
	      uvec newcds = find(cdPropTest == cdswitch(i));
	      uvec oldcds = find(cdCurrTest == cdswitch(i));

	      // SSD for new partition
	      for(int j = 0; j < newcds.size(); j++){
		for(int k = j + 1; k < newcds.size(); k++){
		  newssd += ssd(newcds(j), newcds(k)) * pops(newcds(j)) * 
		    pops(newcds(k));
		}
	      }

	      // SSD for old partition
	      for(int j = 0; j < oldcds.size(); j++){
		for(int k = j + 1; k < oldcds.size(); k++){
		  oldssd += ssd(oldcds(j), oldcds(k)) * pops(oldcds(j)) * 
		    pops(oldcds(k));
		}
	      }

	      // Standard compactness constraint 
	      if(origpartcompact == 0){
		newrpi += newssd / ssdDenom;
		oldrpi += oldssd / ssdDenom;
		distsumssd += (double)(newssd - oldssd) / ssdDenom;
	      }
	      // Compactness constraint relative to original partition
	      if(origpartcompact == 1){
		distsumssd += (double)std::abs((newssd / ssdDenom) - 1);
	      }
	    }

	    compRat = (double)exp((beta * log(2)) * distsumssd);
	     
	    ///////////////////////////////////////////////////////////////////////////
	    // Calculate dissimilarity index for group - for segregation restriction //
	    ///////////////////////////////////////////////////////////////////////////
	    double compRatDissim = 0;
	    diffDissim = 0;
	      
	    // Loop over congressional district pair
	    for(int i = 0; i < cdswitch.size(); i++){
		
	      // Initialize objects
	      int oldpopdist = 0;
	      int newpopdist = 0;
	      int oldpopmin = 0;
	      int newpopmin = 0;
	      uvec newcds = find(cdPropTest == cdswitch(i));
	      uvec oldcds = find(cdCurrTest == cdswitch(i));

	      // Fill objects for old cd assignments
	      for(int j = 0; j < oldcds.size(); j++){		  
		oldpopdist += pops(oldcds(j));
		oldpopmin += grouppop(oldcds(j));
	      }
	      // Fill objects for new cd assignments
	      for(int j = 0; j < newcds.size(); j++){
		newpopdist += pops(newcds(j));
		newpopmin += grouppop(newcds(j));
	      }

	      // Construct proportions
	      double oldminprop = (double)oldpopmin / oldpopdist;
	      double newminprop = (double)newpopmin / newpopdist;

	      // Get difference in dissimilarity
	      diffDissim += (double)(newpopdist * std::abs(newminprop - pAll) - 
				     oldpopdist * std::abs(oldminprop - pAll));
		
	    }

	    compRatDissim = (double)exp((betadiss * log(2)) * 
					denom * diffDissim);

	    ////////////////////////////////////////////////
	    // Calculate penalty for population deviation //
	    ////////////////////////////////////////////////
	    // double compRatPop = 0;
	    // diffPop = 0.0;
	      
	    // // Loop over congressional districts
	    // for(int i = 0; i < cdswitch.size(); i++){
		
	    //   // Population object
	    //   int distpop = 0;
	    //   uvec newcds = find(cdPropTest == cdswitch(i));

	    //   for(int j = 0; j < newcds.size(); j++){
	    // 	distpop += pops(newcds(j));
	    //   }

	    //   // Calculate penalty
	    //   diffPop += (double)std::abs((double)((double)distpop / parity) - 1);
		
	    // }
	      
	    // compRatPop = (double)exp((double)((double)betapop * log(2)) * diffPop);

	    double compRatPop = 0;
	    double diffPop = 0.0;
	    
	    // Initialize psi
	    psiNewPop = 0.0;
	    psiOldPop = 0.0;

	    // Loop over congressional districts
	    for(int i = 0; i < cdswitch.size(); i++){

	      // Population objects
	      int distPopNew = 0;
	      int distPopOld = 0;
	      uvec newcds = find(cdPropTest == cdswitch(i));
	      uvec oldcds = find(cdCurrTest == cdswitch(i));
	      
	      // Get pop of new dists
	      for(int j = 0; j < newcds.size(); j++){
		distPopNew += pops(newcds(j));
	      }
	      // Get pop of old dists
	      for(int j = 0; j < oldcds.size(); j++){
		distPopOld += pops(oldcds(j));
	      }

	      // Calculate penalty
	      psiNewPop += (double)std::abs((double)((double)distPopNew / parity) - 1);
	      psiOldPop += (double)std::abs((double)((double)distPopOld / parity) - 1);

	    }

	    // Get ratio
	    compRatPop = (double)exp((double)((double)betapop * log(2)) * psiNewPop) / 
	      exp((double)((double)betapop * log(2)) * psiOldPop);
	      
	    ////////////////////////////////////////////////////
	    // Calculate penalty for number of units switched //
	    ////////////////////////////////////////////////////
	    double compRatSwitch = 0;
	    diffSwitch = 0;
	      
	    for(int i = 0; i < cdswitch.size(); i++){

	      // Objects to count # similar precs in new cd 
	      int countdists = 0;
	      uvec newcds = find(cdPropTest == cdswitch(i));

	      // Get comparison to original plan, convert to numericvector
	      uvec findOrig = find(getCds == cdswitch(i));
	      Rcpp::NumericVector compareOrig(findOrig.size());
	      for(int j = 0; j < findOrig.size(); j++){
		compareOrig(j) = findOrig(j);
	      }

	      // Loop over objects in newcd, add if in original partition
	      for(int j = 0; j < newcds.size(); j++){
		if(is_true(any(compareOrig == newcds(j))) == TRUE){
		  countdists++;
		}
	      }

	      // Recalculate diffSwitch
	      diffSwitch += (double)countdists / origplan.size();

	    }

	    compRatSwitch = (double)exp((betaswitch * log(2)) * 
					std::abs(diffSwitch - 1));

	    // Set penalties to 1 if not being used
	    if(beta == 0 & annealbeta == 0){
	      compRat = 1;
	    }
	    if(betadiss == 0 & annealbetadiss == 0){
	      compRatDissim = 1;
	    }
	    if(betapop == 0 & annealbetapop == 0){
	      compRatPop = 1;
	    }
	    if(betaswitch == 0 & annealbetaswitch == 0){
	      compRatSwitch = 1;
	    }

	    // Get new mhprob
	    mhprob = (double)mhprob * pow(1 - eprob, c1) / pow(1 - eprob, c2) * 
	      compRat * compRatDissim * compRatPop * compRatSwitch;

	    // Change newCds
	    newCds = Rcpp::clone(cdProp);
	    distParityProp = Rcpp::clone(parThisSwap);

	    // Add switch partition to accunits
	    for(int i = 0; i < switchPartition.size(); i++){
	      accunits.push_back(switchPartition(i));
	    }
	  
	    // Break this loop if good cd
	    breakInner = 1;
	    if(actccs == p){
	      breakOuter = 1;
	    }

	    break; 
	  
	  }

	  if(breakOuter == 1 || breakInner == 1){
	    break;
	  }

	}

      }
           
    }while(breakOuter == 0);

    //////////////////////////////////////////
    // Fifth - Accept with some probability //
    //////////////////////////////////////////
    Rcpp::NumericVector getMin(2);
    getMin(0) = 1;
    getMin(1) = mhprob;
    double accProb = min(getMin);
    vec testkeep = Rcpp::runif(1);
    int decision = 0;
    if(testkeep(0) <= accProb){
      decision++;
    }

    // Propose to change beta if we are annealing
    if(annealbeta == 1){
      beta = changeBeta(betavec, beta, distsumssd, betaweights);
    }
    if(annealbetadiss == 1){
      betadiss = changeBeta(betadissvec, betadiss, diffDissim, betaweights);
    }
    if(annealbetapop == 1){
      // betapop = changeBeta(betapopvec, betapop, diffPop, betaweights);
      if(decision == 1){
	betapop = changeBeta(betapopvec, betapop, psiNewPop, betaweights); // To test
      }
      if(decision == 0){
	betapop = changeBeta(betapopvec, betapop, psiOldPop, betaweights);
      }
    }
    if(annealbetaswitch == 1){
      betaswitch = changeBeta(betaswitchvec, betaswitch, diffSwitch, betaweights);
    }

    /////////////////////////////////////
    // Sixth - clean up, store results //
    /////////////////////////////////////

    // Change congressional district for those switched
    if(decision == 1){
      cds = Rcpp::clone(newCds);
      distParity = Rcpp::clone(distParityProp);
    }

    // Store edges cut
    alCutedge = genCutC(aList, cds);

    // Store previous iteration, and advance the counter by 1
    for(int i = 0; i < units; i++){
      int index = k * units + i;
      cdDist[index] = cds(i);
    }
    
    // Store data about decisions and rejections
    decisionCount(k) = decision;
    numCC(k) = ncc;
    rejSplit(k) = rsplit;
    rejPar(k) = rpar;
    rejAdj(k) = radj;
    drawLam(k) = rlam;
    
    // Store annealed betas and criteria (latter for reweighting)
    betastore(k) = beta;
    diffBetastore(k) = distsumssd;
    betaDisstore(k) = betadiss;
    diffDisstore(k) = diffDissim;
    betaPopstore(k) = betapop;
    // diffPopstore(k) = diffPop;
    betaSwitchstore(k) = betaswitch;
    diffSwitchstore(k) = diffSwitch;
    diffPopstore(k) = psiNewPop; 

    // Increase counter by 1
    Rcpp::Rcout << k << std::endl;

    k++;
    
  }
    
  return Rcpp::List::create(cdDist, numAcc, numCC, rejSplit, rejPar, 
			    rejAdj, drawLam, betastore, betaDisstore, 
			    betaPopstore, betaSwitchstore,
			    diffBetastore, diffDisstore, diffPopstore,
			    diffSwitchstore, decisionCount);

}

