#ifdef INSIDE
#include <stdlib.h>
#include <iostream>
#include <RInside.h>
#endif
#include <Rcpp.h>
#include "nanoflannBits.hpp"


typedef KDTreeVectorOfVectorsAdaptor< point_vov_t, double >  kd_tree_t;
using namespace Rcpp;

class ndIndex{
public:
	ndIndex(SEXP points_){
		//convert to numeric matrix
	    NumericMatrix points(points_);

	   //get the number of rows and columns
		size_t nr = points.nrow(), nc = points.ncol();
		dim = nc;

		//convert to vector of vectors
		data = new point_vov_t(nr);
		for(size_t i=0; i<nr; i++){
			NumericMatrix::Row row = points(i,_);
			(*data)[i].assign(row.begin(), row.end());
		}

		//Construct a kd-tree index:
		mat_index = new kd_tree_t(
				nc /*dim*/,
				*data,
				10 /* max leaf */ );
		mat_index->index->buildIndex();
	}
	~ndIndex(){
		delete data;
		delete mat_index;
	}
	SEXP findNeighbors(SEXP queryPt, int numberOfResults, bool giveData=TRUE){
		NumericVector query_pt(queryPt);

		std::vector<size_t>   ret_indexes(numberOfResults);
    	std::vector<double> outDists(numberOfResults);

    	nanoflann::KNNResultSet<double> resultSet(numberOfResults);

    	resultSet.init(&ret_indexes[0], &outDists[0] );
    	mat_index->index->findNeighbors(resultSet, &query_pt[0], nanoflann::SearchParams(10));

//    	std::cout << "knnSearch(nn="<<numberOfResults<<"): \n";
//    	for (size_t i=0;i<numberOfResults;i++)
//    		std::cout << "ret_index["<<i<<"]=" << ret_indexes[i] << " out_dist_sqr=" << out_dists_sqr[i] << std::endl;
    	//convert the indexing base
    	transform(ret_indexes.begin(),ret_indexes.end(),ret_indexes.begin(),bind2nd( std::plus<double>(), 1.0));
    	//get the actual distance (rather than the square)
    	transform(outDists.begin(),outDists.end(), outDists.begin(), static_cast<double (*)(double)>(std::sqrt));

    	//output
    	if (giveData){
    		//generate the data for output
    		Rcpp::NumericMatrix outData(numberOfResults, dim);
    		for (int i = 0; i<numberOfResults; i++){
    			NumericMatrix::Row row = outData(i,_);
    			for(int j=0; j < dim; j++){
    				outData(i,j) = (*data)[ret_indexes[i]-1][j];
    			}

    		}

			Rcpp::List RResults = Rcpp::List::create(
					Rcpp::Named("Index") = ret_indexes,
					Rcpp::Named("Distance") = outDists,
					Rcpp::Named("Data") = outData
					);
			return RResults;
    	} else {
    		//just output
			Rcpp::List RResults = Rcpp::List::create(
					Rcpp::Named("Index") = ret_indexes,
					Rcpp::Named("Distance") = outDists
					);
    	return RResults;
    	}
	}
	SEXP findRadius(SEXP queryPt, double searchRadius, bool giveData=TRUE){
		NumericVector query_pt(queryPt);

		//need to square the distance
		searchRadius = searchRadius * searchRadius;

		//const num_t search_radius = static_cast<num_t>(0.1);
		std::vector<std::pair<size_t,double> >   ret_matches;

		nanoflann::SearchParams params;
		//params.sorted = false;

		const size_t nMatches = mat_index->index->radiusSearch(&query_pt[0],searchRadius, ret_matches, params);

		//std::cout << "radiusSearch(): radius=" << searchRadius << " -> " << nMatches << " matches\n";


	  	//output
		Rcpp::NumericVector RIndex(nMatches);
		Rcpp::NumericVector RDistance(nMatches);
		for (size_t i=0;i<nMatches;i++){
			RIndex[i] = ret_matches[i].first + 1;
			RDistance[i] = std::sqrt(ret_matches[i].second);
			//std::cout << "idx["<< i << "]=" << ret_matches[i].first << " dist["<< i << "]=" << ret_matches[i].second << std::endl;
		}
		//std::cout << "\n";

		if (giveData){
			//generate the data for output
			Rcpp::NumericMatrix outData(nMatches, dim);
			for (size_t i = 0; i<nMatches; i++){
				NumericMatrix::Row row = outData(i,_);
				for(size_t j=0; j < (*data)[i].size(); j++){
					outData(i,j) = (*data)[RIndex[i]-1][j];
				}
			}
			Rcpp::List RResults = Rcpp::List::create(
					Rcpp::Named("Index") = RIndex,
					Rcpp::Named("Distance") = RDistance,
					Rcpp::Named("Data") = outData
			);
			return RResults;
		} else {
			//just output
			Rcpp::List RResults = Rcpp::List::create(
					Rcpp::Named("Index") = RIndex,
					Rcpp::Named("Distance") = RDistance
			);
			return RResults;
		}

	};

private:
		int dim;
		point_vov_t* data;
		kd_tree_t* mat_index;
};

#ifdef INSIDE
//short test body
int main(int argc, char *argv[]){
        //create an embedded version of R
        RInside R(argc, argv);

        //generate a matrix of 3d points
        SEXP points_;
        //R.parseEval("matrix(runif(4E4),ncol=4)", points_);
        R.parseEval("matrix(c(0,0,3,4,5,8),byrow=TRUE,ncol=2)", points_);
        NumericMatrix points(points_);

    	//create the index
    	ndIndex newObject(points_);

    	// Query point:
    	SEXP query_pt_;
    	//R.parseEval("runif(4)", query_pt_);
    	R.parseEval("c(1,1)", query_pt_);
    	R["resultsSet"] = newObject.findNeighbors(query_pt_, 1, FALSE);
    	R["resultsSet2"] = newObject.findNeighbors(query_pt_, 1, TRUE);
    	R["resultsSet3"] = newObject.findRadius(query_pt_, 4.0, FALSE);
    	R["resultsSet4"] = newObject.findRadius(query_pt_, 4.0, TRUE);


    	R.parseEvalQ("print(resultsSet)");
    	R.parseEvalQ("print(resultsSet2)");
    	R.parseEvalQ("print(resultsSet3)");
    	R.parseEvalQ("print(resultsSet4)");

        exit(0);
}
#else
RCPP_MODULE(ndIndex){
	class_<ndIndex>("ndIndex")

	.constructor<SEXP>()

	.method("findNeighbors", &ndIndex::findNeighbors, "Find the nearest N neighbours.")
	.method("findRadius", &ndIndex::findRadius, "Find all rows within a certain radius.")
	;
}
#endif
