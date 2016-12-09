#include <Rcpp.h>
using namespace Rcpp;

// © Dirk U. Wulff, December 2016


// [[Rcpp::export]]
NumericMatrix tab(std::vector<double> v1, std::vector<double> v2){
  std::map<std::pair<double, double>, int> table;
  std::vector<double>::const_iterator it1 = v1.begin();
  std::vector<double>::const_iterator it2 = v2.begin();
  for(; it1 != v1.end() && it2 != v2.end(); ++it1, ++it2){
    std::pair <double,double> value = std::make_pair(*it1,*it2);
    table[value]++;
    }
  NumericMatrix res_tab(table.size(),3);
  int i = 0;
  std::map<std::pair<double, double>, int>::const_iterator it;
  for(it = table.begin(); it != table.end(); ++it, ++i){
    std::pair <double,double> value = it->first;
    res_tab(i,0) = value.first;
    res_tab(i,1) = value.second;
    res_tab(i,2) = it->second;
    };
  return res_tab;
  }



// [[Rcpp::export]]
NumericMatrix tab_sum(std::vector<double> v1, std::vector<double> v2, std::vector<double> a){
  std::map<std::pair<double, double>, double> table;
  std::vector<double>::const_iterator it1 = v1.begin();
  std::vector<double>::const_iterator it2 = v2.begin();
  std::vector<double>::const_iterator it3 = a.begin();
  for(; it1 != v1.end() && it2 != v2.end(); ++it1, ++it2, ++it3){
    std::pair <double,double> value = std::make_pair(*it1,*it2);
    table[value] += *it3;
    }
  NumericMatrix res_tab(table.size(),3);
  std::map<std::pair<double, double>, double>::const_iterator it;
  int i = 0;
  for(it = table.begin(); it != table.end(); ++it, ++i){
    std::pair <double,double> value = it->first;
    res_tab(i,0) = value.first;
    res_tab(i,1) = value.second;
    res_tab(i,2) = it->second;
    };
  return res_tab;
  }






