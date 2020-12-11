#include <cstdio> 
#include <cstdlib> 
#include <iostream>
#include <vector>
#include <list>

void sortandindices( std::vector<double>& numbers, std::vector<size_t>& indexes);
void dichotomic_grouping( std::vector<std::vector<double> >& p, int ptmax, std::list<std::vector<int> > &part, std::vector<int>& G);
void prctile(double *f, int N, double p, double *c);
