#include <iostream>
#include <boost/lexical_cast.hpp>
#include <boost/random.hpp>
#include "sampling.hpp"

int main(int argc, char** argv) {
 if(argc != 7) {
    std::cerr << "Usage : ./parametrized_sampling.out N k f0 alpha sum_exp _seed" << std::endl;
    std::cerr << "For of p_{ij} = 2*(f_i*f_j)/(f_i+f_j)**{sum_exp}" << std::endl;
    exit(1);
  }

  size_t num_nodes = boost::lexical_cast<size_t>(argv[1]);
  double average_degree = boost::lexical_cast<double>(argv[2]);
  double f0 = boost::lexical_cast<double>(argv[3]);
  double alpha = boost::lexical_cast<double>(argv[4]);
  double sum_exp = boost::lexical_cast<double>(argv[5]);
  uint64_t seed = boost::lexical_cast<uint64_t>(argv[6]);

  boost::random::mt19937 rnd(seed);
  Sampling net(&rnd);
  net.GenerateER( num_nodes, average_degree, &rnd );
  Network* sampled = net.ParametrizedSampling(f0, alpha, sum_exp);
  std::ofstream fout("sampled.edg");
  sampled->Print( fout );
  fout.flush();
  fout.close();

  delete sampled;

  return 0;
}
