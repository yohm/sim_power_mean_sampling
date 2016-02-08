#include <iostream>
#include <boost/lexical_cast.hpp>
#include <boost/random.hpp>
#include "sampling.hpp"

int main(int argc, char** argv) {
 if(argc != 7) {
    std::cerr << "Usage : ./parametrized_sampling.out N k f0 alpha sum_exp _seed" << std::endl;
    std::cerr << "For of p_{ij} = ((f_i^p*f_j^p)/2)^(1/p)" << std::endl;
    exit(1);
  }

  size_t num_nodes = boost::lexical_cast<size_t>(argv[1]);
  double average_degree = boost::lexical_cast<double>(argv[2]);
  double f0 = boost::lexical_cast<double>(argv[3]);
  double alpha = boost::lexical_cast<double>(argv[4]);
  double power = boost::lexical_cast<double>(argv[5]);
  uint64_t seed = boost::lexical_cast<uint64_t>(argv[6]);

  boost::random::mt19937 rnd(seed);
  Sampling net(&rnd);
  net.GenerateER( num_nodes, average_degree, &rnd );
  Network* sampled = net.PowerMeanSampling(f0, alpha, power, true);
  std::ofstream fout("sampled.edg");
  sampled->Print( fout );
  fout.flush();
  fout.close();

  delete sampled;

  return 0;
}
