#include <iostream>
#include <random>
#include "sampling.hpp"

int main(int argc, char** argv) {
 if(argc != 7) {
    std::cerr << "Usage : ./power_mean_sampling.out N k f0 alpha beta _seed" << std::endl;
    std::cerr << "For p_{ij} = ((f_i^p*f_j^p)/2)^(1/p)" << std::endl;
    return 1;
  }

  size_t num_nodes = std::stoul(argv[1]); //  boost::lexical_cast<size_t>(argv[1]);
  double average_degree = std::stod(argv[2]); // boost::lexical_cast<double>(argv[2]);
  double f0 = std::stod(argv[3]); // boost::lexical_cast<double>(argv[3]);
  double alpha = std::stod(argv[4]); // boost::lexical_cast<double>(argv[4]);
  double beta = std::stod(argv[5]);// boost::lexical_cast<double>(argv[5]);
  uint64_t seed = std::stoull(argv[6]); // boost::lexical_cast<uint64_t>(argv[6]);

  std::mt19937 rnd(seed);
  Sampling net(&rnd);
  net.GenerateER( num_nodes, average_degree, &rnd );
  Network* sampled = net.PowerMeanSampling(f0, alpha, beta);
  std::ofstream fout("sampled.edg");
  sampled->Print( fout );
  fout.flush();
  fout.close();

  delete sampled;

  return 0;
}
