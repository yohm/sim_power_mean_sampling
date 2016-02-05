#include <iostream>
#include <fstream>
#include <string>
#include <boost/lexical_cast.hpp>
#include "network.hpp"

int main( int argc, char* argv[]) {
  if(argc != 2) {
    std::cerr << "Usage : ./main.out <edge_file>" << std::endl;
    exit(1);
  }

  Network network;
  std::ifstream fin(argv[1]);
  network.LoadFile( fin );
  // network.Print();
  std::ofstream fout("_output.json");
  fout << "{" << std::endl;
  fout << "  \"NumNodes\": " << network.NumNodes() << ',' << std::endl;
  fout << "  \"NumEdges\": " << network.NumEdges() << ',' << std::endl;
  fout << "  \"AverageDegree\": " << network.AverageDegree() << ',' << std::endl;
  fout << "  \"ClusteringCoefficient\": " << network.ClusteringCoefficient() << ',' << std::endl;
  fout << "  \"AverageEdgeWeight\": " << network.AverageEdgeWeight() << ',' << std::endl;
  fout << "  \"AverageOverlap\": " << network.AverageOverlap() << std::endl;
  fout << "}" << std::endl;

  std::ofstream dd("degree_distribution.dat");
  typedef std::pair<size_t, size_t> Freq;
  BOOST_FOREACH(Freq f, network.DegreeDistribution() ) {
    dd << f.first << ' ' << f.second << std::endl;
  }
  dd.flush();

  typedef std::pair<double, size_t> DoubleSize;
  // double edge_weight_bin_size = 1.0;
  std::ofstream ewd("edge_weight_distribution.dat");
  BOOST_FOREACH(DoubleSize f, network.EdgeWeightDistributionLogBin() ) {
    ewd << f.first << ' ' << f.second << std::endl;
  }
  ewd.flush();

  double strength_bin_size = 1.0;
  std::ofstream sd("strength_distribution.dat");
  BOOST_FOREACH(DoubleSize f, network.StrengthDistribution(strength_bin_size)) {
    sd << f.first << ' ' << f.second << std::endl;
  }
  sd.flush();

  typedef std::pair<size_t, double> SizeDouble;
  std::ofstream cc_d("cc_degree_correlation.dat");
  BOOST_FOREACH(SizeDouble f, network.CC_DegreeCorrelation() ) {
    cc_d << f.first << ' ' << f.second << std::endl;
  }
  cc_d.flush();

  std::ofstream sdc("strength_degree_correlation.dat");
  BOOST_FOREACH(SizeDouble f, network.StrengthDegreeCorrelation() ) {
    sdc << f.first << ' ' << f.second << std::endl;
  }
  sdc.flush();

  std::ofstream ndc("neighbor_degree_correlation.dat");
  BOOST_FOREACH(SizeDouble f, network.NeighborDegreeCorrelation() ) {
    ndc << f.first << ' ' << f.second << std::endl;
  }
  ndc.flush();

  std::ofstream owc("overlap_weight_correlation.dat");
  typedef std::pair<double, double> DD;
  BOOST_FOREACH(DD f, network.OverlapWeightCorrelationLogBin() ) {
    owc << f.first << ' ' << f.second << std::endl;
  }
  owc.flush();

  std::ofstream lrp("link_removal_percolation.dat");
  lrp << "#fraction  weak_link_removal_lcc susceptibility strong_link_removal_lcc susceptibility" << std::endl;
  network.AnalyzeLinkRemovalPercolationVariableAccuracy( 0.02, 0.02, lrp );
  lrp.flush();

  return 0;
}
