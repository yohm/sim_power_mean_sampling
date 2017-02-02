#ifndef SAMPLING_NET_HPP
#define SAMPLING_NET_HPP

#include <iostream>
#include <set>
#include <string>
#include <boost/random.hpp>
#include "network.hpp"

class Sampling : public Network {
public:
  Sampling(boost::random::mt19937* rnd) : m_rnd(rnd) {}
  Network* PowerMeanSampling( double f0, double alpha, double power);
private:
  boost::random::mt19937* const m_rnd;
  double Rand01() {
    boost::random::uniform_01<> uniform;
    return uniform(*m_rnd);
  }
  double RandExp() {
    boost::random::exponential_distribution<double> exp(1.0);
    return exp(*m_rnd);
  }
  double RandWeibull( double exponent, double f0 ) {
    boost::random::weibull_distribution<double> weibull(exponent, f0);
    return weibull(*m_rnd);
  }
  void AssignPreference( std::vector<double>&pref, double f0, double alpha );
  double PowerMean( double fi, double fj, double sum_exp );

  Network* MakeNewNet( const std::set<size_t>& nodes, const std::set<Link>& links );
  std::map<size_t,size_t> CompactIndex( const std::set<size_t>& nodes );
};

#endif

