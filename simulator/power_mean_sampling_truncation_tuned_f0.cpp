#include <iostream>
#include <queue>
#include <boost/lexical_cast.hpp>
#include <boost/random.hpp>
#include "sampling.hpp"

typedef std::deque< std::pair<double, double> > results_t;
typedef boost::random::mt19937 rand_t;

class InputParameters {
public:
  InputParameters() {};
  size_t m_N;
  double m_k0;
  double m_expected_k;
  double m_dk;
  double m_alpha;
  double m_power;
  uint64_t m_seed;
  void Load( int argc, char** argv ) {
    if(argc != 8) {
      std::cerr << "Usage : ./parametrized_sampling.out N k0 expected_k dk alpha mean_power _seed" << std::endl;
      std::cerr << "  For of p_{ij} = ((f_i^p*f_j^p)/2)^(1/p)" << std::endl;
      std::cerr << "  f0 is controlled so that the sampled degree is in expected_k +- dk" << std::endl;
      exit(1);
    }

    m_N = boost::lexical_cast<size_t>(argv[1]);
    m_k0 = boost::lexical_cast<double>(argv[2]);
    m_expected_k = boost::lexical_cast<double>(argv[3]);
    m_dk = boost::lexical_cast<double>(argv[4]);
    m_alpha = boost::lexical_cast<double>(argv[5]);
    m_power = boost::lexical_cast<double>(argv[6]);
    m_seed = boost::lexical_cast<uint64_t>(argv[7]);
  };
};

double FindF0FromResults( const results_t & results, double target ) {
  double xy_sum = 0.0;
  double x_sum = 0.0;
  double y_sum = 0.0;
  double x2_sum = 0.0;
  for( auto result : results ) {
    double x = result.first;
    double y = result.second;
    xy_sum += x*y;
    x_sum += x;
    y_sum += y;
    x2_sum += x*x;
  }

  size_t n = results.size();

  double a = (n*xy_sum - x_sum*y_sum) / (n*x2_sum - x_sum*x_sum);
  double b = (x2_sum*y_sum - xy_sum*x_sum) / (n*x2_sum - x_sum*x_sum);

  std::cerr << "a: " << a << ", b: " << b << std::endl;

  return (target-b) / a;
}

Network* RunSampling( double f0, const InputParameters& input, rand_t& rnd ) {
  Sampling net(&rnd);
  net.GenerateER( input.m_N, input.m_k0, &rnd);

  Network* sampled = net.PowerMeanSampling(f0, input.m_alpha, input.m_power, true);
  return sampled;
}

bool DegreeIsInRange( const Network* net, const InputParameters & input ) {
  double k = net->AverageDegree();
  if( k < input.m_expected_k + input.m_dk && k > input.m_expected_k - input.m_dk ) {
    return true;
  }
  else {
    return false;
  }
}

int main(int argc, char** argv) {

  InputParameters input;
  input.Load(argc, argv);
  rand_t rnd(input.m_seed);

  results_t results;
  results.push_back( std::make_pair(1.0, input.m_k0) );  // heuristic
  results.push_back( std::make_pair(0.0, 0.0) );

  double f0 = FindF0FromResults( results, input.m_expected_k );
  Network* sampled = RunSampling( f0, input, rnd );
  std::cerr << "  f0: " << f0 << ", k: " << sampled->AverageDegree() << std::endl;

  while( DegreeIsInRange( sampled, input ) == false ) {
    results.push_back( std::make_pair( f0, sampled->AverageDegree() ) );
    delete sampled;
    f0 = FindF0FromResults( results, input.m_expected_k );
    sampled = RunSampling( f0, input, rnd );
    std::cerr << "  f0: " << f0 << ", k: " << sampled->AverageDegree() << std::endl;

    if( results.size() > 10 ) {
      results.pop_front();
    }
  }

  std::ofstream fout("sampled.edg");
  sampled->Print( fout );
  fout.flush();
  fout.close();

  delete sampled;

  return 0;
}

