#include <iostream>
#include <queue>
#include <random>
#include "sampling.hpp"

typedef std::deque< std::pair<double, double> > results_t;
typedef std::mt19937 rand_t;

class InputParameters {
public:
  InputParameters() {};
  size_t m_N;
  double m_k0;
  double m_expected_k;
  double m_dk;
  double m_alpha;
  double m_beta;
  uint64_t m_seed;
  void Load( int argc, char** argv ) {
    if(argc != 8) {
      std::cerr << "Usage : ./parametrized_sampling.out N k0 expected_k dk alpha mean_power _seed" << std::endl;
      std::cerr << "  f0 is controlled so that the sampled degree is in expected_k +- dk" << std::endl;
      throw "invalid arguments";
    }

    m_N = std::stoul(argv[1]);
    m_k0 = std::stod(argv[2]);
    m_expected_k = std::stod(argv[3]);
    m_dk = std::stod(argv[4]);
    m_alpha = std::stod(argv[5]);
    m_beta = std::stod(argv[6]);
    m_seed = std::stoull(argv[7]);
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

  Network* sampled = net.PowerMeanSampling(f0, input.m_alpha, input.m_beta);
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

  std::ofstream f0out("f0.dat");
  f0out << f0 << std::endl;
  f0out.flush();
  f0out.close();
  std::ofstream fout("sampled.edg");
  sampled->Print( fout );
  fout.flush();
  fout.close();

  delete sampled;

  return 0;
}

