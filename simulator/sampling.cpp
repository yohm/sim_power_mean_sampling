#include "sampling.hpp"
#include <set>
#include <algorithm>

Network* Sampling::PowerMeanSampling(double f0, double alpha, double power, bool truncate_pf ) {
  std::vector<double> pref( m_nodes.size() );
  AssignPreference( pref, f0, alpha, truncate_pf );

  std::vector<size_t> sampledDegrees( m_nodes.size(), 0 );

  std::set<size_t> sampledNodes;
  std::set<Link> sampledLinks;
  for( const Link& l : m_links ) {
    size_t ni = l.m_node_id1;
    size_t nj = l.m_node_id2;
    double p = PowerMean( pref[ni], pref[nj], power );
    if( Rand01() < p ) {
      sampledLinks.insert( Link(ni,nj,l.m_weight) );
      sampledNodes.insert(ni);
      sampledNodes.insert(nj);
      sampledDegrees[ni] += 1;
      sampledDegrees[nj] += 1;
    }
  }

  std::ofstream fout("node_fitness.dat");
  for( size_t i=0; i < m_nodes.size(); i++ ) {
    fout << i << ' ' << sampledDegrees[i] << ' ' << pref[i] << std::endl;
  }
  fout.close();

  // calculate ki-kj correlation
  std::ofstream kikj_out("ki_kj.dat");
  std::ofstream fifj_out("fi_fj.dat");
  for( const Link& l : sampledLinks ) {
    size_t n1 = l.m_node_id1;
    size_t n2 = l.m_node_id2;
    kikj_out << sampledDegrees[n1] << ' ' << sampledDegrees[n2] << std::endl;
    fifj_out << pref[n1] << ' ' << pref[n2] << std::endl;
  }
  kikj_out.close();
  fifj_out.close();

  Network* net = MakeNewNet( sampledNodes, sampledLinks );
  return net;
}

void Sampling::AssignPreference( std::vector<double>& pref, double f0, double alpha, bool truncate_pf ) {
  for( size_t i=0; i < pref.size(); i++) {
    double x = RandWeibull(alpha, f0);
    while( truncate_pf && x > 1.0 ) {
      x = RandWeibull(alpha, f0);
    }
    pref[i] = x;
  }
}

double Sampling::PowerMean( double fi, double fj, double power ) {
  if( power == 0.0 ) { return std::sqrt( fi*fj ); }
  else {
    return std::pow( 0.5*(std::pow(fi, power) + std::pow(fj, power)), 1.0/power);
  }
}

std::map<size_t,size_t> Sampling::CompactIndex( const std::set<size_t>& nodes ) {
  std::map<size_t, size_t> map;
  size_t index = 0;
  for( size_t n : nodes ) {
    map[n] = index;
    index++;
  }
  return map;
}

Network* Sampling::MakeNewNet(const std::set<size_t>& nodes, const std::set<Link>& links ) {
  std::map<size_t,size_t> indexMap = CompactIndex(nodes);

  Network* net = new Network();

  for( const Link& l : links ) {
    size_t n1 = indexMap[ l.m_node_id1 ];
    size_t n2 = indexMap[ l.m_node_id2 ];
    net->AddLink( n1, n2, l.m_weight );
  }

  return net;
}

