#include <cassert>
#include <queue>
#include <set>
#include <algorithm>
#include <boost/lexical_cast.hpp>
#include "network.hpp"

void Network::LoadFile( std::ifstream& fin ) {
  while( !fin.eof() ) {
    size_t i, j;
    double weight;
    fin >> i;
    if( fin.eof() ) break;
    fin >> j >> weight;
    AddEdge(i, j, weight);
  }

  for( size_t i=0; i < m_nodes.size(); i++) {
    m_nodes[i].m_id = i;
  }
}

void Network::Print( std::ostream& out ) const {
  BOOST_FOREACH(const Link& link, m_links) {
    out << link.m_node_id1 << ' ' << link.m_node_id2 << ' ' << link.m_weight << std::endl;
  }
}

size_t Network::NumEdges() const {
  return m_links.size();
}

double Network::ClusteringCoefficient() const {
  double total = 0.0;
  for( size_t i = 0; i < m_nodes.size(); ++i) {
    total += LocalCC(i);
  }
  return total / m_nodes.size();
}

double Network::LocalCC(size_t idx) const {
  double total = 0.0;
  const Node& node = m_nodes[idx];
  size_t k = node.Degree();
  if( k <= 1 ) { return 0.0; }
  for( size_t i = 0; i < k; ++i) {
    for( size_t j = i+1; j < k; ++j) {
      size_t node1 = node.m_edges[i].m_node_id;
      size_t node2 = node.m_edges[j].m_node_id;
      if( m_nodes[node1].ConnectedTo(node2) ) { total += 1.0; }
    }
  }
  size_t num_pairs = k * (k-1) / 2.0;
  return total / num_pairs;
}

std::map<size_t, size_t> Network::DegreeDistribution() const {
  std::map<size_t, size_t> histo;
  BOOST_FOREACH(const Node& node, m_nodes) {
    if( histo.find(node.Degree()) != histo.end() ) {
      histo[node.Degree()] += 1;
    } else {
      histo[node.Degree()] = 1;
    }
  }
  return histo;
}

double Network::AverageEdgeWeight() const {
  double total = 0.0;
  BOOST_FOREACH(const Link& link, m_links) {
    total += link.m_weight;
  }
  return total / m_links.size();
}

double Network::AverageOverlap() const {
  double total = 0.0;
  for(size_t i = 0; i < m_links.size(); ++i) {
    total += LocalOverlap(i);
  }
  return total / m_links.size();
}

double Network::LocalOverlap(size_t link_id) const {
  size_t i = m_links[link_id].m_node_id1;
  size_t j = m_links[link_id].m_node_id2;

  std::set<size_t> neighbors_i;
  BOOST_FOREACH(const Edge& edge, m_nodes[i].m_edges) {
    neighbors_i.insert( edge.m_node_id );
  }
  size_t num_common = 0;
  BOOST_FOREACH(const Edge& edge, m_nodes[j].m_edges) {
    size_t k = edge.m_node_id;
    if( neighbors_i.find(k) != neighbors_i.end() ) { num_common++; }
  }

  size_t ki = m_nodes[i].Degree();
  size_t kj = m_nodes[j].Degree();
  if( ki == 1 and kj == 1 ) { return 0.0; }
  return static_cast<double>(num_common)/(ki + kj - 2 - num_common);
}

std::map<double, size_t> Network::EdgeWeightDistribution(double bin_size) const {
  std::map<int, size_t> histo;
  BOOST_FOREACH(const Link& link, m_links) {
    int key = static_cast<int>(link.m_weight/bin_size + 0.5);
    if(histo.find(key) == histo.end() ) { histo[key] = 0; }
    histo[key] += 1;
  }

  // convert to double
  std::map<double, size_t> histo_d;
  typedef std::pair<int, size_t> C;
  BOOST_FOREACH(C key_val, histo) {
    double key_d = key_val.first * bin_size;
    histo_d[key_d] = key_val.second;
  }
  return histo_d;
}

std::map<double, size_t> Network::EdgeWeightDistributionLogBin() const {
  std::map<int, size_t> histo;
  auto val_to_binidx = [](double weight)->int {
    if( weight <= 1.0 ) {
      return static_cast<int>( floor( log(weight)/log(2.0) ) );
    } else {
      return static_cast<int>( weight );
    }
  };
  auto binidx_to_val = [](int idx)-> double {
    if( idx <= 0 ) {
      return pow(2.0,idx);
    } else {
      return idx;
    }
  };
  auto binidx_to_binsize = [](int idx)->double {
    if( idx <= 0 ) {
      return pow(2.0,idx);
    } else {
      return 1;
    }
  };
  BOOST_FOREACH(const Link& link, m_links) {
    int key = val_to_binidx(link.m_weight);
    if(histo.find(key) == histo.end() ) { histo[key] = 0; }
    histo[key] += 1;
  }

  // convert to double
  std::map<double, size_t> histo_d;
  typedef std::pair<int, size_t> C;
  BOOST_FOREACH(C key_val, histo) {
    double key_d = binidx_to_val(key_val.first);
    double bin_size = binidx_to_binsize(key_val.first);
    histo_d[key_d] = key_val.second / bin_size; 
  }
  return histo_d;
}

std::map<double, size_t> Network::StrengthDistribution(double bin_size) const {
  std::map<int, size_t> histo;
  BOOST_FOREACH(const Node& node, m_nodes) {
    int key = static_cast<int>( node.Strength()/bin_size + 0.5 );
    if( histo.find(key) == histo.end() ) { histo[key] = 0; }
    histo[key] += 1;
  }

  // convert to double
  std::map<double, size_t> histo_d;
  typedef std::pair<int, size_t> C;
  BOOST_FOREACH(C key_val, histo) {
    double key_d = key_val.first * bin_size;
    histo_d[key_d] = key_val.second;
  }
  return histo_d;
}

std::map<size_t, double> Network::NeighborDegreeCorrelation() const {
  std::map<size_t, size_t> neighbor_degree_counts;
  std::map<size_t, double> neighbor_degree_total;
  for(size_t i = 0; i < m_nodes.size(); ++i) {
    const Node& node = m_nodes[i];
    const size_t k = node.Degree();
    if( k == 0 ) { continue; }
    const double k_nn = AverageNeighborDegree(i);
    if( neighbor_degree_counts.find(k) != neighbor_degree_counts.end() ) {
      neighbor_degree_counts[ k ] += 1;
      neighbor_degree_total[ k ] += k_nn;
    } else {
      neighbor_degree_counts[ k ] = 1;
      neighbor_degree_total[ k ] = k_nn;
    }
  }

  typedef std::pair<size_t, size_t> C;
  std::map<size_t, double> neighbor_degree_average;
  BOOST_FOREACH( C c, neighbor_degree_counts) {
    neighbor_degree_average[c.first] = neighbor_degree_total[ c.first ] / c.second;
  }
  return neighbor_degree_average;
}

double Network::AverageNeighborDegree(size_t i) const {
  const Node& node = m_nodes[i];
  double total = 0.0;
  if( node.m_edges.empty() ) { return 0.0; }
  BOOST_FOREACH(const Edge& e, node.m_edges) {
    size_t neighbor_idx = e.m_node_id;
    total += static_cast<double>( m_nodes[neighbor_idx].Degree() );
  }
  return total / node.m_edges.size();
}

std::map<size_t, double> Network::CC_DegreeCorrelation() const {
  std::map<size_t, size_t> cc_counts;
  std::map<size_t, double> cc_total;
  for(size_t i = 0; i < m_nodes.size(); ++i) {
    const Node& node = m_nodes[i];
    const size_t k = node.Degree();
    if( k == 0 ) { continue; }
    const double lcc = LocalCC(i);
    if( cc_counts.find(k) != cc_counts.end() ) {
      cc_counts[ k ] += 1;
      cc_total[ k ] += lcc;
    } else {
      cc_counts[ k ] = 1;
      cc_total[ k ] = lcc;
    }
  }

  typedef std::pair<size_t, size_t> C;
  std::map<size_t, double> cc_average;
  BOOST_FOREACH( C c, cc_counts) {
    size_t k = c.first;
    cc_average[k] = cc_total[ k ] / c.second;
  }
  return cc_average;
}

std::map<size_t, double> Network::StrengthDegreeCorrelation() const {
  std::map<size_t, size_t> s_counts;
  std::map<size_t, double> s_total;
  BOOST_FOREACH(const Node& node, m_nodes) {
    const size_t k = node.Degree();
    const double s = node.Strength();
    if( s_counts.find(k) == s_counts.end() ) { s_counts[k] = 0; s_total[k] = 0.0;}
    s_counts[k] += 1;
    s_total[k] += s;
  }

  typedef std::pair<size_t, size_t> C;
  std::map<size_t, double> s_average;
  BOOST_FOREACH( C k_counts_pair, s_counts) {
    size_t k = k_counts_pair.first;
    size_t count = k_counts_pair.second;
    s_average[k] = s_total[ k ] / count;
  }
  return s_average;
}

std::map<double, double> Network::OverlapWeightCorrelation(double bin_size) const {
  std::map<int, size_t> counts;
  std::map<int, double> totals;
  for(size_t i = 0; i < m_links.size(); ++i) {
    const Link& link = m_links[i];
    int bin_idx = static_cast<int>(link.m_weight / bin_size + 0.5);
    if( counts.find(bin_idx) == counts.end() ) { counts[bin_idx] = 0; totals[bin_idx] = 0.0; }
    counts[bin_idx] += 1;
    totals[bin_idx] += LocalOverlap(i);
  }

  std::map<double, double> ret;
  typedef std::pair<int, size_t> C;
  BOOST_FOREACH(C bin_idx_count, counts) {
    int bin_idx = bin_idx_count.first;
    ret[ bin_idx * bin_size ] = totals[bin_idx] / counts[bin_idx];
  }
  return ret;
}

std::map<double,double> Network::OverlapWeightCorrelationLogBin() const {
  std::map<int, int> counts;
  std::map<int, double> totals;
  auto val_to_binidx = [](double weight)->int {
    if( weight <= 1.0 ) {
      return static_cast<int>( floor( log(weight)/log(2.0) ) );
    } else {
      return static_cast<int>( weight );
    }
  };
  auto binidx_to_val = [](int idx)-> double {
    if( idx <= 0 ) {
      return pow(2.0,idx);
    } else {
      return idx;
    }
  };
  for(size_t i = 0; i < m_links.size(); ++i) {
    const Link& link = m_links[i];
    int bin_idx = val_to_binidx(link.m_weight);
    if( counts.find(bin_idx) == counts.end() ) { counts[bin_idx] = 0; totals[bin_idx] = 0.0; }
    counts[bin_idx] += 1;
    totals[bin_idx] += LocalOverlap(i);
  }

  std::map<double, double> ret;
  typedef std::pair<int, size_t> C;
  BOOST_FOREACH(C bin_idx_count, counts) {
    int bin_idx = bin_idx_count.first;
    double v = binidx_to_val(bin_idx);
    ret[ v ] = totals[bin_idx] / counts[bin_idx];
  }
  return ret;
}

void Network::SortLinkByWeight() {
  if( !m_is_sorted_by_weight ) {
    std::sort(m_links.begin(), m_links.end(), compareLinkByWeight);
    m_is_sorted_by_weight = true;
  }
}

void Network::AnalyzeLinkRemovalPercolationVariableAccuracy(double df, double d_R, std::ostream& os) {
  typedef std::tuple<double,double,double,double> perc_result_t;  // R_lcc_asc, susc_asc, R_lcc_desc, susc_desc

  auto PercolationAt = [&](double f)->perc_result_t {
    double lcc1, succ1, lcc2, succ2;
    AnalyzeLinkRemovalPercolation(f, true, lcc1, succ1);
    AnalyzeLinkRemovalPercolation(f, false, lcc2, succ2);
    perc_result_t result = std::make_tuple(lcc1, succ1, lcc2, succ2);
    return result;
  };

  std::map<double, perc_result_t> results;
  std::function<void(double,double,double)> investigate_between = [&](double f1, double f2, double d_R)->void {
    if( results.find(f1) == results.end() ) {
      results[f1] = PercolationAt(f1);
    }
    if( results.find(f2) == results.end() ) {
      results[f2] = PercolationAt(f2);
    }
    double dr1 = std::get<0>(results[f1]) - std::get<0>(results[f2]);
    double dr2 = std::get<2>(results[f1]) - std::get<2>(results[f2]);
    double df = std::abs(f2-f1);
    if( df > 0.0001 && (std::abs(dr1) > d_R || std::abs(dr2) > d_R) ) {
      std::cerr << f1 << ' ' << f2 << std::endl;
      std::cerr << dr1 << ' ' << dr2 << ' ' << d_R << std::endl;
      double f12 = (f1+f2)/2.0;
      investigate_between(f1, f12, d_R);
      investigate_between(f12, f2, d_R);
    }
  };

  for( double f = 0.0, f_next = df; f < 1.00; f = f_next, f_next += df ) {
    if( f_next > 1.0 ) { f_next = 1.0; }
    investigate_between(f, f_next, d_R);
  }

  for( const auto& keyval : results ) {
    perc_result_t r = keyval.second;
    double f = keyval.first;
    os << keyval.first << ' '
       << std::get<0>(r) << ' '
       << std::get<1>(r) << ' '
       << std::get<2>(r) << ' '
       << std::get<3>(r) << ' '
       << (1.0-f) * AverageDegree() << std::endl;
  }
}

void Network::AnalyzeLinkRemovalPercolation(double f, bool weak_link_removal, double& r_lcc, double& susceptibility) {
  Network* filtered = MakeFilteredNetwork(f, weak_link_removal);
  filtered->AnalyzePercolation(r_lcc, susceptibility);
  delete filtered;
}

Network* Network::MakeFilteredNetwork(double f, bool weak_link_removal) {
  Network * filtered = new Network(m_nodes.size());
  SortLinkByWeight();
  if( weak_link_removal ) {
    size_t threshold = m_links.size() * f;
    for( size_t i = threshold; i < m_links.size(); i++) {
      Link l = m_links[i];
      filtered->AddEdge( l.m_node_id1, l.m_node_id2, l.m_weight);
    }
  } else {
    size_t threshold = m_links.size() * (1.0 - f);
    for( size_t i = 0; i < threshold; i++) {
      Link l = m_links[i];
      filtered->AddEdge( l.m_node_id1, l.m_node_id2, l.m_weight);
    }
  }
  return filtered;
}

std::vector<size_t> Network::PercolatedClusterSizeDistribution() const {
  std::vector<int> cluster_ids(m_nodes.size(), -1);
  for( size_t i = 0; i < m_nodes.size(); ++i) { // coloring clusters
    if( cluster_ids[i] != -1 ) { continue; }
    SearchConnected(i, i, cluster_ids);
  }

  // get cluster size distribution
  std::map<int, size_t> cluster_sizes;
  BOOST_FOREACH(int id, cluster_ids) {
    if( cluster_sizes.find(id) == cluster_sizes.end() ) { cluster_sizes[id] = 0; }
    cluster_sizes[id] += 1;
  }

  std::vector<size_t> sizes;
  typedef std::pair<int, size_t> IntSize;
  BOOST_FOREACH(IntSize is, cluster_sizes) { sizes.push_back( is.second ); }

  return sizes;
}

void Network::AnalyzePercolation(double& r_lcc, double& susceptibility) const {
  // calculate LCC and susceptibility
  std::vector<size_t> sizes = PercolatedClusterSizeDistribution();
  std::sort(sizes.begin(), sizes.end());
  r_lcc = (1.0 * sizes.back()) / m_nodes.size(); // *(sizes.end()--);
  sizes.pop_back(); // remove the largest cluster
  double total = 0.0;
  BOOST_FOREACH(size_t cluster_size, sizes) {
    total += cluster_size * cluster_size;
  }
  susceptibility = total / m_nodes.size();
}

void Network::SearchConnected(size_t target_node, size_t parent_id, std::vector<int>& cluster_ids) const {
  assert( cluster_ids[target_node] == -1 );
  assert( target_node >= 0 && target_node < m_nodes.size() );

  std::queue<size_t> to_be_searched;
  to_be_searched.push(target_node);
  cluster_ids.at(target_node) = parent_id;
  while( ! to_be_searched.empty() ) {
    size_t searching = to_be_searched.front();
    to_be_searched.pop();
    BOOST_FOREACH(const Edge& edge, m_nodes.at(searching).m_edges) {
      size_t child = edge.m_node_id;
      if( cluster_ids.at(child) == -1 ) {
        cluster_ids.at(child) = parent_id;
        to_be_searched.push(child);
      }
    }
  }
}
