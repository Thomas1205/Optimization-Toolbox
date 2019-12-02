/*** written by Thomas Schoenemann as a private person without employment, October 2009 ***/

#include "submodular_energy_minimization.hh"
#include "vector.hh"

double add_term2(Graph<double,double,double>& graph, uint node1, uint node2,
                 double e00, double e01, double e10, double e11)
{

  //Check for NaN
  assert(e00 == e00);
  assert(e01 == e01);
  assert(e10 == e10);
  assert(e11 == e11);

  assert(e01+e10-e00-e11 >= -1e-6);

  assert(node1 != node2);

  if (e00 == 0.0 && e11 == 0.0) {
    graph.add_edge(node1,node2,e01,e10);
  }
  else {

    graph.add_tweights(node1, e10-e00, 0.0);
    graph.add_tweights(node2, e11-e10, 0.0);

    const double edge_weight = std::max(0.0,e01+e10-e00-e11);

    if (edge_weight > 0.0)
      graph.add_edge(node1,node2,edge_weight,0.0);
  }

  return e00;
}

double add_term3(Graph<double,double,double>& graph, uint node1, uint node2, uint node3,
                 double e000, double e001, double e010, double e011,
                 double e100, double e101, double e110, double e111)
{


  double offs = 0.0;

  double P = e000 + e011 + e101 + e110 - e001 - e010 - e100 - e111;

  if (P >= 0.0) {

    offs += e000;

    //std::cerr << "case P >= 0, P: " << P << std::endl;

    graph.add_tweights(node1,e101-e001,0.0);
    graph.add_tweights(node2,e110-e100,0.0);
    graph.add_tweights(node3,e011-e010,0.0);

    offs += add_term2(graph,node2,node3,0.0, e001+e010 - e000-e011,0.0,0.0);
    offs += add_term2(graph,node1,node3,0.0,0.0, e001+e100 - e000-e101,0.0);
    offs += add_term2(graph,node1,node2,0.0, e010+e100 - e000-e110,0.0,0.0);

    if (P > 0.0) {

      offs -= P;
      int id = graph.add_node();
      graph.add_tweights(id,0,P);
      graph.add_edge(node1,id,P,0.0);
      graph.add_edge(node2,id,P,0.0);
      graph.add_edge(node3,id,P,0.0);
    }
  }
  else {

    std::cerr << "case P < 0" << std::endl;
    std::cerr << "WARNING: untested" << std::endl;

    //TODO("case P < 0");

    offs += e111;

    graph.add_tweights(node1,0.0, e010-e110);
    graph.add_tweights(node2,0.0, e001-e011);
    graph.add_tweights(node3,0.0, e100-e101);

    offs += add_term2(graph,node2,node3, 0.0, 0.0, e101 + e110 - e100 - e111,0.0);
    offs += add_term2(graph,node1,node3, 0.0, e011 + e110 - e010 - e111, 0.0,0.0);
    offs += add_term2(graph,node1,node2, 0.0, 0.0, e011 + e110 - e010 - e111, 0.0);

    offs += P;
    int id = graph.add_node();
    graph.add_tweights(id,-P,0);
    graph.add_edge(node1,id,0.0,-P);
    graph.add_edge(node2,id,0.0,-P);
    graph.add_edge(node3,id,0.0,-P);
  }

  return offs;
}

void add_truncated_linear(Graph<double,double,double>& graph, const std::vector<uint>& participating_nodes,
                          uint cutoff_point, double cutoff_value, bool is_source_term)
{

  assert(cutoff_value >= 0.0);

  uint new_node = graph.add_node();
  if (is_source_term)
    graph.add_tweights(new_node,cutoff_value,0.0);
  else
    graph.add_tweights(new_node,0.0,cutoff_value);

  double single_term = cutoff_value / cutoff_point;

  double forward_cap = (is_source_term) ? single_term : 0.0;
  double backward_cap = (is_source_term) ? 0.0 : single_term;

  for (std::vector<uint>::const_iterator it = participating_nodes.begin();
       it != participating_nodes.end(); it++) {

    graph.add_edge(new_node,*it,forward_cap,backward_cap);
  }
}


void add_concave_increasing_function(Graph<double,double,double>& graph,
                                     const std::vector<uint>& participating_nodes,
                                     const std::map<uint,double>& function, bool is_source_term)
{

  assert(function.find(0) == function.end());

  uint nSegments = uint( function.size() );
  Math1D::NamedVector<double> slope(nSegments,MAKENAME(slope));
  Math1D::NamedVector<uint> breakpoint(nSegments,MAKENAME(breakpoint));

  uint last_breakpoint = 0;
  double last_value = 0.0;

  uint n=0;
  for (std::map<uint,double>::const_iterator it = function.begin(); it != function.end(); it++) {

    uint cur_breakpoint = it->first;
    double cur_value = it->second;

    breakpoint[n] = cur_breakpoint;

    if (cur_value < last_value) {
      INTERNAL_ERROR << " function is not increasing. Exiting..." << std::endl;
      exit(1);
    }

    slope[n] = (cur_value - last_value) / (cur_breakpoint - last_breakpoint);

    if (n > 0 && slope[n] > slope[n-1]) {
      INTERNAL_ERROR << " function is not concave. Exiting..." << std::endl;
    }

    n++;
    last_value = cur_value;
    last_breakpoint = cur_breakpoint;
  }

  if (breakpoint[nSegments-1] > participating_nodes.size()) {
    std::cerr << "WARNING: concave function contains parts that will never be used." << std::endl;
  }

  for (int n=nSegments-1; n > 0; n--) {
    for (int nn=n-1; nn >= 0; nn--)
      slope[nn] -= slope[n];
    //slope[n-1] -= slope[n];
  }

  for (uint n=0; n < nSegments; n++)
    add_truncated_linear(graph,participating_nodes,breakpoint[n],slope[n]*breakpoint[n],is_source_term);
}
