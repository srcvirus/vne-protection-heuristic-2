#ifndef VNE_HEURISTIC_SOLVER_
#define VNE_HEURISTIC_SOLVER_

#include "datastructure.h"
#include "util.h"

class VNEHeuristicSolver {
 public:
  VNEHeuristicSolver() {}
  VNEHeuristicSolver(Graph* physical_topology, Graph* virt_topology,
                     std::vector<std::vector<int>>* location_constraint)
      : physical_topology_(physical_topology),
        virt_topology_(virt_topology),
        location_constraint_(location_constraint) {}
  std::vector<int>& node_map() { return node_map_; }
  std::map<Edge, std::vector<Edge>>& edge_map() { return edge_map_; }
  bool Solve();

 private:
  bool EmbedVLink(const Edge& vlink, std::vector<bool>& v_physical,
                  std::vector<bool>& v_virt);
  Graph* physical_topology_;
  Graph* virt_topology_;
  std::vector<int> node_map_;
  std::map<Edge, std::vector<Edge>> edge_map_;
  std::vector<std::vector<int>>* location_constraint_;
};

#endif  // VNE_HEURISTIC_SOLVER_
