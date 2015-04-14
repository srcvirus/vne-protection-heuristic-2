#include "graph_util.h"
#include "vne_heuristic_solver.h"

bool VNEHeuristicSolver::EmbedVLink(const Edge& vlink,
                                    std::vector<bool>& v_physical,
                                    std::vector<bool>& v_virt) {
  std::vector<int> sources, destinations;
  int u = vlink.u(), v = vlink.v();
  int offset = virt_topology_->node_count();
  Edge shadow_vlink(vlink.u() + offset, vlink.v() + offset, vlink.bandwidth(),
                    vlink.cost());

  if (node_map_[u] == NIL && node_map_[v] == NIL) {
    sources = (*location_constraint_)[u];
    destinations = (*location_constraint_)[v];
  } else if (node_map_[u] == NIL && node_map_[v] != NIL) {
    sources = (*location_constraint_)[u];
    destinations.push_back(node_map_[v]);
  } else if (node_map_[u] != NIL && node_map_[v] == NIL) {
    sources.push_back(node_map_[u]);
    destinations = (*location_constraint_)[v];
  } else {
    sources.push_back(node_map_[u]);
    destinations.push_back(node_map_[v]);
  }
  DEBUG("Embedding primary path\n");
  auto primary_path = ConstrainedSP(this->physical_topology_, v_virt, sources,
                                    destinations, vlink.bandwidth());
  DEBUG("Primary path length = %d\n", primary_path->size());
  if (primary_path->size() <= 0) return false;

  for (auto& x : *primary_path) {
    v_physical[x] = true;
  }
  sources.clear();
  destinations.clear();

  if (node_map_[u + offset] == NIL && node_map_[v + offset] == NIL) {
    sources = (*location_constraint_)[u];
    destinations = (*location_constraint_)[v];
  } else if (node_map_[u + offset] == NIL && node_map_[v + offset] != NIL) {
    sources = (*location_constraint_)[u];
    destinations.push_back(node_map_[v + offset]);
  } else if (node_map_[u + offset] != NIL && node_map_[v + offset] == NIL) {
    sources.push_back(node_map_[u + offset]);
    destinations = (*location_constraint_)[v];
  } else {
    sources.push_back(node_map_[u + offset]);
    destinations.push_back(node_map_[v + offset]);
  }
  DEBUG("Embedding backup path\n");
  auto backup_path = ConstrainedSP(this->physical_topology_, v_physical,
                                   sources, destinations, vlink.bandwidth());
  if (backup_path->size() <= 0) return false;
  for (auto& x : *backup_path) {
    v_virt[x] = true;
  }

  node_map_[u] = primary_path->at(0);
  node_map_[v] = primary_path->at(primary_path->size() - 1);
  node_map_[u + offset] = backup_path->at(0);
  node_map_[v + offset] = backup_path->at(backup_path->size() - 1);
  if (edge_map_.find(vlink) == edge_map_.end()) {
    edge_map_.insert(std::make_pair(vlink, std::vector<Edge>()));
    edge_map_.insert(std::make_pair(shadow_vlink, std::vector<Edge>()));
  }
  for (int i = 1; i < primary_path->size(); ++i) {
    DEBUG("Virt link (%d, %d) mapped to physical link (%d, %d)\n", vlink.u(),
          vlink.v(), primary_path->at(i - 1), primary_path->at(i));
    edge_map_[vlink]
        .push_back(Edge(primary_path->at(i - 1), primary_path->at(i),
                        vlink.bandwidth(), vlink.cost()));
  }
  for (int i = 1; i < backup_path->size(); ++i) {
    DEBUG("Shadow Virt link (%d, %d) mapped to physical link (%d, %d)\n",
          vlink.u(), vlink.v(), backup_path->at(i - 1), backup_path->at(i));
    edge_map_[shadow_vlink]
        .push_back(Edge(backup_path->at(i - 1), backup_path->at(i),
                        shadow_vlink.bandwidth(), shadow_vlink.cost()));
  }
  return true;
}

bool VNEHeuristicSolver::Solve() {
  const int kPnodeCount = physical_topology_->node_count();
  const int kVnodeCount = virt_topology_->node_count();
  node_map_.resize(kVnodeCount * 2, NIL);
  auto edge_list = virt_topology_->ExportAsEdgeList();
  sort(edge_list->begin(), edge_list->end());
  std::vector<bool> v_physical(kPnodeCount, false);
  std::vector<bool> v_virt(kPnodeCount, false);
  for (auto& vlink : *edge_list) {
    DEBUG("Embedding vLink (%d, %d)\n", vlink.u(), vlink.v());
    if (!this->EmbedVLink(vlink, v_physical, v_virt)) return false;
  }
  return true;
}
