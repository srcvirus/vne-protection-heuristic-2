#ifndef _HEURISTIC_H_
#define _HEURISTIC_H_

#include "datastructure.h"
#include "graph_util.h"

#include <algorithm>
#include <math.h>
#include <vector>
#include <utility>

#define PRIMARY 0
#define BACKUP 1

bool vector_size_comparator(const std::pair<int, std::vector<int>>& a,
                            const std::pair<int, std::vector<int>>& b) {
  return a.second.size() < b.second.size();
}

std::unique_ptr<std::vector<std::vector<int>>> CreateInitialNodeMap(
    const Graph* phys_topology,
    const std::vector<std::vector<int>>& location_constraints,
    const std::pair <int, std::vector<int>>& seed) {
  std::unique_ptr<std::vector<std::vector<int>>> node_maps(
      new std::vector<std::vector<int>>(2));
  (*node_maps)[PRIMARY] = std::vector<int>(location_constraints.size(), NIL);
  (*node_maps)[BACKUP] = std::vector<int>(location_constraints.size(), NIL);
  int mapped_node = seed.first;
  (*node_maps)[PRIMARY][mapped_node] = seed.second[PRIMARY];
  (*node_maps)[BACKUP][mapped_node] = seed.second[BACKUP];
  std::vector<std::pair<int, std::vector<int>>> temp_loc_constraints;
  for (int i = 0; i < location_constraints.size(); ++i) {
    temp_loc_constraints.push_back(std::make_pair(i, location_constraints[i]));
  }
  // Sort location constraints according to the number of constraints (low to
  // high)
  std::sort(temp_loc_constraints.begin(), temp_loc_constraints.end(),
            vector_size_comparator);

  std::vector<int> primary, backup;
  // Do the initial node mapping. First map the virtual nodes in the order of
  // the number of possible locations for that virtual node. For a virtual node
  // select the physical node that is not already mapped to some other virtual
  // node or its backup and has minimum overlap in the location_constraint list.
  std::vector<bool> taken(phys_topology->node_count(), false);
  taken[seed.second[PRIMARY]] = true;
  taken[seed.second[BACKUP]] = true;
  for (int i = 0; i < temp_loc_constraints.size(); ++i) {
    int vnode = temp_loc_constraints[i].first;
    if (vnode == mapped_node) continue;
    // printf("Mapping vnode %d\n", vnode);
    int best_candidate = NIL;
    for (int j = 0; j < temp_loc_constraints[i].second.size(); ++j) {
      int candidate = temp_loc_constraints[i].second[j];
      // printf("[Primary] Current candidate: %d\n", candidate);
      if (!taken[candidate]) {
        if (IsFeasibleBetterAssignment(phys_topology, primary, backup,
                                       candidate, best_candidate)) {
          // printf("[Primary] Best candidate updated to: %d\n", candidate);
          best_candidate = candidate;
        }
      }
    }
    if (best_candidate != NIL) {
      taken[best_candidate] = true;
      (*node_maps)[PRIMARY][vnode] = best_candidate;
      primary.push_back(best_candidate);
    }
    best_candidate = NIL;
    for (int j = 0; j < temp_loc_constraints[i].second.size(); ++j) {
      int candidate = temp_loc_constraints[i].second[j];
      // printf("[Backup] Current candidate: %d\n", candidate);
      if (!taken[candidate]) {
        if (IsFeasibleBetterAssignment(phys_topology, backup, primary,
                                       candidate, best_candidate)) {
          // printf("[Backup] Best candidate updated to: %d\n", candidate);
          best_candidate = candidate;
        }
      }
    }
    if (best_candidate != NIL) {
      taken[best_candidate] = true;
      (*node_maps)[BACKUP][vnode] = best_candidate;
      backup.push_back(best_candidate);
    }
  }
  return std::move(node_maps);
}

std::unique_ptr<std::vector<std::vector<int>>> PartitionGraph(
    const Graph* phys_topology, std::vector<std::vector<int>>* node_maps) {
  std::unique_ptr<std::vector<std::vector<int>>> partitions(
      new std::vector<std::vector<int>>(2));
  std::vector<bool> taken(phys_topology->node_count(), false);
  auto& primary = partitions->at(PRIMARY);
  auto& backup = partitions->at(BACKUP);
  std::vector<int> primary_seed, backup_seed;
  for (auto& node : node_maps->at(PRIMARY)) {
    primary.push_back(node);
    primary_seed.push_back(node);
    taken[node] = true;
  }
  for (auto& node : node_maps->at(BACKUP)) {
    backup.push_back(node);
    backup_seed.push_back(node);
    taken[node] = true;
  }
  for (int i = 0; i < phys_topology->node_count(); ++i) {
    if (!taken[i]) {
      if (!IsFeasiblePartition(phys_topology, primary, backup_seed, i)) {
        backup.push_back(i);
      } else if (!IsFeasiblePartition(phys_topology, backup, primary_seed, i)) {
        primary.push_back(i);
      } else {
        int component_decrease_primary =
            NumConnectedComponentsDecrease(phys_topology, primary, i);
        int component_decrease_backup =
            NumConnectedComponentsDecrease(phys_topology, backup, i);
        if (component_decrease_primary > component_decrease_backup) {
          primary.push_back(i);
        } else if (component_decrease_primary < component_decrease_backup) {
          backup.push_back(i);
        } else {
          int cut_edge_primary = NumCutEdges(phys_topology, primary, i);
          int cut_edge_backup = NumCutEdges(phys_topology, backup, i);
          if (cut_edge_primary > cut_edge_backup) {
            primary.push_back(i);
          } else if (cut_edge_primary < cut_edge_backup) {
            backup.push_back(i);
          } else {         
            if (primary.size() > backup.size()) {
              backup.push_back(i);
            } else if (primary.size() < backup.size()) {
                primary.push_back(i);
            } else {
              primary.push_back(i);
              double new_primary_mean_sp = MeanSubsetShortestPathLength(phys_topology, primary, primary_seed);
              primary.pop_back();
              double cur_primary_mean_sp = MeanSubsetShortestPathLength(phys_topology, primary, primary_seed);
              double primary_sp_reduction = new_primary_mean_sp - cur_primary_mean_sp;

              backup.push_back(i);
              double new_backup_mean_sp = MeanSubsetShortestPathLength(phys_topology, backup, backup_seed);
              backup.pop_back();
              double cur_backup_mean_sp = MeanSubsetShortestPathLength(phys_topology, backup, backup_seed);
              double backup_sp_reduction = new_backup_mean_sp - cur_backup_mean_sp;

              if (primary_sp_reduction > backup_sp_reduction) {
                primary.push_back(i);
              } else {
                backup.push_back(i);
              }
            }
          }
        }
      }
    }
  }
  return std::move(partitions);
}

std::vector<std::pair<int, int>> EmbedVLink(
    const Graph* phys_topology, const std::vector<int>& partition, int map_u,
    int map_v, long bw) {
  auto mapped_path = BFS(phys_topology, partition, map_u, map_v, bw);
  std::vector<std::pair<int, int>> edge_map;
  for (int i = 0; i < ((int)mapped_path->size()) - 1; ++i) {
    edge_map.push_back(
        std::pair<int, int>(mapped_path->at(i), mapped_path->at(i + 1)));
  }
  return std::move(edge_map);
}

std::unique_ptr<std::map<std::pair<int, int>, std::vector<std::pair<int, int>>>>
EmbedVN(Graph* phys_topology, const Graph* virt_topology,
        const std::vector<int>& partition, const std::vector<int>& node_map) {
  std::unique_ptr<
      std::map<std::pair<int, int>, std::vector<std::pair<int, int>>>> edge_map(
      new std::map<std::pair<int, int>, std::vector<std::pair<int, int>>>());
  int cost = 0;
  for (int u = 0; u < virt_topology->node_count(); ++u) {
    auto& u_neighbors = virt_topology->adj_list()->at(u);
    for (auto& vend_point : u_neighbors) {
      int v = vend_point.node_id;
      if (u <= v) continue;
      long bw = vend_point.bandwidth;
      auto emap =
          EmbedVLink(phys_topology, partition, node_map[u], node_map[v], bw);
      edge_map->insert(std::make_pair(std::pair<int, int>(u, v), std::move(emap)));
      for (auto& e : emap) {
        phys_topology->reduce_edge_residual_bandwidth(e.first, e.second, bw);
      }
    }
  }
  return std::move(edge_map);
}

long EmbeddingCost(const Graph* phys_topology, const Graph* virt_topology,
                  const VNEmbedding* embedding) {
  long cost = 0;
  for (auto emap_it = embedding->primary_edge_map->begin(); emap_it != embedding->primary_edge_map->end();
        ++emap_it) {
    auto& vlink = emap_it->first;
    auto& plinks = emap_it->second;
    for (auto& e : plinks) {
      cost += phys_topology->get_edge_cost(e.first, e.second) *
                virt_topology->get_edge_bandwidth(vlink.first, vlink.second);
    }
  }
  for (auto semap_it = embedding->backup_edge_map->begin(); semap_it != embedding->backup_edge_map->end();
          ++semap_it) {
    auto& vlink = semap_it->first;
    auto& plinks = semap_it->second;
    for (auto& e : plinks) {
      cost += phys_topology->get_edge_cost(e.first, e.second) *
                virt_topology->get_edge_bandwidth(vlink.first, vlink.second);
    }
  }
  return cost;
}
#endif  // _HEURISTIC_H_
