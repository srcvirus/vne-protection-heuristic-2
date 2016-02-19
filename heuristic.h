#ifndef _HEURISTIC_H_
#define _HEURISTIC_H_

#include "datastructure.h"
#include "graph_util.h"

#include <algorithm>
#include <math.h>
#include <stdlib.h>
#include <vector>
#include <utility>

#define PRIMARY 0
#define BACKUP 1

bool vector_size_comparator(const std::pair<int, std::vector<int>>& a,
                            const std::pair<int, std::vector<int>>& b) {
  return a.second.size() < b.second.size();
}

// Takes a physical network graph phys_topology, a set of location constraints
// location_constraints, and a seed mapping of virtual nodes and returns a
// vector containing exactly two vectors. Each of the vectors in the returned
// vector contains a mapping of the virtual nodes. The 0-th vector represents
// the primary node mapping and the 1-th vector represents the backup node
// mapping, respectively. 
std::unique_ptr<std::vector<std::vector<int>>> CreateInitialNodeMap(
    const Graph* phys_topology,
    const std::vector<std::vector<int>>& location_constraints,
    const std::pair<int, std::vector<int>>& seed) {
  auto tid = pthread_self();
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
  primary.push_back(seed.second[PRIMARY]);
  backup.push_back(seed.second[BACKUP]);
  for (int i = 0; i < temp_loc_constraints.size(); ++i) {
    int vnode = temp_loc_constraints[i].first;
    if (vnode == mapped_node) continue;
    int best_candidate = NIL;
    for (int j = 0; j < temp_loc_constraints[i].second.size(); ++j) {
      int candidate = temp_loc_constraints[i].second[j];
      if (!taken[candidate]) {
        // Try a number of tests on the candidate node and determine if it is a
        // better fit compared to the best_candidate.
        if (IsFeasibleBetterAssignment(phys_topology, primary, backup,
                                       candidate, best_candidate)) {
          best_candidate = candidate;
        }
      }
    }
    if (best_candidate != NIL) {
      taken[best_candidate] = true;
      (*node_maps)[PRIMARY][vnode] = best_candidate;
      primary.push_back(best_candidate);
    }

    // Do the same to find backup node mapping.
    best_candidate = NIL;
    for (int j = 0; j < temp_loc_constraints[i].second.size(); ++j) {
      int candidate = temp_loc_constraints[i].second[j];
      if (!taken[candidate]) {
        if (IsFeasibleBetterAssignment(phys_topology, backup, primary,
                                       candidate, best_candidate)) {
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

// Given a physical topology phys_topology, a set of location constraints for
// virtual nodes location_constraints, vectors primary and backup containing the
// nodes in the primary and backup mapping, respectively, an initial node
// mapping and a seed node mapping, this function re-embeds the virtual nodes
// inside the primary and backup mapping to obtain a more optimized node
// mapping. The return type is a vector of vectors, where each vector contains
// a one to one mapping between a virtual and physical node.
std::unique_ptr<std::vector<std::vector<int>>> ReEmbedVirtualNodes(
    const Graph* phys_topology,
    const std::vector<std::vector<int>>& location_constraints,
    const std::vector<int>& primary, const std::vector<int>& backup,
    const std::vector<std::vector<int>>& initial_node_map,
    const std::pair<int, std::vector<int>>& seed) {
  pthread_t tid = pthread_self();
  std::unique_ptr<std::vector<std::vector<int>>> node_maps(
      new std::vector<std::vector<int>>(2));
  (*node_maps)[PRIMARY] = std::vector<int>(initial_node_map[PRIMARY]);
  (*node_maps)[BACKUP] = std::vector<int>(initial_node_map[BACKUP]);
  int mapped_node = seed.first;
  std::vector<std::pair<int, std::vector<int>>> temp_loc_constraints;
  for (int i = 0; i < location_constraints.size(); ++i) {
    temp_loc_constraints.push_back(std::make_pair(i, location_constraints[i]));
  }
  // Sort location constraints according to the number of constraints (low to
  // high).
  std::sort(temp_loc_constraints.begin(), temp_loc_constraints.end(),
            vector_size_comparator);

  std::vector<int> partition(phys_topology->node_count(), NIL);
  for (auto& n : primary) partition[n] = PRIMARY;
  for (auto& n : backup) partition[n] = BACKUP;

  std::vector<bool> taken(phys_topology->node_count(), false);
  for (auto& m : *node_maps) {
    for (auto& n : m) {
      taken[n] = true;
    }
  }
  for (int i = 0; i < temp_loc_constraints.size(); ++i) {
    int vnode = temp_loc_constraints[i].first;
    int best_candidate = (*node_maps)[PRIMARY][vnode];
    int old_mapping = (*node_maps)[PRIMARY][vnode];
    // For each virtual node iterate over its location constraint set and check
    // if any other candidate inside the current partition can yield a "better"
    // node mapping. Better in this case is computed as the mean shortest path
    // between the physical nodes belonging to a mapping.
    for (int j = 0; j < temp_loc_constraints[i].second.size(); ++j) {
      int candidate = temp_loc_constraints[i].second[j];
      if (!taken[candidate] && partition[candidate] == PRIMARY) {
        double mean_sp_old = MeanSubsetShortestPathLength(
            phys_topology, primary, backup, (*node_maps)[PRIMARY]);
        (*node_maps)[PRIMARY][vnode] = candidate;
        double mean_sp_new = MeanSubsetShortestPathLength(
            phys_topology, primary, backup, (*node_maps)[PRIMARY]);
        (*node_maps)[PRIMARY][vnode] = old_mapping;
        if (mean_sp_new < mean_sp_old) {
          best_candidate = candidate;
        }
      }
    }

    (*node_maps)[PRIMARY][vnode] = best_candidate;
    taken[old_mapping] = false;
    taken[best_candidate] = true;

    // Perform the similar tasks for backup nodes.
    best_candidate = (*node_maps)[BACKUP][vnode];
    old_mapping = (*node_maps)[BACKUP][vnode];
    for (int j = 0; j < temp_loc_constraints[i].second.size(); ++j) {
      int candidate = temp_loc_constraints[i].second[j];
      if (!taken[candidate] && partition[candidate] == BACKUP) {
        double mean_sp_old = MeanSubsetShortestPathLength(
            phys_topology, backup, primary, (*node_maps)[BACKUP]);
        (*node_maps)[BACKUP][vnode] = candidate;
        double mean_sp_new = MeanSubsetShortestPathLength(
            phys_topology, backup, primary, (*node_maps)[BACKUP]);
        (*node_maps)[BACKUP][vnode] = old_mapping;
        if (mean_sp_new < mean_sp_old) {
          best_candidate = candidate;
        }
      }
    }
    (*node_maps)[BACKUP][vnode] = best_candidate;
    taken[old_mapping] = false;
    taken[best_candidate] = true;
  }
  return std::move(node_maps);
}


// Based on a primary and backup node mapping, node_maps[0] and node_maps[1],
// respectively, partitions the physical network phys_topology into two disjoint
// partitions. The return type is a vector of vector, where the 0-th vector
// contains the primary and the 1-th vector contains the backup partitoin,
// respectively.
std::unique_ptr<std::vector<std::vector<int>>> PartitionGraph(
    const Graph* phys_topology, std::vector<std::vector<int>>* node_maps) {
  std::unique_ptr<std::vector<std::vector<int>>> partitions(
      new std::vector<std::vector<int>>(2));
  std::vector<bool> taken(phys_topology->node_count(), false);
  auto& primary = partitions->at(PRIMARY);
  auto& backup = partitions->at(BACKUP);
  std::vector<int> primary_seed, backup_seed;
  // The node mappings are seed for the primary and backup partitions.
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
      // Feasibility Test.
      if (!IsFeasiblePartition(phys_topology, primary, backup_seed, i)) {
        backup.push_back(i);
      } else if (!IsFeasiblePartition(phys_topology, backup, primary_seed, i)) {
        primary.push_back(i);
      } else {
        // Compact Partition Test.
        double primary_sp_reduction = MeanShortestPathReduction(
            phys_topology, primary, backup, primary_seed, i);
        double backup_sp_reduction = MeanShortestPathReduction(
            phys_topology, backup, primary, backup_seed, i);
        if (primary_sp_reduction > backup_sp_reduction) {
          primary.push_back(i);
        } else {
          // Connectivity Contribution Test.
          // First check how much contribution is made in terms of reduction in
          // number of components.
          int component_decrease_primary =
              NumConnectedComponentsDecrease(phys_topology, primary, i);
          int component_decrease_backup =
              NumConnectedComponentsDecrease(phys_topology, backup, i);
          if (component_decrease_primary > component_decrease_backup) {
            primary.push_back(i);
          } else if (component_decrease_primary < component_decrease_backup) {
            backup.push_back(i);
          } else {
            // Then check how many edges from the candidate is incident to the
            // primary and backup partitions, respectively.
            int cut_edge_primary = NumCutEdges(phys_topology, primary, i);
            int cut_edge_backup = NumCutEdges(phys_topology, backup, i);
            if (cut_edge_primary > cut_edge_backup) {
              primary.push_back(i);
            } else if (cut_edge_primary < cut_edge_backup) {
              backup.push_back(i);
            } else {
              // In case of a tie on previous criteria, balance the partitions.
              if (primary.size() > backup.size()) {
                backup.push_back(i);
              } else {
                primary.push_back(i);
              }
            }
          }
        }
      }
    }
  }
  return std::move(partitions);
}

// Given the mapped endpoints of a virtual link map_u, and map_v and a bandwidth
// requirement bw, this function embeds the virtul link on the shortest path
// between map_u and map_v located within partition and having bandwidth bw.
std::vector<std::pair<int, int>> EmbedVLink(const Graph* phys_topology,
                                            const std::vector<int>& partition,
                                            int map_u, int map_v, long bw) {
  auto mapped_path = BFS(phys_topology, partition, map_u, map_v, bw);
  std::vector<std::pair<int, int>> edge_map;
  for (int i = 0; i < ((int)mapped_path->size()) - 1; ++i) {
    edge_map.push_back(
        std::pair<int, int>(mapped_path->at(i), mapped_path->at(i + 1)));
  }
  return std::move(edge_map);
}

// Embed all the virtual links and return the graph embedding.
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
      for (auto& e : emap) {
        phys_topology->reduce_edge_residual_bandwidth(e.first, e.second, bw);
      }
      edge_map->insert(
          std::make_pair(std::pair<int, int>(u, v), std::move(emap)));
    }
  }
  return std::move(edge_map);
}

// Compute the cost of an embedding.
long EmbeddingCost(const Graph* phys_topology, const Graph* virt_topology,
                   const VNEmbedding* embedding) {
  long cost = 0;
  for (auto emap_it = embedding->primary_edge_map->begin();
       emap_it != embedding->primary_edge_map->end(); ++emap_it) {
    auto& vlink = emap_it->first;
    auto& plinks = emap_it->second;
    for (auto& e : plinks) {
      cost += phys_topology->get_edge_cost(e.first, e.second) *
              virt_topology->get_edge_bandwidth(vlink.first, vlink.second);
    }
  }
  for (auto semap_it = embedding->backup_edge_map->begin();
       semap_it != embedding->backup_edge_map->end(); ++semap_it) {
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
