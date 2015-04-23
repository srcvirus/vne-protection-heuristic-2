#include "datastructure.h"
#include "util.h"

#include <algorithm>
#include <queue>

std::unique_ptr<std::vector<int>> BFS(const Graph* graph,
                                      const std::vector<int>& partition,
                                      int src, int dest, long bw) {
  const int kNodeCount = graph->node_count();
  std::vector<bool> visited(kNodeCount, false);
  std::vector<bool> in_partition(kNodeCount, false);
  std::vector<int> pre(kNodeCount, NIL);
  std::queue<int> Q;
  for (auto& node : partition) in_partition[node] = true;
  Q.push(src);
  visited[src] = true;
  while (!Q.empty()) {
    int u = Q.front();
    Q.pop();
    visited[u] = true;
    auto& u_neighbors = graph->adj_list()->at(u);
    for (auto& end_point : u_neighbors) {
      int v = end_point.node_id;
      if (!visited[v] && in_partition[u] && in_partition[v] &&
          end_point.residual_bandwidth >= bw) {
        Q.push(v);
        pre[v] = u;
        if (v == dest) break;
      }
    }
  }
  std::unique_ptr<std::vector<int>> path(new std::vector<int>());
  int node = dest;
  while (pre[node] != NIL) {
    path->push_back(node);
    node = pre[node];
  }
  if (!path->empty()) path->push_back(node);
  std::reverse(path->begin(), path->end());
  return std::move(path);
}

int NumConnectedComponentsDecrease(const Graph* graph,
                                   const std::vector<int>& subset,
                                   int candidate_node) {
  // printf("subset:");
  // for (int a : subset) printf(" %d", a);
  // printf("\n");
  std::vector<bool> is_subset_member(graph->node_count(), false);
  for (auto& node : subset) is_subset_member[node] = true;
  DisjointSet ds(graph->node_count(), subset);
  for (auto& node : subset) {
    auto& adj_list = graph->adj_list()->at(node);
    for (auto& adj_node : adj_list) {
      if (!is_subset_member[adj_node.node_id]) continue;
      ds.Union(node, adj_node.node_id);
    }
  }
  // ds.Print();
  int current_connected_components = ds.CountDistinct();
  ds.AddValidElement(candidate_node);
  auto& adj_list = graph->adj_list()->at(candidate_node);
  for (auto& node : adj_list) {
    if (!is_subset_member[node.node_id]) continue;
    ds.Union(candidate_node, node.node_id);
  }
  // ds.Print();
  int future_connected_components = ds.CountDistinct();
  // printf("current = %d, future = %d\n", current_connected_components,
  //        future_connected_components);
  return current_connected_components - future_connected_components;
}

int NumCutEdges(const Graph* graph, const std::vector<int>& partition,
                int node) {
  int num_cut_edges = 0;
  std::vector<bool> is_partition_member(graph->node_count(), false);
  for (auto& node : partition) is_partition_member[node] = true;
  auto& adj_list = graph->adj_list()->at(node);
  for (auto& n : adj_list) {
    if (is_partition_member[n.node_id]) ++num_cut_edges;
  }
  return num_cut_edges;
}

// Returns true if addition of candidate_a to partition_a does not make it
// impossible for partition_b to become disconnected.
bool IsFeasiblePartition(const Graph* graph,
                         const std::vector<int>& partition_a,
                         const std::vector<int>& partition_b, int candidate_a) {
  std::vector<bool> is_partition_a(graph->node_count(), false);
  is_partition_a[candidate_a] = true;
  for (auto& node : partition_a) is_partition_a[node] = true;
  DisjointSet ds(graph->node_count(), partition_b);
  for (int i = 0; i < graph->node_count(); ++i) {
    auto& adj_list = graph->adj_list()->at(i);
    for (auto& neighbor : adj_list) {
      if (is_partition_a[i] || is_partition_a[neighbor.node_id]) continue;
      // printf("Merging %d and %d\n", i, neighbor.node_id);
      ds.Union(i, neighbor.node_id);
    }
  }
  // ds.Print();
  return ds.CountDistinct() <= 1;
}

bool IsFeasibleBetterAssignment(const Graph* graph, std::vector<int>& partition,
                                const std::vector<int>& forbidden,
                                int candidate, int best_candidate) {
  if (!IsFeasiblePartition(graph, partition, forbidden, candidate)) {
    // printf("Not feasible to assign %d\n", candidate);
    return false;
  }
  if (best_candidate == NIL) return true;
  int component_decrease_current =
      NumConnectedComponentsDecrease(graph, partition, candidate);
  int component_decrease_best =
      NumConnectedComponentsDecrease(graph, partition, best_candidate);
  if (component_decrease_current > component_decrease_best) {
    return true;
  } else if (component_decrease_current == component_decrease_best) {
    int cut_edge_candidate = NumCutEdges(graph, partition, candidate);
    int cut_edge_current = NumCutEdges(graph, partition, best_candidate);
    if (cut_edge_candidate > cut_edge_current) {
      return true;
    }
  }
  return false;
}

std::unique_ptr<std::vector<int>> ShortestPathLengthVector(const Graph* graph, int source, 
                              const std::vector<int>& destinations,
                              const std::vector<int>& subset) {
  std::vector<bool> is_subset(graph->node_count(), false);
  for (auto& node : subset) is_subset[node] = true;
  std::queue<int> Q;
  std::vector<bool> visited(graph->node_count(), false);
  std::vector<int> distance(graph->node_count(), INF);
  Q.push(source);
  distance[source] = 0;
  while(!Q.empty()) {
    int u = Q.front();
    Q.pop();
    visited[u] = true;
    auto& neighbors = graph->adj_list()->at(u);
    for (auto& node : neighbors) {
      if (!visited[node.node_id] && is_subset[node.node_id]) {
        Q.push(node.node_id);
        distance[node.node_id] = distance[u] + 1;
      }
    }
  }
  std::unique_ptr<std::vector<int>> ret_vector(new std::vector<int>());
  for (auto& v : destinations) {
    ret_vector->emplace_back(distance[v]);
  }
  return std::move(ret_vector);
}

double MeanSubsetShortestPathLength(const Graph* graph, 
                                    const std::vector<int>& partition,
                                    const std::vector<int>& mapping_subset) {
  double distance = 0.0;
  int paths = 0;
  for (auto& node : mapping_subset) {
    auto path_vector = ShortestPathLengthVector(graph, node, mapping_subset, partition);
    for (int i = 0; i < path_vector->size(); ++i) {
      distance += path_vector->at(i);
    }
    paths += path_vector->size();
  }
  double mean_path_length = distance / static_cast<double>(paths);
  return mean_path_length;
}


