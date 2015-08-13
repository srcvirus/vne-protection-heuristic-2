#include "datastructure.h"
#include "util.h"

#include <algorithm>
#include <queue>

// Perform a BFS search from src to dest on graph within set partition.
// Discard links with channels < bw. Returns the path from src to dest
// as a pointer to vector.
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
          end_point.residual_channels >= bw) {
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

struct dijkstra_node {
  int u, cost;
  dijkstra_node(int u, int cost): u(u), cost(cost) {}
  bool operator < (const dijkstra_node dnode) const {
    return cost < dnode.cost;
  }
};

std::unique_ptr<std::pair<int, std::unique_ptr<std::vector<int>>>> dwdm_bfs(
    const Graph* graph,
    const std::vector<int>& partition,
    int src, int dest, int ch) {
  const int kNodeCount = graph->node_count();
  std::vector<bool> in_partition(kNodeCount, false);
  for (auto& node : partition) in_partition[node] = true;
  int ch_lo = 0, ch_hi = 6;
  std::unique_ptr<std::pair<int,std::unique_ptr<std::vector<int>>>> dwdm_path(
      new std::pair<int, std::unique_ptr<std::vector<int>>>());
  dwdm_path->first = NIL;
  dwdm_path->second = std::unique_ptr<std::vector<int>>(new std::vector<int>());
  DEBUG("Starting BS\n");
  while(ch_lo < ch_hi) {
    int ch_mid = (ch_lo + ch_hi) / 2;
    std::priority_queue<dijkstra_node> Q;
    std::vector<int> pre(kNodeCount, NIL);
    std::vector<int> d(kNodeCount, INF);
    Q.push(dijkstra_node(src,0));
    d[src] = 0;
    DEBUG("Starting dijkstra\n");
    while(!Q.empty()) {
      dijkstra_node dnode = Q.top();
      int u = dnode.u;
      Q.pop();
      auto& u_neighbors = graph->adj_list()->at(u);
      for (auto& end_point : u_neighbors) {
        int v = end_point.node_id;
        bool has_valid_channel = (end_point.available_channels.find(ch_mid) 
            != end_point.available_channels.end());
        if (has_valid_channel && in_partition[u] && in_partition[v] && 
            end_point.residual_channels >= ch && d[v] > d[u] + end_point.cost) {
          d[v] = d[u] + end_point.cost;
          Q.push(dijkstra_node(v, d[v]));
          pre[v] = u;
        }
      }
    }
    DEBUG("Dijkstra completed\n");
    if (d[dest] != INF) {
      ch_hi = ch_mid;
      int node = dest;
      dwdm_path->first = ch_hi;
      std::vector<int>* path = dwdm_path->second.get();
      path->clear();
      while(pre[node] != NIL) {
        path->push_back(node);
        node = pre[node];
      }
      if (!path->empty()) path->push_back(node);
      std::reverse(path->begin(), path->end());
    } else {
      ch_lo = ch_mid + 1;
    }
  }
  return std::move(dwdm_path);
}

// Returns the reduction in number of components when candidate_node
// is added to the nodes in subset in graph.
int NumConnectedComponentsDecrease(const Graph* graph,
                                   const std::vector<int>& subset,
                                   int candidate_node) {
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
  int current_connected_components = ds.CountDistinct();
  ds.AddValidElement(candidate_node);
  auto& adj_list = graph->adj_list()->at(candidate_node);
  for (auto& node : adj_list) {
    if (!is_subset_member[node.node_id]) continue;
    ds.Union(candidate_node, node.node_id);
  }
  int future_connected_components = ds.CountDistinct();
  return current_connected_components - future_connected_components;
}

// Returns number of links incident to nodes in partition in graph from node.
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
      ds.Union(i, neighbor.node_id);
    }
  }
  return ds.CountDistinct() <= 1;
}

// Given a source node and a set of destinations, this function returns a vector
// containing length of shortest path from node to all destinations in graph
// confined withing nodes in subset and avoiding nodes in forbidden.
std::unique_ptr<std::vector<int>> ShortestPathLengthVector(
    const Graph* graph, int source, const std::vector<int>& destinations,
    const std::vector<int>& subset, const std::vector<int>& forbidden) {
  std::vector<bool> is_subset(graph->node_count(), false);
  std::vector<bool> is_forbidden(graph->node_count(), false);
  for (auto& node : subset) is_subset[node] = true;
  for (auto& node : forbidden) is_forbidden[node] = true;
  std::queue<int> Q;
  std::vector<bool> visited(graph->node_count(), false);
  std::vector<int> distance(graph->node_count(), INF);
  Q.push(source);
  distance[source] = 0;
  while (!Q.empty()) {
    int u = Q.front();
    Q.pop();
    visited[u] = true;
    auto& neighbors = graph->adj_list()->at(u);
    for (auto& node : neighbors) {
      if (!visited[node.node_id] && !is_forbidden[node.node_id]) {
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

// Returns mean shortest path length between vertices in partition while
// avoiding verticesin forbidden under a given virtual network node mapping in
// graph.
double MeanSubsetShortestPathLength(const Graph* graph,
                                    const std::vector<int>& partition,
                                    const std::vector<int>& forbidden,
                                    const std::vector<int>& mapping_subset) {
  double distance = 0.0;
  int paths = 0;
  for (auto& node : mapping_subset) {
    auto path_vector = ShortestPathLengthVector(graph, node, mapping_subset,
                                                partition, forbidden);
    for (int i = 0; i < path_vector->size(); ++i) {
      distance += path_vector->at(i);
    }
    paths += path_vector->size();
  }
  double mean_path_length = distance / static_cast<double>(paths);
  return mean_path_length;
}

// Returns the mean shortest path of the partition U candidate minus the mean
// shortest path of partition U best_candidate, considering the forbidden nodes.
double MeanShortestPathReduction(const Graph* graph,
                                 const std::vector<int>& partition,
                                 const std::vector<int>& forbidden,
                                 int candidate, int best_candidate) {
  std::vector<int> temp_part_a(partition);
  temp_part_a.push_back(candidate);
  double part_a_new_mean_sp =
      MeanSubsetShortestPathLength(graph, temp_part_a, forbidden, temp_part_a);
  temp_part_a.pop_back();
  temp_part_a.push_back(best_candidate);
  double part_a_best_mean_sp =
      MeanSubsetShortestPathLength(graph, temp_part_a, forbidden, temp_part_a);
  return part_a_best_mean_sp - part_a_new_mean_sp;
}

double MeanShortestPathReduction(const Graph* graph,
                                 const std::vector<int>& partition,
                                 const std::vector<int>& forbidden,
                                 const std::vector<int>& destinations,
                                 int candidate) {
  std::vector<int> temp_part_a(partition);
  double part_a_old_mean_sp =
      MeanSubsetShortestPathLength(graph, partition, forbidden, destinations);
  temp_part_a.push_back(candidate);
  double part_a_new_mean_sp =
      MeanSubsetShortestPathLength(graph, partition, forbidden, destinations);
  return part_a_old_mean_sp - part_a_new_mean_sp;
}

// Returns true if it is feasible to assigne candidate to the partition instead
// of best_candidate, considering the forbidden set and the addition yields
// better partition compared to best_candidate. Returns false otherwise.
bool IsFeasibleBetterAssignment(const Graph* graph, std::vector<int>& partition,
                                const std::vector<int>& forbidden,
                                int candidate, int best_candidate) {
  if (!IsFeasiblePartition(graph, partition, forbidden, candidate)) {
    return false;
  }
  if (best_candidate == NIL) return true;
  double shortest_path_reduction = MeanShortestPathReduction(
      graph, partition, forbidden, candidate, best_candidate);
  if (shortest_path_reduction > 0.0) {
    return true;
  } else {
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
      } else {
        return false;
      }
    }
  }
  return false;
}
