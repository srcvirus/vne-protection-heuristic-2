#ifndef DATASTRUCTURE_H_
#define DATASTRUCTURE_H_

#include <list>
#include <set>
#include <string>
#include <sstream>
#include <vector>
#include <math.h>
#include <map>
#include <memory>
#include <stdlib.h>

#define INF 999999999
#define MAXN 1000
#define NIL -1

class Edge {
 public:
  Edge() {}
  Edge(int u, int v, int channels, int cost)
      : u_(u), v_(v), channels_(channels), cost_(cost) {}

  int u() const { return u_; }
  int v() const { return v_; }
  int channels() const { return channels_; }
  int cost() const { return cost_; }

  bool operator<(const Edge& e) const {
    if (channels_ != e.channels_)
      return channels_ < e.channels_;
    else
      return (u_ <= e.u_ || v_ <= e.v_);
  }
  bool operator==(const Edge& e) const { return u_ == e.u_ && v_ == e.v_; }
  std::string GetDebugString() {
    return "[u = " + std::to_string(u_) + ", v = " + std::to_string(v_) +
           ", channels = " + std::to_string(channels_) + "]";
  }

 private:
  int u_, v_;
  int channels_, cost_;
};

// An entry in an adjacent list. An entry contains the node_id of the endpoint.
// The entry contains channels, residual channels, delay and cost of the
// corresponding edge.
struct edge_endpoint {
  int node_id;
  int channels, residual_channels;
  int delay, cost;
  std::set<int> available_channels;
  edge_endpoint(int node_id, int ch, int delay, int cost)
      : node_id(node_id),
        channels(ch),
        delay(delay),
        residual_channels(ch),
        cost(cost) {
    available_channels.clear();
    for (int i = 0; i < channels; ++i) available_channels.insert(i);
  }
  std::string GetDebugString() {
    std::string str_available_channels = "[";
    for (auto ch_it = available_channels.begin();
         ch_it != available_channels.end(); ++ch_it) {
      str_available_channels += " ";
      str_available_channels += std::to_string(*ch_it);
    }
    return "ndoe_id = " + std::to_string(node_id) + ", channels = " +
           std::to_string(channels) + ", delay = " + std::to_string(delay) +
           ", cost = " + std::to_string(cost) + ", available_channels = " +
           str_available_channels;
  }
};

class Graph {
 public:
  Graph() {
    adj_list_ = std::unique_ptr<std::vector<std::vector<edge_endpoint>>>(
        new std::vector<std::vector<edge_endpoint>>);
    node_count_ = edge_count_ = 0;
  }

  Graph(const Graph& graph) {
    adj_list_ = std::unique_ptr<std::vector<std::vector<edge_endpoint>>>(
        new std::vector<std::vector<edge_endpoint>>);
    for (int i = 0; i < graph.node_count(); ++i) {
      auto& adj_list = graph.adj_list()->at(i);
      for (auto& neighbor : adj_list) {
        if (i > neighbor.node_id) {
          this->add_edge(i, neighbor.node_id, neighbor.channels, neighbor.delay,
                         neighbor.cost);
        }
      }
    }
  }

  // Accessor methods.
  int node_count() const { return node_count_; }
  int edge_count() const { return edge_count_; }
  const std::vector<std::vector<edge_endpoint>>* adj_list() const {
    return static_cast<const std::vector<std::vector<edge_endpoint>>*>(
        adj_list_.get());
  }

  // u and v are 0-based identifiers of an edge endpoint. An edge is
  // bi-directional, i.e., calling Graph::add_edge with u = 1, v = 3 will add
  // both (1, 3) and (3, 1) in the graph.
  int add_edge(int u, int v, int ch, int delay, int cost) {
    if (adj_list_->size() < u + 1) adj_list_->resize(u + 1);
    if (adj_list_->size() < v + 1) adj_list_->resize(v + 1);
    adj_list_->at(u).push_back(edge_endpoint(v, ch, delay, cost));
    adj_list_->at(v).push_back(edge_endpoint(u, ch, delay, cost));
    ++edge_count_;
    node_count_ = adj_list_->size();
  }

  void reduce_edge_residual_channels(int u, int v, int ch_delta) {
    auto& u_neighbors = adj_list_->at(u);
    for (auto& end_point : u_neighbors) {
      if (end_point.node_id == v) {
        end_point.residual_channels -= ch_delta;
      }
    }
    auto& v_neighbors = adj_list_->at(v);
    for (auto& end_point : v_neighbors) {
      if (end_point.node_id == u)  {
        end_point.residual_channels -= ch_delta;
      }
    }
  }

  bool remove_edge_available_channel(int u, int v, int ch_index) {
    if (ch_index > get_edge_channels(u, v) - 1) 
      return false;
    auto& u_neighbors = adj_list_->at(u);
    for (auto& end_point : u_neighbors) {
      if (end_point.node_id == v) {
        end_point.available_channels.erase(ch_index);
      }
    }
    auto& v_neighbors = adj_list_->at(v);
    for (auto& end_point : v_neighbors) {
      if (end_point.node_id == u) {
        end_point.available_channels.erase(ch_index);
      }
    }
    return true;
  }

  int get_edge_cost(int u, int v) const {
    auto& neighbors = adj_list_->at(u);
    for (auto& end_point : neighbors) {
      if (end_point.node_id == v) return end_point.cost;
    }
  }

  int get_edge_channels(int u, int v) const {
    auto& neighbors = adj_list_->at(u);
    for (auto& end_point : neighbors) {
      if (end_point.node_id == v) return end_point.channels;
    }
  }

  int get_edge_residual_channels(int u, int v) const {
    auto& neighbors = adj_list_->at(u);
    for (auto& end_point : neighbors) {
      if (end_point.node_id == v) return end_point.residual_channels;
    }
  }

  // Export the graph as a list of edges.
  std::unique_ptr<std::vector<Edge>> ExportAsEdgeList() {
    std::unique_ptr<std::vector<Edge>> edge_list(new std::vector<Edge>());
    for (int u = 0; u < this->node_count(); ++u) {
      auto& neighbors = this->adj_list()->at(u);
      for (auto& end_point : neighbors) {
        int v = end_point.node_id;
        int channels = end_point.channels;
        int cost = end_point.cost;
        if (u < v) continue;
        edge_list->push_back(Edge(u, v, channels, cost));
      }
    }
    return std::move(edge_list);
  }

  std::string GetDebugString() const {
    std::string ret_string = "node_count = " + std::to_string(node_count_);
    ret_string += ", edge_count = " + std::to_string(edge_count_) + "\n";
    for (int i = 0; i < node_count_; ++i) {
      auto& neighbors = adj_list_->at(i);
      ret_string += std::to_string(i) + " --> ";
      for (auto& neighbor : neighbors) {
        ret_string += " (" + neighbor.GetDebugString() + ")";
      }
      ret_string += "\n";
    }
    return ret_string;
  }

 private:
  std::unique_ptr<std::vector<std::vector<edge_endpoint>>> adj_list_;
  int node_count_, edge_count_;
};

class DisjointSet {
 public:
  DisjointSet() {}
  DisjointSet(int n, const std::vector<int>& subset) { Initialize(n, subset); }

  void Union(int x, int y) {
    int x_root = Find(x), y_root = Find(y);
    if (x_root != y_root) {
      if (rank_[x_root] > rank_[y_root])
        set_[y_root] = x_root;
      else if (rank_[x_root] < rank_[y_root])
        set_[x_root] = y_root;
      else {
        set_[y_root] = x_root;
        ++rank_[x_root];
      }
    }
    // this->Print();
  }

  int Find(int x) {
    if (x != set_[x]) set_[x] = Find(set_[x]);
    return set_[x];
  }

  int CountDistinct() {
    int ret = 0;
    std::vector<int> count(n_, 0);
    for (int i = 0; i < set_.size(); ++i) {
      if (valid_mask_[i]) ++count[Find(i)];
    }
    for (auto& x : count)
      if (x != 0) ++ret;
    return ret;
  }

  int n() { return n_; }
  void AddValidElement(int k) { valid_mask_[k] = true; }

  void Print() {
    for (int x = 0; x < set_.size(); ++x) {
      if (valid_mask_[x]) printf("|%d: %d|", x, Find(x));
      // printf("|%d: %d|", x, set_[x]);
    }
    printf("\n");
  }

 private:
  void Initialize(int n, const std::vector<int>& subset) {
    n_ = n;
    set_.resize(n, 0);
    rank_.resize(n, 0);
    valid_mask_.resize(n, false);
    for (int i = 0; i < n; ++i) set_[i] = i;
    for (int i = 0; i < subset.size(); ++i) valid_mask_[subset[i]] = true;
  }
  int n_;
  std::vector<int> set_;
  std::vector<int> rank_;
  std::vector<bool> valid_mask_;
};

struct VNEmbedding {
  std::unique_ptr<std::vector<std::vector<int>>> node_map;
  std::unique_ptr<std::map<std::pair<int, int>,
                           std::pair<int, std::vector<std::pair<int, int>>>>>
      primary_edge_map;
  std::unique_ptr<std::map<std::pair<int, int>,
                           std::pair<int, std::vector<std::pair<int, int>>>>>
      backup_edge_map;
  long cost;
};

struct ThreadParameter {
  int vnode_seed;
  int primary_seed;
  int backup_seed;
  const Graph* phys_topology;
  const Graph* virt_topology;
  const std::vector<std::vector<int>>* location_constraints;
};
#endif  // MIDDLEBOX_PLACEMENT_SRC_DATASTRUCTURE_H_
