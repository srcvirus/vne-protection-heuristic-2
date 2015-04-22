#include "datastructure.h"
#include "heuristic.h"
#include "io.h"
#include "util.h"
#include "vne_heuristic_solver.h"

#include <chrono>
#include <pthread.h>

void* EmbedVNThread(void* args) {
  ThreadParameter* parameter = reinterpret_cast<ThreadParameter*>(args);
  std::unique_ptr<Graph> phys_topology(new Graph(*(parameter->phys_topology)));
  std::unique_ptr<Graph> virt_topology(new Graph(*(parameter->virt_topology)));
  std::unique_ptr<std::vector<std::vector<int>>> location_constraints(
      new std::vector<std::vector<int>>(*(parameter->location_constraints)));
  std::pair<int, std::vector<int>> seed;
  seed.first = parameter->vnode_seed;
  seed.second.emplace_back(parameter->primary_seed);
  seed.second.emplace_back(parameter->backup_seed);
  int cost = 0;
  auto node_map = CreateInitialNodeMap(phys_topology.get(),
                                       *(parameter->location_constraints),
                                       seed);
  // for (int i = 0; i < (*node_map)[PRIMARY].size(); ++i)
  //  printf(" %d", (*node_map)[PRIMARY][i]);
  // printf("\n");
  // for (int i = 0; i < (*node_map)[BACKUP].size(); ++i)
  //   printf(" %d", (*node_map)[BACKUP][i]);
  // printf("\n");
  bool node_map_failed = false;
  for (auto& m : *node_map) {
    for (auto& n : m) {
      if (n == -1) {
        node_map_failed = true;
        break;
      }
    }
    if (node_map_failed) break;
  }
  if (node_map_failed) {
    cost = INF;
  }

  auto& primary_node_map = (*node_map)[PRIMARY];
  auto& backup_node_map = (*node_map)[BACKUP];
  auto partitions = PartitionGraph(phys_topology.get(), node_map.get());
  auto& primary_partition = (*partitions)[PRIMARY];
  auto& backup_partition = (*partitions)[BACKUP];

  bool edge_map_failed = false;
  auto emap = EmbedVN(phys_topology.get(), virt_topology.get(),
                      primary_partition, primary_node_map);
  for (auto emap_it = emap->begin(); emap_it != emap->end(); ++emap_it) {
    if (emap_it->second.empty()) {
      edge_map_failed = true;
      cost = INF;
      break;
    }
  }

  auto semap = EmbedVN(phys_topology.get(), virt_topology.get(),
                       backup_partition, backup_node_map);
  for (auto semap_it = semap->begin(); semap_it != semap->end(); ++semap_it) {
    if (semap_it->second.empty()) {
      edge_map_failed = true;
      cost = INF;
      break;
    }
  }
  VNEmbedding* embedding = new VNEmbedding();
  embedding->node_map = std::move(node_map);
  embedding->primary_edge_map = std::move(emap);
  embedding->backup_edge_map = std::move(semap);
  if (cost != INF) cost = EmbeddingCost(phys_topology.get(), virt_topology.get(), embedding);
  embedding->cost = cost;
  pthread_exit(reinterpret_cast<void*>(embedding));
}

std::unique_ptr<VNEmbedding> ProtectedVNE(const Graph* phys_topology, const Graph* virt_topology,
                  const std::vector<std::vector<int>>& location_constraints) {
  int min_constraint_size = INF, most_constrained_vnode = NIL;
  for (int i = 0; i < location_constraints.size(); ++i) {
    if (location_constraints[i].size() < min_constraint_size) {
      min_constraint_size = location_constraints[i].size();
      most_constrained_vnode = i;
    }
  }
  int n_threads = (min_constraint_size * (min_constraint_size - 1)) / 2;
  std::vector<pthread_t> threads(n_threads);
  int thread_id = 0;
  std::vector<std::unique_ptr<ThreadParameter>> parameters(n_threads);
  for (auto& parameter : parameters) {
    parameter = std::unique_ptr<ThreadParameter>(new ThreadParameter());
  }
  for (int i = 0; i < location_constraints[most_constrained_vnode].size();
       ++i) {
    for (int j = i + 1; j < location_constraints[most_constrained_vnode].size();
         ++j) {
      int primary = location_constraints[most_constrained_vnode][i];
      int backup = location_constraints[most_constrained_vnode][j];
      ThreadParameter* parameter = parameters[thread_id].get();
      parameter->vnode_seed = most_constrained_vnode;
      parameter->primary_seed = primary;
      parameter->backup_seed = backup;
      parameter->phys_topology = phys_topology;
      parameter->virt_topology = virt_topology;
      parameter->location_constraints = &location_constraints;
      pthread_create(&threads[thread_id++], NULL, &EmbedVNThread, parameter);
    }
  }
  std::vector<std::unique_ptr<VNEmbedding>> embeddings(n_threads);
  for (int i = 0; i < n_threads; ++i) {
    void* ret_value;
    pthread_join(threads[i], &ret_value);
    embeddings[i] =
        std::unique_ptr<VNEmbedding>(reinterpret_cast<VNEmbedding*>(ret_value));
  }
  int min_cost_embedding = 0;
  for (int i = 0; i < embeddings.size(); ++i) {
    if (embeddings[i]->cost < embeddings[min_cost_embedding]->cost) {
      min_cost_embedding = i;
    }
  }
  return std::move(embeddings[min_cost_embedding]);
}

int main(int argc, char* argv[]) {
  using std::string;
  auto arg_map = ParseArgs(argc, argv);
  string pn_topology_filename = "";
  string vn_topology_filename = "";
  string location_constraint_filename = "";
  for (auto argument : *arg_map) {
    if (argument.first == "--pn_topology_file") {
      pn_topology_filename = argument.second;
    } else if (argument.first == "--vn_topology_file") {
      vn_topology_filename = argument.second;
    } else if (argument.first == "--location_constraint_file") {
      location_constraint_filename = argument.second;
    } else {
      printf("Invalid command line option: %s\n", argument.first.c_str());
      return 1;
    }
  }

  auto physical_topology =
      InitializeTopologyFromFile(pn_topology_filename.c_str());
  DEBUG(physical_topology->GetDebugString().c_str());
  auto virt_topology = InitializeTopologyFromFile(vn_topology_filename.c_str());
  int offset = virt_topology->node_count();
  DEBUG(virt_topology->GetDebugString().c_str());
  auto location_constraints = InitializeVNLocationsFromFile(
      location_constraint_filename.c_str(), virt_topology->node_count());

  auto solution_start_time = std::chrono::high_resolution_clock::now();

  auto solution_end_time = std::chrono::high_resolution_clock::now();
  auto embedding = ProtectedVNE(physical_topology.get(), virt_topology.get(),
                                    *location_constraints);
  unsigned long long elapsed_time =
      std::chrono::duration_cast<std::chrono::nanoseconds>(
          solution_end_time - solution_start_time).count();
  printf("Solution time: %llu.%llus\n", elapsed_time / 1000000000LL,
         elapsed_time % 1000000000LL);

  // Write node maps to file.
  FILE* nmap_file = fopen((pn_topology_filename + ".nmap").c_str(), "w");
  auto& primary_node_map = (*(embedding->node_map))[PRIMARY];
  for (int i = 0; i < primary_node_map.size(); ++i) {
    fprintf(nmap_file, "Virtual node %d --> physical node %d\n", i,
            primary_node_map[i]);
  }
  fclose(nmap_file);

  // Write shadow nmaps to file.
  FILE* snmap_file = fopen((pn_topology_filename + ".snmap").c_str(), "w");
  auto& backup_node_map = (*(embedding->node_map))[BACKUP];
  for (int i = 0; i < backup_node_map.size(); ++i) {
    fprintf(snmap_file, "Shadow virtual node of %d --> physical node %d\n", i,
            backup_node_map[i]);
  }
  fclose(snmap_file);

  // Write edge maps to file.
  FILE* emap_file = fopen((pn_topology_filename + ".emap").c_str(), "w");
  auto& emap = embedding->primary_edge_map;
  for (auto emap_it = emap->begin(); emap_it != emap->end(); ++emap_it) {
    auto& vlink = emap_it->first;
    auto& plinks = emap_it->second;
    for (auto& e : plinks) {
      fprintf(emap_file, "Virtual edge (%d, %d) --> physical edge (%d, %d)\n",
              vlink.first, vlink.second, e.first, e.second);
    }
  }
  fclose(emap_file);

  // Write shadow edge maps to file.
  FILE* semap_file = fopen((pn_topology_filename + ".semap").c_str(), "w");
  auto& semap = embedding->backup_edge_map;
  for (auto semap_it = semap->begin(); semap_it != semap->end(); ++semap_it) {
    auto& vlink = semap_it->first;
    auto& plinks = semap_it->second;
    for (auto& e : plinks) {
      fprintf(emap_file,
              "Shadow virtual edge of (%d, %d) --> physical edge (%d, %d)\n",
              vlink.first, vlink.second, e.first, e.second);
    }
  }
  fclose(semap_file);

  // Write cost to file.
  FILE* cost_file = fopen((pn_topology_filename + ".cost").c_str(), "w");
  fprintf(cost_file, "%ld\n", embedding->cost);
  fclose(cost_file);

  return 0;
}
