#include "datastructure.h"
#include "heuristic.h"
#include "io.h"
#include "util.h"
#include "vne_heuristic_solver.h"

#include <chrono>
#include <pthread.h>
#include <unistd.h>

// One thread of execution for a protected virtual network embedding. It takes
// a pointer to a ThreadParameter object. The ThreadParameter object contains
// a pointer to a physical network topoogy, a virtual network topology, location
// constraints, and a seed node mapping. It returns the virtual network
// embedding .
void* EmbedVNThread(void* args) {
  // Extract the Thread parameters.
  int tid = pthread_self();
  ThreadParameter* parameter = reinterpret_cast<ThreadParameter*>(args);
  std::unique_ptr<Graph> phys_topology(new Graph(*(parameter->phys_topology)));
  std::unique_ptr<Graph> virt_topology(new Graph(*(parameter->virt_topology)));
  std::unique_ptr<std::vector<std::vector<int>>> location_constraints(
      new std::vector<std::vector<int>>(*(parameter->location_constraints)));
  std::pair<int, std::vector<int>> seed;
  seed.first = parameter->vnode_seed;
  seed.second.emplace_back(parameter->primary_seed);
  seed.second.emplace_back(parameter->backup_seed);
  VNEmbedding* embedding = new VNEmbedding();
  embedding->cost = INF;

  // Find an initial embedding of the virtual nodes.
  auto initial_node_map = CreateInitialNodeMap(
      phys_topology.get(), *(parameter->location_constraints), seed);
  bool node_map_failed = false;
  for (auto& m : *initial_node_map) {
    for (auto& n : m) {
      if (n == -1) {
        node_map_failed = true;
        break;
      }
    }
    if (node_map_failed) break;
  }
  if (node_map_failed) {
    embedding->cost = INF;
    pthread_exit(reinterpret_cast<void*>(embedding));
  }
  DEBUG("[%x]: Initial node mapping done\n", tid);
  // Partition the physical network for primary and backup embedding and then
  // try to re-embed the virtual nodes according to the obtained partition. The
  // idea is to improve the node embedding inside the obtained partitions.
  auto partitions = PartitionGraph(phys_topology.get(), initial_node_map.get());
  auto& primary_partition = (*partitions)[PRIMARY];
  auto& backup_partition = (*partitions)[BACKUP];
  DEBUG("[%x]: Partitioning complete\n", tid);
  auto node_map = ReEmbedVirtualNodes(
      phys_topology.get(), *(parameter->location_constraints),
      primary_partition, backup_partition, *initial_node_map, seed);
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
    pthread_exit(reinterpret_cast<void*>(embedding));
  }
  auto& primary_node_map = (*node_map)[PRIMARY];
  auto& backup_node_map = (*node_map)[BACKUP];
  bool edge_map_failed = false;
  DEBUG("Node mapping complete\n");
  // Primary embedding of virtual links.
  DEBUG("Mapping vlinks\n");
  auto emap = EmbedVN(phys_topology.get(), virt_topology.get(),
                      primary_partition, primary_node_map);
  for (auto emap_it = emap->begin(); emap_it != emap->end(); ++emap_it) {
    if (emap_it->second.second.empty()) {
      DEBUG("Mapping of vlink (%d, %d) failed\n", 
          emap_it->first.first, emap_it->first.second);
      edge_map_failed = true;
      break;
    }
  }
  if (edge_map_failed) {
    embedding->cost = INF;
    pthread_exit(reinterpret_cast<void*>(embedding));
  }
  DEBUG("Vlink mapping complete\n");
  // Backup embedding of virtual links.
  auto semap = EmbedVN(phys_topology.get(), virt_topology.get(),
                       backup_partition, backup_node_map);
  for (auto semap_it = semap->begin(); semap_it != semap->end(); ++semap_it) {
    if (semap_it->second.second.empty()) {
      DEBUG("Mapping of shadow vlink (%d, %d) failed\n", 
          semap_it->first.first, semap_it->first.second);
      edge_map_failed = true;
      break;
    }
  }
  if (edge_map_failed) {
    embedding->cost = INF;
    pthread_exit(reinterpret_cast<void*>(embedding));
  }
  embedding->node_map = std::move(initial_node_map);
  embedding->primary_edge_map = std::move(emap);
  embedding->backup_edge_map = std::move(semap);
  embedding->cost =
      EmbeddingCost(phys_topology.get(), virt_topology.get(), embedding);
  pthread_exit(reinterpret_cast<void*>(embedding));
}

// Takes a virtual network embedding object pointer and writes various
// data to a file whose name is prefixed by filename_prefix.
void WriteSolutionToFile(const std::string& filename_prefix,
                         const VNEmbedding* embedding) {
  printf("Cost = %ld\n", embedding->cost);
  if (embedding->cost == INF) return;
  // Write node maps to file.
  FILE* nmap_file = fopen((filename_prefix + ".nmap").c_str(), "w");
  auto& primary_node_map = (*(embedding->node_map))[PRIMARY];
  for (int i = 0; i < primary_node_map.size(); ++i) {
    fprintf(nmap_file, "Virtual node %d --> physical node %d\n", i,
            primary_node_map[i]);
  }
  fclose(nmap_file);

  // Write shadow nmaps to file.
  FILE* snmap_file = fopen((filename_prefix + ".snmap").c_str(), "w");
  auto& backup_node_map = (*(embedding->node_map))[BACKUP];
  for (int i = 0; i < backup_node_map.size(); ++i) {
    fprintf(snmap_file, "Shadow virtual node of %d --> physical node %d\n", i,
            backup_node_map[i]);
  }
  fclose(snmap_file);

  // Write edge maps to file.
  FILE* emap_file = fopen((filename_prefix + ".emap").c_str(), "w");
  auto& emap = embedding->primary_edge_map;
  for (auto emap_it = emap->begin(); emap_it != emap->end(); ++emap_it) {
    auto& vlink = emap_it->first;
    int channel_index = emap_it->second.first;
    auto& plinks = emap_it->second.second;
    for (auto& e : plinks) {
      fprintf(emap_file,
              "Virtual edge (%d, %d) --> physical edge (%d, %d), channel: %d\n",
              vlink.first, vlink.second, e.first, e.second, channel_index);
    }
  }
  fclose(emap_file);

  // Write shadow edge maps to file.
  FILE* semap_file = fopen((filename_prefix + ".semap").c_str(), "w");
  auto& semap = embedding->backup_edge_map;
  for (auto semap_it = semap->begin(); semap_it != semap->end(); ++semap_it) {
    auto& vlink = semap_it->first;
    int channel_index = semap_it->second.first;
    auto& plinks = semap_it->second.second;
    for (auto& e : plinks) {
      fprintf(semap_file,
              "Shadow virtual edge of (%d, %d) --> physical edge (%d, %d), "
              "channel: %d\n",
              vlink.first, vlink.second, e.first, e.second, channel_index);
    }
  }
  fclose(semap_file);

  // Write cost to file.
  FILE* cost_file = fopen((filename_prefix + ".cost").c_str(), "w");
  fprintf(cost_file, "%ld\n", embedding->cost);
  fclose(cost_file);
}

// This function takes a physical topology, a virtual topology and
// location constraints and returns a pointer to a VNEmbedding
// object. This function creates the all possible initial seed
// mappings and spawns worker threads to compute a whole embedding
// based on the seed node mapping. This function finally aggregates
// results from the worker threads and returns the best obtained
// embedding.
std::unique_ptr<VNEmbedding> ProtectedVNE(
    const Graph* phys_topology, const Graph* virt_topology,
    const std::vector<std::vector<int>>& location_constraints) {
  int min_constraint_size = INF;
  std::vector<int> most_constrained_vnodes;
  for (int i = 0; i < location_constraints.size(); ++i) {
    if (location_constraints[i].size() < min_constraint_size) {
      min_constraint_size = location_constraints[i].size();
    }
  }
  for (int i = 0; i < virt_topology->node_count(); ++i) {
    if (min_constraint_size == location_constraints[i].size()) {
      most_constrained_vnodes.emplace_back(i);
    }
  }
  // int n_threads = most_constrained_vnodes.size() *
  //                ((min_constraint_size * (min_constraint_size - 1)) / 2);
  // Spawn a thread for each of the virtual node that is considered first for
  // mapping.
  int n_threads = 0;
  for (auto& constraint : location_constraints) {
    n_threads += (constraint.size() * (constraint.size() - 1)) / 2;
  }
  std::vector<pthread_t> threads(n_threads);
  int thread_id = 0;
  std::vector<std::unique_ptr<ThreadParameter>> parameters(n_threads);
  for (auto& parameter : parameters) {
    parameter = std::unique_ptr<ThreadParameter>(new ThreadParameter());
  }
  int num_cores = sysconf(_SC_NPROCESSORS_ONLN);
  int current_core = 0;
  // for (auto& most_constrained_vnode : most_constrained_vnodes) {

  // Consider each virtual node as a first node to map and consider all pair
  // of primary, backup mapping for this virtual node in consideration.
  for (int most_constrained_vnode = 0;
       most_constrained_vnode < virt_topology->node_count();
       ++most_constrained_vnode) {
    for (int i = 0; i < location_constraints[most_constrained_vnode].size();
         ++i) {
      for (int j = i + 1;
           j < location_constraints[most_constrained_vnode].size(); ++j) {
        int primary = location_constraints[most_constrained_vnode][i];
        int backup = location_constraints[most_constrained_vnode][j];
        ThreadParameter* parameter = parameters[thread_id].get();
        parameter->vnode_seed = most_constrained_vnode;
        parameter->primary_seed = primary;
        parameter->backup_seed = backup;
        parameter->phys_topology = phys_topology;
        parameter->virt_topology = virt_topology;
        parameter->location_constraints = &location_constraints;
        cpu_set_t cpu_set;
        CPU_ZERO(&cpu_set);
        CPU_SET(current_core, &cpu_set);
        // Create a worker thread and assign it a CPU core in a round robin
        // basis.
        pthread_create(&threads[thread_id], NULL, &EmbedVNThread, parameter);
        pthread_setaffinity_np(threads[thread_id], sizeof(cpu_set), &cpu_set);
        current_core = (current_core + 1) % num_cores;
        ++thread_id;
      }
    }
  }
  std::vector<std::unique_ptr<VNEmbedding>> embeddings(n_threads);
  for (int i = 0; i < n_threads; ++i) {
    void* ret_value;
    pthread_join(threads[i], &ret_value);
    embeddings[i] = std::unique_ptr<VNEmbedding>(
        std::move(reinterpret_cast<VNEmbedding*>(ret_value)));
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
  auto embedding = ProtectedVNE(physical_topology.get(), virt_topology.get(),
                                *location_constraints);
  auto solution_end_time = std::chrono::high_resolution_clock::now();
  unsigned long long elapsed_time =
      std::chrono::duration_cast<std::chrono::nanoseconds>(
          solution_end_time - solution_start_time).count();
  printf("Solution time: %llu.%llus\n", elapsed_time / 1000000000LL,
         elapsed_time % 1000000000LL);
  WriteSolutionToFile(vn_topology_filename, embedding.get());
  return 0;
}
