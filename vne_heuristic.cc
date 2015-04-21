#include "datastructure.h"
#include "heuristic.h"
#include "io.h"
#include "util.h"
#include "vne_heuristic_solver.h"

#include <chrono>

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
  auto node_map = CreateInitialNodeMap(
      physical_topology.get(), *location_constraints.get());
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
    printf("Error: node mapping failed\n");
    return 1;
  }

  auto& primary_node_map = (*node_map)[PRIMARY];
  auto& backup_node_map = (*node_map)[BACKUP];
  auto partitions = PartitionGraph(physical_topology.get(), node_map.get());
  auto& primary_partition = (*partitions)[PRIMARY];
  auto& backup_partition = (*partitions)[BACKUP];
 
  bool edge_map_failed = false;
  auto emap = EmbedVN(physical_topology.get(), virt_topology.get(), primary_partition, primary_node_map);
  for (auto emap_it = emap->begin(); emap_it != emap->end(); ++emap_it) {
    if (emap_it->second.empty()) {
      edge_map_failed = true;
      break;
    }
  }
  if (edge_map_failed) {
    printf("Error: regular edge mapping failed\n");
    return 1;
  }
  auto semap = EmbedVN(physical_topology.get(), virt_topology.get(), backup_partition, backup_node_map);
  for (auto semap_it = semap->begin(); semap_it != semap->end(); ++semap_it) {
    if (semap_it->second.empty()) {
      edge_map_failed = true;
      break;
    }
  }
  if (edge_map_failed) {
    printf("Error: shadow edge map failed\n");
    return 1;
  }

  auto solution_end_time = std::chrono::high_resolution_clock::now();
  unsigned long long elapsed_time = std::chrono::duration_cast<std::chrono::nanoseconds>(
                solution_end_time - solution_start_time).count();
  printf("Solution time: %llu.%llus\n", elapsed_time / 1000000000LL,
                 elapsed_time % 1000000000LL);
  // Write node maps to file.
  FILE* nmap_file = fopen((pn_topology_filename + ".nmap").c_str(), "w");
  for (int i = 0; i < primary_node_map.size(); ++i) {
    fprintf(nmap_file, "Virtual node %d --> physical node %d\n", i, primary_node_map[i]);
  }
  fclose(nmap_file);

  // Write shadow nmaps to file.
  FILE* snmap_file = fopen((pn_topology_filename + ".snmap").c_str(), "w");
  for (int i = 0; i < backup_node_map.size(); ++i) {
    fprintf(snmap_file, "Shadow virtual node of %d --> physical node %d\n", i, backup_node_map[i]);
  }
  fclose(snmap_file);
  
  int cost = 0;
  // Write edge maps to file.
  FILE* emap_file = fopen((pn_topology_filename + ".emap").c_str(), "w");
  for (auto emap_it = emap->begin(); emap_it != emap->end(); ++emap_it) {
    auto& vlink = emap_it->first;
    auto& plinks = emap_it->second;
    for (auto& e : plinks) {
      fprintf(emap_file, "Virtual edge (%d, %d) --> physical edge (%d, %d)\n", vlink.first, vlink.second, e.first, e.second);
      cost += physical_topology->get_edge_cost(e.first, e.second) * virt_topology->get_edge_bandwidth(vlink.first, vlink.second);
    }
  }
  fclose(emap_file);

  // Write shadow edge maps to file.
  FILE* semap_file = fopen((pn_topology_filename + ".semap").c_str(), "w");
  for (auto semap_it = semap->begin(); semap_it != semap->end(); ++semap_it) {
    auto& vlink = semap_it->first;
    auto& plinks = semap_it->second;
    for (auto& e : plinks) {
      fprintf(emap_file, "Shadow virtual edge of (%d, %d) --> physical edge (%d, %d)\n", vlink.first, vlink.second, e.first, e.second);
      cost += physical_topology->get_edge_cost(e.first, e.second) * virt_topology->get_edge_bandwidth(vlink.first, vlink.second);
    }
  }
  fclose(semap_file);

  // Write cost to file.
  FILE* cost_file = fopen((pn_topology_filename + ".cost").c_str(), "w");
  fprintf(cost_file, "%d\n", cost);
  fclose(cost_file);

  return 0;
}
