#ifndef IO_H_
#define IO_H_

#include "datastructure.h"
#include "util.h"

#include <map>
#include <memory>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <string.h>

std::unique_ptr<std::map<std::string, std::string> > ParseArgs(int argc,
                                                               char *argv[]) {
  std::unique_ptr<std::map<std::string, std::string> > arg_map(
      new std::map<std::string, std::string>());
  for (int i = 1; i < argc; ++i) {
    char *key = strtok(argv[i], "=");
    char *value = strtok(NULL, "=");
    DEBUG(" [%s] => [%s]\n", key, value);
    arg_map->insert(std::make_pair(key, value));
  }
  return std::move(arg_map);
}

std::unique_ptr<std::vector<std::vector<std::string> > > ReadCSVFile(
    const char *filename) {
  DEBUG("[Parsing %s]\n", filename);
  FILE *file_ptr = fopen(filename, "r");
  const static int kBufferSize = 1024;
  char line_buffer[kBufferSize];
  std::unique_ptr<std::vector<std::vector<std::string> > > ret_vector(
      new std::vector<std::vector<std::string> >());
  std::vector<std::string> current_line;
  int row_number = 0;
  while (fgets(line_buffer, kBufferSize, file_ptr)) {
    DEBUG("Read %d characters\n", strlen(line_buffer));
    if (strlen(line_buffer) <= 0) continue;
    if (line_buffer[0] == '\n' || line_buffer[0] == '\r') continue;
    current_line.clear();
    char *token = strtok(line_buffer, ",\n\r");
    current_line.push_back(token);
    while ((token = strtok(NULL, ",\n"))) {
      current_line.push_back(token);
    }
    ret_vector->push_back(current_line);
  }
  fclose(file_ptr);
  DEBUG("Parsed %d lines\n", static_cast<int>(ret_vector->size()));
  return std::move(ret_vector);
}

std::unique_ptr<Graph> InitializeTopologyFromFile(const char *filename) {
  int node_count = 0, edge_count = 0;
  auto csv_vector = ReadCSVFile(filename);
  std::unique_ptr<Graph> graph(new Graph());
  for (int i = 0; i < csv_vector->size(); ++i) {
    auto &row = csv_vector->at(i);

    // Each line has the following format:
    // LinkID, SourceID, DestinationID, PeerID, Cost, Bandwidth, Delay.
    int u = atoi(row[1].c_str());
    int v = atoi(row[2].c_str());
    int cost = atoi(row[4].c_str());
    long bw = atol(row[5].c_str());
    int delay = atoi(row[6].c_str());

    DEBUG("Line[%d]: u = %d, v = %d, cost = %d, bw = %ld, delay = %d\n", i, u,
          v, cost, bw, delay);
    graph->add_edge(u, v, bw, delay, cost);
  }
  return std::move(graph);
}

std::unique_ptr<std::vector<std::vector<int> > > InitializeVNLocationsFromFile(
    const char *filename, int num_virtual_nodes) {
  DEBUG("Parsing %s\n", filename);
  auto ret_vector = std::unique_ptr<std::vector<std::vector<int> > >(
      new std::vector<std::vector<int> >(num_virtual_nodes));
  auto csv_vector = ReadCSVFile(filename);
  for (int i = 0; i < csv_vector->size(); ++i) {
    auto &row = csv_vector->at(i);
    int vnode_id = atoi(row[0].c_str());
    for (int j = 1; j < row.size(); ++j) {
      ret_vector->at(vnode_id).push_back(atoi(row[j].c_str()));
    }
  }
  return std::move(ret_vector);
}

#endif  // IO_H_
