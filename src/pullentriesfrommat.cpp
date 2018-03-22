#include <stdio.h>
#include <cstring>
#include <string>
#include <map>
#include <utility>
#include <sstream>
#include <vector>
#include <tclap/CmdLine.h>

#include "coreroutines.h"

#define LINE_SIZE 2048
typedef std::map< std::pair< std::string, std::string >, float> NetMap;
typedef std::map<std::string, int> IndexMap;

bool readnetwork(FILE* fp, NetMap& network, IndexMap& indices) {
  char line[LINE_SIZE];
  int i = 0;
  
  while (fgets(line, LINE_SIZE, fp) != 0) {
    char* a = strtok(line, "\t");
    char* b = strtok(0, "\t");
    char* c = strtok(0, "\t");
    float score = atof(c);
    //printf("%s, %s, %f\n", a, b, score);
    network[std::pair<std::string, std::string>(a, b)] = score;
    if (indices.find(a) == indices.end()) {
      indices[a] = i++;
    }
    if (indices.find(b) == indices.end()) {
      indices[b] = i++;
    }
  }
}


int main(int argc, char** argv)
{
  
  TCLAP::CmdLine cmd("Pull out edges from mat file.", ' ', "0.9");
  TCLAP::ValueArg<std::string> matrixFilename("m", "matrix", "matrix file", true, "", "string");
  cmd.add(matrixFilename);
  TCLAP::ValueArg<std::string> edgesFilename("e", "edges", "wanted edes file", true, "", "string");
  cmd.add(edgesFilename);
  TCLAP::ValueArg<std::string> outFilename("o", "outfile", "Output matrix file", true, "", "string");
  cmd.add(outFilename);
  
  cmd.parse(argc, argv);

  std::string mfilename = matrixFilename.getValue();
  std::ifstream minfile(mfilename);
  std::string line;
  std::getline(minfile, line);
  TIndexMap fullMap, map;
  int numNodes = parseNamesLine(line, fullMap);

  std::string efilename = edgesFilename.getValue();
  FILE* f = fopen(efilename.c_str(), "r");
  if (f == 0) {
    printf("Problem reading network from file %s\n", efilename.c_str());
    return(-1);
  }
  NetMap net;
  IndexMap indices;
  readnetwork(f, net, indices);
  fclose(f);

  // Allocate memory for subnework
  int matrixWidth(indices.size());
  float* mat = 0;
  mat = (float*)malloc(sizeof(float) * matrixWidth * matrixWidth);
  if (mat == 0) {
    printf("Could not allocate matrix buffer of %d bytes.\n",
  	   (int)(matrixWidth * matrixWidth * sizeof(float)));
    exit(-1);
  }

  std::vector<std::string> entries;
  for (const auto &item : indices) {
    entries.push_back(item.first);
  }
  
  readEntriesFromNetwork(minfile, fullMap, entries, map, mat, matrixWidth);

  for (const auto &e : net) {
    const auto edge = e.first;
    const int i1 = map[edge.first];
    const int i2 = map[edge.second];
    const float val = mat[i1 * matrixWidth + i2];
    printf("%s\t%s\t%e\n", edge.first.c_str(), edge.second.c_str(), val);
  }
  
  return 0;
}
