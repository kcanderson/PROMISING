#include "coreroutines.h"
#include "CompleteGraphScorer.h"
#include "FastScorer.h"
#include "PValueModuleScorer.h"
#include <tclap/CmdLine.h>
#include <cstring>

typedef std::map<std::string, TGroups> TDiseases;

bool parseMultipleDiseaseFile(const std::string& filename, TDiseases& diseaseEntries) {
  std::ifstream instream(filename);
  std::string line;
  while (!instream.eof()) {
    TGroups groups;
    std::getline(instream, line);
    std::vector<std::string> split;
    std::stringstream stream(line);
    std::string entry;
    while (std::getline(stream, entry, '\t')) {
      split.push_back(entry);
    }

    for (int i = 1; i < split.size(); ++i) {
      std::stringstream ss(split[i]);
      std::vector<std::string> items;
      while (std::getline(ss, entry, ',')) {
	items.push_back(entry);
      }
      char locusName[8];
      snprintf(locusName, 8, "%d", i);
      //groups.push_back(items);
      groups[locusName] = items;
    }
    if (groups.size() > 0) diseaseEntries[split[0]] = groups;
  }
  return true;
}

// bool parseMultipleDiseaseFile(const std::string& filename, TDiseases& diseaseEntries) {
//   std::ifstream instream(filename);
//   std::string line;
//   TDiseases entries;
//   while (!instream.eof()) {
//     std::getline(instream, line);
//     std::stringstream stream(line);
//     std::string entry;
//     std::string name;
//     stream >> name;
//     std::vector<std::string> genes;
//     TGroups groups;
//     //while (!stream.eof()) {
//     //std::stringstream ss(
//     while (std::getline(stream, entry) {
//       stream >> entry;
//       //entry.erase(entry.find_last_not_of(" \n\r\t")+1);
//       char* dup = strdup(entry.c_str());
//       char* pch = strtok(dup, ",");
//       std::vector<std::string> group;
//       while (pch != 0) {
// 	group.push_back(pch);
// 	pch = strtok(0, ",");
//       }
//       groups.push_back(group);
//     }
//     if (groups.size() > 0) {
//       diseaseEntries[name] = groups;
//     }
//   }
//   return true;
// }

bool readDegreeGroups(const std::string& filename, TIndexMap& geneToIndexMap, std::map<int, int>& nodeDegreeGroup) {
  std::ifstream dinfile(filename);
  TGroups groups;
  
  if (false == readGMT(dinfile, groups)) {
    printf("Could not read degree group file.");
    return false;
  }
  int i = 0;
  for (auto const& g : groups) {
    for (auto const& node : g.second) {
      nodeDegreeGroup[geneToIndexMap[node]] = i;
    }
    i++;
  }

  dinfile.close();
  return true;
}

int main(int argc, char** argv) {
  try {

    // Command-line parsing
    TCLAP::CmdLine cmd("Prioritization of candidate genes in disjoint sets.", ' ', "0.9");
    TCLAP::ValueArg<std::string> netFilename("m", "matrix", "Adjacency matrix", true, "", "string");
    cmd.add(netFilename);
    TCLAP::ValueArg<std::string> groupsFilename("g", "groups", "Groups file", true, "", "string");
    cmd.add(groupsFilename);
    TCLAP::ValueArg<int> pvalIterations("p", "pval", "Pvalue iterations", false, -1, "int");
    cmd.add(pvalIterations);
    TCLAP::ValueArg<std::string> outFilename("o", "outfile", "Output summary file", false, "", "string");
    cmd.add(outFilename);
    TCLAP::ValueArg<int> cGroupSize("s", "size", "Complete graph group size", false, 3, "int");
    cmd.add(cGroupSize);

    TCLAP::ValueArg<std::string> method("z", "method", "Scoring method, complete or fast", false, "complete", "string");
    cmd.add(method);

    TCLAP::ValueArg<std::string> degree("d", "degree_groups", "Degree groups for node permutations", false, "", "string");
    cmd.add(degree);

    TCLAP::SwitchArg multiple("l", "multiple", "Do multiple runs", false);
    cmd.add(multiple);
    
    // Read in command-line options.
    cmd.parse(argc, argv);

    //// Read in groups from file
    // std::map<int, int> nodeDegreeGroups;
    // std::string gfilename = groupsFilename.getValue();
    // std::ifstream ginfile(gfilename);
    // TGroups groups;
    // readGroups(ginfile, groups);

    // // Flatten out groups
    // std::vector<std::string> entries;
    // flattenGroups<std::string>(groups, entries);
    // const int numEntriesInGroups = entries.size();

    int cgraphSize = cGroupSize.getValue();
    if (cgraphSize != 3 && cgraphSize != 4) {
      printf("Inappropriate complete graph size, %d. Only accept 3 or 4.", cgraphSize);
      exit(-1);
    }

    // Do we calculate a p-value?
    int pIterations = pvalIterations.getValue();

    // Get network filename
    std::string nfilename = netFilename.getValue();
    std::ifstream ninfile(nfilename);
    std::string line;

    // // Grab degree groups
    // std::map<int, int> nodeDegreeGroups;
    // std::string dFileName = degree.getValue();
    // if (dFileName == "") {
    //   for (auto const& i : fullMap) {
    // 	nodeDegreeGroups[i.second] = 0;
    //   }
    // } else {
    //   if (false == readDegreeGroups(dFileName, fullMap, nodeDegreeGroups)) {
    // 	printf("Problem reading node degree groups.");
    // 	return(-1);
    //   }
    // }
    
    // Grab network names
    std::getline(ninfile, line);
    TIndexMap fullMap, map;
    int numNodes = parseNamesLine(line, fullMap);

    float* mat(0);
    IModuleScorer* moduleScorer(0);
    int matrixWidth(0);

    // Make module scorer
    if (cgraphSize == 4) {
      moduleScorer = new CompleteGraphScorer4();
    }
    else {
      if (method.getValue() == "complete") {
	moduleScorer = new CompleteGraphScorer();
      } else {
	moduleScorer = new FastScorer();
      }
    }
    
    if (!multiple.getValue()) {

      // Read in groups from file
      std::string gfilename = groupsFilename.getValue();
      std::ifstream ginfile(gfilename);
      TGroups groups;
      readGMT(ginfile, groups);

      // Flatten out groups
      std::vector<std::string> entries;
      flattenGroups<std::string>(groups, entries);
      const int numEntriesInGroups = entries.size();

      // Grab degree groups
      std::map<int, int> nodeDegreeGroups;
      std::string dFileName = degree.getValue();
      if (dFileName == "") {
	for (auto const& i : fullMap) {
	  nodeDegreeGroups[i.second] = 0;
	}
      } else {
	if (false == readDegreeGroups(dFileName, fullMap, nodeDegreeGroups)) {
	  printf("Problem reading node degree groups.");
	  return(-1);
	}
      }
      
      if (pIterations == -1) {
	// Just score the nodes. No p-value calculation.
	//moduleScorer = new CompleteGraphScorer4();
	matrixWidth = numEntriesInGroups;
      
	// Allocate memory for subnework
	mat = (float*)malloc(sizeof(float) * matrixWidth * matrixWidth);
	if (mat == 0) {
	  printf("Could not allocate matrix buffer of %d bytes.\n",
		 (int)(numEntriesInGroups * numEntriesInGroups * sizeof(float)));
	  exit(-1);
	}

	std::cout << "Pulling out necessary subnetwork from file." << std::endl;
	// Pull in values from full matrix
	readEntriesFromNetwork(ninfile, fullMap, entries, map, mat, matrixWidth);
      }
      else {
	// We need to calculate empirical p-values.
	// PValueModulesScorer will take ownership of the complete graph scorer.
	moduleScorer = new PValueModuleScorer(pIterations, moduleScorer, nodeDegreeGroups);
	if (moduleScorer == 0) {
	  std::cerr << "ModuleScorer is nil. Something went wrong." << std::endl;
	  exit(-1);
	}
      
	matrixWidth = numNodes;
	map = fullMap;
	std::cout << "Reading entire network into memory. This may take a while." << std::endl;

	// Allocate a big block of memory. This could easily fail.
	mat = (float*)malloc(sizeof(float) * matrixWidth * matrixWidth);
	if (mat == 0) {
	  std::cerr << "Could not allocate matrix buffer of ";
	  std::cerr << (int)(matrixWidth * matrixWidth * sizeof(float)) << " bytes." << std::endl;
	  exit(-1);
	}

	// Read the network from the stream
	readEntireNetwork(ninfile, numNodes, mat);
	std::cout << "Completed reading network of " << sizeof(float) * matrixWidth * matrixWidth / (1024*1024*1024.0) << "GB." << std::endl;
      
      }

      TIndicesGroups igroups;
      mapGroupsToIndices<std::string>(groups, map, igroups);
        
      // Score genes based on strong modules
      std::cout << "Scoring genes." << std::endl;
      TScoreMap scores;
      TIndices inds;
      for (auto const& g : igroups) {
	for (auto i : g.second) {
	  inds.push_back(i);
	}
      }
    
      moduleScorer->ScoreModule(mat, matrixWidth, igroups, inds, scores);

      TReverseIndexMap rmap2;
      for (auto const &p : map) {
	rmap2[p.second] = p.first;
      }
      moduleScorer->BriefSummary(scores, rmap2, std::cout);

      // Write summary to file.
      std::string outFile = outFilename.getValue();
      if (outFile != "") {
	std::cout << std::endl << "Writing summary to output file " << outFile << std::endl;
	std::ofstream outStream(outFile, std::ofstream::out);
	moduleScorer->LongSummary(scores, rmap2, igroups, outStream);
      }
    }
    else {
      // Deal with multiple thingies.
      std::string gfilename = groupsFilename.getValue();
      TDiseases diseases;
      if (parseMultipleDiseaseFile(gfilename, diseases)) {
	// Allocate a big block of memory. This could easily fail.
	matrixWidth = numNodes;
	map = fullMap;
	std::cout << "Reading entire network into memory. This may take a while." << std::endl;
	mat = (float*)malloc(sizeof(float) * matrixWidth * matrixWidth);
	if (mat == 0) {
	  std::cerr << "Could not allocate matrix buffer of ";
	  std::cerr << (int)(matrixWidth * matrixWidth * sizeof(float)) << " bytes." << std::endl;
	  exit(-1);
	}
	// Read the network from the stream
	readEntireNetwork(ninfile, numNodes, mat);
	std::cout << "Completed reading network of " << sizeof(float) * matrixWidth * matrixWidth / (1024*1024*1024.0) << "GB." << std::endl;

	std::map<int, int> nodeDegreeGroups;
	if (pIterations > 0) {
	  // We need to calculate empirical p-values.
	  // PValueModulesScorer will take ownership of the complete graph scorer.
	  std::string dFileName = degree.getValue();
	  if (dFileName == "") {
	    for (auto const& i : fullMap) {
	      nodeDegreeGroups[i.second] = 0;
	    }
	  } else {
	    if (false == readDegreeGroups(dFileName, fullMap, nodeDegreeGroups)) {
	      printf("Problem reading node degree groups.");
	      return(-1);
	    }
	  }
	  moduleScorer = new PValueModuleScorer(pIterations, moduleScorer, nodeDegreeGroups);
	  if (moduleScorer == 0) {
	    std::cerr << "ModuleScorer is nil. Something went wrong." << std::endl;
	    exit(-1);
	  }
	}

	std::string outFile = outFilename.getValue();
	std::ofstream* outStream = 0;
	if (outFile != "") {
	  outStream = new std::ofstream(outFile, std::ofstream::out);
	}
	for (const auto &disease : diseases) {
	  TIndicesGroups igroups;
	  mapGroupsToIndices<std::string>(disease.second, map, igroups);
	  TIndices inds;
	  for (auto const& g : igroups) {
	    for (auto i : g.second) {
	      inds.push_back(i);
	    }
	  }
	  TScoreMap scores;
      	  moduleScorer->ScoreModule(mat, matrixWidth, igroups, inds, scores);
	  
	  TReverseIndexMap rmap2;
	  for (auto const &p : map) {
	    rmap2[p.second] = p.first;
	  }
	  printf("%s\n", disease.first.c_str());
	  moduleScorer->BriefSummary(scores, rmap2, std::cout);
	  printf("\n\n");
	  if (outFile != "") {
	    *outStream << "-------------------" << std::endl;
	    *outStream << disease.first << std::endl;
	    *outStream << "-------------------" << std::endl;
	    moduleScorer->LongSummary(scores, rmap2, igroups, *outStream);
	  }
	}
	if (outStream) delete outStream;
      }
      else {
	printf("Problem reading entries from file.\n");
      }
    }
    

    // Delete heap stuffs.
    free(mat);
    delete moduleScorer;
      
  } catch (TCLAP::ArgException &e) {
    std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;
  }

}
