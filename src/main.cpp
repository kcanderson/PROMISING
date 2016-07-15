#include "coreroutines.h"
#include "CompleteGraphScorer.h"
#include "PValueModuleScorer.h"
#include <tclap/CmdLine.h>

int main(int argc, char** argv) {
  try {

    // Command-line parsing
    TCLAP::CmdLine cmd("Command description message", ' ', "0.9");
    TCLAP::ValueArg<std::string> netFilename("n", "network", "Network file", true, "", "string");
    cmd.add(netFilename);
    TCLAP::ValueArg<std::string> groupsFilename("g", "groups", "Groups file", true, "", "string");
    cmd.add(groupsFilename);
    TCLAP::ValueArg<int> pvalIterations("p", "pval", "Pvalue iterations", false, -1, "int");
    cmd.add(pvalIterations);
    TCLAP::ValueArg<std::string> outFilename("o", "outfile", "Output summary file", false, "", "string");
    cmd.add(outFilename);
    // Read in command-line options.
    cmd.parse(argc, argv);

    // Read in groups from file
    std::string gfilename = groupsFilename.getValue();
    std::ifstream ginfile(gfilename);
    TGroups groups;
    readGroups(ginfile, groups);

    // Flatten out groups
    std::vector<std::string> entries;
    flattenGroups<std::string>(groups, entries);
    const int numEntriesInGroups = entries.size();

    // Do we calculate a p-value?
    int pIterations = pvalIterations.getValue();
    
    // Get network filename
    std::string nfilename = netFilename.getValue();
    std::ifstream ninfile(nfilename);
    std::string line;
    
    // Grab network names
    std::getline(ninfile, line);
    TIndexMap fullMap, map;
    int numNodes = parseNamesLine(line, fullMap);

    
    float* mat(0);
    int matrixWidth(0);
    IModuleScorer* moduleScorer(0);
    
    if (pIterations == -1) {
      // Just score the nodes. No p-value calculation.
      moduleScorer = new CompleteGraphScorer4();
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
    else if (pIterations > 0) {
      // We need to calculate empirical p-values.
      // PValueModulesScorer will take ownership of the complete graph scorer.
      moduleScorer = new PValueModuleScorer(pIterations, new CompleteGraphScorer4());
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
    if (mat == 0) {
      std::cerr << "Problem reading matrix from file." << std::endl;
      exit(-1);
    }
    if (moduleScorer == 0) {
      std::cerr << "ModuleScorer is nil. Something went wrong." << std::endl;
      exit(-1);
    }

    TIndicesGroups igroups;
    mapGroupsToIndices<std::string>(groups, map, igroups);
        
    // Score genes based on strong modules
    std::cout << "Scoring genes." << std::endl;
    TScoreMap scores;
    moduleScorer->ScoreModule(mat, matrixWidth, igroups, scores);

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

    // Delete heap stuffs.
    free(mat);
    delete moduleScorer;
    
  } catch (TCLAP::ArgException &e) {
    std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;
  }

}
