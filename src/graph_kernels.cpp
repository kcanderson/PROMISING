#include <stdio.h>
#include <cstring>
#include <string>
#include <map>
#include <utility>
#include <sstream>
#include <vector>
#include <tclap/CmdLine.h>

typedef std::map< std::pair< std::string, std::string >, float> NetMap;
typedef std::map<std::string, int> IndexMap;

extern "C" {
    // LU decomoposition of a general matrix
    void sgetrf_(int* M, int *N, float* A, int* lda, int* IPIV, int* INFO);

    // generate inverse of a matrix given its LU decomposition
    void sgetri_(int* N, float* A, int* lda, int* IPIV, float* WORK, int* lwork, int* INFO);
}

void inverse(float* A, int N)
{
    int *IPIV = new int[N+1];
    int LWORK = N*N;
    float *WORK = new float[LWORK];
    int INFO;

    sgetrf_(&N,&N,A,&N,IPIV,&INFO);
    sgetri_(&N,A,&N,IPIV,WORK,&LWORK,&INFO);

    delete [] IPIV;
    delete [] WORK;
}

void laplacian(float* a, int n) {
  float* ptr(a);
  for (int i = 0; i < n; ++i, ptr += n) {
    float acc = 0.0;
    for (int j = 0; j < n; ++j) {
      acc += ptr[j];
      ptr[j] = -ptr[j];
    }
    ptr[i] = acc;
  }
}

void scalarmultiply(float scalar, float* a, int n) {
  for (int i = 0; i < n*n; ++i) {
    a[i] = a[i] * scalar;
  }
}

void regularizedlaplacian(float alpha, float* a, int n) {
  laplacian(a, n);
  scalarmultiply(alpha, a, n);
  // Add identity
  float* ptr(a);
  for (int i = 0; i < n; ++i, ptr += n) {
    ptr[i] = ptr[i] + 1.0f;
  }
  inverse(a, n);
}

#define LINE_SIZE 2048

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

bool populatematrix(const NetMap& net, const IndexMap& indices, float* const mat, int n) {
  for (auto const &v : net) {
    IndexMap::const_iterator it;
    it = indices.find(v.first.first);
    if (it == indices.end()) {
      printf("Index out of range\n");
    } else {
      const int i1 = it->second;
      it = indices.find(v.first.second);
      if (it == indices.end()) {
	printf("Index out of range\n");
      } else {
	const int i2 = it->second;
	mat[i1 * n + i2] = v.second;
	// Assume undirected: symmetric
	mat[i2 * n + i1] = v.second;
      }
    }
  }
}

void writematrixtofile(const float* mat, int n, FILE* f) {
  
  for (int i = 0; i < n; ++i, mat += n) {
    std::stringstream sstream;
    
    for (int j = 0; j < n; ++j) {
      //if (j != 0) v = v + "\t";
      //v = v + std::to_string(mat[j]);
      if (j != 0) sstream << "\t";
      sstream << mat[j];
    }
    fprintf(f, "%s\n", sstream.str().c_str());
  }
}

int main(int argc, char** argv)
{
  
  TCLAP::CmdLine cmd("Calc graph kernel of network.", ' ', "0.9");
  TCLAP::ValueArg<std::string> netFilename("n", "network", "Network file", true, "", "string");
  cmd.add(netFilename);
  TCLAP::ValueArg<std::string> outFilename("o", "outfile", "Output matrix file", true, "", "string");
  cmd.add(outFilename);
  TCLAP::ValueArg<float> alpha("a", "alpha", "Alpha parameter", false, 0.01f, "float");
  cmd.add(alpha);
  TCLAP::ValueArg<std::string> method("k", "kernel", "Graph kernel", false, "rl", "string");
  cmd.add(method);
  
  cmd.parse(argc, argv);
  
  std::string nfilename = netFilename.getValue();
  
  FILE* f = fopen(nfilename.c_str(), "r");
  if (f == 0) {
    printf("Problem reading network from file %s\n", nfilename.c_str());
    return(-1);
  }
  NetMap net;
  IndexMap indices;
  readnetwork(f, net, indices);
  fclose(f);
  
  IndexMap::const_iterator it = indices.begin();
  // for (; it != indices.end(); it++) {
  //   printf("%s, %i\n", it->first.c_str(), it->second);
  // }

  const int n = indices.size();
  float* mat = (float*)malloc(sizeof(float) * n * n);
  if (mat == 0) {
    printf("Problem allocating adjacency matrix of %i bytes.", sizeof(float) * n * n);
    return(-1);
  }
  // Initialize to zeros.
  for (int i = 0; i < n * n; ++i) {
    mat[i] = 0.0f;
  }
  populatematrix(net, indices, mat, n);


  const std::string m = method.getValue();

  if (m == "amat") {
    printf("Adjacency matrix\n");
  } else {
    printf("Regularized Laplacian\n");
    float a = alpha.getValue();
    regularizedlaplacian(a, mat, n);
  }
  std::string ofilename = outFilename.getValue();
  f = fopen(ofilename.c_str(), "w");
  if (f == 0) {
    printf("Bad output file name: %s\n", ofilename.c_str());
    return(-1);
  }
  
  std::map<int, std::string> names;
  for (auto const &v : indices) {
    names[v.second] = v.first;
  }
  std::stringstream ss;
  for (int i = 0; i < n; ++i) {
    if (i != 0) ss << "\t";
    ss << names[i];
  }
  fprintf(f, "%s\n", ss.str().c_str());
  writematrixtofile(mat, n, f);
  fclose(f);
  // for (int i = 0; i < n; ++i) {
  //   for (int j = 0; j < n; ++j) {
  //     printf("%f ", mat[i*n + j]);
  //   }
  //   printf("\n");
  // }
  
  free(mat);
  return 0;
}
