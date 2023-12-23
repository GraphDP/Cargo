#include <iostream>
#include <fstream>
#include <sstream>
#include <unordered_map>
#include <set>
#include <vector>
#include <map>
#include <algorithm>
#include <cstdio>
#include <stdlib.h>
#include <math.h>
#include <random>
#include <time.h>

using namespace std;

FILE *FileOpen(string filename, const char *mode) {
    FILE *fp;

    if ((fp = fopen(filename.c_str(), mode)) == NULL) {
        cout << "cannot open " << filename << endl;
        exit(-1);
    }
    return fp;
}

void Sample_Node(string filename, set<int>&N, int num_n){
    set<int>Nodes;
    ifstream file;
    file.open(filename, ios::in);
    string line;
    srand((unsigned) time(NULL));
    if (file.is_open()) {
        cout << "success open" << endl;
        while (getline(file, line)) {
            if (line[0] == '#')
                continue;
            istringstream iss(line);
            int firstNode, secNode;
            iss >> firstNode >> secNode;
            Nodes.insert(firstNode);
            Nodes.insert(secNode);
        }
    }
    vector<int>vec_node(Nodes.begin(),Nodes.end());
    std::random_device rd;
    std::mt19937 g(rd());
    shuffle(vec_node.begin(),vec_node.end(),g);

    cout<<"Node size: "<<vec_node.size()<<endl;
    cout<<"Sampled node size: "<<num_n<<endl;
    for(int i=0;i<num_n;i++){
        N.insert(vec_node[i]);
    }
    file.close();
}
void ReadFile_Generate(string filename, string file1, set<int>Nodes)
{
    ifstream file;
    file.open(filename, ios::in);
    string line;
    srand((unsigned) time(NULL));
    if (file.is_open()) {
        cout << "success open" << endl;
        while (getline(file, line)) {
            if (line[0] == '#')
                continue;
            istringstream iss(line);
            int firstNode, secNode;
            iss >> firstNode >> secNode;

            if (Nodes.count(firstNode)>0) {
                if (Nodes.count(secNode)>0) {
                    FILE *fp;
                    fp = FileOpen(file1, "a+");
                    fprintf(fp, "%d %d\n", firstNode, secNode);
                    fclose(fp);
                }
            }
        }
    }
    else{
        cout << "Open Failure!" << endl;
    }
    file.close();
}

int main(int argc, char *argv[]) {

    string filename="/data/facebook.txt";
    string file1="/data/facebook2000.txt";

    set<int>Nodes;
    int num_n=2000;
    Sample_Node(filename, Nodes,num_n);
    std::ofstream file_1(file1, std::ios::out | std::ios::trunc);
    if(!file_1) {
        // Handle the error, e.g., the file might not exist or other reasons.
        return 1;
    }
    file_1.close();

    cout<<"******************"<<endl;
    cout<<"dataset: "<<filename<<endl;
    ReadFile_Generate(filename, file1, Nodes);
    return 0;
}