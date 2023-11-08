#include <iostream>
#include <fstream>
#include <sstream>
#include <unordered_map>
#include <set>
#include <vector>
#include <algorithm>
#include <cstdio>
#include <stdlib.h>
#include <math.h>
#include <random>
#include <time.h>
#include "mt19937ar.h"

#include "cryptopp/integer.h"
#include "cryptopp/osrng.h"
#include "cryptopp/nbtheory.h"
#include "cryptopp/modarith.h"
#include "include/stats.hpp"

using namespace std;
using namespace CryptoPP;
// Initialization of statslib

stats::rand_engine_t engine(1776);

int get_NodeNum(string filename){
    set<int>node;
    ifstream file;
    file.open(filename, ios::in);
    string line;
    if (file.is_open())
    {
        cout<<"success open"<<endl;
        while (getline(file, line)) {
            if (line[0] == '#') {
                continue;
            }
            istringstream iss(line);
            int firstNode, secNode;
            int firstNode_map,secNode_map;
            iss >> firstNode >> secNode;
            node.insert(firstNode);
            node.insert(secNode);
        }
    }
    file.close();
    return node.size();
}

void ReadFile(string filename, map<int,int> & deg, vector<vector<int>>& graph, int & d_max, int n)
{
    graph.resize(n, vector<int>(n));
    for(int i=0;i<n;i++) {
        deg[i]=0;
        for (int j = 0; j < n; j++) {
            graph[i][j] = 0;
        }
    }
    map<int, int> nodeID;
    int flag=0;
    ifstream file;
    file.open(filename, ios::in);
    string line;
    if (file.is_open())
    {
        cout<<"success open"<<endl;
        while (getline(file, line)) {
            if (line[0] == '#') {
                continue;
            }
            istringstream iss(line);
            int firstNode, secNode;
            int firstNode_map,secNode_map;
            iss >> firstNode >> secNode;
            if(nodeID.find(firstNode)==nodeID.end()){
                nodeID[firstNode]=flag;
                flag++;
            }
            if(nodeID.find(secNode)==nodeID.end()){
                nodeID[secNode]=flag;
                flag++;
            }
            firstNode_map=nodeID[firstNode];
            secNode_map=nodeID[secNode];

            // store graph by adjacent list, undirected graph
            if(graph[firstNode_map][secNode_map]==0){
                graph[firstNode_map][secNode_map]=graph[secNode_map][firstNode_map]=1;
                deg[firstNode_map]++;
                deg[secNode_map]++;
            }
        }
    }else
        cout << "Open Failure!" << endl;
    d_max=0;
    for (const auto &pair : deg) {
        if (pair.second > d_max) {
            d_max = pair.second;
        }
    }
    file.close();
}

FILE *FileOpen(string filename, const char *mode) {
    FILE *fp;

    if ((fp = fopen(filename.c_str(), mode)) == NULL) {
        cout << "cannot open " << filename << endl;
        exit(-1);
    }
    return fp;
}

int TrueTriangle(vector<vector<int>> graph,int num_user){
    int sum=0;
    int temp1,temp2,temp3;
    for(int i=0;i<num_user;i++){
        for(int j=i+1;j<num_user;j++) {
            for (int k = j + 1; k < num_user; k++) {
                if (graph[i][j] * graph[i][k] * graph[j][k] == 1) {
                    sum++;
                    //cout<<i<<" "<<j<<" "<<k<<endl;
                }
            }
        }
    }
    cout<<endl;
    return sum;
}

bool cmp(const pair<int, double> a, const pair<int, double> b) {
    return a.second>b.second;
}

void Local_Project(vector<vector<int>> & graph, int parameter ,map<int,int> deg, int n) {
    for(int i=0;i<n;i++){
        if(deg[i]>parameter){
            int num_delete = deg[i] - parameter;
            //map<int,int> deg_similarity;
            vector< std::pair<int, double> > deg_similarity;
            int deg_ID=deg[i];
            for(int j=0;j<n;j++){
                if(graph[i][j]==1){
                    int deg_j=deg[j];
                    double score= abs(deg_ID-deg_j)/double(deg_ID);
                    deg_similarity.push_back(make_pair(j,score));
                }
            }
            sort(deg_similarity.begin(),deg_similarity.end(),cmp);
            for(int k=0;k<num_delete;k++){
                int delete_ID=deg_similarity[k].first;
                graph[i][delete_ID]=0;
            }
        }
    }
}

int getRand(int min, int max) {
    return ( rand() % (max - min + 1) ) + min ;
}

void SecureTriangle(vector<vector<int>>  graph, int num_user, double &T1, double &T2){
    double x,y,z,o,p,q,w;
    double x1,y1,z1,o1,p1,q1,w1;
    double x2,y2,z2,o2,p2,q2,w2;
    double e,f,g;
    double e1,f1,g1;
    double e2,f2,g2;
    double u1,u2;
    int a,b,c;
    double a1,b1,c1;
    double a2,b2,c2;
    for(int i=0;i<num_user;i++)
        for(int j=i+1;j<num_user;j++)
            for(int k=j+1;k<num_user;k++) {
                a=graph[i][j];
                b=graph[i][k];
                c=graph[j][k];

                //cout<<"data: "<<i<<" "<<j<<" "<<k<<endl;
                //cout<<"temp: "<<a<<" "<<b<<" "<<c<<endl;
                a1 = getRand(1, num_user);
                a2 = a - a1;
                b1 = getRand(1, num_user);
                b2 = b - b1;
                c1 = getRand(1, num_user);
                c2 = c - c1;
                x = getRand(1, num_user);
                y = getRand(1, num_user);
                z = getRand(1, num_user);
                w = x * y * z;
                o = x * y;
                p = x * z;
                q = y * z;
                x1 = getRand(1, num_user);
                x2 = x - x1;
                y1 = getRand(1, num_user);
                y2 = y - y1;
                z1 = getRand(1, num_user);
                z2 = z - z1;
                w1 = getRand(1, num_user);
                w2 = w - w1;
                o1 = getRand(1, num_user);
                o2 = o - o1;
                p1 = getRand(1, num_user);
                p2 = p - p1;
                q1 = getRand(1, num_user);
                q2 = q - q1;

                //Server S1
                e1 = a1 - x1;
                f1 = b1 - y1;
                g1 = c1 - z1;
                //Server S2
                e2 = a2 - x2;
                f2 = b2 - y2;
                g2 = c2 - z2;
                //exchange and compute
                e = e1 + e2;
                f = f1 + f2;
                g = g1 + g2;
                //server S1
                u1 = w1 + o1 * g + p1 * f + q1 * e + x1 * f * g + y1 * e * g + z1 * e * f;
                //server S2
                u2 = w2 + o2 * g + p2 * f + q2 * e + x2 * f * g + y2 * e * g + z2 * e * f + e * f * g;
                T1 += u1;
                T2 += u2;
            }
}

double gamma_noise(double eps, int num_user, double d_ma_prime){
    double k=1.0/num_user;
    double theta=d_ma_prime/eps;
    double lap_noise=0.0;

    std::random_device rd;
    std::mt19937 gen(rd());

    double gamma_1 = stats::rgamma(k,theta,engine);
    double gamma_2 = stats::rgamma(k,theta,engine);
    lap_noise=lap_noise+gamma_1-gamma_2;

    return lap_noise;
}

double Noisy_Triangle(double T1, double T2, double eps, int num_user, double d_max_prime){
    double noise_1=0.0;
    double noise_2=0.0;
    double noise;
    for(int i=0;i<num_user;i++){
        noise=gamma_noise(eps, num_user, d_max_prime);
        double temp_1=getRand(0, num_user);
        double temp_2=noise-temp_1;
        noise_1+=temp_1;
        noise_2+=temp_2;
    }
    return T1+T2+noise_1+noise_2;
}

void Noisydegree(map<int,int> deg, double & d_max_prime, double eps){
    double max=0.0;
    for (const auto &pair: deg) {
        double noisy_deg=pair.second+stats::rlaplace(0.0, 2.0/eps);
        if(noisy_deg>max)
            max=noisy_deg;
    }
    d_max_prime=max;
}

double Compute_l2loss(int trueTriangle, int noisyTriangle){
    double l2loss=pow(trueTriangle-noisyTriangle,2);
    return l2loss;
}

double Compute_RE(int trueTriangle, int noisyTriangle){
    double re=abs(trueTriangle-noisyTriangle)/double(trueTriangle);
    return re;
}

void printG(map<int,set<int>> a_list,int num_user){
    for(int i=0;i<num_user;i++){
        cout<<"-----: "<<i<<" :"<<endl;
        for(const auto& neighbor_ID : a_list[i]){
            cout<<neighbor_ID<<" ";
        }
        cout<<endl;
    }

}
int main(int argc, char *argv[]) {
    map<int,int> deg;
    map<int,set<int>> a_list;
    map<int,int> local_triangle;

    double d_max_prime;
    int d_max=0;
    int max_ID=0;
    int parameter=0;
    int true_triangle=0;
    int pro_triangle=0;
    double noisy_triangle=0;
    double eps=0;
    double alpha=0.1;
    double eps_deg=alpha*eps;
    double eps_tri=(1-alpha)*eps;
    int num_user;
    double l2loss=0.0; //l2 loss
    double re=0.0;    // relative error
    double l2loss_pro=0.0; //l2 loss
    double re_pro=0.0;    // relative error
    int num_iter=30;
    double T1,T2;
    clock_t start, end;
    vector<vector<int>> graph;

    double epsilons[6]={0.5,1,1.5,2,2.5,3};
    string filename="wiki-Vote.txt";
    cout << "******************" << endl;
    cout << "dataset: " << filename << endl;
    num_user = get_NodeNum(filename);
    ReadFile(filename, deg, graph, d_max, num_user);

    cout << "number of nodes: " << num_user << endl;
    cout<<"true maximum degree: "<<d_max<<endl;

    true_triangle = TrueTriangle(graph, num_user);
    cout << "number of triangles: " << true_triangle << endl;

    for (int i = 0; i < 6; i++) {
        map<int, set<int>> pro_a_list;
        map<int, int> pro_local_triangle;
        map<int, int> noisy_local_triangle;

        eps = epsilons[i];
        cout << "-------------------" << endl;
        cout << "epsilon: " << eps << endl;
        alpha = 0.1;
        eps_deg = alpha * eps;
        eps_tri = (1 - alpha) * eps;
        double temp_time;
        double sum_d_max_prime = 0;

        double sum_time = 0.0;
        double sum_noisy_triangle = 0.0;
        start = clock();
        for (int j = 0; j < num_iter; j++) {
            Noisydegree(deg, d_max_prime, eps_deg);
            sum_d_max_prime = sum_d_max_prime + d_max_prime;
        }
        end = clock();
        sum_time += double(end - start) / CLOCKS_PER_SEC;

        start = clock();
        d_max_prime = round(sum_d_max_prime / num_iter);
        //pro_a_list.insert(a_list.begin(),a_list.end());
        parameter = d_max_prime;
        cout << "d_max_prime: " << d_max_prime << endl;
        Local_Project(graph, parameter, deg, num_user);
        pro_triangle = TrueTriangle(graph, num_user);
        cout << "pro_triangle: " << pro_triangle << endl;

        l2loss_pro = Compute_l2loss(true_triangle, pro_triangle);
        re_pro = Compute_RE(true_triangle, pro_triangle);
        cout << "l2loss of projection: " << re_pro << " " << " re of projection: " << re << endl;
        end = clock();
        temp_time = double(end - start) / CLOCKS_PER_SEC;
        sum_time = sum_time + temp_time;

        // securely compute the secret shares of triangles
        start = clock();
        T1 = T2 = 0;
        SecureTriangle(graph, num_user, T1, T2);

        end = clock();
        temp_time = double(end - start) / CLOCKS_PER_SEC;
        cout << "time of secure computation: " << temp_time << endl;
        sum_time = sum_time + temp_time;
        //cout<<"test: "<<T1+T2<<endl;

        //perturb using distributed noise
        start = clock();
        for (int j = 0; j < num_iter; j++) {
            noisy_triangle = Noisy_Triangle(T1, T2, eps, num_user, d_max_prime);
            sum_noisy_triangle += noisy_triangle;
        }
        end = clock();
        temp_time = double(end - start) / CLOCKS_PER_SEC;
        temp_time = temp_time / num_iter;
        sum_time = sum_time + temp_time;

        noisy_triangle = round(sum_noisy_triangle / num_iter);
        cout << "noisy_triangle: " << noisy_triangle << endl;

        l2loss = Compute_l2loss(true_triangle, noisy_triangle);
        re = Compute_RE(true_triangle, noisy_triangle);
        cout << "l2loss: " << l2loss << " " << " re: " << re << endl;
        cout << "running time: " << sum_time << "s" << endl;
    return 0;
}
