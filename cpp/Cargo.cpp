#include <iostream>
#include <fstream>
#include <sstream>
#include <unordered_map>
#include <set>
#include <vector>
#include <algorithm>
#include<cstdio>
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

AutoSeededRandomPool rng;

FILE *FileOpen(string filename, const char *mode) {
    FILE *fp;

    if ((fp = fopen(filename.c_str(), mode)) == NULL) {
        cout << "cannot open " << filename << endl;
        exit(-1);
    }
    return fp;
}

Integer getRandomPrime(int bits) {
    PrimeAndGenerator prime(1,rng,bits);
    return prime.Prime();
}

Integer getRandomInRange(int a, Integer p) {
    Integer random;
    do {
        random.Randomize(rng, a, p.ConvertToLong() - 1);
    } while (random.IsZero());
    return random;
}

Integer powm(Integer a, Integer b, Integer mod) {
    ModularArithmetic ma(mod);
    return ma.Exponentiate(a, b);
}

Integer invmod(Integer a, Integer p) {
    ModularArithmetic ma(p);
    return ma.MultiplicativeInverse(a);
}

int psica(vector<int> adj_1, vector<int> adj_2)
{
    set<int>s_1(adj_1.begin(), adj_1.end());
    set<int>s_2(adj_2.begin(),adj_2.end());

    int bits = 32;
    Integer p = getRandomPrime(bits);  //modulus
    Integer g = getRandomInRange(3, p);  //Radix
    // user_i's pairkeys
    Integer x1 = getRandomInRange(3, p);
    Integer Y1 = powm(g, x1, p);
    //user_j's pairkeys
    Integer x2 = getRandomInRange(3, p);
    Integer Y2 = powm(g, x2, p);

    Integer r1 = getRandomInRange(3, p);
    Integer r2 = getRandomInRange(3, p);
    Integer c11 = powm(g, r1, p);
    Integer c12 = powm(g, r2, p);

    Integer temp1=powm(Y1, r1, p);
    Integer temp2=powm(Y2, r2, p);
    Integer temp3=powm(c11, x1, p);
    Integer temp4=powm(Y2, r2, p);
    Integer temp5=invmod(temp3, p);

    int num_psica=0;

    for(set<int>::iterator it_1 = s_1.begin();it_1!=s_1.end();it_1++) {
        for (set<int>::iterator it_2 = s_2.begin(); it_2 != s_2.end(); it_2++) {
            // encode message1
            Integer message_1 = ( temp1* *it_1) % p;  // user_1
            Integer message_12 = (temp2 * message_1) % p; //user_2
            Integer message_2 = (message_12 * temp5) % p; //user_1 removes user_1's mask

            //encode message2
            Integer message_22 = ( temp4 * *it_2) % p;  // user_2
            int result = message_2.Compare(message_22); // common or not

            if (result == 0) {
                num_psica++;
                continue;
            }
        }
    }
    return num_psica;
}

void ReadFile(string filename, map<int,int> & deg, map<int,vector<int>> & a_list, int & d_max, int & max_ID)
{
    ifstream file;
    file.open(filename, ios::in);
    string line;
    max_ID=0;
    if (file.is_open())
    {
        cout<<"success open"<<endl;
        while (getline(file, line)) {
            if (line[0] == '#') {
                continue;
            }
            istringstream iss(line);
            int firstNode, secNode;
            iss >> firstNode >> secNode;
            // store graph by adjacent list, undirected graph
            if (max_ID<firstNode){
                max_ID=firstNode;
            }
            else if(max_ID<secNode){
                max_ID=secNode;
            }
            if (a_list.find(firstNode) != a_list.end()) {
                a_list[firstNode].push_back(secNode);
                deg[firstNode]++;
            } else {
                a_list[firstNode] = vector<int>();
                a_list[firstNode].push_back(secNode);
                deg[firstNode] = 1;
            }
            if (a_list.find(secNode) != a_list.end()) {  //directed --> undirected
                a_list[secNode].push_back(firstNode);
                deg[secNode]++;
            } else {
                a_list[secNode] = vector<int>();
                a_list[secNode].push_back(firstNode);
                deg[secNode] = 1;
            }
        }
    }else{
        cout << "Open Failure!" << endl;
    }

    d_max=0;
    for (const auto &pair : deg) {
        if (pair.second > d_max) {
            d_max = pair.second;
        }
    }
}

int TrueTriangle(map<int,vector<int>> a_list, map<int,int> & local_triangle, string filename)
{
    string outdir;
    string outfilename;
    outdir=filename.substr(0,3);
    outfilename=outdir+"_trueTriangle.txt";
    FILE *fp;

    fp = FileOpen(outfilename, "w+");
    int num_triangle=0;
    for(const auto &pair : a_list){
        int ID_1=pair.first;
        vector<int> adj_1=pair.second;
        int num_ID_1=0;

        for(vector<int>::iterator it = adj_1.begin();it!=adj_1.end();it++){
            int ID_2=*it;
            if(ID_1<ID_2){
                vector<int> adj_2=a_list[ID_2];
                int temp= psica(adj_1,adj_2);
                num_ID_1=num_ID_1+temp;
            }
        }
        local_triangle[ID_1]=num_ID_1;
        fprintf(fp,"%d %d\n",ID_1,local_triangle[ID_1]);
        num_triangle=num_triangle+local_triangle[ID_1];
    }
    fclose(fp);
    return num_triangle/3;
}

int ProTriangle(map<int,vector<int>> pro_a_list,map<int,vector<int>> a_list, map<int,int> & pro_local_triangle, string filename,map<int,int> deg, double d_max_prime)
{
    string outdir;
    string outfilename;
    outdir=filename.substr(0,3);
    outfilename=outdir+"_proTriangle.txt";
    FILE *fp;

    fp = FileOpen(outfilename, "w+");
    int num_triangle=0;
    for(const auto &pair : pro_a_list){
        int ID_1=pair.first;
        vector<int> adj_1;
        if (deg[ID_1]<=d_max_prime)
            adj_1=a_list[ID_1];
        else
            adj_1=pair.second;
        int num_ID_1=0;
        for(vector<int>::iterator it = adj_1.begin();it!=adj_1.end();it++){
            int ID_2=*it;
            if(ID_1<ID_2){
                vector<int> adj_2;
                if (deg[ID_2]<=d_max_prime)
                    adj_2=a_list[ID_2];
                else
                    adj_2=pro_a_list[ID_2];
                int temp= psica(adj_1,adj_2);
                num_ID_1=num_ID_1+temp;
            }
        }
        pro_local_triangle[ID_1]=num_ID_1;

        fprintf(fp,"%d %d\n",ID_1,pro_local_triangle[ID_1]);
        num_triangle=num_triangle+pro_local_triangle[ID_1];
    }
    fclose(fp);
    return num_triangle/3;
}

bool cmp(const pair<int, double> a, const pair<int, double> b) {
    return a.second>b.second;
}

void LocalProject (map<int,vector<int>> & pro_a_list, int parameter,int max_ID,map<int,int> deg, string filename){
    int flag = 1;
    for (const auto &pair: pro_a_list) {
        int ID = pair.first;
        if (pro_a_list[ID].size() < parameter) {
            int num_add = parameter - pro_a_list[ID].size();
            pro_a_list[ID].insert(pro_a_list[ID].end(), num_add, max_ID + flag);
            flag++;
        } else {
            int num_delete = pro_a_list[ID].size() - parameter;
            vector< std::pair<int, double> > deg_similarity;
            int deg_ID=deg[ID];

            for(int i=0;i<pro_a_list[ID].size();i++){
                int neighbor_ID=pro_a_list[ID][i];
                int deg_i=deg[neighbor_ID];
                double score= abs(deg_ID-deg_i)/double(deg_ID);
                deg_similarity.push_back(make_pair(neighbor_ID,score));
            }
            sort(deg_similarity.begin(),deg_similarity.end(),cmp);
            for(int j=0;j<num_delete;j++){
                int delete_ID=deg_similarity[j].first;
                std::vector<int>::iterator pos;
                pos = find(pro_a_list[ID].begin(), pro_a_list[ID].end(), delete_ID);
                if (pos != pro_a_list[ID].end()){
                    pro_a_list[ID].erase(pos);
                }
            }
        }
    }
}

void RandomProject(map<int,vector<int>> & pro_a_list, int parameter,int max_ID,map<int,int> deg, string filename) {
    int flag = 1;
    for (const auto &pair: pro_a_list) {
        int ID = pair.first;
        if (pro_a_list[ID].size() < parameter) {
            int num_add = parameter - pro_a_list[ID].size();
            pro_a_list[ID].insert(pro_a_list[ID].end(), num_add, max_ID + flag);
            flag++;
        } else {
            int num_delete = pro_a_list[ID].size() - parameter;
            //map<int,int> deg_similarity;
            vector< std::pair<int, double> > deg_similarity;
            int deg_ID=deg[ID];

            for(int j=0;j<num_delete;j++){
                int random_ID=getRand(0, pro_a_list[ID].size()-1);
                std::vector<int>::iterator pos;
                pos = pro_a_list[ID].begin()+random_ID;
                if (pos != pro_a_list[ID].end())
                    pro_a_list[ID].erase(pos);
            }
        }
    }
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

void Local_Perturb(map<int,int> local_triangle, map<int,int> & noisy_local_triangle, double eps, int num_user, double d_max_prime){

    for (const auto &pair: local_triangle) {
        int user_ID=pair.first;
        int true_triangle=pair.second;
        double noise=gamma_noise(eps, num_user, d_max_prime);
        double noisy_triangle=true_triangle+noise;
        if (noisy_triangle<0)
            noisy_triangle=0;
        noisy_local_triangle[user_ID]=round(noisy_triangle);
    }
}

int getRand(int min, int max) {
    return ( rand() % (max - min + 1) ) + min ;
}

void SS_Agg_int(map<int,int> noisy_local_triangle, int & noisy_trianle){

    //secret sharing
    int sum_triangle=0;
    int sum_triangle_1=0;
    int sum_triangle_2=0;

    for (const auto &pair: noisy_local_triangle) {
        int local_triangle=pair.second;
        int split_1=getRand(0, local_triangle);
        int split_2=local_triangle-split_1;
        sum_triangle_1=sum_triangle_1+split_1;
        sum_triangle_2=sum_triangle_2+split_2;
    }
    //Aggregation
    sum_triangle=sum_triangle_1+sum_triangle_2;
    noisy_trianle=sum_triangle/3;
}

void SS_Agg(map<int,int> noisy_local_triangle, int & noisy_trianle){
    vector<Integer>server_1;
    vector<Integer>server_2;
    //secret sharing
    for (const auto &pair: noisy_local_triangle) {
        double local_triangle=pair.second;
        AutoSeededRandomPool rng;
        Integer split_1=rng.GenerateWord32(0, int(local_triangle));
        Integer split_2=local_triangle-split_1;
        server_1.push_back(split_1);
        server_2.push_back(split_2);
    }
    //Aggregation
    Integer sum_triangle=0;
    for(int i=0;i<server_1.size();i++){
        sum_triangle=sum_triangle+server_1[i]+server_2[i];
    }
    noisy_trianle=sum_triangle.ConvertToLong()/3;
}

int test_SS_Agg(map<int,int> noisy_local_triangle){
    int sum_triangle=0;
    for (const auto &pair: noisy_local_triangle) {
        sum_triangle=sum_triangle+pair.second;
    }
    return sum_triangle/3;
}

double Compute_l2loss(int trueTriangle, int noisyTriangle){
    double l2loss=pow(trueTriangle-noisyTriangle,2);
    return l2loss;
}

double Compute_RE(int trueTriangle, int noisyTriangle){
    double re=abs(trueTriangle-noisyTriangle)/double(trueTriangle);
    return re;
}

int main(int argc, char *argv[]) {
    map<int,int> deg;
    map<int,vector<int>> a_list;
    map<int,int> local_triangle;

    double d_max_prime;
    int d_max=0;
    int max_ID=0;
    int parameter=0;
    int true_triangle=0;
    int pro_triangle=0;
    int noisy_triangle=0;
    double eps=0;
    double alpha=0.1;
    double eps_deg=alpha*eps;
    double eps_tri=(1-alpha)*eps;
    int num_user;
    double l2loss=0.0; //l2 loss
    double re=0.0;    // relative error
    int num_iter=1;
    clock_t start,end;

    double epsilons[6]={0.5,1,1.5,2,2.5,3};

    string filename=argv[1];  //input graph dataset
    cout<<"******************"<<endl;
    cout<<"dataset: "<<filename<<endl;
    ReadFile(filename, deg, a_list, d_max,max_ID);
    num_user=deg.size();
    true_triangle=TrueTriangle(a_list,local_triangle,filename);
    cout<<"number of triangles: "<<true_triangle<<endl;

    for(int i=0;i<1;i++){
        map<int,vector<int>> pro_a_list;
        map<int,int> pro_local_triangle;
        map<int,int> noisy_local_triangle;

        eps=epsilons[i];
        cout<<"-------------------"<<endl;
        cout<<"epsilon: "<<eps<<endl;
        alpha=0.1;
        eps_deg=alpha*eps;
        eps_tri=(1-alpha)*eps;
        double sum_time=0.0;
        double sum_noisy_triangle=0.0;
        start=clock();
        for(int j=0;j<num_iter;j++){

            // noisy maximum degree
            Noisydegree(deg, d_max_prime, eps_deg);
            parameter=round(d_max_prime);
            cout<<"d_max_prime: "<<d_max_prime<<endl;

            // local projection
            pro_a_list.insert(a_list.begin(),a_list.end());
            LocalProject(pro_a_list, parameter,max_ID,deg,filename); 

            // local triangle counting
            pro_triangle=ProTriangle(pro_a_list,a_list,pro_local_triangle,filename, deg, d_max_prime);
            cout<<"pro_triangle: "<<pro_triangle<<endl;
            noisy_local_triangle.insert(pro_local_triangle.begin(),pro_local_triangle.end());

            // local perturbation with distributed noise
            Local_Perturb(pro_local_triangle, noisy_local_triangle, eps_tri, num_user, d_max_prime);

            // secret sharing and aggregation
            SS_Agg(noisy_local_triangle, noisy_triangle);
            sum_noisy_triangle=sum_noisy_triangle+noisy_triangle;
        }
        end = clock();
        sum_time+=double(end-start)/CLOCKS_PER_SEC;
        sum_time=sum_time/num_iter;
        noisy_triangle=sum_noisy_triangle/num_iter;

        l2loss= Compute_l2loss(true_triangle,noisy_triangle);
        re= Compute_RE(true_triangle,noisy_triangle);
        cout<<"l2loss: "<<l2loss<<" "<<" re: "<<re<<endl;
        cout<<"running time: "<<sum_time<<"s"<<endl;
    }

    return 0;
}
