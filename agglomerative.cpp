//
//  agglomerative.cpp
//  agglomerative
//
//  Created by Bahadır on 8.01.2018.
//  Copyright © 2018 Bahadır. All rights reserved.
//

#include <stdio.h>
#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>
#include <string>
#include <limits.h>
#include "strutils.h"
using namespace std;


// For the euclidian distance calculation
// Create an array of these then calculate the distances between them.
struct point
{
    int id;
    double rtt1, rtt2, rtt3, rtt4, rtt5, q1, q2, q3, q4, q5;
    
    point()
    :id(0), rtt1(0), rtt2(0), rtt3(0), rtt4(0), rtt5(0), q1(0), q2(0), q3(0), q4(0), q5(0)
    {}
    
    point(int id_no, double num1, double num2, double num3, double num4, double num5, double num6,double num7, double num8, double num9, double num10)
    :id(id_no), rtt1(num1), rtt2(num2), rtt3(num3), rtt4(num4), rtt5(num5), q1(num6), q2(num7), q3(num8), q4(num9), q5(num10)
    {}
    
};


// After distance matrix is constructed, all we need for clustering is the id of the point
struct simplePoint {
    int id;
    simplePoint* right;
    
    simplePoint(int i): id(i), right(NULL) {}
};


// Cluster nodes
struct clusterNode {
    int clusterID;
    clusterNode* down;
    clusterNode* up;
    simplePoint* right;
    clusterNode(int id): clusterID(id), down(NULL), up(NULL), right(NULL) {}
    
};


vector<double> separatedByComma(string line){
    string numStr = "";
    vector<double> vec(11);
    int index = 0;
    for (int i = 0; i < line.length(); i++){
        if (line[i] != ',')
            numStr += line[i];
        else{
            vec[index++] = atof(numStr);
            numStr = "";
        }
        
        if (i == line.length()-1){
            vec[index++] = atof(numStr);
        }
        
    }
    
    return vec;
}


// Calculates Euclidian distance between two point
double euclid_dist(const point &p1, const point &p2){
    double distance = 0;
    distance += pow(abs(p1.rtt1-p2.rtt1), 2);
    distance += pow(abs(p1.rtt2-p2.rtt2), 2);
    distance += pow(abs(p1.rtt3-p2.rtt3), 2);
    distance += pow(abs(p1.rtt4-p2.rtt4), 2);
    distance += pow(abs(p1.rtt5-p2.rtt5), 2);
    distance += pow(abs(p1.q1-p2.q1), 2);
    distance += pow(abs(p1.q2-p2.q2), 2);
    distance += pow(abs(p1.q3-p2.q3), 2);
    distance += pow(abs(p1.q4-p2.q4), 2);
    distance += pow(abs(p1.q5-p2.q5), 2);
    
    return pow(distance, 0.5);
}


// Alternative distance formula
double alter_dist(const point &p1, const point &p2){
    
    double up = abs(p1.rtt1-p2.rtt1) + abs(p1.rtt2-p2.rtt2) + abs(p1.rtt3-p2.rtt3) + abs(p1.rtt4-p2.rtt4) + abs(p1.rtt5-p2.rtt5) + abs(p1.q1-p2.q1) + abs(p1.q2-p2.q2) + abs(p1.q3-p2.q3) + abs(p1.q4-p2.q4) + abs(p1.q5-p2.q5);
    
    double down = abs(p1.rtt1+p2.rtt1) + abs(p1.rtt2+p2.rtt2) + abs(p1.rtt3+p2.rtt3) + abs(p1.rtt4+p2.rtt4) + abs(p1.rtt5+p2.rtt5) + abs(p1.q1+p2.q1) + abs(p1.q2+p2.q2) + abs(p1.q3+p2.q3) + abs(p1.q4+p2.q4) + abs(p1.q5+p2.q5);
    
    double d = up / down;
    return d;
}


// Enter i,j get the index in the distance .
long int index(int i, int j){
    long int sum;
    sum = (19999*19998)/2 - (19999-i)*(19998-i)/2;
    return sum + j;
}


// create "size * clusterNode" with one simplePoint. simplePoint id's 0,1, .... , size-1;
clusterNode* constructor(int size){
    clusterNode* head;
    head = new clusterNode(0);
    head->right = new simplePoint(0);
    
    clusterNode* temp = head;
    
    for(int i = 1; i < size; i++){
        temp->down = new clusterNode(i);
        temp->down->up = temp;
        temp = temp->down;
        temp->right = new simplePoint(i);
    }
    return head;
}


// Free the memory
void destructor(clusterNode* &head){
    clusterNode* temp = head;
    while(head != NULL){
        
        simplePoint * p = head->right;
        head->right = NULL;
        simplePoint * tp = p;
        
        while(p != NULL){
            p = p->right;
            delete tp;
            tp = p;
        }
        
        head = head->down;
        delete temp;
        temp = head;
        
    }
}


// Calculate the complete(farthest) distance between two cluster node.
// Use a lookup table. Update the table each time of clustering.
double completeDistance(const clusterNode* c1, const clusterNode* c2, const double* dm){
    simplePoint* p1 = c1->right;
    simplePoint* p2 = c2->right;
    double max = __DBL_MIN__;
    
    while (p1 != NULL) {
        while (p2 != NULL) {
            double temp = dm[index(p1->id, p2->id)];
            if (temp > max)
                max = temp;
            p2 = p2->right;
        }
        p1 = p1->right;
    }
    return max;
}


// Calculate the single(closest) distance between two cluster node.
// Use a lookup table. Update the table each time of clustering.
double singleDistance(const clusterNode* c1, const clusterNode* c2, const double* dm){
    simplePoint* p1 = c1->right;
    simplePoint* p2 = c2->right;
    double min =__DBL_MAX__;
    
    while (p1 != NULL){
        while (p2 != NULL){
            double temp = dm[index(p1->id, p2->id)];
            if(temp < min)
                min = temp;
            p2 = p2->right;
        }
        p1 = p1->right;
    }
    return min;
}

// input: two cluster pointer
// output: merge two cluster pointer into one cluster
void merge(clusterNode* head, clusterNode* c1, clusterNode* c2, double* cm, const double* dm, const bool &t){
    // Merge c1 and c2
    simplePoint* p1 = c1->right;
    while(p1->right != NULL)
        p1 = p1->right;
    p1->right = c2->right;
    c2->right = NULL;
    
    c2->up->down = c2->down;
    if(c2->down != NULL)
        c2->down->up = c2->up;
    delete c2;
    
    // Update clusterMatrix.
    clusterNode* temp = head;
    if(temp->clusterID == c1->clusterID)
        temp = temp->down;
    
    // start to compare merged node to other nodes
    while (temp != NULL) {
        double d;
        if(t)
            d = completeDistance(c1, temp, dm);
        else
            d = singleDistance(c1,temp,dm); // Write
        
        cm[index((int)min(c1->clusterID,temp->clusterID), (int)max(c1->clusterID,temp->clusterID))] = d;
        
        temp = temp->down;
        if(temp != NULL && temp->clusterID == c1->clusterID)
            temp = temp->down;
    }
}



// Start to cluster
void agglomerative(clusterNode* &head, const double* dm, double* cm, int c_s, const bool &t){
    //int counter = 1;
    for (int i = 0; i < 20000-c_s; i++){
        clusterNode* first = head;
        clusterNode* second = head->down;
        
        // check distances between 'CLUSTERS'
        double min = __DBL_MAX__;
        
        clusterNode* c1 = NULL;
        clusterNode* c2 = NULL;
        // put flag(point to) to clusterNodes which are going to be merged.
        
        // Compare distances between clusters
        while(first->down != NULL){
            while(second != NULL){
                // Get distance from the clusterDistance - prev: completeDistance(first, second, dm);
                double d = cm[index(first->clusterID, second->clusterID)];
                
                if (d < min){
                    min = d;
                    //cout << "new min: " << min << endl;
                    c1 = first;
                    c2 = second;
                }
                second = second->down;
            }
            first = first->down;
            second = first->down;
        }
        //cout << counter++ << " - link " << c1 << " and " << c2 << endl;
        merge(head, c1, c2, cm, dm, t);
    }
}


//int arr[5] = {30,25,10,10,15};
// Find the server with min time
int minServer(point p){
    int arr[5];
    arr[0] = p.rtt1 + p.q1 + 30;
    arr[1] = p.rtt2 + p.q2 + 25;
    arr[2] = p.rtt3 + p.q3 + 10;
    arr[3] = p.rtt4 + p.q4 + 10;
    arr[4] = p.rtt5 + p.q5 + 15;
    
    // find max in array of five
    
    int index = 0;
    double min = arr[0];
    for(int i = 1; i < 5; i++){
        if(arr[i]<min){
            min = arr[i];
            index = i;
        }
    }
    
    return index;
}


// Get the max value in a vector
int getMaxInt(vector<int> &v){
    int max = v[0];
    for(int i=1;i<v.size();i++){
        if(v[i]>max)
            max=v[i];
    }
    return max;
}


// Get the min value in a vector
int getMinInt(vector<int> &v){
    int min = v[0];
    for(int i=1;i<v.size();i++){
        if(v[i]<min)
            min=v[i];
    }
    return min;
}


// Labels the clusters according to their time spent.
// vec --> pointVec in the main
vector<int> label_clusters(clusterNode* head, const vector<point> &vec){
    vector<int> labels(5);
    int c = 0;
    int i;
    while(head != NULL){
        vector<int> labelCounter(5,0); // C1, C2, C3, C4, C5
        simplePoint* p = head->right;
        while(p != NULL){
            i = minServer(vec[p->id]);
            labelCounter[i]++;
            p = p->right;
        }
        // Find index
        int currentLabel = 0;
        int m = getMaxInt(labelCounter);
        for(int k = 0; k < 5; k++){
            cout << k << ":" << labelCounter[k] << " | ";
            if(labelCounter[k] == m)
                currentLabel = k;
        }
        cout << endl;
        labels[c++] = currentLabel;
        head = head->down;
    }
    
    return labels;
}

double getMin(const vector<double> &v){
    double min = __DBL_MAX__;
    for(int i = 0; i < v.size(); i++){
        if(v[i] < min)
            min = v[i];
    }
    return min;
}

double getMax(const vector<double> &v){
    double max = __DBL_MIN__;
    for(int i = 0; i < v.size(); i++){
        if(v[i] < max)
            max = v[i];
    }
    return max;
}

void preprocess(vector<double> &vec){
    for(int i = 0; i < vec.size(); i++){
        double max = getMax(vec);
        double min = getMin(vec);
        vec[i] = (vec[i] - min) / (max - min);
    }
}


///////////////////////////////////////////////////////////////

int main(){
    ifstream file;
    file.open("RTT_Queue.csv");
    
    vector<point> pointVec(20000);
    double* distanceMatrix = new double[10000*20000];
    double* clusterMatrix= new double[10000*20000];
    
    string line;
    getline(file, line);
    
    
    // PREPROCESSING - OPTIONAL
    bool is_preprocess = true;
    
    while(getline(file, line)){
        vector<double> commaSepVec = separatedByComma(line);
        
        //Preprocessing
        if(is_preprocess){
            preprocess(commaSepVec);
        }
        
        point p(commaSepVec[0], commaSepVec[1], commaSepVec[2], commaSepVec[3], commaSepVec[4], commaSepVec[5], commaSepVec[6], commaSepVec[7], commaSepVec[8], commaSepVec[9], commaSepVec[10]);

        pointVec[p.id] = p;
    }
    
    // 1 --> EUCLID
    // 0 --> ALTERNATE
    bool DISTANCE = 0;
    
    // Create distance matrix.
    for (int i = 0; i < 20000; i++){
        for (int j = i+1; j < 20000; j++){
            if(DISTANCE){
                distanceMatrix[index(i, j)] = euclid_dist(pointVec[i], pointVec[j]);
                clusterMatrix[index(i, j)] = euclid_dist(pointVec[i], pointVec[j]);
            }
            else{
                distanceMatrix[index(i, j)] = alter_dist(pointVec[i], pointVec[j]);
                clusterMatrix[index(i, j)] = alter_dist(pointVec[i], pointVec[j]);
            }
        }
    }
    
    cout << "Distance matrix is ready!" << endl;
    
    
    
    // Construct the clusterNodes with initial nodes.
    
    clusterNode* head = constructor(20000);
    
    // 1 --> COMPLETE
    // 0 --> SINGLE
    agglomerative(head, distanceMatrix, clusterMatrix, 5, 0);
    
    vector<int> labels = label_clusters(head, pointVec);
    
    for(int i = 0; i < 5; i++){
        printf("%d", labels[i]);
    }
    
    
    
    destructor(head);
    delete[] distanceMatrix;
    delete [] clusterMatrix;
    return 0;
}

