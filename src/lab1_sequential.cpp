extern "C" {
    #include "lab1_sequential.h"
}

#include <vector>
#include <tuple>
#include <stdlib.h>
#include <time.h>
#include <limits>
#include <cmath>
#include <iostream>

void kmeans_sequential(int N, int K, int* data_points, int** data_point_cluster, 
int** centroids, int* num_iterations) {
    std::vector<int>* centroidSaveDB = new std::vector<int>();
    std::vector<double> currentCentroidsDouble(3*K);
    
    // Initialize centroids
    srand(time(0));
    int numRandomCentInit = 0;
    while(numRandomCentInit<K) {
        int pickIndex = rand()%N;
        int pickedX = *(data_points+(3*pickIndex));
        int pickedY = *(data_points+(3*pickIndex)+1);
        int pickedZ = *(data_points+(3*pickIndex)+2);
        
        // check no centroid is repeated
        bool matched = false;
        for (int i=0; i<numRandomCentInit; i++) {
            int iThCurrCentX = currentCentroidsDouble[(3*i)];
            int iThCurrCentY = currentCentroidsDouble[(3*i)+1];
            int iThCurrCentZ = currentCentroidsDouble[(3*i)+2];
            if (pickedX==iThCurrCentX & pickedY==iThCurrCentY & pickedZ==iThCurrCentZ) {
                matched = true;
                break;
            }
        }

        // pickable centroid found
        if (!matched) {
            currentCentroidsDouble[3*numRandomCentInit] = pickedX;
            currentCentroidsDouble[3*numRandomCentInit+1] = pickedY;
            currentCentroidsDouble[3*numRandomCentInit+2] = pickedZ;
            numRandomCentInit++;
        }
    }

    // successfully pick random points as centroids (Forgy)
    centroidSaveDB->insert(centroidSaveDB->end(), currentCentroidsDouble.begin(), currentCentroidsDouble.end());

    // update centroids
    *num_iterations = -1;
    double error = 100000.0;

    std::vector<int> partitionEntries(N);
    std::vector<int> numOfPointsInPartition(K);
    std::vector<double> prevCentroidsDouble(3*K);
    bool zeroThIteration = true;

    while (error>1.0) {
        if (!zeroThIteration) {
            // compute centroids
            for (int i=0; i<K; i++) {
                // initialize memory of next set of centroids this will contain sum of points of partition temporarily
                currentCentroidsDouble[(3*i)] = 0;
                currentCentroidsDouble[(3*i)+1] = 0;
                currentCentroidsDouble[(3*i)+2] = 0;
                numOfPointsInPartition[i] = 0;
            }
            
            // prepare sums and nums for avg
            for (int i=0; i<N; i++) {
                int partitionId = partitionEntries[i];
                int pickedX = *(data_points+(3*i));
                int pickedY = *(data_points+(3*i)+1);
                int pickedZ = *(data_points+(3*i)+2);
                currentCentroidsDouble[(3*partitionId)] += pickedX;
                currentCentroidsDouble[(3*partitionId)+1] += pickedY;
                currentCentroidsDouble[(3*partitionId)+2] += pickedZ;
                numOfPointsInPartition[partitionId] += 1;
            }

            // calculate avgs i.e. centroids
            for (int i=0; i<K; i++) {
                int iThCurrCentX = currentCentroidsDouble[(3*i)];
                int iThCurrCentY = currentCentroidsDouble[(3*i)+1];
                int iThCurrCentZ = currentCentroidsDouble[(3*i)+2];
                double denominator = double(numOfPointsInPartition[i]);
                currentCentroidsDouble[(3*i)] = round(double(iThCurrCentX)/denominator);
                currentCentroidsDouble[(3*i)+1] = round(double(iThCurrCentY)/denominator);
                currentCentroidsDouble[(3*i)+2] = round(double(iThCurrCentZ)/denominator);
            }

            // complete update of centroids
            centroidSaveDB->insert(centroidSaveDB->end(), currentCentroidsDouble.begin(), currentCentroidsDouble.end());
        }

        // assign partition
        for (int i=0; i<N; i++) {
            double minDistance = std::numeric_limits<double>::max();
            int pickedX = *(data_points+(3*i));
            int pickedY = *(data_points+(3*i)+1);
            int pickedZ = *(data_points+(3*i)+2);
            for (int j=0; j<K; j++) {
                double jThCurrCentX = currentCentroidsDouble[(3*j)];
                double jThCurrCentY = currentCentroidsDouble[(3*j)+1];
                double jThCurrCentZ = currentCentroidsDouble[(3*j)+2];
                double currCentDistance = pow(pickedX-jThCurrCentX,2)+pow(pickedY-jThCurrCentY,2)+pow(pickedZ-jThCurrCentZ,2);
                if (currCentDistance<minDistance) {
                    minDistance = currCentDistance;
                    partitionEntries[i] = j;
                }
            }
        }
        
        // error computation
        if (!zeroThIteration) {
            error = 0.0;
            for (int i=0; i<K; i++) {
                double iThCurrCentX = currentCentroidsDouble[(3*i)];
                double iThCurrCentY = currentCentroidsDouble[(3*i)+1];
                double iThCurrCentZ = currentCentroidsDouble[(3*i)+2];
                double iThPrevCentX = prevCentroidsDouble[(3*i)];
                double iThPrevCentY = prevCentroidsDouble[(3*i)+1];
                double iThPrevCentZ = prevCentroidsDouble[(3*i)+2];
                error += pow(iThCurrCentX-iThPrevCentX, 2) + pow(iThCurrCentY-iThPrevCentY, 2) + pow(iThCurrCentZ-iThPrevCentZ, 2);
            }
        } else {
            zeroThIteration = false;
        }
        // copy current centroid, they are prev for next iteration
        prevCentroidsDouble = currentCentroidsDouble;
        *num_iterations += 1;
    }

    int* clusterPartIdList = (int*) malloc(sizeof(int)*N*4);
    for (int i=0; i<N; i++) {
        int partitionId = partitionEntries[i];
        int pickedX = *(data_points+(3*i));
        int pickedY = *(data_points+(3*i)+1);
        int pickedZ = *(data_points+(3*i)+2);
        *(clusterPartIdList+(4*i)) = pickedX;
        *(clusterPartIdList+(4*i+1)) = pickedY;
        *(clusterPartIdList+(4*i)+2) = pickedZ;
        *(clusterPartIdList+(4*i)+3) = partitionId;
    }
    *centroids = centroidSaveDB->data();
    *data_point_cluster = clusterPartIdList;
}