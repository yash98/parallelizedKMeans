extern "C" {
    #include "lab1_sequential.h"
}

#include <vector>
#include <tuple>
#include <stdlib.h>
#include <time.h>
#include <limits>
#include <cmath>

void kmeans_sequential(int N, int K, int* data_points, int** data_point_cluster, 
int** centroids, int* num_iterations) {
    std::vector<int*>* centroidListHeads = new std::vector<int*>();
    int* currentCentroids = (int*) malloc(sizeof(int)*K*3);
    
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
            int iThCurrCentX = *(currentCentroids+(3*i));
            int iThCurrCentY = *(currentCentroids+(3*i)+1);
            int iThCurrCentZ = *(currentCentroids+(3*i)+2);
            if (pickedX==iThCurrCentX & pickedY==iThCurrCentY & pickedZ==iThCurrCentZ) {
                matched = true;
                break;
            }
        }

        // pickable centroid found
        if (!matched) {
            *(currentCentroids+(3*numRandomCentInit)) = pickedX;
            *(currentCentroids+(3*numRandomCentInit)+1) = pickedY;
            *(currentCentroids+(3*numRandomCentInit)+2) = pickedZ;
            numRandomCentInit++;
        }
    }

    // successfully pick random points as centroids (Forgy)
    centroidListHeads->push_back(currentCentroids);

    // update centroids
    *num_iterations = -1;
    int maxError = 4;

    std::vector<int> partitionEntries(N);
    std::vector<std::tuple<int, int, int>> sumOfPointsInPartition(K);
    std::vector<int> numOfPointsInPartition(K);
    bool zeroThIteration = true;

    while (maxError>3) {
        if (!zeroThIteration) {
            // compute centroids
            int* nextCentroids = (int*) malloc(sizeof(int)*K*3);
            for (int i=0; i<K; i++) {
                // initialize memory of next set of centroids this will contain sum of points of partition temporarily
                *(nextCentroids+(3*i)) = 0;
                *(nextCentroids+(3*i)+1) = 0;
                *(nextCentroids+(3*i)+2) = 0;
                numOfPointsInPartition[i] = 0;
            }
            
            // prepare sums and nums for avg
            for (int i=0; i<N; i++) {
                int partitionId = partitionEntries[i];
                int pickedX = *(data_points+(3*i));
                int pickedY = *(data_points+(3*i)+1);
                int pickedZ = *(data_points+(3*i)+2);
                *(nextCentroids+(3*partitionId)) += pickedX;
                *(nextCentroids+(3*partitionId)+1) += pickedY;
                *(nextCentroids+(3*partitionId)+2) += pickedZ;
                numOfPointsInPartition[partitionId] += 1;
            }

            // calculate avgs i.e. centroids
            for (int i=0; i<K; i++) {
                int iThNextCentX = *(nextCentroids+(3*i));
                int iThNextCentY = *(nextCentroids+(3*i)+1);
                int iThNextCentZ = *(nextCentroids+(3*i)+2);
                float denominator = float(numOfPointsInPartition[i]);
                iThNextCentX = round(float(iThNextCentX)/denominator);
                iThNextCentY = round(float(iThNextCentY)/denominator);
                iThNextCentZ = round(float(iThNextCentZ)/denominator);
            }

            // complete update of centroids
            centroidListHeads->push_back(nextCentroids);
            currentCentroids = nextCentroids;
        } else {
            // skip above on 0th iteration as we have to use init centroids
            zeroThIteration = false;
        }

        // assign partition
        for (int i=0; i<N; i++) {
            float minDistance = std::numeric_limits<float>::max();
            int pickedX = *(data_points+(3*i));
            int pickedY = *(data_points+(3*i)+1);
            int pickedZ = *(data_points+(3*i)+2);
            for (int j=0; j<K; j++) {
                int jThCurrCentX = *(currentCentroids+(3*j));
                int jThCurrCentY = *(currentCentroids+(3*j)+1);
                int jThCurrCentZ = *(currentCentroids+(3*j)+2);
                float currCentDistance = sqrt(powf(pickedX-jThCurrCentX,2)+powf(pickedY-jThCurrCentY,2)+powf(pickedZ-jThCurrCentZ,2));
                if (currCentDistance<minDistance) {
                    minDistance = currCentDistance;
                    partitionEntries[i] = j;
                }
            }
        }
        
        // error computation
        int* prevCentroids = currentCentroids-K*3;
        maxError = 3;
        for (int i=0; i<K; i++) {
            int iThCurrCentX = *(currentCentroids+(3*i));
            int iThCurrCentY = *(currentCentroids+(3*i)+1);
            int iThCurrCentZ = *(currentCentroids+(3*i)+2);
            int iThPrevCentX = *(prevCentroids+(3*i));
            int iThPrevCentY = *(prevCentroids+(3*i)+1);
            int iThPrevCentZ = *(prevCentroids+(3*i)+2);
            int eachCentroidError = abs(iThCurrCentX-iThPrevCentX) + abs(iThCurrCentY-iThPrevCentY) + abs(iThCurrCentZ-iThPrevCentZ);
            if (eachCentroidError>maxError) {
                maxError = eachCentroidError;
            }
        }
        *num_iterations += 1;
    }

    int* clusterPartIdList = (int*) malloc(sizeof(int)*K*4);
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
    centroids = centroidListHeads->data();
    *data_point_cluster = clusterPartIdList;
}