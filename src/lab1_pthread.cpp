extern "C" {
    #include "lab1_pthread.h"
}

#include <vector>
#include <tuple>
#include <stdlib.h>
#include <time.h>
#include <limits>
#include <cmath>
#include <iostream>
#include <pthread.h>

std::vector<int>* numOfPointsInPartition;
std::vector<double>*currentCentroidsDouble;
std::vector<int>* partitionEntries;
int* clusterPartIdList;
int* data_pointsG;
int Ng;
int Kg;
int threadNumG;
int* tidArr;

void* meanOfPartition(void* tid) {
    int* id = (int*) tid;
    for (int i=0; i<Ng; i++) {
        int partitionId = (*partitionEntries)[i];
        if (partitionId%threadNumG==*id) {
            int pickedX = *(data_pointsG+(3*i));
            int pickedY = *(data_pointsG+(3*i)+1);
            int pickedZ = *(data_pointsG+(3*i)+2);
            (*currentCentroidsDouble)[(3*partitionId)] += pickedX;
            (*currentCentroidsDouble)[(3*partitionId)+1] += pickedY;
            (*currentCentroidsDouble)[(3*partitionId)+2] += pickedZ;
            (*numOfPointsInPartition)[partitionId] += 1;
        }
    }
    return NULL;
}

void* assignPartition(void* tid) {
    int* id = (int*) tid;
    int istart = (Ng/threadNumG)*(*id);
    int iend;
    if (*id==threadNumG-1){
        iend = Ng;
    } else {
        iend = (Ng/threadNumG)*((*id)+1);
    }
    for (int i=istart; i<iend; i++) {
        double minDistance = std::numeric_limits<double>::max();
        int pickedX = *(data_pointsG+(3*i));
        int pickedY = *(data_pointsG+(3*i)+1);
        int pickedZ = *(data_pointsG+(3*i)+2);
        for (int j=0; j<Kg; j++) {
            double jThCurrCentX = (*currentCentroidsDouble)[(3*j)];
            double jThCurrCentY = (*currentCentroidsDouble)[(3*j)+1];
            double jThCurrCentZ = (*currentCentroidsDouble)[(3*j)+2];
            double currCentDistance = pow(pickedX-jThCurrCentX,2)+pow(pickedY-jThCurrCentY,2)+pow(pickedZ-jThCurrCentZ,2);
            if (currCentDistance<minDistance) {
                minDistance = currCentDistance;
                (*partitionEntries)[i] = j;
            }
        }
    }
    return NULL;
}

void kmeans_pthread(int num_threads, int N, int K, int* data_points, int** data_point_cluster, 
int** centroids, int* num_iterations) {
    // pthread init and global init
    pthread_t multThreads[num_threads];
    tidArr = new int[num_threads];
    data_pointsG = data_points;
    Ng = N;
    Kg = K;
    threadNumG = num_threads;

    std::vector<int>* centroidSaveDB = new std::vector<int>();
    currentCentroidsDouble = new std::vector<double>(3*K);
    
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
            int iThCurrCentX = (*currentCentroidsDouble)[(3*i)];
            int iThCurrCentY = (*currentCentroidsDouble)[(3*i)+1];
            int iThCurrCentZ = (*currentCentroidsDouble)[(3*i)+2];
            if (pickedX==iThCurrCentX & pickedY==iThCurrCentY & pickedZ==iThCurrCentZ) {
                matched = true;
                break;
            }
        }

        // pickable centroid found
        if (!matched) {
            (*currentCentroidsDouble)[3*numRandomCentInit] = pickedX;
            (*currentCentroidsDouble)[3*numRandomCentInit+1] = pickedY;
            (*currentCentroidsDouble)[3*numRandomCentInit+2] = pickedZ;
            numRandomCentInit++;
        }
    }

    // successfully pick random points as centroids (Forgy)
    centroidSaveDB->insert(centroidSaveDB->end(), currentCentroidsDouble->begin(), currentCentroidsDouble->end());

    // update centroids
    *num_iterations = -1;
    double error = std::numeric_limits<double>::max();

    partitionEntries = new std::vector<int>(N);
    numOfPointsInPartition = new std::vector<int>(K);
    std::vector<double>* prevCentroidsDouble;
    bool zeroThIteration = true;

    while (error>1.0) {
        if (!zeroThIteration) {
            // compute centroids
            currentCentroidsDouble = new std::vector<double>(3*K);
            for (int i=0; i<K; i++) {
                // initialize memory of next set of centroids this will contain sum of points of partition temporarily
                (*numOfPointsInPartition)[i] = 0;
            }
            
            // prepare sums and nums for avg
            // threading
            for (int i=0; i<num_threads; i++) {
                tidArr[i] = i;
                pthread_create(&multThreads[i], NULL, meanOfPartition, &tidArr[i]);
            }
            // synchronize
            for (int i = 0; i < num_threads; i ++) {
                pthread_join(multThreads[i], NULL);
            }

            // calculate avgs i.e. centroids
            for (int i=0; i<K; i++) {
                int iThCurrCentX = (*currentCentroidsDouble)[(3*i)];
                int iThCurrCentY = (*currentCentroidsDouble)[(3*i)+1];
                int iThCurrCentZ = (*currentCentroidsDouble)[(3*i)+2];
                double denominator = double((*numOfPointsInPartition)[i]);
                (*currentCentroidsDouble)[(3*i)] = round(double(iThCurrCentX)/denominator);
                (*currentCentroidsDouble)[(3*i)+1] = round(double(iThCurrCentY)/denominator);
                (*currentCentroidsDouble)[(3*i)+2] = round(double(iThCurrCentZ)/denominator);
            }

            // complete update of centroids
            centroidSaveDB->insert(centroidSaveDB->end(), currentCentroidsDouble->begin(), currentCentroidsDouble->end());
        }

        // assign partition
        // threading
        for (int i=0; i<num_threads; i++) {
            tidArr[i] = i;
            pthread_create(&multThreads[i], NULL, assignPartition, &tidArr[i]);
        }
        // synchronize
        for (int i = 0; i < num_threads; i ++) {
            pthread_join(multThreads[i], NULL);
        }
        
        // error computation
        if (!zeroThIteration) {
            error = 0.0;
            for (int i=0; i<K; i++) {
                double iThCurrCentX = (*currentCentroidsDouble)[(3*i)];
                double iThCurrCentY = (*currentCentroidsDouble)[(3*i)+1];
                double iThCurrCentZ = (*currentCentroidsDouble)[(3*i)+2];
                double iThPrevCentX = (*prevCentroidsDouble)[(3*i)];
                double iThPrevCentY = (*prevCentroidsDouble)[(3*i)+1];
                double iThPrevCentZ = (*prevCentroidsDouble)[(3*i)+2];
                error += pow(iThCurrCentX-iThPrevCentX, 2) + pow(iThCurrCentY-iThPrevCentY, 2) + pow(iThCurrCentZ-iThPrevCentZ, 2);
            }
            prevCentroidsDouble->clear();
            delete prevCentroidsDouble;
        } else {
            zeroThIteration = false;
        }
        // copy current centroid, they are prev for next iteration
        prevCentroidsDouble = currentCentroidsDouble;
        *num_iterations += 1;
    }

    clusterPartIdList = (int*) malloc(sizeof(int)*N*4);
    for (int i=0; i<N; i++) {
        int partitionId = (*partitionEntries)[i];
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
    delete[] tidArr;
}