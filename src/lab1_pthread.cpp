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
#include <unistd.h>

std::vector<int>* numOfPointsInPartition;
std::vector<double>*currentCentroidsDouble;
std::vector<int>* partitionEntries;
int* clusterPartIdList;
int* data_pointsG;
int Ng;
int Kg;
int threadNumG;
int* tidArr;
int cacheLineSizeBlock;

void* meanOfPartition(void* tid) {
    int* id = (int*) tid;
    for (int i=0; i<Ng; i++) {
        int partitionId = (*partitionEntries)[i];
        if (partitionId%threadNumG==*id) {
            int pickedX = *(data_pointsG+(3*i));
            int pickedY = *(data_pointsG+(3*i)+1);
            int pickedZ = *(data_pointsG+(3*i)+2);
            (*currentCentroidsDouble)[(partitionId*cacheLineSizeBlock)] += pickedX;
            (*currentCentroidsDouble)[(partitionId*cacheLineSizeBlock)+1] += pickedY;
            (*currentCentroidsDouble)[(partitionId*cacheLineSizeBlock)+2] += pickedZ;
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
            double jThCurrCentX = (*currentCentroidsDouble)[(j*cacheLineSizeBlock)];
            double jThCurrCentY = (*currentCentroidsDouble)[(j*cacheLineSizeBlock)+1];
            double jThCurrCentZ = (*currentCentroidsDouble)[(j*cacheLineSizeBlock)+2];
            double currCentDistance = pow(pickedX-jThCurrCentX,2)+pow(pickedY-jThCurrCentY,2)+pow(pickedZ-jThCurrCentZ,2);
            if (currCentDistance<minDistance) {
                minDistance = currCentDistance;
                (*partitionEntries)[i] = j;
            }
        }
    }
    return NULL;
}

void* writeClusterList(void* tid) {
    int* id = (int*) tid;
    int istart = (Ng/threadNumG)*(*id);
    int iend;
    if (*id==threadNumG-1){
        iend = Ng;
    } else {
        iend = (Ng/threadNumG)*((*id)+1);
    }
    for (int i=istart; i<iend; i++) {
        int partitionId = (*partitionEntries)[i];
        int pickedX = *(data_pointsG+(3*i));
        int pickedY = *(data_pointsG+(3*i)+1);
        int pickedZ = *(data_pointsG+(3*i)+2);
        *(clusterPartIdList+(4*i)) = pickedX;
        *(clusterPartIdList+(4*i+1)) = pickedY;
        *(clusterPartIdList+(4*i)+2) = pickedZ;
        *(clusterPartIdList+(4*i)+3) = partitionId;
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

    cacheLineSizeBlock = (int) ceil(double(sysconf(_SC_LEVEL1_DCACHE_LINESIZE))/double(sizeof(double)));

    std::vector<int>* centroidSaveDB = new std::vector<int>();
    centroidSaveDB->reserve(3*K*201);
    currentCentroidsDouble = new std::vector<double>(3*K);
    currentCentroidsDouble = new std::vector<double>(K*cacheLineSizeBlock);
    
    // Initialize centroids
    // srand(time(0));
    srand(0);
    int numRandomCentInit = 0;
    while(numRandomCentInit<K) {
        int pickIndex = rand()%N;
        int pickedX = *(data_points+(3*pickIndex));
        int pickedY = *(data_points+(3*pickIndex)+1);
        int pickedZ = *(data_points+(3*pickIndex)+2);
        
        // check no centroid is repeated
        bool matched = false;
        for (int i=0; i<numRandomCentInit; i++) {
            int iThCurrCentX = (*currentCentroidsDouble)[(i*cacheLineSizeBlock)];
            int iThCurrCentY = (*currentCentroidsDouble)[(i*cacheLineSizeBlock)+1];
            int iThCurrCentZ = (*currentCentroidsDouble)[(i*cacheLineSizeBlock)+2];
            if (pickedX==iThCurrCentX & pickedY==iThCurrCentY & pickedZ==iThCurrCentZ) {
                matched = true;
                break;
            }
        }

        // pickable centroid found
        if (!matched) {
            (*currentCentroidsDouble)[numRandomCentInit*cacheLineSizeBlock] = pickedX;
            (*currentCentroidsDouble)[numRandomCentInit*cacheLineSizeBlock+1] = pickedY;
            (*currentCentroidsDouble)[numRandomCentInit*cacheLineSizeBlock+2] = pickedZ;
            numRandomCentInit++;
        }
    }

    // successfully pick random points as centroids (Forgy)
    for (int i=0; i<K; i++) {
        centroidSaveDB->insert(centroidSaveDB->end(), currentCentroidsDouble->begin()+i*cacheLineSizeBlock, 
        currentCentroidsDouble->begin()+i*cacheLineSizeBlock+3);
    }

    // update centroids
    *num_iterations = -1;
    double error = std::numeric_limits<double>::max();

    partitionEntries = new std::vector<int>(N);
    numOfPointsInPartition = new std::vector<int>(K);
    std::vector<double>* prevCentroidsDouble;
    bool zeroThIteration = true;

    while (error>1.0 && (*num_iterations<200)) {
        if (!zeroThIteration) {
            // compute centroids
            currentCentroidsDouble = new std::vector<double>(K*cacheLineSizeBlock);
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
                int iThCurrCentX = (*currentCentroidsDouble)[(i*cacheLineSizeBlock)];
                int iThCurrCentY = (*currentCentroidsDouble)[(i*cacheLineSizeBlock)+1];
                int iThCurrCentZ = (*currentCentroidsDouble)[(i*cacheLineSizeBlock)+2];
                double denominator = double((*numOfPointsInPartition)[i]);
                (*currentCentroidsDouble)[(i*cacheLineSizeBlock)] = round(double(iThCurrCentX)/denominator);
                (*currentCentroidsDouble)[(i*cacheLineSizeBlock)+1] = round(double(iThCurrCentY)/denominator);
                (*currentCentroidsDouble)[(i*cacheLineSizeBlock)+2] = round(double(iThCurrCentZ)/denominator);
            }

            // complete update of centroids
            for (int i=0; i<K; i++) {
                centroidSaveDB->insert(centroidSaveDB->end(), currentCentroidsDouble->begin()+i*cacheLineSizeBlock, 
                currentCentroidsDouble->begin()+i*cacheLineSizeBlock+3);
            }
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
                double iThCurrCentX = (*currentCentroidsDouble)[(i*cacheLineSizeBlock)];
                double iThCurrCentY = (*currentCentroidsDouble)[(i*cacheLineSizeBlock)+1];
                double iThCurrCentZ = (*currentCentroidsDouble)[(i*cacheLineSizeBlock)+2];
                double iThPrevCentX = (*prevCentroidsDouble)[(i*cacheLineSizeBlock)];
                double iThPrevCentY = (*prevCentroidsDouble)[(i*cacheLineSizeBlock)+1];
                double iThPrevCentZ = (*prevCentroidsDouble)[(i*cacheLineSizeBlock)+2];
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
    prevCentroidsDouble->clear();
    delete prevCentroidsDouble;
    delete currentCentroidsDouble;

    clusterPartIdList = (int*) malloc(sizeof(int)*N*4);
    // write data point with clusters
    // threading
    for (int i=0; i<num_threads; i++) {
        tidArr[i] = i;
        pthread_create(&multThreads[i], NULL, writeClusterList, &tidArr[i]);
    }
    // synchronize
    for (int i = 0; i < num_threads; i ++) {
        pthread_join(multThreads[i], NULL);
    }

    *centroids = centroidSaveDB->data();
    *data_point_cluster = clusterPartIdList;
    delete[] tidArr;
}