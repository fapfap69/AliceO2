// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file Clusterer2.cxx
/// \brief Implementation of the HMPID cluster finder
#include <algorithm>
#include "FairLogger.h" // for LOG
#include "Framework/Logger.h"
#include "HMPIDBase/Param.h"
#include "DataFormatsHMP/Cluster.h"
#include "HMPIDReconstruction/Clusterer2.h"
#include "SimulationDataFormat/MCCompLabel.h"
#include "SimulationDataFormat/MCTruthContainer.h"
#include <TStopwatch.h>

using namespace o2::hmpid;

void Clusterer2::Dig2Clu(gsl::span<const o2::hmpid::Digit> digs, std::vector<o2::hmpid::Cluster>& clus, float* pUserCut, bool isUnfold)
{
  // Finds all clusters for a given digits list provided not empty. Currently digits list is a list of all digits for a single chamber.
  // Puts all found clusters in separate lists, one per clusters.
  // Arguments: pDigAll     - list of digits for all chambers
  //            pCluAll     - list of clusters for all chambers
  //            isTryUnfold - flag to choose between CoG and Mathieson fitting
  //  Returns: none
 
  mapIndex_t pIdxMap;
  mapCharge_t pQMap;
  
  // reset the maps with -1 (empty pad)
  std::fill(&pIdxMap[0][0][0] , &pIdxMap[0][0][0] + (sizeof(pIdxMap) / sizeof(size_t)), -1);
  std::fill(&pQMap[0][0][0] , &pQMap[0][0][0] + (sizeof(pQMap) / sizeof(int16_t)), -1);

  std::vector<DigCoord> pDigCoord;
  int padChX = 0, padChY = 0, module = 0;

  for (size_t iDigIdx = 0; iDigIdx < digs.size(); iDigIdx++) { // loop to fill all the data structures.
    o2::hmpid::Digit::pad2Absolute(digs[iDigIdx].getPadID(), &module, &padChX, &padChY);
    pIdxMap[module][padChX][padChY] = iDigIdx;
    pQMap[module][padChX][padChY] = digs[iDigIdx].mQ;
    pDigCoord.push_back({module, padChX, padChY});
  }
  
  DigCoord localMax;
  int count = 0;
  int prev = clus.size();
  for (size_t iDigIdx = 0; iDigIdx < digs.size(); iDigIdx++) { // loop to calculate clusters
    if(pIdxMap[pDigCoord.at(iDigIdx).ch][pDigCoord.at(iDigIdx).x][pDigCoord.at(iDigIdx).y] == -1) {
      continue;
    }
    std::vector<size_t> theCluster;
    std::vector<DigCoord> theClusterCoord;
    regionGrowing(pQMap, pIdxMap, pDigCoord.at(iDigIdx), theCluster, localMax);
    Cluster clu;
    clu.setCh(pDigCoord.at(iDigIdx).ch);
    for(int i=0; i<theCluster.size(); i++) {
      clu.digAdd(&digs[theCluster.at(i)], pDigCoord.at(theCluster.at(i)).x, pDigCoord.at(theCluster.at(i)).y);
      theClusterCoord.push_back(pDigCoord.at(theCluster.at(i)));
    }
    clu.solve2(&clus, pUserCut, isUnfold, pQMap, theClusterCoord); // sol  
    count++;
  }
  std::cout << "Dig2Clu() : digits " << digs.size() << " clusters =" << count << " prev clust=" << prev << " actual clusters=" << clus.size() << std::endl;

  return;
}


// Funzione per eseguire l'algoritmo di region growing sull'immagine
void Clusterer2::regionGrowing(const mapCharge_t &pQMap, mapIndex_t &pIdxMap, DigCoord seed, std::vector<size_t> &theCluster, DigCoord &localMax)
{
    // Inizializza una coda per i punti adiacenti
    std::queue<DigCoord> queue;

    // Aggiungi il punto di seme iniziale alla coda e imposta come visitato
    queue.push(seed);
    theCluster.push_back(pIdxMap[seed.ch][seed.x][seed.y]);
    localMax = seed; // assume il massimo locale

    // Finchè la coda non è vuota
    while (!queue.empty()) {
        // Prendi un punto dalla coda
        DigCoord current = queue.front();
        queue.pop();

        // Verifica i punti adiacenti al punto corrente
        int neighborCh = current.ch;
        int neighborX = current.x;
        int neighborY = current.y;

        // Punto adiacente a sinistra
        neighborX = current.x - 1;
        if (neighborX >= 0 && pIdxMap[neighborCh][neighborX][neighborY] != -1) {
            DigCoord neighbor = {neighborCh, neighborX, neighborY};
            queue.push(neighbor);
            theCluster.push_back(pIdxMap[neighborCh][neighborX][neighborY]);
            pIdxMap[neighborCh][neighborX][neighborY] = -1;
        }

        // Punto adiacente a destra
        neighborX = current.x + 1;
        if (neighborX < numPadX && pIdxMap[neighborCh][neighborX][neighborY] != -1) {
            DigCoord neighbor = {neighborCh, neighborX, neighborY};
            queue.push(neighbor);
            theCluster.push_back(pIdxMap[neighborCh][neighborX][neighborY]);
            pIdxMap[neighborCh][neighborX][neighborY] = -1;
        }

        // Punto adiacente in alto
        neighborX = current.x;
        neighborY = current.y - 1;
        if (neighborY >= 0 && pIdxMap[neighborCh][neighborX][neighborY] != -1) {
            DigCoord neighbor = {neighborCh, neighborX, neighborY};
            queue.push(neighbor);
            theCluster.push_back(pIdxMap[neighborCh][neighborX][neighborY]);
            pIdxMap[neighborCh][neighborX][neighborY] = -1;
        }

        // Punto adiacente in basso
        neighborY = current.y + 1;
        if (neighborY < numPadY && pIdxMap[neighborCh][neighborX][neighborY] != -1) {
            DigCoord neighbor = {neighborCh, neighborX, neighborY};
            queue.push(neighbor);
            theCluster.push_back(pIdxMap[neighborCh][neighborX][neighborY]);
            pIdxMap[neighborCh][neighborX][neighborY] = -1;
        }
    }
    std::cout << "regionGrowing() : seed " << seed.ch << ":" << seed.x << "," << seed.y << " Cluster :" << theCluster.size() << std::endl;
  return;  

} 


//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

