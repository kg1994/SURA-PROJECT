
#include "viterbi.h"
#include <iostream>
#include <algorithm>
using namespace std;




std::vector<int> ViterbiPath(std::vector<double> init, std::vector<vector<double> > trans, std::vector<vector<double> > obs, double *delta,std::vector<string> m_chordnames,vector<string> scale_dict) {

    


    int nState = init.size();
 
    cerr << trans.size() << " " <<trans[0].size() <<endl;
    int nFrame = obs.size(); 
    //cerr << "VALUES - - "  << nState << " " << nFrame << endl;
    
    // check for consistency
    if (trans[0].size() != nState || trans.size() != nState || obs[0].size() != nState) {
        cerr << "ERROR: matrix sizes inconsistent." << endl;
    }
      
   // vector<string> scale_dict = chords[index];
    // F#m- scale_dict.push_back("F#m");scale_dict.push_back("G#dim");scale_dict.push_back("A");scale_dict.push_back("Bm");scale_dict.push_back("C#m");scale_dict.push_back("D");scale_dict.push_back("E");
    //scale_dict.push_back("F");scale_dict.push_back("Gm");scale_dict.push_back("Am");scale_dict.push_back("Bb");scale_dict.push_back("C");scale_dict.push_back("Dm");scale_dict.push_back("Edim");
   // E MAJOR scale_dict.push_back("E");scale_dict.push_back("F#m");scale_dict.push_back("G#m");scale_dict.push_back("A");scale_dict.push_back("B");scale_dict.push_back("C#m");scale_dict.push_back("D#dim");
  // DMAJOR -  scale_dict.push_back("D");scale_dict.push_back("Em");scale_dict.push_back("F#m");scale_dict.push_back("G");scale_dict.push_back("A");scale_dict.push_back("Bm");scale_dict.push_back("C#dim");

   // A MAJOR- scale_dict.push_back("A");scale_dict.push_back("Bm");scale_dict.push_back("C#m");scale_dict.push_back("D");scale_dict.push_back("E");scale_dict.push_back("F#m");scale_dict.push_back("G#dim");
    //scale_dict.push_back("C");scale_dict.push_back("Dm");scale_dict.push_back("Em");scale_dict.push_back("F");scale_dict.push_back("G");scale_dict.push_back("Am");scale_dict.push_back("Bdim");

    vector<int> factor(m_chordnames.size());
    for(int i=0;i<m_chordnames.size();i++)
    {
        if(find(scale_dict.begin(), scale_dict.end(), m_chordnames[i]) != scale_dict.end())
            {factor[i] = 5; }
        else
            factor[i] = 1;

    }

    // double maxvalue=0;
    // double minvalue=10000000000;

    // vector<vector<double> > delta; // "matrix" of conditional probabilities    
    vector<vector<int> > psi; //  "matrix" of remembered indices of the best transitions
    vector<int> path = vector<int>(nFrame, nState-1); // the final output path (current assignment arbitrary, makes sense only for Chordino, where nChord-1 is the "no chord" label)
    vector<double> scale = vector<double>(nFrame, 0); // remembers by how much the vectors in delta are scaled.
    
    double deltasum = 0;
    
    



    /* initialise first frame */
    // delta.push_back(init);    */
    for (int iState = 0; iState < nState; ++iState) {
        delta[iState] = init[iState] * obs[0][iState];
        deltasum += delta[iState];
    }

    
    for (int iState = 0; iState < nState; ++iState) delta[iState] /= deltasum; // normalise (scale)
    scale.push_back(1.0/deltasum);
    psi.push_back(vector<int>(nState,0));
    
    /* rest of the forward step */
    for (int iFrame = 1; iFrame < nFrame; ++iFrame) {
        // delta.push_back(vector<double>(nState,0));
        deltasum = 0;
        psi.push_back(vector<int>(nState,0));

        /* every state wants to know which previous state suits him best */
        for (int jState = 0; jState < nState; ++jState) {            
            int bestState = nState - 1;
            double bestValue = 0;
            if (obs[iFrame][jState] > 0) {
                for (int iState = 0; iState < nState; ++iState) {
                    double currentValue = delta[(iFrame-1) * nState + iState] * trans[iState][jState]*factor[iState];


            //         if(currentValue>maxvalue) maxvalue= currentValue;
            // if(currentValue<minvalue && currentValue!=0) minvalue = currentValue;
                    // if((factor[iState]!=1 && (bestValue-currentValue)<=0.000000000000001 && ((bestValue-currentValue)>=0) ))
                    //     cerr << " ENTER \n";

                    //between 0.0001 and 0.00001 (factor[iState]!=1 && (bestValue-currentValue)<=0.00002 && ((bestValue-currentValue)>=0) )
                    if ( currentValue > bestValue ) {
                        bestValue = currentValue;
                        bestState = iState;
                    }
                }
            }
            
            // cerr << bestState <<" ::: " << bestValue << endl ;
            delta[iFrame * nState + jState] = bestValue * obs[iFrame][jState];
            deltasum += delta[iFrame * nState + jState];
            psi[iFrame][jState] = bestState;
        }
        if (deltasum > 0) {
            for (int iState = 0; iState < nState; ++iState) {            
                delta[iFrame * nState + iState] /= deltasum; // normalise (scale)
            }
            scale.push_back(1.0/deltasum);
        } else {
            for (int iState = 0; iState < nState; ++iState) {            
                delta[iFrame * nState + iState] = 1.0/nState;
            }
            scale.push_back(1.0);
        }
        
    }
    //cerr << maxvalue << " -------- " << minvalue << endl;
    /* initialise backward step */
    int bestValue = 0;
    for (int iState = 0; iState < nState; ++iState) {
        double currentValue = delta[(nFrame-1) * nState + iState];
        if (currentValue > path[nFrame-1]) {
            cerr << " entering the loop \n";
            bestValue = currentValue;            
            path[nFrame-1] = iState;
        }
    }
    //cout << path[nFrame-1] << endl;
    /* rest of backward step */
    for (int iFrame = nFrame-2; iFrame > -1; --iFrame) {
       // cerr << "  _--------_ " << iFrame+1 << " " << path[iFrame+1] <<  endl;

        path[iFrame] = psi[iFrame+1][path[iFrame+1]];
        //cerr << path[iFrame] << endl;
    }    


    return path;
}
