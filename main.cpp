#include <iostream>
#include <cmath>
#include <fstream>
#include <ctime>

// particle parameters
const double dparticleL{2 * (3.985714285714286e-5)};            // diameter of large particles
const double dparticleS{dparticleL / 10};                       // diameter of small particles
const double dparticleM{(dparticleL + dparticleS) / 10};        // median diameter
const double dNozzle{8.37e-4};                                  // diamter of nozzle
const double zHeight{6.2e-4};                                   // zHeight of nozzle
const int scale = 1;                                            // scaling factor
const double xmin{(dparticleL / 2) * scale};                    // xmin for location
const double xmax{(dNozzle - xmin) * scale};                    // xmax for location
const double ymin{(zHeight / 10) * scale};                      // ymin for location
const double ymax{2 * dNozzle * scale};                         // ymax for location

const double areaNozzle{(xmax - xmin) * (ymax - ymin)};
const int nParticles{500};
double randPos;                                                 // new random position
double rState[nParticles][3];                                   // state matrix

int main() {
    std::cout << "xmin = " << xmin << " xmax = " << xmax << " ymin = " << ymin << " ymax = " << ymax << std::endl;

    srand(time(NULL));      // seed for random number
    // initialize random positions in state matrix
    for (int i=0; i<nParticles; i++)
        for (int j=0; j<3; ++j){
            randPos = (double) rand() / RAND_MAX;
            if (j == 0){
                rstate[i][j] = randPos*(xmax - xmin + 1) + xmin;
            }
            else if ( j == 1 ){
                rstate[i][j] = randPos*(ymax - ymin + 1) + ymin;
            }
            else if ( j == 2 ){
                rstate[i][j] = randPos*()
            }

    }
    return 0;
}
