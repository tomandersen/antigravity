//
//  main.cpp
//  antigravity
//
//  Created by Tom Andersen on 2019-08-21.
//  Copyright © 2019 Tom Andersen. All rights reserved.
//

//#include <iostream>
//
//int main(int argc, const char * argv[]) {
//    // insert code here...
//    std::cout << "Hello, Worldjjjfff!\n";
//    return 0;
//}
//
#include <cstring>

#include <boost/lambda/lambda.hpp>
#include <iostream>
#include <iterator>
#include <algorithm>
#include <boost/program_options.hpp>
#include <exception>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/tensor.hpp>
//#include <boost/multiprecision/float128.hpp>
//using namespace boost::multiprecision;

typedef long double bigFloat;
//typedef float128 bigFloat;



using namespace boost::numeric::ublas;
using std::cout;



// units are SI
const bigFloat khbar = 1.0545718e-34; // 1.0545718 × 10-34 m2 kg / s google search
const bigFloat kBigG = 6.67430e-11;     // https://en.wikipedia.org/wiki/Gravitational_constant
const bigFloat kc = 299792458.0; // m/s
//const bigFloat kMeterPerNM = 1e-9; 1e9
const bigFloat km_p = 1.6726219e-27; // kilograms

namespace po = boost::program_options;
// Make struct/class for settings
class Settings
{
    public:
        bigFloat lattice;
    
        // inspection point === 'i'
        bigFloat ix;
        bigFloat iy;
        bigFloat iz;

        bigFloat angularMomentum;
        bigFloat mass; // in kg
    
        // number of grid points
        long nX;
        long nY;
        long nZ;
};
// make functions to evaluate Kerr metric at given x, y, z
static void playWithTensor(void);
static tensor<bigFloat> calcMetric(const Settings& settings);
static vector<bigFloat> calcGravityForceBertschinger(const Settings& settings);
static tensor<bigFloat> minkowskiMetric(void);
static void addMetricContribution(tensor<bigFloat>& metric, bigFloat x, bigFloat y, bigFloat z, bigFloat a, bigFloat mass);
static bigFloat calcLittleR(bigFloat aX, bigFloat aY, bigFloat aZ, bigFloat as);
static tensor<bigFloat> calcDifferentialMetric(const Settings& settings, bigFloat delX, bigFloat delY, bigFloat delZ);
static vector<bigFloat> calculateNewtonianGravity(const Settings& settings);
static vector<bigFloat> crossProd3d(vector<bigFloat> a, vector<bigFloat> b);

// util to read an item
static bigFloat readItem(po::variables_map &vm, const char* itemName, bigFloat defaultVal)
{
    using std::cout;
    if (vm.count(itemName)) {
        cout << itemName << " set to " << vm[itemName].as<bigFloat>() << " \n";
        return vm[itemName].as<bigFloat>();
    }
    cout << itemName << " is default " << defaultVal << " .\n";
    return defaultVal;
}
// read commands
static int inputVariables(int ac, char** av, Settings& outSettings)
{
    using std::cerr;
    using std::cout;
    using std::endl;
    using std::exception;
    try {
        
        po::options_description desc("Allowed options");
        desc.add_options()
        ("help", "This program calculates the metric underneath an array of spinning atoms, assuming that each atom is some kind of spinning Kerr like system\n The atoms are spinning in the z direction, so for max effect locate the observation point ")
        ("lattice", po::value<bigFloat>(), "atomic lattice spacing, meters")
        ("ix", po::value<bigFloat>(), "x distance from crystal face, meters")
        ("iy", po::value<bigFloat>(), "y distance from crystal face, meters")
        ("iz", po::value<bigFloat>(), "z distance from crystal face, meters")
        ("angularMomentum", po::value<bigFloat>(), "spin per atom, SI units  - eg enter hbar 1.0545718 × 10-34 m2 kg")
        ("mass", po::value<bigFloat>(), "mass per atom, kg - eg 1.6726219e-27")
        ("nX", po::value<bigFloat>(), "number of atoms in sample, X eg 1000")
        ("nY", po::value<bigFloat>(), "number of atoms in sample, Y eg 20")
        ("nZ", po::value<bigFloat>(), "number of atoms in sample, Z eg 1000");

        po::variables_map vm;
        po::store(po::parse_command_line(ac, av, desc), vm);
        po::notify(vm);
        
        if (vm.count("help")) {
            cout << desc << "\n";
            return 0;
        }
        
        outSettings.lattice = readItem(vm, "lattice", 1e-10); // 1e-10 for atomic lattice, 1000 for earth
        outSettings.ix = readItem(vm, "ix", 0);
        outSettings.iy = readItem(vm, "iy", (1.0L + 2.74525499e-11L)*khbar/(km_p*kc));//khbar/(km_p*kc)); //); // 1e-10 or near there for atomic lattice, 6378000 for earth
        outSettings.iz = readItem(vm, "iz", 1e-26);
        outSettings.angularMomentum = readItem(vm, "angularMomentum", khbar); // proton is khbar Earth is 7.07e33
        outSettings.mass = readItem(vm, "mass", km_p); // km_p proton // 5.972e24  kg for earth
        outSettings.nX = readItem(vm, "nX", 1);
        outSettings.nY = readItem(vm, "nY", 1);
        outSettings.nZ = readItem(vm, "nZ", 1);
    }
    catch(exception& e) {
        cerr << "error: " << e.what() << "\n";
        return 1;
    }
    catch(...) {
        cerr << "Exception of unknown type!\n";
        return 1;
    }
    return 0;
}

// Section on Kerr metric calculations:



// Following  notation,
// Matt Visser: The Kerr spacetime: A brief introduction
int main(int ac, char** av)
{
    cout << " our floating point type is " << 8*sizeof(bigFloat) << "bits\n";
    
    using namespace boost::lambda;
    using std::cout;
    Settings settings;
    inputVariables(ac, av, settings);
    
    bigFloat a = settings.angularMomentum/(settings.mass*kc);
    cout << "Ring is at " << a << " metres \n";

    tensor<bigFloat> metricPertubations = calcMetric(settings);
    cout << "\n\nFinal Metric pertubation is\n" << metricPertubations << std::endl;

    // now do it without spin:
//    settings.angularMomentum = 0.0;
//    tensor<bigFloat> metricPertubationsNoSpin = calcMetric(settings);
//    cout << "\n\n No Spin Metric pertubation is \n" << metricPertubationsNoSpin << std::endl;
//    tensor<bigFloat> diff = (metricPertubations - metricPertubationsNoSpin)/metricPertubationsNoSpin;
//    cout << "\n\n Diff in Metric pertubation is \n" << diff << std::endl;
    
    // as a test, etc.
    vector<bigFloat> hField = calcGravityForceBertschinger(settings);
    cout << "\nhField is " << hField << "\n\n";
    vector<bigFloat> velocity (3); // velocity needs to be in m/s for my way of calculating.
    velocity[0] = 0;
    velocity[1] = -kc;
    velocity[2] = 0;
    vector<bigFloat> gravimagneticAccel = crossProd3d(velocity, hField);
    cout << "\n gravimagneticAccel is " << gravimagneticAccel << "\n\n";
    // now do the classical newton gravity calculation:
    vector<bigFloat> grav = calculateNewtonianGravity(settings);
    cout << "\ngrav is " << grav << "\n\n";

}

// Metric calculation: returns metric pertubation...
static tensor<bigFloat> calcMetric(const Settings& settings) {
    tensor<bigFloat> metricPertubation{4,4}; // all zeros
    // run loop
    long halfNX = settings.nX/2;
    long halfNZ = settings.nZ/2;
    bigFloat lattice = settings.lattice;

    // calculate a in dimensionless units based on real units.
    bigFloat mass = settings.mass; // perhaps a setting
    bigFloat a = settings.angularMomentum/(mass*kc); // a is J/Mc https://en.wikipedia.org/wiki/Kerr_metric

    bigFloat ix = settings.ix;
    bigFloat iy = settings.iy;
    bigFloat iz = settings.iz;
    // the chunk of matter is centered on x = 0, z = 0, and goes back along neg Y axis

   // net location is described in coordinates using the atom as the base, so we use stuff like netY = iy - aY:
   for (long ny = 0; ny < settings.nY; ny++){
        bigFloat aY = -lattice*ny; // aY === 'atom Y' - position of atom
        bigFloat netY = iy - aY;
        for (long nx = -halfNX; nx <= halfNX; nx++){
            bigFloat aX = lattice*nx;
            bigFloat netX = ix - aX;
            for (long nz = -halfNZ; nz <= halfNZ; nz++){
                bigFloat aZ = lattice*nz;
                bigFloat netZ = iz - aZ;
                
                // pertubations to the metric add linearly
                addMetricContribution(metricPertubation, netX, netY, netZ, a, mass);
            }
        }
    }
//    tensor<bigFloat> minkowski = minkowskiMetric();
//    tensor<bigFloat> metric = minkowski + metricPertubation;
    
    return metricPertubation;
}

static vector<bigFloat> crossProd3d(vector<bigFloat> a, vector<bigFloat> b)
{
 //A x B = (a2b3  -   a3b2,     a3b1   -   a1b3,     a1b2   -   a2b1)
    vector<bigFloat> out(3);
    out[0] = a[1]*b[2] - a[2]*b[1];
    out[1] = a[2]*b[0] - a[0]*b[2];
    out[2] = a[0]*b[1] - a[1]*b[0];
    return out;
}


// a has units of metres, mass in kg, the is SI
static void addMetricContribution(tensor<bigFloat>& metric, bigFloat x, bigFloat y, bigFloat z, bigFloat a, bigFloat mass)
{
    // calculate all the metric add ons from a Kerr solution located at aX, aY, aZ, at point xi, yi, zi
    // use equation 43 in Visser.
    bigFloat r = calcLittleR(x, y, z, a);
    bigFloat r2 = r*r;
    bigFloat z2 = z*z;
    bigFloat a2 = a*a;
    bigFloat factor = 2.0*kBigG*mass*r2*r/(kc*kc)/(r2*r2 + a2*z2);
    

    //cout << "\nfactor is " << factor;
    //bigFloat derivative = kBigG*mass*(1/(r + 0.001*r) - 1/(r - 0.001*r))/(0.002*r);
    //cout << "\nderivative is " << derivative;
   //cout << "\nfactor field is " << factor/(2.0*r);. 
    
    // A note on units and including 'c' in the metric
    // I include c in the metric so that calculations like getting g from
    // Gravitation in the weak-field limit
    // equation 15 Bertschinger, Edmund 2000
    // works directly. In other words the units on the metric are as follows:
    // h_00 -- m^2/t^2 (so Del h_00 is the grav force in m/s^2)
    // h_i0, h_0i m/s - so that v cross curl w is a force m/s^2 (i = 1, 2, 3)
    // h_ij are unitless, so that equation 15 second term has m/s^2 after all the velocity multiplying and derivatives
    

    // calc little ls equation 44 VISSER, with c added in for the oth component.
    bigFloat l0 = kc;
    bigFloat l1 = (r*x + a*y)/(r2 + a2);
    bigFloat l2 = (r*y - a*x)/(r2 + a2);
    bigFloat l3 = (z/r);

    // we have 4 components on diagonal, 6 off diagonal
    bigFloat m00 = factor*l0*l0;
    bigFloat m01 = factor*l0*l1;
    bigFloat m02 = factor*l0*l2;
    bigFloat m03 = factor*l0*l3;
    
    bigFloat m10 = m01;
    bigFloat m11 = factor*l1*l1;
    bigFloat m12 = factor*l1*l2;
    bigFloat m13 = factor*l1*l3;
    
    bigFloat m20 = m02;
    bigFloat m21 = m12;
    bigFloat m22 = factor*l2*l2;
    bigFloat m23 = factor*l2*l3;
    
    bigFloat m30 = m03;
    bigFloat m31 = m13;
    bigFloat m32 = m23;
    bigFloat m33 = factor*l3*l3;
    
    metric.at(0,0) += m00;
    metric.at(0,1) += m01;
    metric.at(0,2) += m02;
    metric.at(0,3) += m03;
    
    metric.at(1,0) += m10;
    metric.at(1,1) += m11;
    metric.at(1,2) += m12;
    metric.at(1,3) += m13;
    
    metric.at(2,0) += m20;
    metric.at(2,1) += m21;
    metric.at(2,2) += m22;
    metric.at(2,3) += m23;
    
    metric.at(3,0) += m30;
    metric.at(3,1) += m31;
    metric.at(3,2) += m32;
    metric.at(3,3) += m33;
    static int count = 1;
    if (count++ % 1000000000 == 0)
        cout << "Metric pertubation is now is\n" << metric << std::endl;
}

static bigFloat calcLittleR(bigFloat netX, bigFloat netY, bigFloat netZ, bigFloat as)
{
    // as per equation 47 in The Kerr spacetime: A brief introduction - Matt Visser
    bigFloat rSQ = netX*netX + netY*netY + netZ*netZ;
    bigFloat rSQMaSQ = rSQ - as*as;
    bigFloat smallRoot = sqrt(rSQMaSQ*rSQMaSQ + 4*as*as*netZ*netZ);
    bigFloat r = sqrt((rSQMaSQ + smallRoot)/2.0);
    
    // improve using expansion when in tight
    bigFloat delta = 4*as*as*netZ*netZ/(rSQMaSQ*rSQMaSQ);
    if (rSQMaSQ < 0.0 && delta < 0.01) {
        bigFloat checkr = sqrt(-4*as*as*netZ*netZ/rSQMaSQ/4.0);
        //cout << "\ncheck r is " << checkr << " r is " << r << " error is " << (checkr - r)/(checkr + r) << " \n";
        return checkr;
    }
    
    

    return r;
}

static vector<bigFloat> calcGravityForceBertschinger(const Settings& settings)
{
    // we want to get the gravity
    // Equation 11 Visser:
    // -2phi = h_00
    // And Equation 15:
    // g is the gradient of phi + the time derivative of the first column, w, since our source is static, that second part is zero
    // so we need to get some gradients:
    // use the lattice spacing to come up with a step to do the first derivative around.
    bigFloat del = fmin(settings.lattice*0.001, sqrt(settings.ix*settings.ix + settings.iy*settings.iy + settings.iz*settings.iz)*0.002);
    
    // we need the distance from the r = 0 disc. The closer we are to the disk the smaller step we need to take
    bigFloat a = settings.angularMomentum/(settings.mass*kc); // a is J/Mc https://en.wikipedia.org/wiki/Kerr_metric
    // This calculates the distance to a ring located at 0,0,0 with us where we are, which is the ring at 0.0.
    // If you are only working with one ring that's fine. Otherwise you may be near another ring, and we will fail at this point.
    bigFloat r = calcLittleR(settings.ix, settings.iy, settings.iz, a);
    if (del > 0.00000001*r) {
        del = 0.00000001*r;
    }
    
    tensor<bigFloat> differentialX = calcDifferentialMetric(settings, del, 0, 0);
    //cout << "Metric diffX is:\n" << differentialX << std::endl;
    // Note on signs: h_00 is minus 2 Phi (eqn 11 Visser), and grav field is minus del Phi, so final sign is positive
    cout << "\ng in X direction is -DelPhi, so -h_00/2: " << differentialX.at(0,0)/2.0;
    
    tensor<bigFloat> differentialY = calcDifferentialMetric(settings, 0, del, 0);
    //cout << "\nMetric diffY is:\n" << differentialY << std::endl;
    cout << "\ng in Y direction is -DelPhi, so -h_00/2: " << differentialY.at(0,0)/2.0;///8.394046e-10;
    
    tensor<bigFloat> differentialZ = calcDifferentialMetric(settings, 0, 0, del);
    //cout << "\nMetric diffZ is:\n" << differentialZ << std::endl;
    cout << "\ng in Z direction is -DelPhi, so -h_00/2: " << differentialZ.at(0,0)/2.0;
    
    bigFloat xComp = differentialY.at(0,3) - differentialZ.at(0,2);
    bigFloat yComp = differentialZ.at(0,1) - differentialX.at(0,3);
    bigFloat zComp = differentialX.at(0,2) - differentialY.at(0,1);
    vector<bigFloat> hField (3); // https://en.wikipedia.org/wiki/Gravitoelectromagnetism gives 1.012×10−14 Hz I am a factor of 4 more. That's oK. I'm using Bertschinger. https://en.wikipedia.org/wiki/Talk:Gravitoelectromagnetism
    hField[0] = xComp;
    hField[1] = yComp;
    hField[2] = zComp;
    return hField;
}



static vector<bigFloat> calculateNewtonianGravity(const Settings& settings){
    long halfNX = settings.nX/2;
    long halfNZ = settings.nZ/2;
    bigFloat lattice = settings.lattice;
    
    // calculate a in dimensionless units based on real units.
    bigFloat mass = settings.mass; // perhaps a setting
    bigFloat fieldX = 0;
    bigFloat fieldY = 0;
    bigFloat fieldZ = 0;

    bigFloat ix = settings.ix;
    bigFloat iy = settings.iy;
    bigFloat iz = settings.iz;

    for (long ny = 0; ny < settings.nY; ny++){
        bigFloat aY = -lattice*ny; // aY === 'atom Y' - position of atom
        bigFloat netY = iy - aY;
        for (long nx = -halfNX; nx <= halfNX; nx++){
            bigFloat aX = lattice*nx;
            bigFloat netX = ix - aX;
            for (long nz = -halfNZ; nz <= halfNZ; nz++){
                bigFloat aZ = lattice*nz;
                bigFloat netZ = iz - aZ;

                bigFloat distSq = netX*netX + netY*netY + netZ*netZ;
                bigFloat dist = sqrt(distSq);
                
                // Grav field is minus -- https://en.wikipedia.org/wiki/Gravitational_field
                // netX/dist is the vector component in X direction.
                fieldX += -netX/dist*kBigG*mass/(distSq);
                fieldY += -netY/dist*kBigG*mass/(distSq);
                fieldZ += -netZ/dist*kBigG*mass/(distSq);
            }
        }
    }
    
    bigFloat totalNumParticles = settings.nX*settings.nY*settings.nZ;
    bigFloat totalMass = mass*totalNumParticles;
    bigFloat approxDistance = sqrt(settings.ix*settings.ix + settings.iy*settings.iy + settings.iz*settings.iz);
    approxDistance += lattice*(settings.nY - 1)/2.0;
    bigFloat sillyApprox = kBigG*totalMass/(approxDistance*approxDistance);
    cout << "\nApprox is: " << sillyApprox << "\n";
    
    vector<bigFloat> grav (3);
    grav[0] = fieldX;
    grav[1] = fieldY;
    grav[2] = fieldZ;
    return grav;
}

static tensor<bigFloat> calcDifferentialMetric(const Settings& settings, bigFloat delX, bigFloat delY, bigFloat delZ) {
    Settings settingsPlus = settings;
    settingsPlus.ix += delX;
    settingsPlus.iy += delY;
    settingsPlus.iz += delZ;
    tensor<bigFloat> metricPlus = calcMetric(settingsPlus);
    Settings settingsMinus = settings;
    settingsMinus.ix -= delX;
    settingsMinus.iy -= delY;
    settingsMinus.iz -= delZ;
    tensor<bigFloat> metricMinus = calcMetric(settingsMinus);
    
    tensor<bigFloat> differential = (metricPlus - metricMinus)/(2.0*sqrt(delX*delX + delY*delY + delZ*delZ));
    //tensor<bigFloat> differential = (metricPlus - metricMinus)/(2.0*(delX*delX + delY*delY + delZ*delZ));
    return differential;
}


static tensor<bigFloat> minkowskiMetric(void){
    tensor<bigFloat> metric{4,4};
    metric.at(0,0) = -1.0;
    metric.at(1,1) = 1.0;
    metric.at(2,2) = 1.0;
    metric.at(3,3) = 1.0;
    return metric;
}


static void playWithTensor(void)
{
    tensor<bigFloat> metric{4,4}; // see https://www.boost.org/doc/libs/1_70_0/libs/numeric/ublas/doc/tensor.html#tensor_1
    
    // initialize the tensor as per equation 43 with only the flat metric
    // I assume linearized metric. The expression for each atom in the lattice is exact,
    // but the field is weak at the measure point, ix, iy, iz
    // so we can add the second halfs linearly it seems
    // the Minkowski metric has signature - + + +
    metric.at(0,0) = -1.0;
    metric.at(1,1) = 1.0;
    metric.at(2,2) = 1.0;
    metric.at(3,3) = 1.0;
    cout << "Metric is initally\n" << metric << std::endl;

    // https://github.com/BoostGSoC18/tensor/wiki/Documentation
    vector<bigFloat> la (4);
    // see https://www.boost.org/doc/libs/1_70_0/libs/numeric/ublas/doc/tensor.html#tensor_1
    la[0] = -1;
    la[1] = 1;
    la[2] = 4;
    la[3] = 6;
    
    vector<bigFloat> lb (4);
    // see https://www.boost.org/doc/libs/1_70_0/libs/numeric/ublas/doc/tensor.html#tensor_1
    lb[0] = 1;
    lb[1] = 2;
    lb[2] = 3;
    lb[3] = 4;

    tensor<bigFloat> vec{4,1};
    vec[0] = -0.1;
    vec[1] = 0.2;
    vec[2] = 0.3;
    vec[3] = 0.4;
    auto Z = outer_prod(vec,vec);
    tensor<bigFloat> zTensor(Z); // see https://www.boost.org/doc/libs/1_70_0/libs/numeric/ublas/doc/tensor.html#tensor_1
//    zTensor = Z;
    cout << "Metric add is\n" << zTensor << std::endl;
    //    typedef std::istream_iterator<bigFloat> in;
//
//    std::for_each(in(std::cin), in(),
//        std::cout << "*3\n----\n"  << (_1 * 3) << "\n\n"
//    );
    
}
