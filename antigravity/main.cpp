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

using namespace boost::numeric::ublas;
using std::cout;

// units are SI
const double khbar = 1.0545718e-34; // 1.0545718 × 10-34 m2 kg / s google search
const double kBigG = 6.67430e-11;     // https://en.wikipedia.org/wiki/Gravitational_constant
const double kc = 299792458.0; // m/s
//const double kMeterPerNM = 1e-9; 1e9
const double km_p = 1.6726219e-27; // kilograms

namespace po = boost::program_options;
// Make struct/class for settings
class Settings
{
    public:
        double lattice;
    
        // inspection point === 'i'
        double ix;
        double iy;
        double iz;

        double angularMomentum;
        double mass; // in kg
    
        // number of grid points
        long nX;
        long nY;
        long nZ;
};
// make functions to evaluate Kerr metric at given x, y, z
static void playWithTensor(void);
static tensor<double> calcMetric(const Settings& settings);
static void calcNewtonianGravity(const Settings& settings);
static tensor<double> minkowskiMetric(void);
static void addMetricContribution(tensor<double>& metric, double x, double y, double z, double a, double mass);
static double calcLittleR(double aX, double aY, double aZ, double as);
static tensor<double> calcDifferentialMetric(const Settings& settings, double delX, double delY, double delZ);
static vector<double> calculateForceMagnitude(const Settings& settings);

// util to read an item
static double readItem(po::variables_map &vm, const char* itemName, double defaultVal)
{
    using std::cout;
    if (vm.count(itemName)) {
        cout << itemName << " set to " << vm[itemName].as<double>() << " \n";
        return vm[itemName].as<double>();
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
        ("lattice", po::value<double>(), "atomic lattice spacing, meters")
        ("ix", po::value<double>(), "x distance from crystal face, meters")
        ("iy", po::value<double>(), "y distance from crystal face, meters")
        ("iz", po::value<double>(), "z distance from crystal face, meters")
        ("angularMomentum", po::value<double>(), "spin per atom, SI units  - eg enter hbar 1.0545718 × 10-34 m2 kg")
        ("mass", po::value<double>(), "mass per atom, kg - eg 1.6726219e-27")
        ("nX", po::value<double>(), "number of atoms in sample, X eg 1000")
        ("nY", po::value<double>(), "number of atoms in sample, Y eg 20")
        ("nZ", po::value<double>(), "number of atoms in sample, Z eg 1000");

        po::variables_map vm;
        po::store(po::parse_command_line(ac, av, desc), vm);
        po::notify(vm);
        
        if (vm.count("help")) {
            cout << desc << "\n";
            return 0;
        }
        
        outSettings.lattice = readItem(vm, "lattice", 1000);
        outSettings.ix = readItem(vm, "ix", 0.0);
        outSettings.iy = readItem(vm, "iy", 6.371e6);
        outSettings.iz = readItem(vm, "iz", 0.0);
        outSettings.angularMomentum = readItem(vm, "angularMomentum", 7.07e33);
        outSettings.mass = readItem(vm, "mass", 5.972e24); // km_p proton // 5.972 × 10^24 kg for earth
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
    using namespace boost::lambda;
    using std::cout;
    Settings settings;
    inputVariables(ac, av, settings);
    
    tensor<double> metricPertubations = calcMetric(settings);
    cout << "\n\nFinal Metric pertubation is\n" << metricPertubations << std::endl;

    // now do it without spin:
    settings.angularMomentum = 0.0;
    tensor<double> metricPertubationsNoSpin = calcMetric(settings);
    cout << "\n\n No Spin Metric pertubation is \n" << metricPertubationsNoSpin << std::endl;
    tensor<double> diff = (metricPertubations - metricPertubationsNoSpin)/metricPertubationsNoSpin;
    cout << "\n\n Diff in Metric pertubation is \n" << diff << std::endl;
    
    // as a test, etc.
    calcNewtonianGravity(settings);
}

// Metric calculation: returns metric pertubation...
static tensor<double> calcMetric(const Settings& settings) {
    tensor<double> metricPertubation{4,4}; // all zeros
    // run loop
    long halfNX = settings.nX/2;
    long halfNZ = settings.nZ/2;
    double lattice = settings.lattice;

    // calculate a in dimensionless units based on real units.
    double mass = settings.mass; // perhaps a setting
    double a = settings.angularMomentum/(mass*kc); // a is J/Mc https://en.wikipedia.org/wiki/Kerr_metric
    
    double ix = settings.ix;
    double iy = settings.iy;
    double iz = settings.iz;
    // the chunk of matter is centered on x = 0, z = 0, and goes back along neg Y axis

   // net location is described in coordinates using the atom as the base, so we use stuff like netY = iy - aY:
   for (long ny = 0; ny < settings.nY; ny++){
        double aY = -lattice*ny; // aY === 'atom Y' - position of atom
        double netY = iy - aY;
        for (long nx = -halfNX; nx <= halfNX; nx++){
            double aX = lattice*nx;
            double netX = ix - aX;
            for (long nz = -halfNZ; nz <= halfNZ; nz++){
                double aZ = lattice*nz;
                double netZ = iz - aZ;
                
                // pertubations to the metric add linearly
                addMetricContribution(metricPertubation, netX, netY, netZ, a, mass);
            }
        }
    }
//    tensor<double> minkowski = minkowskiMetric();
//    tensor<double> metric = minkowski + metricPertubation;
    
    return metricPertubation;
}




static void addMetricContribution(tensor<double>& metric, double x, double y, double z, double a, double mass)
{
    // calculate all the metric add ons from a Kerr solution located at aX, aY, aZ, at point xi, yi, zi
    // use equation 43 in Visser.
    double r = calcLittleR(x, y, z, a);
    double r2 = r*r;
    double z2 = z*z;
    double a2 = a*a;
    double factor = 2.0*kBigG*mass*r2*r/(kc*kc)/(r2*r2 + a2*z2);
    //cout << "\nfactor is " << factor;
    //double derivative = kBigG*mass*(1/(r + 0.001*r) - 1/(r - 0.001*r))/(0.002*r);
    //cout << "\nderivative is " << derivative;
   //cout << "\nfactor field is " << factor/(2.0*r);

    // calc little ls equation 44
    double l0 = kc;
    double l1 = (r*x + a*y)/(r2 + a2);
    double l2 = (r*y - a*x)/(r2 + a2);
    double l3 = (z/r);

    // we have 4 components on diagonal, 6 off diagonal
    double m00 = factor*l0*l0;
    double m01 = factor*l0*l1;
    double m02 = factor*l0*l2;
    double m03 = factor*l0*l3;
    
    double m10 = m01;
    double m11 = factor*l1*l1;
    double m12 = factor*l1*l2;
    double m13 = factor*l1*l3;
    
    double m20 = m02;
    double m21 = m12;
    double m22 = factor*l2*l2;
    double m23 = factor*l2*l3;
    
    double m30 = m03;
    double m31 = m13;
    double m32 = m23;
    double m33 = factor*l3*l3;
    
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

static double calcLittleR(double netX, double netY, double netZ, double as)
{
    // as per equation 47 in The Kerr spacetime: A brief introduction - Matt Visser
    double rSQ = netX*netX + netY*netY + netZ*netZ;
    double rSQMaSQ = rSQ - as*as;
    double smallRoot = sqrt(rSQMaSQ*rSQMaSQ + 4*as*as*netZ*netZ);
    double r = sqrt((rSQMaSQ + smallRoot)/2.0);
    return r;
}

static void calcNewtonianGravity(const Settings& settings)
{
    // we want to get the gravity
    // Equation 11 Visser:
    // -2phi = h_00
    // And Equation 15:
    // g is the gradient of phi + the time derivative of the first column, w, since our source is static, that second part is zero
    // so we need to get some gradients:
    // use the lattice spacing to come up with a step to do the first derivative around.
    double del = settings.lattice*0.001;
    tensor<double> differentialX = calcDifferentialMetric(settings, del, 0, 0);
    //cout << "Metric diffX is:\n" << differentialX << std::endl;
    // Note on signs: h_00 is minus 2 Phi (eqn 11 Visser), and grav field is minus del Phi, so final sign is positive
    cout << "\ng in X direction is -DelPhi, so -h_00/2: " << differentialX.at(0,0)/2.0;
    
    tensor<double> differentialY = calcDifferentialMetric(settings, 0, del, 0);
    //cout << "\nMetric diffY is:\n" << differentialY << std::endl;
    cout << "\ng in Y direction is -DelPhi, so -h_00/2: " << differentialY.at(0,0)/2.0;///8.394046e-10;
    
    tensor<double> differentialZ = calcDifferentialMetric(settings, 0, 0, del);
    //cout << "Metric diffZ is:\n" << differentialZ << std::endl;
    cout << "\ng in Z direction is -DelPhi, so -h_00/2: " << differentialZ.at(0,0)/2.0;
    
    // now do the classical newton gravity calculation:
    vector<double> grav = calculateForceMagnitude(settings);
    cout << "\ngrav is " << grav << "\n\n";
}

static vector<double> calculateForceMagnitude(const Settings& settings){
    long halfNX = settings.nX/2;
    long halfNZ = settings.nZ/2;
    double lattice = settings.lattice;
    
    // calculate a in dimensionless units based on real units.
    double mass = settings.mass; // perhaps a setting
    double fieldX = 0;
    double fieldY = 0;
    double fieldZ = 0;

    double ix = settings.ix;
    double iy = settings.iy;
    double iz = settings.iz;

    for (long ny = 0; ny < settings.nY; ny++){
        double aY = -lattice*ny; // aY === 'atom Y' - position of atom
        double netY = iy - aY;
        for (long nx = -halfNX; nx <= halfNX; nx++){
            double aX = lattice*nx;
            double netX = ix - aX;
            for (long nz = -halfNZ; nz <= halfNZ; nz++){
                double aZ = lattice*nz;
                double netZ = iz - aZ;

                double distSq = netX*netX + netY*netY + netZ*netZ;
                double dist = sqrt(distSq);
                
                // Grav field is minus -- https://en.wikipedia.org/wiki/Gravitational_field
                // netX/dist is the vector component in X direction.
                fieldX += -netX/dist*kBigG*mass/(distSq);
                fieldY += -netY/dist*kBigG*mass/(distSq);
                fieldZ += -netZ/dist*kBigG*mass/(distSq);
            }
        }
    }
    
    double totalNumParticles = settings.nX*settings.nY*settings.nZ;
    double totalMass = mass*totalNumParticles;
    double approxDistance = sqrt(settings.ix*settings.ix + settings.iy*settings.iy + settings.iz*settings.iz);
    approxDistance += lattice*(settings.nY - 1)/2.0;
    double sillyApprox = kBigG*totalMass/(approxDistance*approxDistance);
    cout << "\nApprox is: " << sillyApprox << "\n";
    
    vector<double> grav (3);
    grav[0] = fieldX;
    grav[1] = fieldY;
    grav[2] = fieldZ;
    return grav;
}

static tensor<double> calcDifferentialMetric(const Settings& settings, double delX, double delY, double delZ) {
    Settings settingsPlus = settings;
    settingsPlus.ix += delX;
    settingsPlus.iy += delY;
    settingsPlus.iz += delZ;
    tensor<double> metricPlus = calcMetric(settingsPlus);
    Settings settingsMinus = settings;
    settingsMinus.ix -= delX;
    settingsMinus.iy -= delY;
    settingsMinus.iz -= delZ;
    tensor<double> metricMinus = calcMetric(settingsMinus);
    
    tensor<double> differential = (metricPlus - metricMinus)/(2.0*sqrt(delX*delX + delY*delY + delZ*delZ));
    //tensor<double> differential = (metricPlus - metricMinus)/(2.0*(delX*delX + delY*delY + delZ*delZ));
    return differential;
}


static tensor<double> minkowskiMetric(void){
    tensor<double> metric{4,4};
    metric.at(0,0) = -1.0;
    metric.at(1,1) = 1.0;
    metric.at(2,2) = 1.0;
    metric.at(3,3) = 1.0;
    return metric;
}


static void playWithTensor(void)
{
    tensor<double> metric{4,4}; // see https://www.boost.org/doc/libs/1_70_0/libs/numeric/ublas/doc/tensor.html#tensor_1
    
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
    vector<double> la (4);
    // see https://www.boost.org/doc/libs/1_70_0/libs/numeric/ublas/doc/tensor.html#tensor_1
    la[0] = -1;
    la[1] = 1;
    la[2] = 4;
    la[3] = 6;
    
    vector<double> lb (4);
    // see https://www.boost.org/doc/libs/1_70_0/libs/numeric/ublas/doc/tensor.html#tensor_1
    lb[0] = 1;
    lb[1] = 2;
    lb[2] = 3;
    lb[3] = 4;

    tensor<double> vec{4,1};
    vec[0] = -0.1;
    vec[1] = 0.2;
    vec[2] = 0.3;
    vec[3] = 0.4;
    auto Z = outer_prod(vec,vec);
    tensor<double> zTensor(Z); // see https://www.boost.org/doc/libs/1_70_0/libs/numeric/ublas/doc/tensor.html#tensor_1
//    zTensor = Z;
    cout << "Metric add is\n" << zTensor << std::endl;
    //    typedef std::istream_iterator<double> in;
//
//    std::for_each(in(std::cin), in(),
//        std::cout << "*3\n----\n"  << (_1 * 3) << "\n\n"
//    );
    
}
