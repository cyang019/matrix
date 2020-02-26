#include "matrix.h"
#include "gtest/gtest.h"
#include <complex>
#include <iostream>

namespace {
    using MatrixCD = matrix::Matrix<std::complex<double>>; 
    using cxdbl = std::complex<double>;
    using namespace std;

    TEST(TestMatrix, Lstsq){
      MatrixCD L = {
        {
          cxdbl(500.0507882919981, 0.0), cxdbl(-0.3124355472758323, -500000.0), cxdbl(-1.897304479874008, -500000.0), cxdbl(6.01144125944056e-06, 0.0), cxdbl(-0.3124355472758323, 500000.0), cxdbl(-0.050781088867451814, 0.0), cxdbl(-3.608470504652113e-09, 0.0), cxdbl(-1.927862737394476e-07, 0.0), cxdbl(-1.897304479874008, 500000.0), cxdbl(-3.608470504652113e-09, 0.0), cxdbl(-500.00000720312994, 0.0), cxdbl(7.696423535276798e-11, 0.0), cxdbl(6.01144125944056e-06, 0.0), cxdbl(-1.927862737394476e-07, 0.0), cxdbl(7.696423535276798e-11, 0.0), cxdbl(-7.320009932344586e-13, 0.0), 
        },
        {
          cxdbl(-0.3124355472758323, -500000.0), cxdbl(750.0492261214687, 400000000.0), cxdbl(6.029709186910314e-06, 0.0), cxdbl(-1.8992027334806911, -500000.0), cxdbl(-0.0007810888681724938, 0.0), cxdbl(0.3124355472709935, 500000.0), cxdbl(-2.969919640254635e-09, 0.0), cxdbl(3.609214576170885e-09, 0.0), cxdbl(3.598273912084254e-09, 0.0), cxdbl(-1.8992027334896961, 500000.0), cxdbl(2.7165163571595614e-11, 0.0), cxdbl(499.999992789663, 0.0), cxdbl(-2.9609130712630664e-09, 0.0), cxdbl(8.38017099709243e-06, 0.0), cxdbl(-1.1362098572002174e-14, 0.0), cxdbl(-2.2326373502262248e-11, 0.0), 
        },
        {
          cxdbl(-1.897304479874008, -500000.0), cxdbl(6.029709186910314e-06, 0.0), cxdbl(1000500.0507738852, 400000000.0), cxdbl(0.3124354475252178, 500000.0), cxdbl(3.598273912084254e-09, 0.0), cxdbl(-1.927952803084399e-07, 0.0), cxdbl(-0.3124355472198245, 500000.0), cxdbl(0.050781088867482685, 0.0), cxdbl(-7.203130143060898e-06, 0.0), cxdbl(2.7401965359177662e-11, 0.0), cxdbl(1.8973044798830145, 500000.0), cxdbl(-3.609396083902739e-09, 0.0), cxdbl(2.7318930920922638e-11, 0.0), cxdbl(-7.320354496385455e-13, 0.0), cxdbl(-8.380170807336724e-06, 0.0), cxdbl(1.927862737402832e-07, 0.0), 
        },
        {
          cxdbl(6.01144125944056e-06, 0.0), cxdbl(-1.8992027334806911, -500000.0), cxdbl(0.3124354475252178, 500000.0), cxdbl(1000750.0492117009, 800000000.0), cxdbl(-2.9609130712630664e-09, 0.0), cxdbl(-3.5975320276926647e-09, 0.0), cxdbl(0.0007810888682033656, 0.0), cxdbl(0.31243554732699635, 500000.0), cxdbl(2.7318930920922638e-11, 0.0), cxdbl(-7.210336878269498e-06, 0.0), cxdbl(3.5973426489666474e-09, 0.0), cxdbl(1.899202733480691, 500000.0), cxdbl(-1.1327642163009527e-14, 0.0), cxdbl(2.7318930921279657e-11, 0.0), cxdbl(2.9609130720972735e-09, 0.0), cxdbl(-6.011441070061833e-06, 0.0), 
        },
        {
          cxdbl(-0.3124355472758323, 500000.0), cxdbl(-0.0007810888681724938, 0.0), cxdbl(3.598273912084254e-09, 0.0), cxdbl(-2.9609130712630664e-09, 0.0), cxdbl(750.0492261214687, -400000000.0), cxdbl(0.3124355472709935, -500000.0), cxdbl(-1.8992027334896961, -500000.0), cxdbl(8.38017099709243e-06, 0.0), cxdbl(6.029709186910314e-06, 0.0), cxdbl(-2.969919640254635e-09, 0.0), cxdbl(2.7165163571595614e-11, 0.0), cxdbl(-1.1362098572002174e-14, 0.0), cxdbl(-1.8992027334806911, 500000.0), cxdbl(3.609214576170885e-09, 0.0), cxdbl(499.999992789663, 0.0), cxdbl(-2.2326373502262248e-11, 0.0), 
        },
        {
          cxdbl(-0.050781088867451814, 0.0), cxdbl(0.3124355472709935, 500000.0), cxdbl(-1.927952803084399e-07, 0.0), cxdbl(-3.5975320276926647e-09, 0.0), cxdbl(0.3124355472709935, -500000.0), cxdbl(500.05078829199834, 0.0), cxdbl(8.405645620914765e-06, 0.0), cxdbl(-1.8973044798830154, -500000.0), cxdbl(-1.927952803084399e-07, 0.0), cxdbl(8.405645620914765e-06, 0.0), cxdbl(-7.32070010877866e-13, 0.0), cxdbl(2.7638767153216658e-11, 0.0), cxdbl(-3.5975320276926647e-09, 0.0), cxdbl(-1.8973044798830154, 500000.0), cxdbl(2.7638767153216658e-11, 0.0), cxdbl(-500.00000720313017, 0.0), 
        },
        {
          cxdbl(-3.608470504652113e-09, 0.0), cxdbl(-2.969919640254635e-09, 0.0), cxdbl(-0.3124355472198245, 500000.0), cxdbl(0.0007810888682033656, 0.0), cxdbl(-1.8992027334896961, -500000.0), cxdbl(8.405645620914765e-06, 0.0), cxdbl(1000750.0492117002, 0.0), cxdbl(-0.31243564704631216, 500000.0), cxdbl(2.7401965359177662e-11, 0.0), cxdbl(-1.1396659811346025e-14, 0.0), cxdbl(-8.405645423666798e-06, 0.0), cxdbl(2.969919641092091e-09, 0.0), cxdbl(-7.210336878269498e-06, 0.0), cxdbl(2.7401965373938054e-11, 0.0), cxdbl(1.8992027334896966, 500000.0), cxdbl(3.6082732566853943e-09, 0.0), 
        },
        {
          cxdbl(-1.927862737394476e-07, 0.0), cxdbl(3.609214576170885e-09, 0.0), cxdbl(0.050781088867482685, 0.0), cxdbl(0.31243554732699635, 500000.0), cxdbl(8.38017099709243e-06, 0.0), cxdbl(-1.8973044798830154, -500000.0), cxdbl(-0.31243564704631216, 500000.0), cxdbl(1000500.0507738855, 400000000.0), cxdbl(-7.320354496385455e-13, 0.0), cxdbl(2.7401965373938054e-11, 0.0), cxdbl(1.9279528030927736e-07, 0.0), cxdbl(-6.029708989660571e-06, 0.0), cxdbl(2.7318930921279657e-11, 0.0), cxdbl(-7.203130143060903e-06, 0.0), cxdbl(-3.5984632890568127e-09, 0.0), cxdbl(1.8973044798740089, 500000.0), 
        },
        {
          cxdbl(-1.897304479874008, 500000.0), cxdbl(3.598273912084254e-09, 0.0), cxdbl(-7.203130143060898e-06, 0.0), cxdbl(2.7318930920922638e-11, 0.0), cxdbl(6.029709186910314e-06, 0.0), cxdbl(-1.927952803084399e-07, 0.0), cxdbl(2.7401965359177662e-11, 0.0), cxdbl(-7.320354496385455e-13, 0.0), cxdbl(1000500.0507738852, -400000000.0), cxdbl(-0.3124355472198245, -500000.0), cxdbl(1.8973044798830145, -500000.0), cxdbl(-8.380170807336724e-06, 0.0), cxdbl(0.3124354475252178, -500000.0), cxdbl(0.050781088867482685, 0.0), cxdbl(-3.609396083902739e-09, 0.0), cxdbl(1.927862737402832e-07, 0.0), 
        },
        {
          cxdbl(-3.608470504652113e-09, 0.0), cxdbl(-1.8992027334896961, 500000.0), cxdbl(2.7401965359177662e-11, 0.0), cxdbl(-7.210336878269498e-06, 0.0), cxdbl(-2.969919640254635e-09, 0.0), cxdbl(8.405645620914765e-06, 0.0), cxdbl(-1.1396659811346025e-14, 0.0), cxdbl(2.7401965373938054e-11, 0.0), cxdbl(-0.3124355472198245, -500000.0), cxdbl(1000750.0492117002, 0.0), cxdbl(-8.405645423666798e-06, 0.0), cxdbl(1.8992027334896966, -500000.0), cxdbl(0.0007810888682033656, 0.0), cxdbl(-0.31243564704631216, -500000.0), cxdbl(2.969919641092091e-09, 0.0), cxdbl(3.6082732566853943e-09, 0.0), 
        },
        {
          cxdbl(-500.00000720312994, 0.0), cxdbl(2.7165163571595614e-11, 0.0), cxdbl(1.8973044798830145, 500000.0), cxdbl(3.5973426489666474e-09, 0.0), cxdbl(2.7165163571595614e-11, 0.0), cxdbl(-7.32070010877866e-13, 0.0), cxdbl(-8.405645423666798e-06, 0.0), cxdbl(1.9279528030927736e-07, 0.0), cxdbl(1.8973044798830145, -500000.0), cxdbl(-8.405645423666798e-06, 0.0), cxdbl(500.0507882919982, 0.0), cxdbl(0.3124355472334609, 500000.0), cxdbl(3.5973426489666474e-09, 0.0), cxdbl(1.9279528030927736e-07, 0.0), cxdbl(0.3124355472334609, -500000.0), cxdbl(-0.05078108886751355, 0.0), 
        },
        {
          cxdbl(7.696423535276798e-11, 0.0), cxdbl(499.999992789663, 0.0), cxdbl(-3.609396083902739e-09, 0.0), cxdbl(1.899202733480691, 500000.0), cxdbl(-1.1362098572002174e-14, 0.0), cxdbl(2.7638767153216658e-11, 0.0), cxdbl(2.969919641092091e-09, 0.0), cxdbl(-6.029708989660571e-06, 0.0), cxdbl(-8.380170807336724e-06, 0.0), cxdbl(1.8992027334896966, -500000.0), cxdbl(0.3124355472334609, 500000.0), cxdbl(750.0492261214686, 400000000.0), cxdbl(2.9609130720972735e-09, 0.0), cxdbl(-3.5984632890568127e-09, 0.0), cxdbl(-0.00078108886823424, 0.0), cxdbl(-0.3124355473380639, -500000.0), 
        },
        {
          cxdbl(6.01144125944056e-06, 0.0), cxdbl(-2.9609130712630664e-09, 0.0), cxdbl(2.7318930920922638e-11, 0.0), cxdbl(-1.1327642163009527e-14, 0.0), cxdbl(-1.8992027334806911, 500000.0), cxdbl(-3.5975320276926647e-09, 0.0), cxdbl(-7.210336878269498e-06, 0.0), cxdbl(2.7318930921279657e-11, 0.0), cxdbl(0.3124354475252178, -500000.0), cxdbl(0.0007810888682033656, 0.0), cxdbl(3.5973426489666474e-09, 0.0), cxdbl(2.9609130720972735e-09, 0.0), cxdbl(1000750.0492117009, -800000000.0), cxdbl(0.31243554732699635, -500000.0), cxdbl(1.899202733480691, -500000.0), cxdbl(-6.011441070061833e-06, 0.0), 
        },
        {
          cxdbl(-1.927862737394476e-07, 0.0), cxdbl(8.38017099709243e-06, 0.0), cxdbl(-7.320354496385455e-13, 0.0), cxdbl(2.7318930921279657e-11, 0.0), cxdbl(3.609214576170885e-09, 0.0), cxdbl(-1.8973044798830154, 500000.0), cxdbl(2.7401965373938054e-11, 0.0), cxdbl(-7.203130143060903e-06, 0.0), cxdbl(0.050781088867482685, 0.0), cxdbl(-0.31243564704631216, -500000.0), cxdbl(1.9279528030927736e-07, 0.0), cxdbl(-3.5984632890568127e-09, 0.0), cxdbl(0.31243554732699635, -500000.0), cxdbl(1000500.0507738855, -400000000.0), cxdbl(-6.029708989660571e-06, 0.0), cxdbl(1.8973044798740089, -500000.0), 
        },
        {
          cxdbl(7.696423535276798e-11, 0.0), cxdbl(-1.1362098572002174e-14, 0.0), cxdbl(-8.380170807336724e-06, 0.0), cxdbl(2.9609130720972735e-09, 0.0), cxdbl(499.999992789663, 0.0), cxdbl(2.7638767153216658e-11, 0.0), cxdbl(1.8992027334896966, 500000.0), cxdbl(-3.5984632890568127e-09, 0.0), cxdbl(-3.609396083902739e-09, 0.0), cxdbl(2.969919641092091e-09, 0.0), cxdbl(0.3124355472334609, -500000.0), cxdbl(-0.00078108886823424, 0.0), cxdbl(1.899202733480691, -500000.0), cxdbl(-6.029708989660571e-06, 0.0), cxdbl(750.0492261214686, -400000000.0), cxdbl(-0.3124355473380639, 500000.0), 
        },
        {
          cxdbl(-7.320009932344586e-13, 0.0), cxdbl(-2.2326373502262248e-11, 0.0), cxdbl(1.927862737402832e-07, 0.0), cxdbl(-6.011441070061833e-06, 0.0), cxdbl(-2.2326373502262248e-11, 0.0), cxdbl(-500.00000720313017, 0.0), cxdbl(3.6082732566853943e-09, 0.0), cxdbl(1.8973044798740089, 500000.0), cxdbl(1.927862737402832e-07, 0.0), cxdbl(3.6082732566853943e-09, 0.0), cxdbl(-0.05078108886751355, 0.0), cxdbl(-0.3124355473380639, -500000.0), cxdbl(-6.011441070061833e-06, 0.0), cxdbl(1.8973044798740089, -500000.0), cxdbl(-0.3124355473380639, 500000.0), cxdbl(500.0507882919984, 0.0), 
        }
      };
      MatrixCD gamma_rho0_super = {
        {
          cxdbl(-15.778991805212845, 0.0), 
        },
        {
          cxdbl(3.7811030673545742e-06, 0.0), 
        },
        {
          cxdbl(-5.990504709861602e-05, 0.0), 
        },
        {
          cxdbl(2.27084151437427e-10, 0.0), 
        },
        {
          cxdbl(3.7811030673406965e-06, 0.0), 
        },
        {
          cxdbl(-15.78201668892785, 0.0), 
        },
        {
          cxdbl(2.2831399407883732e-10, 0.0), 
        },
        {
          cxdbl(-5.9916549312110007e-05, 0.0), 
        },
        {
          cxdbl(-5.9905047098670805e-05, 0.0), 
        },
        {
          cxdbl(2.2831399407883732e-10, 0.0), 
        },
        {
          cxdbl(15.77898700596819, 0.0), 
        },
        {
          cxdbl(3.793101180246672e-06, 0.0), 
        },
        {
          cxdbl(2.27084151437427e-10, 0.0), 
        },
        {
          cxdbl(-5.9916549312097626e-05, 0.0), 
        },
        {
          cxdbl(3.793101180246672e-06, 0.0), 
        },
        {
          cxdbl(15.782021488172493, 0.0), 
        }
      };

      auto [rho_eq, status] = lstsq(L, gamma_rho0_super);
      auto calculated_rhs = L * rho_eq;
      cout << "calculated rhs\n" << calculated_rhs << endl;
      cout << "rhs" << gamma_rho0_super << endl;

      MatrixCD desired_rho_eq = {
        {
          cxdbl(-0.030875055742808813, 2.9518099639110554e-11), 
        },
        {
          cxdbl(-3.798444301407042e-05, 1.5170429889616723e-09), 
        },
        {
          cxdbl(-3.919893940999409e-05, -9.941460438025703e-08), 
        },
        {
          cxdbl(1.5179601483929651e-09, 1.2822446385072126e-10), 
        },
        {
          cxdbl(-3.7984443014296486e-05, -1.517224693785075e-09), 
        },
        {
          cxdbl(-0.00048740181672443845, 1.0219991508927755e-10), 
        },
        {
          cxdbl(-9.800205485118254e-08, 1.2135860550272787e-06), 
        },
        {
          cxdbl(-3.9206465896110466e-05, -9.943340175500611e-08), 
        },
        {
          cxdbl(-3.919893940976975e-05, 9.941460381767523e-08), 
        },
        {
          cxdbl(-9.80022358636542e-08, -1.2135860550246832e-06), 
        },
        {
          cxdbl(0.0004843911986832152, 2.974301956984995e-11), 
        },
        {
          cxdbl(3.7991969547203735e-05, -1.5170427664939994e-09), 
        },
        {
          cxdbl(1.5179601475486182e-09, -1.2822469031258052e-10), 
        },
        {
          cxdbl(-3.920646589588768e-05, 9.943340051341728e-08), 
        },
        {
          cxdbl(3.799196954743034e-05, 1.517225150878216e-09), 
        },
        {
          cxdbl(0.030878066351274228, 1.0269651033665897e-10), 
        }
      };
      cout << "Differences:\n" << rho_eq - desired_rho_eq << endl;
      ASSERT_TRUE(allclose(rho_eq, desired_rho_eq, 1.0e-5, 1.0e-1));
    }
}


