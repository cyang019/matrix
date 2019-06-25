namespace matrix { inline namespace v1 {
  namespace exponential {
    constexpr double c(int i)
    {
      switch(i){
        case 7: return 9.92063492063492e-06;
        case 11: return 9.94131285136576e-11;
        case 15: return 2.22819456055356e-16;
        case 19: return 1.69079293431187e-22;
        case 23: return 5.48340615485360e-29;
        case 27: return 8.82996160201868e-36;
        default: return 0.0;
      }
    }

    constexpr double b(int m, int i)
    {
      switch(m){
        case 3:
          switch(i){
            case 0: return 120.0;
            case 1: return 60.0;
            case 2: return 12.0;
            case 3: return 1.0;
            default: return 0;
          }
          break;
        case 5:
          switch(i){
            case 0: return 30240.0;
            case 1: return 15120.0;
            case 2: return 3360.0;
            case 3: return 420.0;
            case 4: return 30.0;
            case 5: return 1.0;
            default: return 0;
          }
          break;
        case 7:
          switch(i){
            case 0: return 17297280.0;
            case 1: return 8648640.0;
            case 2: return 1995840.0;
            case 3: return 277200.0;
            case 4: return 25200.0;
            case 5: return 1512.0;
            case 6: return 56.0;
            case 7: return 1.0;
            default: return 0;
          }
          break;
        case 9:
          switch(i){
            case 0: return 17643225600.0;
            case 1: return 8821612800.0;
            case 2: return 2075673600.0;
            case 3: return 302702400.0;
            case 4: return 30270240.0;
            case 5: return 2162160.0;
            case 6: return 110880.0;
            case 7: return 3960.0;
            case 8: return 90.0;
            case 9: return 1.0;
            default: return 0;
          }
          break;
        case 11:
          switch(i){
            case 0: return 281585880576 * 100.0;
            case 1: return 140792940288 * 100.0;
            case 2: return 3352212864 * 1000.0;
            case 3: return 502831929600.0;
            case 4: return 52929676800.0;
            case 5: return 4116752640.0;
            case 6: return 242161920.0;
            case 7: return 10810800.0;
            case 8: return 360360.0;
            case 9: return 8580.0;
            case 10: return 132.0;
            case 11: return 1.0;
            default: return 0;
          }
          break;
        case 13:
          switch(i){
            case 0: return 6476475253248 * 10000.0;
            case 1: return 3238237626624 * 10000.0;
            case 2: return 77717703038976 * 100.0;
            case 3: return 11873537964288 * 100.0;
            case 4: return 129060195264 * 1000.0;
            case 5: return 105594705216 * 100.0;
            case 6: return 670442572800.0;
            case 7: return 33522128640.0;
            case 8: return 1323241920.0;
            case 9: return 40840800.0;
            case 10: return 960960.0;
            case 11: return 16380.0;
            case 12: return 182.0;
            case 13: return 1.0;
            default: return 0;
          }
          break;
        default: 
          return 0.0;
      }
    } // b()

    constexpr double theta_3 = 1.495585217958292e-2;
    constexpr double theta_5 = 2.539398330063230e-1;
    constexpr double theta_7 = 9.504178996162932e-1;
    constexpr double theta_9 = 2.097847961257068e0;
    constexpr double theta_13 = 4.25;

    constexpr double q3_upper = 1.0e0;
    constexpr double q5_upper = 1.3e0;
    constexpr double q7_upper = 2.6e0;
    constexpr double q9_upper = 8.2e0;
    constexpr double q13_upper = 7.1e1;
  } // namespace exponential
} // namespace v1
} // namespace matrix
