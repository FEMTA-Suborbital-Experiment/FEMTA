/**
 * @Authors Noah
 * @Compiler_Version Release 1.0
 * @Licence GPLv3
 * 
 * @Description Example with all features in release 1.0
 **/


define enter pad;        // rocket sitting on pad
define leave ascent;     // rocket increasing in altitude
define leave appogee;    // rocket at the edge of atmosphere
define leave descent;    // rocket falling back to earth


ad15_gnd 5Hz {
  
  if (A0 < 10Kpa) {
    set broadcom 17 neg;
    set broadcom 27 pos;
    set broadcom 27 pos after 250ms;    // we take power off so as to avoid burning the solenoid
    leave pad;
    enter ascent after 250ms;
  }
  
  if (A1 > 10torr | reverses, singular) {
    connect broadcom 17 to source;
    leave appogee;
    enter descent;
  }
  
  [calibrate | A0, poly, V, kPa |    // we expect 11.3v at 1atm | V->kPa
     -0.00000805518095144006,        // degrees are from high to low
      0.000179617025688186,
     -0.00149263466053365,
      0.0137083951784765,
      3.84928878749102,
     -0.389869854380195
   ]
  
  [calibrate   | A0, poly, raw, V | 0.0001874287, -0.0009012627]
  [conversions | A0, raw, V, kPa  |                            ]
  [print       | gray             |                            ]
}

ds32 1Hz {
  [calibrate   | Time, poly, raw, s | 0.0009765625, 0.0]
  [conversions | Time, raw, s       |                  ]
  [print       | gray               |                  ]
}
