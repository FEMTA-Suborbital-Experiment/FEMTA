
define enter start;
define enter waiting;
define leave apogee;
define leave done;
define leave chilly;

ds32 1Hz {
  
  if (State waiting, apogee | Time > 0s | single) {
    set pin 27 neg;
    enter done;
  }
  
  if (State waiting | Time > 0.125min | single, reverses) {
    set pin 17 neg;
    set pin 27 pos after 2s;
    leave start;
    enter apogee after 4s;
  }
  
  if (Temperature < 10raw) {
    enter chilly;
  }
  
  /* intentional errors - should produce meaningful messages
  if (Temperature > 90raw) {}                                    // empty if body
  if (State nonexistant | Time > 0raw) { enter done };           // invalid state
  if (State waiting, nonexistant | Time > 0raw) {enter done};    // invalid state
  if (Nonexistant > 90raw) { enter done };                       // invalid target
  if (Temperature < 90nil) { enter done };                       // invalid unit
  */
  
  [calibrate   | Time, poly, raw, s | 0.0009765625, 0.0];
  [conversions | Time, raw, s, ms   |                  ];
  //[print       | gray               |                  ];
}

ad15_gnd 1Hz {
  
  if (State apogee, done | A0 > 0.553mV | forever) {
    set pin 27 pos after 10000ms;
    leave done;
  }
  
  [calibrate   | A0, poly, raw, mV | 0.0001874287, -0.0009012627];
  [conversions | A0, raw, mV, V    |                            ];
  [print       | gray              |                            ];

  [smooth      | A0 | 0.1 ];
}
