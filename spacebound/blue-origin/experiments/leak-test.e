
define enter start;
define leave descent;
define leave rising;

Sensor ad15_vdd 5Hz {

  where A0 is HFE;
    
  [calibrate | HFE, poly, V, kPa |
    -0.00000805518095144006,
     0.000179617025688186,
    -0.00149263466053365,
     0.0137083951784765,
     3.84928878749102,
    -0.389869854380195
   ];
  
  [calibrate   | HFE, poly, raw, V | 0.00041234314, -0.00198277794];
  [conversions | HFE, raw, V, kPa  |                              ];
  [print       | blue, HFE         | 3                            ];

  where A1 is prop1;
  where A2 is prop2;
  
  [calibrate | prop1, poly, V, kPa |
    -0.00000805518095144006,
     0.000179617025688186,
    -0.00149263466053365,
     0.0137083951784765,
     3.84928878749102,
    -0.389869854380195
   ];
  
  [calibrate | prop2, poly, V, kPa |
    -0.00000805518095144006,
     0.000179617025688186,
    -0.00149263466053365,
     0.0137083951784765,
     3.84928878749102,
    -0.389869854380195
   ];

  [calibrate   | prop1, poly, raw, V | 0.00041234314, -0.00198277794];
  [calibrate   | prop2, poly, raw, V | 0.00041234314, -0.00198277794];
  [conversions | prop1, raw, V, kPa  |                              ];
  [conversions | prop2, raw, V, kPa  |                              ];
  [print       | salmon, prop1       | 3                            ];
  [print       | orange, prop2       | 3                            ];

}

Sensor ad15_sda 20Hz {
  
  
  /*
  if (State start | HFE > 0kPa) {

    // open solenoid 2
    // set pin 20 pos;
    // set pin 21 neg;
    // set pin 20 neg after 350ms;

    // open the vent valve
    // set pin 17 pos;
    // set pin  4 neg;
    // set pin 17 neg after 350ms;
    
    // open solenoid 1
    set pin 26 pos;
    set pin 19 neg;
    set pin 26 neg after 350ms;
    
    leave start;
    enter descent after 1s;
  }
  
  if (State descent | HFE < 1kPa) {
    
    // close the vent valve
    // set pin  4 pos;
    // set pin 17 neg;
    // set pin  4 neg after 350ms;
    
    // close solenoid 2
    // set pin 21 pos;
    // set pin 20 neg;
    // set pin 21 neg after 350ms;
    
    leave descent;
    enter rising after 60s;
  }
  
  if (State rising | HFE > 0kPa) {
    
    // open the vent valve
    set pin 17 pos;
    set pin  4 neg;
    set pin 17 neg after 350ms;
    
    // open solenoid 2
    set pin 20 pos;
    set pin 21 neg;
    set pin 20 neg after 350ms;
    
    leave rising;
  }
  */
  
}

Sensor ds32 1Hz {
  [calibrate   | Time, poly, raw, s | 0.0009765625, 0.0]
  [conversions | Time, raw, s       |                  ]
}
