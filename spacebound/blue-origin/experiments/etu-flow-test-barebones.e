/*
** Flow Test Fall 2021
**      author(s): Jake Valdez, Joe Kawiecki
**        date: September 10, 2021
** description:
**     Flow test involving propellant sending units, solenoids,
**     4 pressure sensors, and the collection chamber.
*/


// Define loiter states to idle and transition between states, only dependent on
// conditionals including sensor values and the current state
define leave preflight
define leave launch
define leave ascent
define leave microgravity
define leave descent
define leave postflight


//Notes/Tips:
// LargeSol1 -> Open:12  Close:16
// LargeSol2 -> Open:20  Close:21
// SmallSol1 -> Open:14  Close:15
// See flight.e doc for solenoid movements
// Open  solenoid: Open pin +  ->  Close pin -  ->  Open pin + after delay
// Close solenoid: Close pin +  ->  Open pin -  ->  Close pin - after delay


// Following sensors from flight.e, on Evan's branch
Sensor ds32 1Hz {  // Real Time Clock
    [calibrate   | Time, poly, raw, s | 0.0009765625, 0.0];
    [conversions | Time, raw, s       |                  ];
}

Sensor slf3_left 100Hz {  // Flow
    [calibrate   | Flow, poly, raw, ml/min | 1.0, 0.0];
    [conversions | Flow, raw, ml/min       |         ];
}

Sensor slf3_right 100Hz {  // Flow
    [calibrate   | Flow, poly, raw, ml/min | 1.0, 0.0];
    [conversions | Flow, raw, ml/min       |         ];
}

Sensor hscm 100Hz {  // (Atmospheric) Pressure

}

Sensor mcp9 100Hz {  // Temperature

}


Sensor ad15_vdd 10Hz {   // Why is this at 1Hz? How did Jake decide that?
  where A0 is HFE1;
  where A1 is prop;
  where A2 is CC;
  where A3 is HFE2;
  // A3 and A0 varies depending on where pressure sensors connected, verify with
  // PCBv3 connections

  [calibrate   | HFE1, poly, raw, V | 0.00041234314, -0.00198277794];
  [conversions | HFE1, raw, V, kPa  |                              ];
  [calibrate   | HFE2, poly, raw, V | 0.00041234314, -0.00198277794];
  [conversions | HFE2, raw, V, kPa  |                              ];

  [calibrate | prop, poly, V, kPa |
    -0.00000805518095144006,
     0.000179617025688186,
    -0.00149263466053365,
     0.0137083951784765,
     3.84928878749102,
    -0.389869854380195
   ];

  [calibrate | CC, poly, V, kPa |
    -0.00000805518095144006,
     0.000179617025688186,
    -0.00149263466053365,
     0.0137083951784765,
     3.84928878749102,
    -0.389869854380195
   ];

   [print       | blue, HFE1         | 3                            ];
   [print       | green, CC       | 3                            ];


  if (State preflight| CC > 0kPa) {
    // Close all solenoids at start
    set pin 16 pos;
    set pin 12 neg;
    set pin 16 neg after 350ms;

    set pin 21 pos;
    set pin 20 neg;
    set pin 21 neg after 350ms;

    set pin 15 pos;
    set pin 14 neg;
    set pin 15 neg after 350ms;

    leave prefligt;
    enter launch after 1200s; 		//Time for vacuum chamber to set up
  }


  // t = 0s + 20min
  if (State launch) {
    // Open the collection chamber (small solenoid) on launch
    set pin 14 pos;
    set pin 15 neg;
    set pin 14 neg after 350ms;

    leave launch;
    enter ascent after 8s;
  }

  //t = 8s + 20min
  if (State ascent | CC < 0.5kPa) {
    leave ascent;
    enter meco after 133s;
  }


  // Action state: Open big sol valves at meco, close small solenoid
  // t = 141s + 20min
  if (State meco) {
    
    // Close the vent valve
    set pin 15 pos;
    set pin 14 neg;
    set pin 15 neg after 350ms;

    // Open flow valve after small delay
    set pin 12 pos;
    set pin 16 neg;
    set pin 12 neg after 5s;

    set pin 20 pos;
    set pin 21 neg;
    set pin 20 neg after 5s;

    leave meco;
    enter microgravity after 40s;
  }


  // t = 181s + 20min
  if (State microgravity) {

    leave microgravity;
    enter descent after 165s;
  }


  // t=346s + 20min
  if (State descent) {
    // Close all solenoids
    set pin 16 pos;
    set pin 12 neg;
    set pin 16 neg after 350ms;

    set pin 21 pos;
    set pin 20 neg;
    set pin 21 neg after 350ms;

    set pin 15 pos;
    set pin 14 neg;
    set pin 15 neg after 350ms;

    leave descent;
    enter postflight after 104s;
  }


  // t=450s + 20min
  if (State postflight) {

    leave postflight;
  }
}


// Thermistor sensor
Sensor ad15_gnd 1Hz {
  //TODO: determine the usual thermistor order (HFE, Propellant1, propellant2 etc.)
  where A0 is thermistor1 // for ..
  where A1 is thermistor2 // for ..
  where A2 is thermistor3 // for ..
  where A3 is thermistor4 // for ..

  // From adc_therm-test.e under experiments
  [conversions | A0, raw, V, C    |                              ];
  [conversions | A1, raw, V, C    |                              ];
  [conversions | A2, raw, V, C    |                              ];
  [conversions | A3, raw, V, C    |                              ];

  [calibrate   | A0, poly, raw, V | 0.00041234314, -0.00198277794];
  [calibrate   | A1, poly, raw, V | 0.00041234314, -0.00198277794];
  [calibrate   | A2, poly, raw, V | 0.00041234314, -0.00198277794];
  [calibrate   | A3, poly, raw, V | 0.00041234314, -0.00198277794];

  [calibrate   | A0, poly, V, C   | -0.0000003, -0.000004, 0.0023, -0.1077, 1.9093];
  [calibrate   | A1, poly, V, C   | -0.0000003, -0.000004, 0.0023, -0.1077, 1.9093];
  [calibrate   | A2, poly, V, C   | -0.0000003, -0.000004, 0.0023, -0.1077, 1.9093];
  [calibrate   | A3, poly, V, C   | -0.0000003, -0.000004, 0.0023, -0.1077, 1.9093];

}
