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
define leave ascent
define leave microgravity
define leave descent
define leave postflight

// Define action states to enter and immediately leave the state when a
// conditional is met. Only dependent on state only
define enter start
define leave microgravity_end
define leave launch
define leave cleanup

// Notes and questions for state pairs
// Define state transitions at the end of action states?
//


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


Sensor ad15_vdd 1Hz {   // Why is this at 1Hz? How did Jake decide that?
  where A0 is HFE1
  where A1 is prop
  where A2 is CC
  where A3 is HFE2
  // A3 and A0 varies depending on where pressure sensors connected, verify with
  // PCBv3 connections

  [calibrate   | HFE, poly, raw, V | 0.00041234314, -0.00198277794];
  [conversions | HFE, raw, V, kPa  |                              ];

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
   [print       | red, prop          | 3                            ];
   [print       | green, CC       | 3                            ];
   [print       | yellow,HFE2         | 3                            ];

  if (State start | CC > 0kPa) {
  // Enter phases of launch, I added
    enter launch
    leave launch after 8s

    //enter ascent after 8s
    //leave ascent after 181s

    //enter meco after 141s
    //leave meco after 181s

    //enter microgravity after 181s
    //leave microgravity after 346s

    //enter descent after 346s
    //leave descent after 450s

    //enter postflight after 450s
    //leave postflight after 500s

    leave start
  }

  if (State launch | CC > 0kPa) {
    // Open the vent valve on launch
    set pin 14 pos;
    set pin 15 neg;
    set pin 14 neg after 350ms;

    leave launch;
    enter ascent after 8s;
    enter meco after 141s;
  }

  // From Mark Hartigan's full-flow.e
  // Don't think we do anything in ascent, right?
  // Loiter state with pressure sensor condition
  if (State ascent | CC < 0.5kPa) {
    // Close the vent valve
    // set pin 15 pos;
    // set pin 14 neg;
    // set pin 15 neg after 350ms;

    // Open solenoid after 60 seconds
    // set pin 12 pos after 60000ms;
    // set pin 21 neg after 60000ms;
    // set pin 12 neg after 60350ms;

    //enter test after 30s;  // end test after a half minute

    leave ascent;
  }


  // Conditional?
  // Don't think we do anything in meco either right?
  // Action state: Open and close valves at meco, leave meco immediately
  // t ~ 141s, so enter microgravity at microgravity_time-141
  if (State meco) {

    enter microgravity after 40s
    leave meco;
  }


  // Conditional?
  if (State microgravity) {
    // This will evaluate 165 times
    // Close the vent valve
    set pin 15 pos;
    set pin 14 neg;
    set pin 15 neg after 350ms;

    // Open flow/run valves at start of microgravity
    set pin 12 pos;
    set pin 16 neg;
    set pin 12 neg after 350ms;

    set pin 21 pos;
    set pin 20 neg;
    set pin 21 neg after 350ms;

    leave microgravity;
    // Should this evaluate as 346s from t=0, or as 165s from microgravity/meco?
    enter descent after 346s;
  }


  // Conditional?
  if (State descent) {
    // Close flow/run valves at end of microgravity (=start of descent)
    // For operating solenoids, is it the same syntax for both open and close?
    set pin 16 pos;
    set pin 12 neg;
    set pin 16 neg after 350ms;

    set pin 21 pos;
    set pin 20 neg;
    set pin 21 neg after 350ms;

    leave descent;
    enter postflight after 450s;
  }


  // Conditional?
  // Don't think we do anything in postflight right?
  if (State postflight) {

    leave postflight;
    enter cleanup after 500s;
  }

  if (State cleanup | CC > 0kPa) {
    // Close all solenoids
    set pin 15 pos;
    set pin 14 neg;
    set pin 15 neg after 350ms;

    set pin 12 pos;
    set pin 16 neg;
    set pin 12 neg after 350ms;

    set pin 21 pos;
    set pin 20 neg;
    set pin 21 neg after 350ms;

    leave cleanup;
  }
}



// Thermistor sensor
Sensor ad15_gnd 1Hz {
  //TODO: determine the usual thermistor order (HFE, Propellant1, propellant2 etc.)
  where A0 is thermistor1 // for ..
  where A1 is thermistor2 // for ..
  where A2 is thermistor3 // for ..
  where A3 is thermistor4 // for ..

}
