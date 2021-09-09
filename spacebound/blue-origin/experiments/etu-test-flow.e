/*
** Flow Test Fall 2021
**      author(s): Jake Valdez, Joe Kawiecki
**        date: September 10, 2021
** description:
**     Flow test involving propellant sending units, solenoids,
**     4 pressure sensors, and the collection chamber.
*/

define enter start
define leave launch
define leave ascent
define leave meco
define leave microgravity
define leave descent
define leave postflight
//define leave completed
define leave cleanup


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
  where A0 is HFE
  where A1 is prop
  where A2 is CC
  where A3 is A3        // What is A3 supposed to be?

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

   [print       | blue, HFE         | 3                            ];
   [print       | red, prop          | 3                            ];
   [print       | maroon, CC       | 3                            ];
   [print       | pink, A3         | 3                            ];

  if (State start | CC > 0kPa) {
  // Enter phases of launch, I added 
    //enter launch
    //leave launch after 8s

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
    enter cleanup after 505s
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
  if (State ascent | collection < 0.5kPa) {
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
    enter microgravity after 181s;
  }


  // Conditional?
  // Don't think we do anything in meco either right?
  if (State meco) {
    
    leave meco;
    enter microgravity after 181s;
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
