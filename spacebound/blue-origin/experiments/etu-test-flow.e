/*
** Flow Test Fall 2021
**      author(s): Jake Valdez, Joe Kawiecki
**        date: September 7, 2021
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
define leave completed
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


Sensor ad15_vdd 1Hz {
  where A0 is HFE
  where A1 is prop
  where A2 is CC
  where A3 is

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
  // Enter phases of launch
    enter launch
    leave launch after 8s

    enter ascent after 8s
    leave ascent after 181s

    enter meco after 141s
    leave meco after 181s

    enter microgravity after 181s
    leave microgravity after 346s

    enter descent after 346s
    leave descent after 450s

    enter postflight after 450s
    leave postflight after 500s

    leave start
    enter cleanup after 505s
  }

  // From Mark Hartigan's full-flow.e
  if (State ascent | collection < 0.5kPa) {
    // close the vent valve
    set pin 15 pos;
    set pin 14 neg;
    set pin 15 neg after 350ms;

    // open solenoid after 60 seconds
    set pin 12 pos after 60000ms;
    set pin 21 neg after 60000ms;
    set pin 12 neg after 60350ms;

    leave ascent;
    enter test after 30s;  // end test after a half minute
  }

  if (State meco) {
    // This will evaluate 40 times
  }

  if (State microgravity) {
    // This will evaluate 165 times
  }

  if (State cleanup | CC > 0kPa) {
    // Close solenoids


  }
}
