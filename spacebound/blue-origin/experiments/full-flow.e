/*
** Flow Test Summer 2021
**      author: Mark Hartigan
**        date: June 9, 2021
** description:
**     Flow test involving propellant sending units, solenoids,
**     3 pressure sensors, and the collection chamber.
*/

define enter start;
define leave ascent;
define leave test;
define leave descent;


Sensor ds32 1Hz {
  [calibrate   | Time, poly, raw, s | 0.0009765625, 0.0];
  [conversions | Time, raw, s       |                  ];
}

Sensor slf3_left 5Hz {

  [calibrate   | Flow, poly, raw, mLpmin | 1, 0];
  [conversions | Flow, raw, mLpmin       |     ];
  [print       | white, Flow             | 3   ];

}

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

  where A1 is prop;
  where A2 is collection;

  [calibrate | prop, poly, V, kPa |
    -0.00000805518095144006,
     0.000179617025688186,
    -0.00149263466053365,
     0.0137083951784765,
     3.84928878749102,
    -0.389869854380195
   ];

  [calibrate | collection, poly, V, kPa |
    -0.00000805518095144006,
     0.000179617025688186,
    -0.00149263466053365,
     0.0137083951784765,
     3.84928878749102,
    -0.389869854380195
   ];

  [calibrate   | prop, poly, raw, V       | 0.00041234314, -0.00198277794];
  [calibrate   | collection, poly, raw, V | 0.00041234314, -0.00198277794];
  [conversions | prop, raw, V, kPa        |                              ];
  [conversions | collection, raw, V, kPa  |                              ];
  [print       | salmon, prop             | 3                            ];
  [print       | mint, collection         | 3                            ];

  if (State start | collection > 0kPa) {

    // close solenoid 2
    // set pin 21 pos;
    // set pin 20 neg;
    // set pin 21 neg after 350ms;

    // open the vent valve
    set pin 14 pos;
    set pin 15 neg;
    set pin 14 neg after 350ms;

    // close solenoid
    set pin 21 pos;
    set pin 12 neg;
    set pin 21 neg after 350ms;

    leave start;
    enter ascent after 1s;
  }

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

  if (State test | collection > 0kPa) {
    // close solenoid
    set pin 21 pos;
    set pin 12 neg;
    set pin 21 neg after 350ms;

    leave test;
    enter descent after 1s;
  }
}
