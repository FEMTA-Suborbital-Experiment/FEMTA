/*
** Flow Test Fall 2021
**      author(s): Jake Valdez, Joe Kawiecki
**        date: October 6th, 2021
** description:
**     Flow test involving propellant sending units, solenoids,
**     4 pressure sensors, 4 thermistors, and the collection chamber.
*/


// Define loiter states to idle and transition between states, only dependent on
// conditionals including sensor values and the current state
define enter preflight
define leave launch
define leave ascent
define leave meco
define leave microgravity
define leave descent
define leave postflight

// Reference Note:
/*
  Our reference frame is based on facing the collection chamber from
  the propellant tank assembly.
*/


//Notes/Tips:
// LargeSol1 -> Open:16  Close:12 Left Solenoid
// LargeSol2 -> Open:20  Close:21 Right Solenoid
// SmallSol1 -> Open:14  Close:15
// SmallSol2 -> Open:18  Close:23
// SmallSol3 -> Open:24  Close:25
// See flight.e doc for solenoid movements
// Open  solenoid: Open pin  +  ->  Open pin  - after delay of 20ms
// Close solenoid: Close pin +  ->  Close pin - after delay of 20ms


// following sensors from flight.e, on Evan's branch
Sensor ds32 10Hz {  // Real Time Clock
    [calibrate   | Time, poly, raw, s | 0.0009765625, 0.0];
    [conversions | Time, raw, s       |                  ];
}

/*
Sensor slf3_left 100Hz {  // Flow
    [calibrate   | Flow, poly, raw, mLpmin | 1.0, 0.0];
    [conversions | Flow, raw, mLpmin       |         ];
}

Sensor slf3_right 100Hz {  // Flow
    [calibrate   | Flow, poly, raw, mLpmin | 1.0, 0.0];
    [conversions | Flow, raw, mLpmin       |         ];
}
*/

/*
Sensor hscm 100Hz {  // (Atmospheric) Pressure

}

Sensor mcp9 100Hz {  // Temperature

}
*/

Sensor ad15_vdd 10Hz {

  // See reference notes
  where A0 is propL; // Pressure sensor for the left propellant tank
  where A1 is propR; // Ditto for right propellant tank
  where A2 is CC; // Ditto for HFE
  where A3 is HFE; // Ditto for Collection Chamber

  // From PCBv3 ADC Calibration under Electronics/Calibrations
  // Calibrations for ADC channels: A0, A1, A2, A3
  [calibrate   | A0, poly, raw, V | 0.00061138, -0.00240832];
  [calibrate   | A1, poly, raw, V | 0.00061158, -0.00045386];
  [calibrate   | A2, poly, raw, V | 0.00061155, -0.00029600];
  [calibrate   | A3, poly, raw, V | 0.00061156, -0.00031852];


  // Linear Calibration Curves for propellant tanks

  // 9/22/2021 Calibration curves from
  // Calibrations/Pressure Sensors/Pressure Sensor Recalibration

  // DMP331 10581670
  // A0
  [calibrate | propL, poly, V, kPa | 3.95405772, -0.01466985];

  // DMP331 10581671
  // A1
  [calibrate | propR, poly, V, kPa | 3.94825760, -1.31862661];

  // DMP331 10581672
  // A2
  [calibrate | CC, poly, V, kPa | 3.94928477, 0.04101970];

  // DMP331 11110585
  // A3
  [calibrate | HFE, poly, V, kPa | 3.9580327, -0.03188930];

  [conversions | A0, raw, V, kPa  |                              ];
  [conversions | A1, raw, V, kPa  |                              ];
  [conversions | A2, raw, V, kPa  |                              ];
  [conversions | A3, raw, V, kPa  |                              ];

  [print       | red, HFE         | 3                            ];
  [print       | green, propL     | 3                            ];
  [print       | white, propR     | 3                            ];
  [print       | orange, CC       | 3                            ];


  if (State preflight | CC > 0kPa) {

    // Close all flow solenoids at start
    set pin 12 pos;
    set pin 12 neg after 20ms;

    set pin 21 pos;
    set pin 21 neg after 20ms;

    // Open vent valve
    set pin 14 neg;
    set pin 15 pos;
    set pin 14 pos after 20ms;
    set pin 14 neg after 40ms;
    set pin 15 neg after 40ms;

    leave preflight;
    enter launch after 10s; 		//Time for vacuum chamber to set up
  }


  // t = 0s + 6min
  if (State launch | CC > 0kPa) {
    // Open the collection chamber (small solenoid) on launch
    //set pin 14 pos;
    //set pin 15 neg;
    //set pin 14 neg after 350ms;


    leave launch;
    enter ascent after 600s;
  }

  //t = 8s + 20min
  if (State ascent | CC < 1.5kPa) {
    leave ascent;
    enter meco after 5s;
  }


  // Action state: Open big sol valves at meco, close small solenoid
  // t = 141s + 20min
  if (State meco | CC > 0kPa) {

    // Close the vent valve
    set pin 14 pos;
    set pin 15 neg;
    set pin 14 neg after 20ms;

    // Open two flow valves at the same time after small delay
    set pin 16 pos;
    set pin 16 neg after 20ms;

    set pin 20 pos;
    set pin 20 neg after 20ms;

    leave meco;
    enter microgravity after 10s;
  }


  // t = 181s + 20min
  if (State microgravity | CC > 0kPa) {

    leave microgravity;
    enter descent after 600s;			//Observe time
  }


  // t=346s + 20min
  if (State descent | CC > 0kPa) {
    // Open for testing  all solenoids
    set pin 16 pos;
    set pin 16 neg after 20ms;

    set pin 20 pos;
    set pin 20 neg after 20ms;

    set pin 14 neg;
    set pin 15 pos;
    set pin 14 pos after 20ms;
    set pin 14 neg after 40ms;
    set pin 15 neg after 40ms;

    leave descent;
    enter postflight after 600s; 		//Pump back out
  }


  // t=450s + 20min
  if (State postflight | CC > 0kPa) {

    leave postflight;
  }
}


/*
// Thermistor sensor
Sensor ad15_gnd 10Hz {
  where A0 is thermistor1 // for Left Propellant Tank
  where A1 is thermistor2 // for Right Propellant Tank
  where A2 is thermistor3 // for HFE Thermistor 1 (Back)
  where A3 is thermistor4 // for HFE Thermistor 2 (Top)

  // From PCBv3 ADC Calibration under Electronics/Calibrations
  // ADC Channel calibrations
    [calibrate   | A0, poly, raw, V | 0.000187760043219144,
                                      0.0000186149452909308
    ];
    [calibrate   | A1, poly, raw, V | 0.000187759867294351,
                                      0.0000207491490691839
    ];
    [calibrate   | A2, poly, raw, V | 0.000187761014729653,
                                      0.0000193254277389214
    ];
    [calibrate   | A3, poly, raw, V | 0.000187750489295807,
                                      0.000013227022497931
    ];


  // Thermistor calibrations taken from "Voltage to Temperature Calibrations",
  // which is under Calibrations/Thermistors
    [calibrate   | thermistor1, poly, V, C   | 0.0000000077,
                                               -0.0000009762,
                                                0.0000469620,
                                               -0.0011800366,
                                                0.00181440736
    ];

    [calibrate   | thermistor2, poly, V, C   | 26939605.0709972,
                                               -84025.306425667,
                                               65.792989370689
    ];


    [calibrate   | thermistor3, poly, V, C   | 29992593.4712639,
                                               -89247.9017336437,
                                               67.2483882781364
    ];

    [calibrate   | thermistor4, poly, V, C   | 29664639.0926322,
                                               -88686.8609090078,
                                               67.2298598866639
    ];

    [conversions | thermistor1, raw, V, C    |  ];
    [conversions | thermistor2, raw, V, C    |  ];
    [conversions | thermistor3, raw, V, C    |  ];
    [conversions | thermistor4, raw, V, C    |  ];

    // options are Flow, Temperature, or Bubble
    [print | white, thermistor1 | 1 ];
    [print | white, thermistor2 | 1 ];
    [print | white, thermistor3 | 1 ];
    [print | white, thermistor4 | 1 ];
}
*/
