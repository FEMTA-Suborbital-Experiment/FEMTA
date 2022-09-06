/*
** Flow Test Fall 2021
**      author(s): Jake Valdez, Joe Kawiecki
**        date: October 22, 2021
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
// LargeSol1 -> Open:12  Close:16 Left Solenoid
// LargeSol2 -> Open:20  Close:21 Right Solenoid
// SmallSol1 -> Open:14  Close:15
// SmallSol2 -> Open:18  Close:23
// SmallSol3 -> Open:24  Close:25
// See flight.e doc for solenoid movements
// Open  solenoid:  Open pin +  ->  Close pin -  ->  Open pin  - after delay
// Close solenoid: Close pin +  ->  Open pin -   ->  Close pin - after delay


// Real time clock
Sensor ds32 10Hz {
    [calibrate   | Time, poly, raw, s | 0.0009765625, 0.0];
    [conversions | Time, raw, s       |                  ];
}

/*
// Left flow sensor
// Note: currently you can only use one flow sensor at a time
Sensor slf3_left 10Hz {  // Flow
    // Calibrate and convert Flow from raw to mL/min
    [calibrate   | Flow, poly, raw, mLpmin | 1.0, 0.0];
    [conversions | Flow, raw, mLpmin       |         ];

    // slf3 contains an internal calibration, so just use the raw temperature
    // from slf3
    // see slf3.c in blue-origin/sensors
    [calibrate   | Temperature, poly, raw, C | 1.0, 0.0];
    [conversions | Temperature, raw, C       |         ];

    // Similarly, detect if bubble in flow, like IR Flow sensor
    // No further calibrations
    [calibrate   | Bubble, poly, raw, V | 1.0, 0.0];
    [conversions | Bubble, raw, V       |         ];

    [print       | pink, Flow        | 3 ];
    [print       | pink, Temperature | 3 ];
    [print       | pink, Bubble      | 3 ];
}
*/

// Right flow sensor
// Note: currently you can only use one flow sensor at a time
Sensor slf3_right 10Hz {
    // Calibrate and convert Flow from raw to mL/min
    [calibrate   | Flow, poly, raw, mLpmin | 1.0, 0.0];
    [conversions | Flow, raw, mLpmin       |         ];

    // slf3 contains an internal calibration, so just use the raw temperature
    // from slf3
    // see slf3.c in blue-origin/sensors
    [calibrate   | Temperature, poly, raw, C | 1.0, 0.0];
    [conversions | Temperature, raw, C       |         ];

    // Similarly, detect if bubble in flow, like IR Flow sensor
    // No further calibrations
    [calibrate   | Bubble, poly, raw, V | 1.0, 0.0];
    [conversions | Bubble, raw, V       |         ];

    [calibrate   | HighFlow, poly, raw, V | 1.0, 0.0];
    [conversions | HighFlow, raw, V       |         ];

    [calibrate   | ExponentialSmoothing, poly, raw, V | 1.0, 0.0];
    [conversions | ExponentialSmoothing, raw, V       |         ];

    [print       | cyan, Flow        | 3 ];
    [print       | cyan, Temperature | 3 ];
    [print       | salmon, Bubble      | 3 ];
    [print       | yellow, HighFlow    | 3 ];
}

/* HSCM Sensor not defined on I2CBus
// TODO: Verify working syntax and test these at some point
Sensor hscm 100Hz {  // (Atmospheric) Pressure

    // HSCM sensor on FEMTA returns values in the range 0-15 psi. However,
    // Exp only supports kPa, psi, and torr. Convert to kPa by using
    // coefficients of polynomial
    // 1 psi = 6.89475728 kPa
    [calibrate   | Pressure, poly, raw, kPa | 6.89475728, 0.0];
    [conversions | Pressure, raw, kPa       |         ];

    // HSCM Sensor also includes a way to measure temperature
    // There is Temperature Conversion in hscm.c; no need to convert.
    [calibrate   | Temperature, poly, raw, mLpmin | 1.0, 0.0];
    [conversions | Temperature, raw, mLpmin       |         ];

    // Optional Temperature conversion that uses Exp poly to do the conversion,
    // rather existing conversion on hscm.c
    // First coefficient: 200/2047
    // Second coefficient: -50 offset
    // [calibrate   | Temperature, poly, raw, mLpmin | 0.09770395701, -50.0];
    // [conversions | Temperature, raw, mLpmin       |         ];

    [print       | mint, Pressure         | 3                         ];
    [print       | mint, Temperature         | 3                         ];
}
*/

/*
// atmospheric temperature sensor
Sensor mcp9 10Hz {
    [calibrate   | Temperature, poly, raw, C | 1.0, 0.0];
    [conversions | Temperature, raw, C       |         ];

    [print       | blue, Temperature | 3 ];
}
*/

/*
// IR Flow Sensors
// Note: plugging these sensors in PCBv3 causes the i2C bus to "scramble"
Sensor ad15_sda 10Hz {
  where A0 is flow_left;
  where A1 is flow_right;

  // IR flow sensor definitions
  // If IR sensors detect a flow, there should be voltage
  [calibrate | A0, poly, raw, V | 1, 0];
  [calibrate | A1, poly, raw, V | 1, 0];

  [conversions | A0, raw, V | ];
  [conversions | A1, raw, V | ];
  
  [print       | yellow, flow_left | 3];
  [print       | yellow, flow_right | 3];
}
*/

// Pressure Transducers
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
  if (State ascent | CC < 4.0kPa) {
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

    leave microgravity after 165s;
    enter descent after 600s;			//Observe time
  }


  // t=346s + 20min
  if (State descent | CC > 0kPa) {
    // Open for testing  all solenoids
    set pin 16 pos;
    set pin 16 neg after 20ms;

    set pin 20 pos;
    set pin 20 neg after 20ms;

    // 11/17/2021 Keep the vent valve closed.
    //set pin 14 neg;
    //set pin 15 pos;
    //set pin 14 pos after 20ms;
    //set pin 14 neg after 40ms;
    //set pin 15 neg after 40ms;

    leave descent;
    enter postflight after 600s; 		//Pump back out
  }


  // t=450s + 20min
  if (State postflight | CC > 0kPa) {

    leave postflight;
  }
}



/*
// TODO: Revisit thermistor calibrations after thermistor rewiring and
// possible recalibration

// Thermistor sensor
Sensor ad15_gnd 10Hz {
  where A0 is thermistor1; // for Left Propellant Tank    (Blue Tape 1)
  where A1 is thermistor2; // for Right Propellant Tank   (Blue Tape 2)
  where A2 is thermistor3; // for HFE Thermistor 1 (Back) (Blue Tape 3)
  where A3 is thermistor4; // for HFE Thermistor 2 (Top)  (Blue Tape 4)

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
    [calibrate   | thermistor1, poly, V, C   | 28893702.7627466,
                                             -87623.1450618716,
                                              67.0367407808325

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
