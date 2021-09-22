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
// SmallSol2 -> Open:18  Close:23
// SmallSol3 -> Open:24  Close:25
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


Sensor ad15_vdd 10Hz {

  // Facing the collection chamber assembly from the propellant tank assembly..
  where A0 is propL; // Pressure sensor for the left propellant tank
  where A1 is propR; // Ditto for right propellant tank
  where A2 is HFE; // Ditto for HFE
  where A3 is CC; // Ditto for Collection Chamber
  // A3 and A0 varies depending on where pressure sensors connected, verify with
  // PCBv3 connections

  // Calibrations for ADC channels: A0, A1, A2, A3
  // 9/22/21: Arbitrary calibration curves
  [calibrate   | A0, poly, raw, V | 0.00041234314, -0.00198277794];
  [conversions | A0, raw, V, kPa  |                              ];
  [calibrate   | A1, poly, raw, V | 0.00041234314, -0.00198277794];
  [conversions | A1, raw, V, kPa  |                              ];
  [calibrate   | A2, poly, raw, V | 0.00041234314, -0.00198277794];
  [conversions | A2, raw, V, kPa  |                              ];
  [calibrate   | A3, poly, raw, V | 0.00041234314, -0.00198277794];
  [conversions | A3, raw, V, kPa  |                              ];


  // Linear Calibration Curves for propellant tanks

  // 9/22/2021 Calibration curves from
  // Calibrations/Pressure Sensors/Pressure Sensor Recalibration
  // DMP331 10581670
  [calibrate | propL, poly, V, kPa | 3.95, -0.0147];

  // Ditto for propR
  // DMP331 10581671
  [calibrate | propR, poly, V, kPa | 3.95, -1.32];

   // Ditto for HFE
   // DMP331 10581672
  [calibrate | HFE, poly, V, kPa | 3.95, 0.041];

    // Ditto for CC
    // DMP331 10581673
  [calibrate | CC, poly, V, kPa | 3.96, -0.0319];

   [print       | blue, propL         | 3                            ];
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

  // From thermistor.e under experiments in Evan branch
  // ADC Channel calibrations
  // 9/22/2021: Arbitrary calibration values
    [calibrate   | A0, poly, raw, V | 1, 0 ];  //placeholder calibration curve, get from electronics, simbox electronics, or simbox software
    [calibrate   | A1, poly, raw, V | 1, 0 ];
    [calibrate   | A2, poly, raw, V | 1, 0 ];
    [calibrate   | A3, poly, raw, V | 1, 0 ];


  // Thermistor calibrations taken from "Calibration Attempt 1 - Thermistors",
  // which is under Calibrations/Thermistors
  // From Eric Mesina's excel sheet, the coefficients correspond to the
  // coefficients of the polynomial slope under "Relevant Values and Graphs"
  // Should we use the coefficients from "Relevant Values and Graphs" or the
  // coefficients from "Linear Curve Approxs"?
  // Also, the calibration assumes voltage in the x-axis and temperature in the
  // y, but in the excel sheet, these two variables are switched.
    [calibrate   | thermistor1, poly, V, C   | -0.0000003,
                                               -0.000004,
                                                0.0023,
                                               -0.1077,
                                                1.9093
    ];

    [calibrate   | thermistor2, poly, V, C   | -0.00002,
                                                0.0023,
                                               -0.106,
                                                2.003
    ];


    [calibrate   | thermistor3, poly, V, C   | -0.00002,
                                                0.0024,
                                               -0.1074,
                                                1.9844
    ];

    [calibrate   | thermistor4, poly, V, C   | -0.00002,
                                                0.0025,
                                               -0.1009,
                                                1.828
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
