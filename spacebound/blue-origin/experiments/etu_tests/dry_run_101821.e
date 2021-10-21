/*
** Flow Test Fall 2021
**      author(s): Jake Valdez, Joe Kawiecki, Evan Rittner
**        date: October 18th, 2021
** description:
**     Flow test involving propellant sending units, solenoids,
**     4 pressure sensors, 4 thermistors, and the collection chamber.
*/


// Define loiter states to idle and transition between states, only dependent on
// conditionals including sensor values and the current state
define enter test
define leave posttest

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

Sensor ad15_sda 10Hz {
  where A0 is flow_left;
  where A1 is flow_right;

  // IR flow sensor definitions
  // If IR sensors detect a flow, there should be voltage
  [calibrate | A0, poly, raw, V | 1, 0];
  [calibrate | A1, poly, raw, V | 1, 0];

  [conversions | A0, raw, V | 		];
  [conversions | A1, raw, V | 		];
}


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


  // Action state: Open big sol valves at meco, close small solenoid
  // t = 141s + 20min

  /*
  if (State test | CC > 0kPa) {

    // Close the vent valve
    set pin 15 pos;
    set pin 14 neg;
    set pin 15 neg after 20ms;

    leave test;
    enter posttest after 10s;
  }
  */

  /*
  // t=450s + 20min
  if (State posttest | CC > 0kPa) {
  }
  */
}

