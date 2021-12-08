/*
** MCP9 Test Fall 2021
**      author(s): Jake Valdez
**        date: October 28, 2021
** description:
**     Test IR and Digital flow sensors only to conveniently measure
**     flow rate, temperature, and bubble
*/

define enter start


Sensor ds32 10Hz {
    [calibrate   | Time, poly, raw, s | 0.0009765625, 0.0];
    [conversions | Time, raw, s       |                  ];
}

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

    [print       | pink, Flow         | 3                         ];
    [print       | pink, Temperature     | 3                            ];
    [print       | pink, Bubble     | 3                            ];
}


// Comments for implementation above also concern "right" flow sensor
Sensor slf3_right 100Hz {  // Flow
    [calibrate   | Flow, poly, raw, mLpmin | 1.0, 0.0];
    [conversions | Flow, raw, mLpmin       |         ];

    [calibrate   | Temperature, poly, raw, C | 1.0, 0.0];
    [conversions | Temperature, raw, C       |         ];

    [calibrate   | Bubble, poly, raw, V | 1.0, 0.0];
    [conversions | Bubble, raw, V       |         ];

    [print       | cyan, Flow         | 3                         ];
    [print       | cyan, Temperature     | 3                            ];
    [print       | cyan, Bubble     | 3                            ];
}

/*
// IR Flow Sensors
Sensor ad15_sda 10Hz {
  where A0 is flow_left;
  where A1 is flow_right;

  // IR flow sensor definitions
  // If IR sensors detect a flow, there should be voltage
  [calibrate | A0, poly, raw, V | 1, 0];
  [calibrate | A1, poly, raw, V | 1, 0];

  [conversions | A0, raw, V | 		];
  [conversions | A1, raw, V | 		];
  [print       | mint, flow_left         | 3                         ];
  [print       | mint, flow_right     | 3                            ];
}
*/

// Pressure Transducers
Sensor ad15_vdd 10Hz {

  // See reference notes
  where A0 is propL; // Pressure sensor for the left propellant tank
  where A1 is propR; // Ditto for right propellant tank
  where A2 is CC; // Ditto for HFE
  //where A3 is HFE; // Ditto for Collection Chamber

  // From PCBv3 ADC Calibration under Electronics/Calibrations
  // Calibrations for ADC channels: A0, A1, A2, A3
  [calibrate   | A0, poly, raw, V | 0.00061138, -0.00240832];
  [calibrate   | A1, poly, raw, V | 0.00061158, -0.00045386];
  [calibrate   | A2, poly, raw, V | 0.00061155, -0.00029600];
  //[calibrate   | A3, poly, raw, V | 0.00061156, -0.00031852];


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
  //[calibrate | HFE, poly, V, kPa | 3.9580327, -0.03188930];

  [conversions | A0, raw, V, kPa  |                              ];
  [conversions | A1, raw, V, kPa  |                              ];
  [conversions | A2, raw, V, kPa  |                              ];
  //[conversions | A3, raw, V, kPa  |                              ];

  //[print       | red, HFE         | 3                            ];
  [print       | green, propL     | 3                            ];
  [print       | white, propR     | 3                            ];
  [print       | orange, CC       | 3                            ];

}

