/*
** All Analog to Digital Converter Test
**      author: Evan Rittner
**        date: June 29, 2021
** description:
**     Experiment to verify/debug the ADCs. Each of the four ADCs will
**     output readings from all four of their channels. This can be used
**     to enusre the ADCs work, and to see which pins connect to which ADC.
**     Outputs are printed in volts (as opposed to a converted unit like kPa).
*/


define enter test;


Sensor ds32 1Hz {
    [calibrate   | Time, poly, raw, s | 0.0009765625, 0.0];
    [conversions | Time, raw, s       |                  ];
}

Sensor ad15_gnd 1Hz {
    [conversions | A0, raw, V       |                              ];
    [conversions | A1, raw, V       |                              ];
    [conversions | A2, raw, V       |                              ];
    [conversions | A3, raw, V       |                              ];

    [calibrate   | A0, poly, raw, V | 0.00041234314, -0.00198277794];
    [calibrate   | A1, poly, raw, V | 0.00041234314, -0.00198277794];
    [calibrate   | A2, poly, raw, V | 0.00041234314, -0.00198277794];
    [calibrate   | A3, poly, raw, V | 0.00041234314, -0.00198277794];
    
    [print       | red, A0          | 3                            ];
    [print       | maroon, A1       | 3                            ];
    [print       | salmon, A2       | 3                            ];
    [print       | pink, A3         | 3                            ];
}

Sensor ad15_vdd 1Hz {
    [conversions | A0, raw, V       |                              ];
    [conversions | A1, raw, V       |                              ];
    [conversions | A2, raw, V       |                              ];
    [conversions | A3, raw, V       |                              ];

    [calibrate   | A0, poly, raw, V | 0.00041234314, -0.00198277794];
    [calibrate   | A1, poly, raw, V | 0.00041234314, -0.00198277794];
    [calibrate   | A2, poly, raw, V | 0.00041234314, -0.00198277794];
    [calibrate   | A3, poly, raw, V | 0.00041234314, -0.00198277794];
    
    [print       | blue, A0         | 3                            ];
    [print       | cyan, A1         | 3                            ];
    [print       | purple, A2       | 3                            ];
    [print       | brown, A3        | 3                            ];
}

Sensor ad15_sda 1Hz {
    [conversions | A0, raw, V       |                              ];
    [conversions | A1, raw, V       |                              ];
    [conversions | A2, raw, V       |                              ];
    [conversions | A3, raw, V       |                              ];

    [calibrate   | A0, poly, raw, V | 0.00041234314, -0.00198277794];
    [calibrate   | A1, poly, raw, V | 0.00041234314, -0.00198277794];
    [calibrate   | A2, poly, raw, V | 0.00041234314, -0.00198277794];
    [calibrate   | A3, poly, raw, V | 0.00041234314, -0.00198277794];
    
    [print       | green, A0        | 3                            ];
    [print       | mint, A1         | 3                            ];
    [print       | yellow, A2       | 3                            ];
    [print       | orange, A3       | 3                            ];
}

Sensor ad15_scl 1Hz {
    [conversions | A0, raw, V       |                              ];
    [conversions | A1, raw, V       |                              ];
    [conversions | A2, raw, V       |                              ];
    [conversions | A3, raw, V       |                              ];

    [calibrate   | A0, poly, raw, V | 0.00041234314, -0.00198277794];
    [calibrate   | A1, poly, raw, V | 0.00041234314, -0.00198277794];
    [calibrate   | A2, poly, raw, V | 0.00041234314, -0.00198277794];
    [calibrate   | A3, poly, raw, V | 0.00041234314, -0.00198277794];
    
    [print       | magenta, A0      | 3                            ];
    [print       | gray, A1         | 3                            ];
    [print       | red, A2          | 3                            ];
    [print       | pink, A3         | 3                            ];
}