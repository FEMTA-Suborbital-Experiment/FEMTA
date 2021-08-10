/*
** ADC Thermistor Test
**      author: Evan Rittner
**        date: August 10, 2021
** description:
**     Experiment to verify/debug the ADCs. The ADC ad15_gnd will read
**     the inputted analog data, convert it to temperature values, and 
**     print them to the terminal.
*/


define enter test;


Sensor ds32 1Hz {
    [calibrate   | Time, poly, raw, s | 0.0009765625, 0.0];
    [conversions | Time, raw, s       |                  ];
}

Sensor ad15_gnd 1Hz {
    [conversions | A0, raw, V, C    |                              ];
    [conversions | A1, raw, V, C    |                              ];
    [conversions | A2, raw, V, C    |                              ];
    [conversions | A3, raw, V, C    |                              ];

    [calibrate   | A0, poly, raw, V | 0.00041234314, -0.00198277794];
    [calibrate   | A1, poly, raw, V | 0.00041234314, -0.00198277794];
    [calibrate   | A2, poly, raw, V | 0.00041234314, -0.00198277794];
    [calibrate   | A3, poly, raw, V | 0.00041234314, -0.00198277794];

    [calibrate   | A0, poly, V, C   | -0.0000003, -0.000004, 0.0023, -0.1077, 1.9093];
    [calibrate   | A1, poly, V, C   | -0.0000003, -0.000004, 0.0023, -0.1077, 1.9093];
    [calibrate   | A2, poly, V, C   | -0.0000003, -0.000004, 0.0023, -0.1077, 1.9093];
    [calibrate   | A3, poly, V, C   | -0.0000003, -0.000004, 0.0023, -0.1077, 1.9093];

    [print       | red, A0          | 3                            ];
    [print       | maroon, A1       | 3                            ];
    [print       | salmon, A2       | 3                            ];
    [print       | pink, A3         | 3                            ];
}
