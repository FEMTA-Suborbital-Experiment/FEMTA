/*
 * Honeywell HSCMNNN015PA2A3 Test
 *
 * Description: Test of board mounted pressure sensor to confirm functionality
 *
 * Author: Mark Hartigan
 * Date: June 10, 2021
 */

define enter start;

Sensor hscm 1Hz {

  [calibrate   | Temperature, poly, raw, C | 1, 0];
  [conversions | Temperature, raw, C       |     ];

  [print | salmon, Pressure  | 0];
  [print | mint, Temperature | 3];
}

Sensor ad15_vdd 1Hz {

  [calibrate | A0, poly, V, kPa |
   -0.00000805518095144006,
    0.000179617025688186,
   -0.00149263466053365,
    0.0137083951784765,
    3.84928878749102,
   -0.389869854380195
  ];

  [calibrate   | A0, poly, raw, V | 0.00041234314, -0.00198277794];
  [conversions | A0, raw, V, kPa  |                              ];
  [print       | blue, A0         | 3                            ];
}
