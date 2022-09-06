/*
 * Vishay NTCASCWE3 Test
 * 
 * Description: Implementation test of FEMTA's thermistor 
 * to determine functionality. NOTE - NOT A CALIBRATION TEST.
 *
 * Author: Vishal Ravi, Gouri Bellad, Mark Hartigan, Evan Rittner
 * Date: February 25 2021 (updated July 16 2021)
 */
 
define enter start;

Sensor ad15_gnd 1Hz { //vcc is placeholder adc. thermistors may be connected to different adc

    where A0 is thermistor1; //thermistors may be connected to different ports A0-A3 for this adc or they might be spread out among other adcs
    where A1 is thermistor2;
    where A2 is thermistor3;
    where A3 is thermistor4;
    
    [calibrate   | thermistor1, poly, raw, V | 1, 0 ];  //placeholder calibration curve, get from electronics, simbox electronics, or simbox software
    [calibrate   | thermistor1, poly, V, C   | -0.0000003,
                                               -0.000004,
                                                0.0023,
                                               -0.1077,
                                                1.9093
    ];
 
    [calibrate   | thermistor2, poly, raw, V | 1, 0 ];
    [calibrate   | thermistor2, poly, V, C   | -0.00002,
                                                0.0023,
                                               -0.106,
                                                2.003
    ];

    [calibrate   | thermistor3, poly, raw, V | 1, 0 ];
    [calibrate   | thermistor3, poly, V, C   | -0.00002,
                                                0.0024,
                                               -0.1074,
                                                1.9844
    ];
 
    [calibrate   | thermistor4, poly, raw, V | 1, 0 ];
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


Sensor ad15_sda 1Hz {

    where A0 is thermistor5; //thermistors may be connected to different ports A0-A3 for this adc or they might be spread out among other adcs
    where A1 is thermistor6;
    where A2 is thermistor7;
    where A3 is thermistor8;
 
    [calibrate   | thermistor5, poly, raw, V | 1, 0 ;
    [calibrate   | thermistor5, poly, V, C   | -0.00002,
                                                0.0025,
                                               -0.1003,
                                                1.8398
    ];

    [calibrate   | thermistor6, poly, raw, V | 1, 0 ];
    [calibrate   | thermistor6, poly, V, C   | -0.0001,
                                                0.0022,
                                               -0.1076,
                                                2.0743
    ];
    
    [calibrate   | thermistor7, poly, raw, V | 1, 0 ];
    [calibrate   | thermistor7, poly, V, C   | -0.00002,
                                                0.0025,
                                               -0.1031,
                                                1.8698
    ];
    
    [calibrate   | thermistor8, poly, raw, V | 1, 0 ];
    [calibrate   | thermistor8, poly, V, C   | -0.00002,
                                                0.0024,
                                               -0.106,
                                                1.9608
    ];
    
    [conversions | thermistor5, raw, V, C    |  ];
    [conversions | thermistor6, raw, V, C    |  ];
    [conversions | thermistor7, raw, V, C    |  ];
    [conversions | thermistor8, raw, V, C    |  ];
 
    // options are Flow, Temperature, or Bubble
    [print | red, thermistor5 | 1 ];
    [print | red, thermistor6 | 1 ];
    [print | red, thermistor7 | 1 ];
    [print | red, thermistor8 | 1 ];
}


Sensor ad15_vdd 1Hz {

    where A0 is thermistor9; //thermistors may be connected to different ports A0-A3 for this adc or they might be spread out among other adcs
    where A1 is thermistor10;
    where A2 is thermistor11;
    where A3 is thermistor12;

    [calibrate   | thermistor9, poly, raw, V | 1, 0 ];
    [calibrate   | thermistor9, poly, V, C   | -0.00002,
                                                0.0023,
                                               -0.1048,
                                                1.9794
    ];
 
    [calibrate   | thermistor10, poly, raw, V | 1, 0 ];
    [calibrate   | thermistor10, poly, V, C   | -0.00002,
                                                 0.0023,
                                                -0.1065,
                                                 2.0187
    ];
 
    [calibrate   | thermistor11, poly, raw, V | 1, 0 ];
    [calibrate   | thermistor11, poly, V, C   | -0.00002,
                                                 0.0022,
                                                -0.1075,
                                                 2.0676
    ];
 
    [calibrate   | thermistor12, poly, raw, V | 1, 0 ];
    [calibrate   | thermistor12, poly, V, C   | -0.00002,
                                                 0.0023,
                                                -0.1009,
                                                 1.9147
    ];
    
    [conversions | thermistor9, raw, V, C     |  ];
    [conversions | thermistor10, raw, V, C    |  ];
    [conversions | thermistor11, raw, V, C    |  ];
    [conversions | thermistor12, raw, V, C    |  ];

    // options are Flow, Temperature, or Bubble
    [print | red, thermistor9 | 1 ];
    [print | red, thermistor10 | 1 ];
    [print | red, thermistor11 | 1 ];
    [print | red, thermistor12 | 1 ];
}


Sensor ds32 1Hz { 
    [calibrate   | Time, poly, raw, s | 0.0009765625, 0.0];
    [conversions | Time, raw, s       |                  ];
}
