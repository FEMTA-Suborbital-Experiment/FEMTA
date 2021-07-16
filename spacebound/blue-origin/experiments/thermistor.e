/*
 * Vishay NTCASCWE3 Test
 * 
 * Description: Implementation test of FEMTA's thermistor 
 * to determine functionality. NOTE - NOT A CALIBRATION TEST.
 *
 * Author: Vishal Ravi, Gouri Bellad, Mark Hartigan
 * Date: February 25 2021
 */
 
define enter start;

Sensor ad15_vcc 1Hz { //vcc is placeholder adc. thermistors may be connected to different adc

 where A0 is thermistor1; //thermistors may be connected to different ports A0-A3 for this adc or they might be spread out among other adcs
 where A1 is thermistor2;
 where A2 is thermistor3;
 where A3 is thermistor4;
 
[calibrate   | thermistor1, poly, raw, V | 1 0 ];  //placeholder calibration curve, get from electronics, simbox electronics, or simbox software
[calibrate   | thermistor1, poly, V, C   | -0.0000003
                                            -0.000004
                                             0.0023
                                            -0.1077
                                            1.9093
];
[conversions | thermistor1, raw, V, C    |  ];
 
[calibrate   | thermistor2, poly, raw, V | 1 0 ];
[calibrate   | thermistor2, poly, V, C   | -0.00002
                                             0.0023
                                            -0.106
                                            2.003
];
[conversions | thermistor2, raw, V, C    |  ];

[calibrate   | thermistor3, poly, raw, V | 1 0 ]
[conversions | thermistor3, raw, V, C    | -.00002
                                             .0024
                                            -.1074
                                            1.9844
] 
[conversions | thermistor3, raw, V, C
 
[calibrate   | thermistor4, poly, raw, V | 1 0 ]
[conversions | thermistor4, raw, V, C    | -.00002
                                            .0025
                                           -.1009
                                           1.828
]


 // options are Flow, Temperature, or Bubble
   [print | white, thermistor1 | 1 ];
   [print | white, thermistor2 | 1 ];
   [print | white, thermistor3 | 1 ];
   [print | white, thermistor4 | 1 ];

}



Sensor ad15_sda 1 Hz {

 where A0 is thermistor5 //thermistors may be connected to different ports A0-A3 for this adc or they might be spread out among other adcs
 where A1 is thermistor6
 Where A2 is thermistor7
 where A3 is thermistor8
 
 [calibrate   | thermistor5, poly, raw, V | 1 0 ]
 [conversions | thermistor5, raw, V, C    |-.00002
                                            .0025
                                           -.1003
                                           1.8398
 ]  

 [calibrate   | thermistor6, poly, raw, V | 1 0 ]
 [conversions | thermistor6, raw, V, C    | -.0001
                                             .0022
                                            -.1076
                                            2.0743
 ]  
 
 [calibrate   | thermistor7, poly, raw, V | 1 0 ]
 [conversions | thermistor7, raw, V, C    | -.00002
                                             .0025
                                            -.1031
                                            1.8698
 ]  

 [calibrate   | thermistor8, poly, raw, V | 1 0 ]
 [conversions | thermistor8, raw, V, C    | -.00002
                                             .0024
                                            -.106
                                            1.9608
 ]  
 
  // options are Flow, Temperature, or Bubble
  [print | red, thermistor5 | 1 ];
  [print | red, thermistor6 | 1 ];
  [print | red, thermistor7 | 1 ];
  [print | red, thermistor8 | 1 ];
}

Sensor ad15_vdd 1H {
where A0 is thermistor9 //thermistors may be connected to different ports A0-A3 for this adc or they might be spread out among other adcs
where A1 is thermistor10
where A2 is thermistor11
where A3 is thermistor12

 [calibrate   | thermistor9, poly, raw, V | 1 0 ]
 [conversions | thermistor9, raw, V, C    | -.00002
                                             .0023
                                            -.1048
                                            1.9794
 ]
 
  [calibrate   | thermistor10, poly, raw, V | 1 0 ]
 [conversions | thermistor10, raw, V, C    | -.00002
                                              .0023
                                             -.1065
                                             2.0187
 ]
 
  [calibrate   | thermistor11, poly, raw, V | 1 0 ]
 [conversions | thermistor11, raw, V, C    | -.00002
                                              .0022
                                             -.1075
                                             2.0676
 ]
 
 
  [calibrate   | thermistor12, poly, raw, V | 1 0 ]
 [conversions | thermistor12, raw, V, C    | -.00002
                                              .0023
                                             -.1009
                                             1.9147
 ]


  // options are Flow, Temperature, or Bubble
  [print | red, thermistor9 | 1 ];
  [print | red, thermistor10 | 1 ];
  [print | red, thermistor11 | 1 ];
  [print | red, thermistor12 | 1 ];
}
 ds32 1Hz { 
  [calibrate   | Time, poly, raw, s | 0.0009765625, 0.0]
  [conversions | Time, raw, s       |                  ]
}

