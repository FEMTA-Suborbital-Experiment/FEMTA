/*
** Collection Chamber Solenoid Test
**      author: Evan Rittner
**        date: June 23, 2021
** description:
**     Experiment to verify the collection chamber solenoid.
**     One pressure sensor records continuously; the
**     solenoid opens at a certain point.
*/


define enter test;
define leave complete;


Sensor ds32 1Hz {
  [calibrate   | Time, poly, raw, s | 0.0009765625, 0.0];
  [conversions | Time, raw, s       |                  ];
}


Sensor ad15_vdd 1Hz {
    
    where A0 is CC;
    
    [calibrate | CC, poly, V, kPa |
        -0.00000805518095144006,
         0.000179617025688186,
        -0.00149263466053365,
         0.0137083951784765,
         3.84928878749102,
        -0.389869854380195
    ];

    [calibrate   | CC, poly, raw, V | 0.00041234314, -0.00198277794];
    [conversions | CC, raw, V, kPa  |                              ];
    [print       | mint, CC         | 3                            ];
    
    
    if (State test | CC > 0kPa) {
        // Open solenoid after delay
        set pin 16 pos after 120s; // Change time until solenoid opens here and next line
        set pin 12 neg after 120s;
        set pin 16 neg after 350ms;
        
        leave test;
        enter complete after 1s;
    }
}
