define enter start;

Sensor test 1Hz {
// [calibrate   | random, poly, raw, C | 2.0, 2.0]
//  [conversions | random, raw, C       |       ]
  [print       | yellow, random      |       ]
}


//Sensor test 1Hz {
//  [calibrate   | random, poly, raw, C | 0.000000005, 1]
//  [conversions | random, raw, C       |       ]
//  [smooth      | random               | 0.01  ];
//  [print       | orange, random       |       ]
//}

//Sensor test 1Hz {
//  [calibrate   | saw_tooth, poly, raw, C | 1, 0]
//  [conversions | saw_tooth, raw, C       |       ]
//  [smooth      | saw_tooth               | 0.8  ]
//  [print       | orange, saw_tooth       |       ]
//}

Sensor test 1Hz {
//  [calibrate   | randsmooth, poly, raw | 3, 1  ]
//  [conversions | randsmooth, raw       |       ]
  [print       | orange, randsmooth    |       ]
}

