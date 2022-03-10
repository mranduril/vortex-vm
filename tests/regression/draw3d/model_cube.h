#pragma once

#include "model.h"

const model_t model_cube = {
    true,  // depth
    false, // color
    true,  // tex
    { // vertices
        {-6.337301, 0.000000, 24.177938, 24.949747, 0xffffffff, 0.000000, 0.000000},
        {6.337301, 0.000000, 24.177938, 24.949747, 0xffffffff, 1.000000, 0.000000},
        {6.337301, 11.949731, 18.974358, 20.000000, 0xffffffff, 1.000000, 1.000000},
        {-6.337301, 11.949731, 18.974358, 20.000000, 0xffffffff, 0.000000, 1.000000},
        {-6.337301, -11.949731, 18.974358, 20.000000, 0xffffffff, 1.000000, 1.000000},
        {6.337301, -11.949731, 18.974358, 20.000000, 0xffffffff, 0.000000, 1.000000},
        {6.337301, 0.000000, 13.770777, 15.050253, 0xffffffff, 0.000000, 0.000000},
        {-6.337301, 0.000000, 13.770777, 15.050253, 0xffffffff, 1.000000, 0.000000},
        {6.337301, 11.949731, 18.974358, 20.000000, 0xffffffff, 1.000000, 1.000000},
        {-6.337301, 11.949731, 18.974358, 20.000000, 0xffffffff, 0.000000, 1.000000},
        {-6.337301, -11.949731, 18.974358, 20.000000, 0xffffffff, 1.000000, 1.000000},
        {6.337301, -11.949731, 18.974358, 20.000000, 0xffffffff, 0.000000, 1.000000},
        {6.337301, 0.000000, 13.770777, 15.050253, 0xffffffff, 0.000000, 0.000000},
        {-6.337301, 0.000000, 13.770777, 15.050253, 0xffffffff, 1.000000, 0.000000},
        {-6.337301, 0.000000, 24.177938, 24.949747, 0xffffffff, 0.000000, 0.000000},
        {6.337301, 0.000000, 24.177938, 24.949747, 0xffffffff, 1.000000, 0.000000},
        {6.337301, 11.949731, 18.974358, 20.000000, 0xffffffff, 1.000000, 1.000000},
        {-6.337301, 11.949731, 18.974358, 20.000000, 0xffffffff, 0.000000, 1.000000},
        {-6.337301, -11.949731, 18.974358, 20.000000, 0xffffffff, 1.000000, 1.000000},
        {6.337301, -11.949731, 18.974358, 20.000000, 0xffffffff, 0.000000, 1.000000},
        {-6.337301, 0.000000, 24.177938, 24.949747, 0xffffffff, 0.000000, 0.000000},
        {6.337301, 0.000000, 24.177938, 24.949747, 0xffffffff, 1.000000, 0.000000},
        {6.337301, 11.949731, 18.974358, 20.000000, 0xffffffff, 1.000000, 1.000000},
        {-6.337301, 11.949731, 18.974358, 20.000000, 0xffffffff, 0.000000, 1.000000},
        {-6.337301, -11.949731, 18.974358, 20.000000, 0xffffffff, 1.000000, 1.000000},
        {6.337301, -11.949731, 18.974358, 20.000000, 0xffffffff, 0.000000, 1.000000},
        {6.337301, 0.000000, 13.770777, 15.050253, 0xffffffff, 0.000000, 0.000000},
        {-6.337301, 0.000000, 13.770777, 15.050253, 0xffffffff, 1.000000, 0.000000},
        {6.337301, 0.000000, 24.177938, 24.949747, 0xffffffff, 1.000000, 0.000000},
        {6.337301, 11.949731, 18.974358, 20.000000, 0xffffffff, 1.000000, 1.000000},
        {-6.337301, 11.949731, 18.974358, 20.000000, 0xffffffff, 0.000000, 1.000000},
        {-6.337301, -11.949731, 18.974358, 20.000000, 0xffffffff, 1.000000, 1.000000},
        {6.337301, -11.949731, 18.974358, 20.000000, 0xffffffff, 0.000000, 1.000000},
        {6.337301, 0.000000, 13.770777, 15.050253, 0xffffffff, 0.000000, 0.000000}
    }, 
    { // primitives
        {2, 1, 3},
        {3, 1, 0},
        {1, 2, 0},    
        {0, 2, 3},
        {4, 0, 5},
        {5, 0, 1},    
        {1, 5, 0},
        {0, 5, 4},
        {3, 0, 7},    
        {7, 0, 4},
        {5, 4, 1},
        {1, 4, 0}
    },
    "fire.png"
};