//
//  Integration.h
//  UniDecLibrary
//
//  Created by Austin  Carr on 3/8/22.
//

#ifndef Integration_h
#define Integration_h

#include <stdio.h>
int integrate_dd(double* kernel_x, double* kernel_y, int kernellen, double* data_x,
                 double* data_y, int datalen, double** kernel_x_new, double** kernel_y_new);



#endif /* Integration_h */
