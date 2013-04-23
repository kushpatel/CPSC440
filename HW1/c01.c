#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>


main()
{
  int ints[20];
  double a[2][2];

  double b[4];

  //
  // the output will go to the screen and c.13
  //
  cprini("stdout","c.13");
  fill_integer_array(&ints, 20);
  cprinf("ints = ", ints, 20);  

  //
  //          |  1.0     3.0  |
  //     a =  |               |
  //          |  2.0     4.0  |
  //
  a[0][0] = 1;
  a[1][0] = 2;
  a[0][1] = 3;
  a[1][1] = 4;

  //
  //     print the 2x2 matrix a as a chunk of memory
  //     note the order: first row, then second row
  //
  cprin2("a = ", &a[0][0], 4);

  //
  //     print the second row of a
  //
  cprin2("2nd row of a = ", &a[1][0], 2);
  
  cprina("string = ", "big and small", 3);

  treat_as_matrix(&b, 2, 2);
  cprin2("b = ", b, 4);

  int j = 1;
  cprinf("j = ", &j, 1);
  change_int_dumb(j);
  cprinf("j = ", &j, 1);
  change_int(&j);
  cprinf("j = ", &j, 1);

}

void change_int(int* a) {
  *a = 2;
}

void change_int_dumb(int a) {
  a = 2;
}


void fill_integer_array(int* ks, int n) {

  int i;
  int* p = ks;
 
  for (i = 0; i < n; i++) {
    *p = i + 11;
    p = p + 1;
  }

}


//
//  n = number of rows, m = number of columns
//
void treat_as_matrix(double *a, int n, int m) {
  int i, j;
  for (i = 0; i < n; i++) {
    for (j = 0; j < m; j++) {
      //
      // this means a[i][j] = i + 3*j;
      // however a[i][j] will not compile of course...
      //
      *(a+i*m+j) = i + 3*j;
    }
  }
}
