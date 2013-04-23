#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>


/*
*
***************************************************************************
*
*        This is the end of the debugging code, and the beginning
*        of the printing functions proper. The code below contains
*        five user-callable functions: cprini, cprin, cprin2, cprinf,
*        cprin_start_stop. Their functions are as follows.
*
*  cprini - communicates to the function cprin8 (NOT a user-callable
*        function) the names of the two files on which the output 
*        will be written. Please note that the name "stdout" is 
*        treated differently from all other names: specifying it 
*        will cause the output to be written on the standard output
*        - screen, most probably.
*  cprin - prints a real *4 (of the type "float") array
*  cprin2 - prints a double precision (of the type "double") array
*  cprinf - prints an integer *4 array
*
***************************************************************************
*
*/


void cprini(str17,str27)
char *str17,*str27;
{
  int *adp,*ap,itype,i1,i2,mes,afp,n,acp;
/*
*        . . . Initialize
*/
        itype=11;
        cprin8(mes,ap,afp,adp,acp,n,itype,str17,str27,i1,i2);
}



void cprina(mes,acp,n)
char *mes,*acp;
int n;
{
    int *adp,*ap,itype,str17,str27,i1,i2,afp;
/*
*        Print character data
*/
        itype=4;
        cprin8(mes,ap,afp,adp,acp,n,itype,str17,str27,i1,i2);

}


void cprinf(mes,afp,n)
char *mes;
int n;
double *afp;
{
    int *adp,*ap,itype,str17,str27,i1,i2,acp;
/*
*        Print integer data
*/
        itype=2;
        cprin8(mes,ap,afp,adp,acp,n,itype,str17,str27,i1,i2);

}



void cprin(mes,ap,n)
char *mes;
int n;
double *ap;
{
    int *adp,*afp,itype,str17,str27,i1,i2,acp;
/*
*        Print single precision data
*/
        itype=1;
        cprin8(mes,ap,afp,adp,acp,n,itype,str17,str27,i1,i2);
}



void cprin2(mes,adp,n)
char *mes;
int n;
double *adp;
{
    int *ap,*afp,itype,str17,str27,i1,i2,acp;
/*
*        Print double precision data
*/
        itype=3;
        cprin8(mes,ap,afp,adp,acp,n,itype,str17,str27,i1,i2);
}



void cprin_start_stop(i1,i2)
int i1,i2;
{
  int *adp,*ap,itype,mes,afp,n,acp;
  char *str17,*str27;
/*
*        . . . stop/resume
*/
        itype=21;
        cprin8(mes,ap,afp,adp,acp,n,itype,str17,str27,i1,i2);

}


void cprin8(mes,ap,afp,adp,acp,n,itype,str17,str27,i1,i2)
     char *mes,*str17,*str27,*acp;
int n,*afp,itype,i1,i2;
float *ap;
double *adp;
{
    static int ifprint1,ifprint2;
    static FILE *st1,*st2,*dumfile;
/*
*        If this is the initialization call - open the 
*        files str17, str27
*/

//    if (strcmp(str17, "stdout") == 0) printf("A");

    if(itype==11) {
        ifprint1=0;
        ifprint2=0;


        if(!str17) goto jump1;
        if ( strcmp(str17, "stdout") != 0 ) st1=fopen(str17,"w");
        if ( strcmp(str17, "stdout") == 0 ) st1=stdout;

    jump1:

        if(!str27) goto jump2;
        if ( strcmp(str27, "stdout") != 0 ) st2=fopen(str27,"w");
        if ( strcmp(str27, "stdout") == 0 ) st2=stdout;
    jump2:

        if(str17) ifprint1=1;
        if(str27) ifprint2=1;
        return;
    }
/*
*        If this is the "stop/resume" call - stop/resume printing 
*/
    if(itype==21) {

        

        if(i1==0) ifprint1=0;
        if(i1!=0) ifprint1=1;

        if(i2==0) ifprint2=0;
        if(i2!=0) ifprint2=1;

        return;
    }


        if(ifprint1 !=0) cprin7(mes,ap,afp,adp,acp,n,itype,st1);
        if(ifprint2 !=0) cprin7(mes,ap,afp,adp,acp,n,itype,st2);
        return;
}



void cprin7(mes,ap,afp,adp,acp,n,itype,str)
char *mes,*acp;
FILE *str;
int n,*afp,itype;
float *ap;
double *adp;
{
    int i,c;
/*
*            Process the message
*/
      if(mes) fprintf(str,"%s\n",mes);

/*
*          Process the double precision data to be printed
*/
      if(itype ==3) {
          for(i=0; i<n; i=i+1) {
              fprintf(str,"  %11.4le", adp[i]);
              if(i%6==5 || i==n-1) fprintf(str,"\n");
          }
          return;
      }

/*
*          Process the integer data to be printed
*/
    if(itype ==2) {
        for(i=0; i<n; i=i+1) {
            fprintf(str," %7d", afp[i]);
            if(i%10==9 || i==n-1) fprintf(str,"\n");
        }
        return;
    }

/*
*          Process the single precision data to be printed
*/
    if(itype ==1) {
        for(i=0; i<n; i=i+1) {
            fprintf(str,"  %11.4le", *ap++);
            if(i%6==5 || i==n-1) fprintf(str,"\n");
        }
        return;
    }

/*
*          Process the character data to be printed
*/

    if(itype==4) {
        for(i=0; i<n; i++) {
            fprintf(str,"%c", acp[i]);
	    if(i%60==59 || i==n-1) fprintf(str,"\n");
        }
        return;
    }

}



















