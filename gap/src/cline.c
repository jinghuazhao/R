/* Last modified Dec 13, 2001 */

#define MFLEN 80
#define MNLEN 16
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

int get_flag(int argc, char **argv, char *flag, int vtype, void *val) {
   int found, i, j, len, iv, negate=0;
   long lv;
   float fv;
   double dv;
   int *ival;
   long *lval;
   float *fval;
   double *dval;
   char fs[MFLEN], fss[MNLEN];
   char *argi;
   
   ival = (int *) val;
   lval = (long *) val;
   fval = (float *) val;
   dval = (double *) val;
   len = strlen(flag);
   if (!len) return 0;
   for (i=1; i<argc; i++) {
      argi = argv[i];
      if (argi && argi[0] == '-') {
         if (argi[1]=='n' && argi[2]=='o') {
	   argi += 3;
	   negate = 1;
	 }
	 else {
	   argi++;
	   negate = 0;
	 }
         found = 1;
         for (j=0; j<len; j++) found = found && (argi[j] == flag[j]);
      } 
      else {
         found = 0;
      }
      if (found) {
	argv[i] = (char *) 0; /* Delete flag */
         if (negate && vtype) return -1; 
         strcpy(fs, argi+len);
	 /* If no value in flag, look at next argument */
	 if (vtype>0 && strlen(fs)==0 && (i+1)<argc) {
	   argi = argv[i+1];
	   if (argi[0]!='-') {
	     strcpy(fs, argi);
	     argv[i+1] = (char *) 0;
	   }
	 } 
	 if (vtype==0) {
	   if (strlen(fs)>1) return -1;
	   if (strlen(fs)==0) {
	     *(int*)val = negate? 0 : 1;
	     return i;
	   }
	   if (fs[0] ==  '+') {
	     *(int*)val = negate? 0 : 1;
	     return i;
	   }
	   if (fs[0] ==  '-') {
	     *(int*)val = negate? 1 : 0;
	     return i;
	   }
	   return -1;
	 }
	 else if (vtype==4) {
	   if (!val) val = (char *) malloc(1+strlen(fs));
	   strcpy((char*)val, fs);
	   return i;
	 }
	 else {
	   j = 0;
	   do {
	     for (i=0; fs[j] && fs[j]!=':'; i++, j++) fss[i] = fs[j];
	     fss[i] = 0;
	     switch (vtype){
	     case 1: 
               if (sscanf(fss, "%d", &iv)) {
		 *(ival++) = iv;
                  return i;
	       } else {
                  return -1;
	       }
	     case 5:
	       if (sscanf(fss, "%ld", &lv)) {
		 *(lval++) = lv;
		 return i;
	       }
	       else {
		 return -1;
	       }
	     case 2:
               if (sscanf(fss, "%f", &fv)) {
                  *(fval++) = fv;
                  return i;
	       } else {
                  return -1;
	       }
            case 3:
               if (sscanf(fss, "%lf", &dv)) {
                  *(dval++) = dv;
                  return i;
	       } else {
                  return -1;
	       }
	     }
	   } while (fs[j++]);
	 }
      }
   }
   return 0; /* Flag not found */
}

int get_arg(int argc, char** argv, char* val) {
   int i;
   char* argi;
   for (i=1; i<argc; i++) {
      argi = argv[i];
      if (argi && argi[0] != '-') {
	argv[i] = (char *) 0;
	if (argi[0]=='.' && !argi[1]) {
	  *val = (char) 0;
	  return 0; /* Argument missing */
	}
	else {
	  strcpy(val, argi);
	  return i;
	}
      }
   }
   *val = (char) 0;
   return 0; /* Argument not found */
}

char* unrec(int argc, char** argv) {
  int i;
  for (i=1; i<argc; i++) {
    if (argv[i]) return argv[i];
  }
  return (char *) 0;
}
