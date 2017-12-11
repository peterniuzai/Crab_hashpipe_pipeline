//gcc fwrite_test.c -o fwrute_test
#include <stdio.h>

void main(){
FILE * fp;
int a=1234567;
fp = fopen("fwrite_test.txt","w");
fwrite("123456789a",sizeof(char),10,fp);
fwrite(&a, sizeof(a),1,fp);
fclose(fp);
}
