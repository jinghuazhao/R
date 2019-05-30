#include <stdio.h>
//  https://www.geeksforgeeks.org/snprintf-c-library/
int main()
{
    char buffer[50];
    char* s = "geeksforgeeks";
//  Counting and buffering the character using snprintf
    int j = snprintf(buffer, 6, "%s.%s\n", s);
//  Print the buffered string and character count
    printf("%s, character count = %04d, sizes = %03d/%02d\n", buffer, j, sizeof(buffer), sizeof(j));
    return 0;
}
