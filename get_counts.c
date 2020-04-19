#include <stdio.h>
#include <stdlib.h>

int main()
{
    int c;
    while ((c = getchar()) != EOF)
    {
        if (c != 'a')
        {
            putchar(c);
        }
    }

    char *genome = argv[1];

    printf("X\n");

    return 0;
}