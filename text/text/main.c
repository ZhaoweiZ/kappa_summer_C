//
//  main.c
//  text
//
//  Created by Zhaowei Zhang on 6/5/15.
//  Copyright (c) 2015 Zhaowei Zhang. All rights reserved.
//

#include <stdio.h>

int main(int argc, const char * argv[])
{
    
}

double check():
{
    double array_1[10];
    double array_2[10];
    for(int i = 0; i<10;++i)
    {
        array_1[i] = i;
        array_2[i] = 10*i;
    }

    for(int k = 0; k<10; ++k)
    {
        if(array_1[k] > 5)
        {
            return (array_2)[k];
        }
    }
}
