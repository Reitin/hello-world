#include <stdio.h>
#include <stdlib.h>

int main(int argc, char *argv[])
{
    FILE *in;
    FILE *out;

    float *time;
    float *value;
    
    float t, v;
    
    int Nch;
    
    
//time [s], value
    printf("Enter number of spectral channels Nch=");
    scanf("%d", &Nch);
    
    time = new float[Nch];
    value = new float[Nch];
    printf ("Ok!\n");
    
    in = fopen(argv[1], "r+");
    if (in == NULL) {printf("Error opening file: %s\n", argv[1]); exit(1);}
    
    out = fopen(argv[2], "wb");
    if (out == NULL) {printf("Error creating file: %s\n", argv[2]); exit(1);}
    
    printf("Files opened!\n");
	while(!feof(in))
	{
        for (int i=0; i < Nch; i++)
        {
            fscanf(in, "%f %f", &t, &v);
            time[i] = t/86400.0;
            value[i] = v;
            printf("Time: %f, Value: %f\n", time[i], value[i]);
        }
		if (value[0]!=0)
		{
            float p = time[0];
			fwrite(&p, sizeof(float), 1, out);
			fwrite(value, sizeof(float), (Nch-1), out);
		}
	}

    fclose(in);
    fclose(out);

}
