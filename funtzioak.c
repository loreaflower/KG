#include "cargar-triangulo.h"
#include <stdbool.h>
#include <math.h>
#include <stdio.h>

void print_matrizea(char *str, double *m)
{
    int i;

    printf("%s\n",str);
    for (i = 0;i<4;i++)
    printf("%lf, %lf, %lf, %lf\n",m[i*4],m[i*4+1],m[i*4+2], m[i*4+3]);
}

void mxp(punto *pptr, double m[16], punto p)
{
    pptr->x = p.x*m[0] + p.y*m[1] + p.z*m[2] + m[3];
    pptr->y = p.x*m[4] + p.y*m[5] + p.z*m[6] + m[7];
    pptr->z = p.x*m[8] + p.y*m[9] + p.z*m[10] + m[11];
 
    pptr->u = p.u;
    pptr->v = p.v;
}


void interpolatu_barizentroarekin(float t, punto *p1, punto *p2, punto *emaitzaptr)
{
    emaitzaptr->x= t*p1->x+(1-t)*p2->x;
    emaitzaptr->y= t*p1->y+(1-t)*p2->y;
    emaitzaptr->z= t*p1->z+(1-t)*p2->z;    
    emaitzaptr->u= t*p1->u+(1-t)*p2->u;  
    emaitzaptr->v= t*p1->v+(1-t)*p2->v; 
}



void matrize_biderketa(double *emptr, double *ezkmptr, double *eskmptr)
{
    
    int i,j,k;
    double emaitza;
    for(i=0; i<4; i++){
        for(j=0;j<4;j++){
            emaitza=0.0;
            for(k=0;k<4;k++){
                emaitza+=ezkmptr[i*4+k]*eskmptr[k*4+j];
            }
            emptr[i*4+j]=emaitza;
        }
    }
    /*
    print_matrizea("Eskuineko matrizea", eskmptr);
    print_matrizea("Ezkerreko matrizea", ezkmptr);
    print_matrizea("Emaitzaren matrizea", emptr);
    */
}


void biderketa_bektoriala(double *vem, double *v1, double *v2)
{
    vem[0] = v1[1]*v2[2]-v2[1]*v1[2];
    vem[1] = -(v1[0]*v2[2]-v2[0]*v1[2]);
    vem[2] = v1[0]*v2[1]-v2[0]*v1[1];
}

double biderketa_eskalarra(double *v1, double *v2, int tam)
{
    int i;
    double batura;

    batura=0;
    for(i=0; i<tam; i++)
    {
        batura+= v1[i]*v2[i];
    }
    return batura;
}




void kalkulatu_normala(hiruki *tri)
{
    double norma;
    double v1[3],v2[3],vbiderketa[3];


    v1[0]=tri->p2.x-tri->p1.x;
    v1[1]=tri->p2.y-tri->p1.y;
    v1[2]=tri->p2.z-tri->p1.z;

    v2[0]=tri->p3.x-tri->p1.x;
    v2[1]=tri->p3.y-tri->p1.y;
    v2[2]=tri->p3.z-tri->p1.z;

    biderketa_bektoriala(vbiderketa,v1,v2);

    norma =sqrt(pow(vbiderketa[0],2)+pow(vbiderketa[1],2)+pow(vbiderketa[2],2));

    if(norma>0)
    {
        tri->normala[0]=vbiderketa[0]/norma;
        tri->normala[1]=vbiderketa[1]/norma;
        tri->normala[2]=vbiderketa[2]/norma;
    }
    else
    {
        tri->normala[0]=-1;
        tri->normala[1]=-1;
        tri->normala[2]=-1;
    }
  

}

void mesa_kalkulatu(double *emptr, double *mobj)
{
    double xyz[3];
    double e[3];

    e[0]=mobj[3];
    e[1]=mobj[7];
    e[2]=mobj[11];
    
    //x bektorea
    xyz[0]=mobj[0];
    xyz[1]=mobj[4];
    xyz[2]= mobj[8];

    emptr[0]= xyz[0];
    emptr[1]= xyz[1];
    emptr[2]= xyz[2];
    emptr[3]= -biderketa_eskalarra(xyz, e,3);

    //y bektorea
    xyz[0]=mobj[1];
    xyz[1]=mobj[5];
    xyz[2]= mobj[9];

    emptr[4]= xyz[0];
    emptr[5]= xyz[1];
    emptr[6]= xyz[2];
    emptr[7]= -biderketa_eskalarra(xyz, e,3);

    //z bektorea
    xyz[0]=mobj[2];
    xyz[1]=mobj[6];
    xyz[2]=mobj[10];

    emptr[8]= xyz[0];
    emptr[9]= xyz[1];
    emptr[10]= xyz[2];
    emptr[11]= -biderketa_eskalarra(xyz, e,3);

    emptr[12]= 0;
    emptr[13]= 0;
    emptr[14]= 0;
    emptr[15]= 1;
}
