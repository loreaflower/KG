//	Program developed by
//	
//	Informatika Fakultatea
//	Euskal Herriko Unibertsitatea
//	http://www.ehu.eus/if
//
// to compile it: gcc dibujar-triangulos-y-objetos.c -lGL -lGLU -lglut -lm
//                gcc *.c -lGL -lGLU -lglut -lm guztiak batera egiteko


//TODO NORMALAK ONDO FUNTZIONATZEN DUELA AZTERTU


#include <GL/glut.h>
#include <stdio.h>
#include <string.h>
#include "cargar-triangulo.h"
#include <math.h>
#include <stdbool.h>
#include "funtzioak.h"

typedef struct mlist
{
    double m[16];
    double mesa[16];
    struct mlist *hptr;
} mlist;
    
typedef struct triobj
{
    hiruki *triptr;
    int num_triangles;
    mlist *mptr;
    struct triobj *hptr;
    unsigned char *kolorea;
} triobj;

// testuraren informazioa
// información de textura

extern int load_ppm(char *file, unsigned char **bufferptr, int *dimxptr, int * dimyptr);
unsigned char *bufferra;
int dimx,dimy;

int indexx;
hiruki *triangulosptr;
triobj *foptr; //first object pointer
triobj *sel_ptr; //selected object
triobj *sel_cam_ptr; //selected camera
triobj *fcptr; //first camera
int denak;
int lineak;
int objektuak;
char aldaketa;
int ald_lokala;
int ald_kamera;
int paralelo;
int back_culling;
int normala;
int dimentsioak;

double mlag[16];
double mjatorrira[16];
double mesamobj[16];


char fitxiz[100];


unsigned char * color_textura(float u, float v)
{
    int desplazamendua;
    int xind, yind;
    char * lag;

    //u eta v [0,1] artean egoteko
    if(u<0) u=0;
    if(u>1) u=1;
    if(v<0) v=0;
    if(v>1) v=1;
    
    //desplazamendua lortu
    xind = u*(dimx-1); 
    xind=xind % dimx;
    yind = v*(dimy-1); 
    yind=yind % dimy;
    yind = dimy-yind-1;
    desplazamendua = yind*dimx + xind;



    lag = (unsigned char *)bufferra;  // pixel on the left and top
    return(lag+(3*desplazamendua));
}



void  dibujar_linea_z(int linea,float c1x, float c1z, float c1u,float c1v,float c2x,float c2z,float c2u,float c2v)
{
    float xkoord,zkoord;
    float u,v;
    float t,deltat;
    unsigned char r,g,b;
    unsigned char *colorv;

    glBegin( GL_POINTS );
    if((c2x -c1x)!=0)
        deltat = 1.0/(c2x-c1x);

    for (t=1;t>=0;t-=deltat)
    {
        //xkoord, zkoor, u eta v interpolatu
        xkoord = t*c1x+(1-t)*c2x;
        zkoord = t*c1z+(1-t)*c2z;
        u =  t*c1u+(1-t)*c2u;
        v=  t*c1v+(1-t)*c2v;
        colorv=  color_textura(u, v); 
        r= colorv[0];
        g=colorv[1];
        b=colorv[2];    
        glColor3ub(r,g,b);
        glVertex3f(xkoord, linea, zkoord );
        
    }
    glEnd();
}


void marraztu_normala(hiruki *tptr, punto p1)
{
    punto p4, p5;
    p4.x=tptr->p1.x+tptr->normala[0]*40;
    p4.y=tptr->p1.y+tptr->normala[1]*40;
    p4.z=tptr->p1.z+tptr->normala[2]*40;

    mxp(&p5, mesamobj, p4);

    glBegin(GL_LINES);
    glVertex3d(p1.x, p1.y, p1.z);
    glVertex3d(p5.x, p5.y, p5.z);
    glEnd();
}

int triangelua_marraztu(triobj *optr,hiruki *tptr)
{
    punto e, esae, zkam, esazkam;
    double akenb[3], vxn;
    if(paralelo==0)
    {
        e.x= sel_cam_ptr->mptr->m[3];
        e.y= sel_cam_ptr->mptr->m[7];
        e.z= sel_cam_ptr->mptr->m[11];
        e.u=0;
        e.v=0;

        mxp(&esae, optr->mptr->mesa, e);

        akenb[0] = esae.x - tptr->p1.x;
        akenb[1] = esae.y - tptr->p1.y;
        akenb[2] = esae.z - tptr->p1.z;

    }
    else{
        zkam.x= sel_cam_ptr->mptr->m[2];
        zkam.y= sel_cam_ptr->mptr->m[6];
        zkam.z= sel_cam_ptr->mptr->m[10];
        zkam.u=0;
        zkam.v=0;

        mxp(&esazkam, optr->mptr->mesa, zkam);
        akenb[0] = esazkam.x;
        akenb[1] = esazkam.y;
        akenb[2] = esazkam.z;
        
    }
        vxn =biderketa_eskalarra(akenb, tptr->normala,3);

    if(vxn<0) return 0;
    else return 1;


}

void dibujar_triangulo(triobj *optr, int i, double *m)
{
    hiruki *tptr;

    punto *pgoiptr, *pbeheptr, *perdiptr;
    float x1,h1,z1,u1,v1,x2,h2,z2,u2,v2,x3,h3,z3,u3,v3;
    float c1x,c1z,c1u,c1v,c2x,c2z,c2u,c2v;
    int linea;
    float cambio1,cambio1z,cambio1u,cambio1v,cambio2,cambio2z,cambio2u,cambio2v;
    punto p1,p2,p3,eb1,eb2; 
    float t, delta12; //pgoiptr eta perdptr artean dagoen lerro batetik besterako diferentzia kalkulatzeko
    float s, delta13; //pgoiptr eta pbeheptr artean dagoen lerro batetik besterako diferentzia kalkulatzeko
    float q, delta23; //perdptr eta pbeheptr artean dagoen lerro batetik besterako diferentzia kalkulatzeko
    int aurrea;

    if (i >= optr->num_triangles) return;
    //taula horretako i triangelua hartzen du eta triangeluak p1,p2 eta p3 ditu
    tptr = optr->triptr+i;
    mxp(&p1,m,tptr->p1);
    mxp(&p2,m,tptr->p2);
    mxp(&p3,m,tptr->p3);

    
    aurrea = triangelua_marraztu(optr, tptr);
    if(aurrea==1 || back_culling==0)
    {
        if (lineak == 1) //triangeluaren lineak marrazteko
        {
            if(aurrea==1) 
            {   
                if(optr->kolorea!=0)
                    glColor3ub(optr->kolorea[0],optr->kolorea[1], optr->kolorea[2]);
                else
                    glColor3d(255.0,255.0,255.0);
            }
            else glColor3d(255.0,0.0,0.0);

            if(normala==1)
                marraztu_normala(tptr, p1);

            glBegin(GL_POLYGON);
            glVertex3d(p1.x, p1.y, p1.z);
            glVertex3d(p2.x, p2.y, p2.z);
            glVertex3d(p3.x, p3.y, p3.z);
            glEnd();
            return;
        }
        else 
        {
            //Triangeluko puntuak ordenatu
            if(p1.y>p2.y)
            {
                pgoiptr = &p1;
                pbeheptr = &p2;
            }
            else
            {
                pgoiptr = &p2;
                pbeheptr = &p1;
            }

            if(p3.y > pgoiptr->y)
            {
                perdiptr = pgoiptr;
                pgoiptr = &p3;
            }
            else 
            {
                if(p3.y<pbeheptr->y)
                {
                    perdiptr = pbeheptr;
                    pbeheptr = &p3;
                }
                else
                {
                    perdiptr = &p3;
                }
            }



            //interpolaziorako puntuen arteko diferentzien alderantzizkoa lortu
            if(pgoiptr->y!=perdiptr->y) delta12 = 1/(pgoiptr->y - perdiptr->y); 
            else delta12 = 2; //puntuak altura berdinean badaude =2 esleitu begizta behin bakarrik egieko (lerroa marraztu)


            if(perdiptr->y!=pbeheptr->y) delta23 = 1/(perdiptr->y - pbeheptr->y);
            else delta23 = 2;


            if(pgoiptr->y!=pbeheptr->y) delta13 = 1/(pgoiptr->y - pbeheptr->y);
            else delta13 = 2;

        

            for(t=1,s=1; t>=0; t-=delta12, s-=delta13)
            {
                interpolatu_barizentroarekin(t, pgoiptr, perdiptr, &eb1);
                interpolatu_barizentroarekin(s, pgoiptr, pbeheptr, &eb2);

                if(eb1.x<eb2.x) //ezkerrerago zein dagoen erabaki
                    dibujar_linea_z(eb1.y,eb1.x,eb1.z,eb1.u,eb1.v,eb2.x,eb2.z,eb2.u,eb2.v);
                else
                    dibujar_linea_z(eb2.y,eb2.x,eb2.z,eb2.u,eb2.v,eb1.x,eb1.z,eb1.u,eb1.v);
            }
            

            for(q=1; q>=0; q-=delta23, s-=delta13)
            {
                interpolatu_barizentroarekin(q, perdiptr, pbeheptr, &eb1);
                interpolatu_barizentroarekin(s, pgoiptr, pbeheptr, &eb2);   

                if(eb1.x<eb2.x)
                    dibujar_linea_z(eb1.y,eb1.x,eb1.z,eb1.u,eb1.v,eb2.x,eb2.z,eb2.u,eb2.v);
                else
                    dibujar_linea_z(eb2.y,eb2.x,eb2.z,eb2.u,eb2.v,eb1.x,eb1.z,eb1.u,eb1.v);    
            }
            
        } 
    }  
}





static void marraztu(void)
{
    float u,v;
    int i,j;
    triobj *auxptr;
  
    /*
    unsigned char* colorv;
    unsigned char r,g,b;
    */

    // marrazteko objektuak behar dira
    // no se puede dibujar sin objetos
    if (foptr ==0) return;

    // clear viewport...
    if (objektuak == 1) glClear( GL_DEPTH_BUFFER_BIT|GL_COLOR_BUFFER_BIT );
        else 
        {
        if (denak == 0) glClear( GL_DEPTH_BUFFER_BIT|GL_COLOR_BUFFER_BIT );
        }

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(-500.0, 500.0, -500.0, 500.0,0.,1000.0);
 
    
    if (objektuak == 1)
        {
        if (denak == 1)
            {
            for (auxptr =foptr; auxptr != 0; auxptr = auxptr->hptr)
                {
                    matrize_biderketa(mesamobj, sel_cam_ptr->mptr->mesa,auxptr->mptr->m);
                for (i =0; i < auxptr->num_triangles; i++)
                    {
                    dibujar_triangulo(auxptr,i, mesamobj);
                    }
                }
            }
        else
            {
                matrize_biderketa(mesamobj, sel_cam_ptr->mptr->mesa,sel_ptr->mptr->m);
            for (i =0; i < sel_ptr->num_triangles; i++)
                {
                dibujar_triangulo(sel_ptr,i,mesamobj);
                }
            }
        }
    else
        {
        matrize_biderketa(mesamobj, sel_cam_ptr->mptr->mesa,sel_ptr->mptr->m);
        dibujar_triangulo(sel_ptr,indexx, mesamobj);
        }
    glFlush();
}



void analisi_aldaketa(double x, double y, double z)
{
    double alpha = 0.05;
    int i;
    mlist *em1ptr= (mlist *)malloc(sizeof(mlist));
    mlist *em2ptr= (mlist *)malloc(sizeof(mlist));


    //1. Mugitu kamera jatorrira
    for (i=0; i<16; i++) mjatorrira[i] =0;
         mjatorrira[0] = 1.0;
         mjatorrira[5] = 1.0;
         mjatorrira[10] = 1.0;
         mjatorrira[15] = 1.0;

    if(ald_kamera==1){
        mjatorrira[3]=-sel_ptr->mptr->m[3];
        mjatorrira[7]=-sel_ptr->mptr->m[7];
        mjatorrira[11]=-sel_ptr->mptr->m[11];
    }else{
        mjatorrira[3]=-sel_cam_ptr->mptr->m[3];
        mjatorrira[7]=-sel_cam_ptr->mptr->m[7];
        mjatorrira[11]=-sel_cam_ptr->mptr->m[11];
    }

    //2. biratu adierazitako ardatz baten inguruan
    mlag[0]=cos(alpha)+(1-cos(alpha))*x*x;        mlag[1]=(1-cos(alpha))*x*y-z*sin(alpha);      mlag[2]=(1-cos(alpha))*x*z+y*sin(alpha);      mlag[3]=0;
    mlag[4]=(1-cos(alpha))*x*y+z*sin(alpha);      mlag[5]=cos(alpha)+(1-cos(alpha))*y*y;        mlag[6]=(1-cos(alpha))*y*z-x*sin(alpha);      mlag[7]=0;
    mlag[8]=(1-cos(alpha))*x*z-y*sin(alpha);      mlag[9]=(1-cos(alpha))*y*z+x*sin(alpha);      mlag[10]=cos(alpha)+(1-cos(alpha))*z*z;       mlag[11]=0;
    mlag[12]=0;                                   mlag[13]=0;                                   mlag[14]=0;                                   mlag[15]=1;

    matrize_biderketa(em1ptr->m, mjatorrira, mlag);

    //3. Eraman kamera hasierako lekura
    mjatorrira[3]=-mjatorrira[3];
    mjatorrira[7]=-mjatorrira[7];
    mjatorrira[11]=-mjatorrira[11];
    
    matrize_biderketa(em2ptr->m,em1ptr->m,mjatorrira);

    if(ald_kamera==1){
        em2ptr->hptr = sel_cam_ptr->mptr;
        sel_cam_ptr->mptr= em2ptr;
    }else{
        em2ptr->hptr = sel_ptr->mptr;
        sel_ptr->mptr= em2ptr;
    }
    
}

void lortu_cam_e(double *e)
{
    e[0]= sel_cam_ptr->mptr->m[3];
    e[1]= sel_cam_ptr->mptr->m[7];
    e[2]= sel_cam_ptr->mptr->m[11];
    e[3]=1;
}

void lortu_obj_at(double *at)
{
    at[0]= sel_ptr->mptr->m[3];
    at[1]= sel_ptr->mptr->m[7];
    at[2]= sel_ptr->mptr->m[11];
    at[3]=1;
}

void lortu_vup(double *vup)
{
    vup[0]=0;
    vup[1]=1;
    vup[2]=0;
    vup[3]=0;
}

void birkolokatu(double *posberria, double *e, double *vup, double *at)
{ 
    double xc[3]; 
    double yc[3]; 
    double zc[3]; 
    double mods; 
    double lag[3];

    lag[0] =e[0]-at[0]; 
    lag[1] =e[1]-at[1]; 
    lag[2] =e[2]-at[2]; 

    mods= sqrt(pow(lag[0], 2)+pow(lag[1], 2)+pow(lag[2], 2)); 

    if (mods==0) return; 

    zc[0]= (lag[0])/mods; 
    zc[1]= (lag[1])/mods; 
    zc[2]= (lag[2])/mods; 

    biderketa_bektoriala(lag, vup, zc); 

    mods= sqrt(pow(lag[0], 2)+pow(lag[1], 2)+pow(lag[2], 2)); 

    if (mods==0) return; 

    xc[0]= (lag[0])/mods; 
    xc[1]= (lag[1])/mods; 
    xc[2]= (lag[2])/mods; 

    biderketa_bektoriala(yc, zc, xc); 

    posberria[0] = xc[0];   posberria[1] = yc[0];   posberria[2] = zc[0];   posberria[3] = e[0]; 
    posberria[4] = xc[1];   posberria[5] = yc[1];   posberria[6] = zc[1];   posberria[7] = e[1];  
    posberria[8] = xc[2];   posberria[9] = yc[2];   posberria[10] = zc[2];  posberria[11] = e[2]; 
    posberria[12] = 0;      posberria[13] = 0;      posberria[14] = 0;      posberria[15] = 1; 
}

void analisi_kokaketa()
{
    double e[4], at[4],vup[4];
    mlist *emptr;
    lortu_cam_e(e);
    lortu_obj_at(at);
    lortu_vup(vup);

    birkolokatu(emptr->m,e,at,vup);
    mesa_kalkulatu(emptr->mesa,emptr->m);

    emptr->hptr=0;
    sel_cam_ptr->mptr->hptr=emptr;
}

void read_from_file(char *fitx, int camera) //camera==0 bada objektu gisa gordeko dugu, camera==1 bada kamera gisa
{
int i,retval;
triobj *optr;
//kamera posizio egokian kokatzeko beharrezko aldagaiak
double e[4];
double vup[4];
double at[4];

    //printf("%s fitxategitik datuak hartzera\n",fitx);
    optr = (triobj *)malloc(sizeof(triobj));
    optr->kolorea=0;
    //retval = cargar_triangulos(fitx, &(optr->num_triangles), &(optr->triptr));
    retval = cargar_triangulos_color(fitx,&(optr->num_triangles), &(optr->triptr), &(optr->kolorea));
    if (retval == -1) 
         {
         printf("%s fitxategitik datuak hartzerakoan arazoak izan ditut\n    Problemas al leer\n",fitxiz);
         free(optr);
         }
    else
    {
        
         triangulosptr = optr->triptr;
         //printf("objektuaren matrizea...\n");
         optr->mptr = (mlist *)malloc(sizeof(mlist));

         for (i=0; i<16; i++) optr->mptr->m[i] =0;
         optr->mptr->m[0] = 1.0;
         optr->mptr->m[5] = 1.0;
         optr->mptr->m[10] = 1.0;
         optr->mptr->m[15] = 1.0;


         for (i=0; i<16; i++) optr->mptr->mesa[i] =0;
         optr->mptr->mesa[0] = 1.0;
         optr->mptr->mesa[5] = 1.0;
         optr->mptr->mesa[10] = 1.0;
         optr->mptr->mesa[15] = 1.0;

        optr->mptr->hptr = 0;


         //printf("objektu zerrendara doa informazioa...\n");
        if(camera==0)
        {
            optr->hptr = foptr;
            foptr = optr;
            sel_ptr = optr;
         }else{
            e[0]=0; e[1]=0; e[2]=-700; e[3]=1;

            lortu_vup(vup);


//TODO VUP[3]=0 EDO VUP[3]=0?



            //vup[0]=0; vup[1]=1; vup[2]=0; vup[3]=0;

            at[0]=0; at[1]=0; at[2]=0; at[3]=0;

            birkolokatu(optr->mptr->m, e, vup, at);
            mesa_kalkulatu(optr->mptr->mesa, optr->mptr->m);

            optr -> hptr = fcptr;
            fcptr = optr;
            sel_cam_ptr = optr;
         }
    }
     printf("datuak irakurrita\nLecura finalizada\n");
}

void x_aldaketa(int dir)
{
    int a;
    double alpha = 0.05;
    mlist *emptr= (mlist *)malloc(sizeof(mlist));
    if(dir) a=1.5;
    else a=-1.5;
    if(aldaketa == 'r') //rotazioa
    {
        mlag[0]=1;      mlag[1]=0;             mlag[2]=0;                 mlag[3]=0;
        mlag[4]=0;      mlag[5]=cos(alpha);    mlag[6]=-sin(alpha)*a;     mlag[7]=0;
        mlag[8]=0;      mlag[9]=sin(alpha)*a;  mlag[10]=cos(alpha);       mlag[11]=0;
        mlag[12]=0;     mlag[13]=0;            mlag[14]=0;                mlag[15]=1;
    }
    else //traslazioa
    {
        mlag[0]=1;      mlag[1]=0;             mlag[2]=0;                 mlag[3]=4*a;
        mlag[4]=0;      mlag[5]=1;             mlag[6]=0;                 mlag[7]=0;
        mlag[8]=0;      mlag[9]=0;             mlag[10]=1;                mlag[11]=0;
        mlag[12]=0;     mlag[13]=0;            mlag[14]=0;                mlag[15]=1;
    }

    if(ald_kamera)
    {
        if(aldaketa=='r'){
            if(ald_lokala)
            {
                matrize_biderketa(emptr->m, sel_cam_ptr->mptr->m, mlag);
                emptr->hptr = sel_cam_ptr->mptr;
                sel_cam_ptr->mptr= emptr;

                mesa_kalkulatu(sel_cam_ptr->mptr->mesa, sel_cam_ptr->mptr->m);
            }
            //else{analisi_aldaketa(sel_ptr->mptr->m[1],sel_cam_ptr->mptr->m[5],sel_ptr->mptr->m[9]);}

            
        }
    }
    else
    {
        if(ald_lokala){ matrize_biderketa(emptr->m, sel_ptr->mptr->m, mlag);}
        else{ matrize_biderketa(emptr->m, mlag, sel_ptr->mptr->m);}

        emptr->hptr = sel_ptr->mptr;
        sel_ptr->mptr= emptr;
        mesa_kalkulatu(sel_ptr->mptr->mesa, sel_ptr->mptr->m);
    }
   
}


void y_aldaketa(int dir)
{
    mlist *emptr= (mlist *)malloc(sizeof(mlist));
    int a;
    double alpha = 0.05;


    if(dir) a=1.5;
    else a=-1.5;
    if(aldaketa == 'r') //rotazioa
    {

        mlag[0]=cos(alpha);       mlag[1]=0;        mlag[2]=sin(alpha)*a;      mlag[3]=0;
        mlag[4]=0;                mlag[5]=1;        mlag[6]=0;                 mlag[7]=0;
        mlag[8]=-sin(alpha)*a;    mlag[9]=0;        mlag[10]=cos(alpha);       mlag[11]=0;
        mlag[12]=0;               mlag[13]=0;       mlag[14]=0;                mlag[15]=1;
    }
    else //traslazioa
    {
        mlag[0]=1;      mlag[1]=0;             mlag[2]=0;                 mlag[3]=0;
        mlag[4]=0;      mlag[5]=1;             mlag[6]=0;                 mlag[7]=4*a;
        mlag[8]=0;      mlag[9]=0;             mlag[10]=1;                mlag[11]=0;
        mlag[12]=0;     mlag[13]=0;            mlag[14]=0;                mlag[15]=1;
    }



    if(ald_kamera==1)
    {
        if(aldaketa=='r'){
            if(ald_lokala==1)
            {
                matrize_biderketa(emptr->m, sel_cam_ptr->mptr->m, mlag);
                emptr->hptr = sel_cam_ptr->mptr;
                sel_cam_ptr->mptr= emptr;

                mesa_kalkulatu(sel_cam_ptr->mptr->mesa, sel_cam_ptr->mptr->m);
            }
            //else{analisi_aldaketa(sel_ptr->mptr->m[1],sel_cam_ptr->mptr->m[5],sel_ptr->mptr->m[9]);}

            
        }
    }
    else
    {
        if(ald_lokala){ matrize_biderketa(emptr->m, sel_ptr->mptr->m, mlag);}
        else{ matrize_biderketa(emptr->m, mlag, sel_ptr->mptr->m);}

        emptr->hptr = sel_ptr->mptr;
        sel_ptr->mptr= emptr;
        mesa_kalkulatu(sel_ptr->mptr->mesa, sel_ptr->mptr->m);
    }
    
}



void z_aldaketa(int dir)
{
    mlist *emptr= (mlist *)malloc(sizeof(mlist));
    int a;
    double alpha = 0.05;
    if(dir) a=1.5;
    else a=-1.5;
    if(aldaketa == 'r') //rotazioa
    {
        mlag[0]=cos(alpha);       mlag[1]=-sin(alpha)*a;   mlag[2]=0;     mlag[3]=0;
        mlag[4]=sin(alpha)*a;     mlag[5]=cos(alpha);      mlag[6]=0;     mlag[7]=0;
        mlag[8]=0;                mlag[9]=0;               mlag[10]=1;    mlag[11]=0;
        mlag[12]=0;               mlag[13]=0;              mlag[14]=0;    mlag[15]=1;
    }
    else //traslazioa
    {
        mlag[0]=1;      mlag[1]=0;             mlag[2]=0;                 mlag[3]=0;
        mlag[4]=0;      mlag[5]=1;             mlag[6]=0;                 mlag[7]=0;
        mlag[8]=0;      mlag[9]=0;             mlag[10]=1;                mlag[11]=4*a;
        mlag[12]=0;     mlag[13]=0;            mlag[14]=0;                mlag[15]=1;
    }

   if(ald_kamera)
    {
  
        if(ald_lokala)
        {
            
            matrize_biderketa(emptr->m, sel_cam_ptr->mptr->m, mlag);
            emptr->hptr = sel_cam_ptr->mptr;
            sel_cam_ptr->mptr= emptr;
            mesa_kalkulatu(sel_cam_ptr->mptr->mesa, sel_cam_ptr->mptr->m);
        }
        //else{analisi_aldaketa(sel_ptr->mptr->m[1],sel_cam_ptr->mptr->m[5],sel_ptr->mptr->m[9]);}

            

    }
    else
    {
        if(ald_lokala){ matrize_biderketa(emptr->m, sel_ptr->mptr->m, mlag);}
        else{ matrize_biderketa(emptr->m, mlag, sel_ptr->mptr->m);}

        emptr->hptr = sel_ptr->mptr;
        sel_ptr->mptr= emptr;
        mesa_kalkulatu(sel_ptr->mptr->mesa, sel_ptr->mptr->m);
    }
}

void kamera_obj_ald()
{
    triobj *lag;
    triobj *folag;


    lag = sel_ptr;
    folag = foptr;

    sel_ptr = sel_cam_ptr;
    foptr = fcptr;

    sel_cam_ptr = lag;
    fcptr = folag;
}

void undo()
{
    if(ald_kamera)
    {
        if(sel_cam_ptr->mptr->hptr!=0)
                sel_cam_ptr->mptr = sel_cam_ptr->mptr->hptr;
    }
    else
        {
        if(sel_ptr->mptr->hptr!=0)
                sel_ptr->mptr = sel_ptr->mptr->hptr;
        }
    
}

void eskalatu(int size)
{
    mlist *emptr= (mlist *)malloc(sizeof(mlist));
    double a;
    if(size==0) a=1/1.05;
    else a=1.05;
    mlag[0]=a;      mlag[1]=0;             mlag[2]=0;                 mlag[3]=0;
    mlag[4]=0;      mlag[5]=a;             mlag[6]=0;                 mlag[7]=0;
    mlag[8]=0;      mlag[9]=0;             mlag[10]=a;                mlag[11]=0;
    mlag[12]=0;     mlag[13]=0;            mlag[14]=0;                mlag[15]=1;


    matrize_biderketa(emptr->m, mlag, sel_ptr->mptr->m);

    emptr->hptr = sel_ptr->mptr;
    sel_ptr->mptr= emptr;
}



// This function will be called whenever the user pushes one key
static void teklatua (unsigned char key, int x, int y)
{
int retval;
int i;
FILE *obj_file;

switch(key)
	{
	case 13: 
	        if (foptr != 0)  // objekturik ez badago ezer ez du egin behar
	                         // si no hay objeto que no haga nada
	            {
	            indexx ++;  // azkena bada lehenengoa bihurtu
		                // pero si es el último? hay que controlarlo!
		    if (indexx == sel_ptr->num_triangles) 
		        {
		        indexx = 0;
		        if ((denak == 1) && (objektuak == 0))
		            {
		            glClear( GL_DEPTH_BUFFER_BIT|GL_COLOR_BUFFER_BIT );
		            glFlush();
		            }
		        }
		    }
		break;
	case 'd':
		if (denak == 1) denak = 0;
		    else denak = 1;
		break;
	case 'o':
		if (objektuak == 1) objektuak = 0;
		    else objektuak = 1;
		break;
	case 'l':
		if (lineak == 1) lineak = 0;
		    else lineak = 1;
		break;
	case 't':
	        aldaketa = 't';
            printf("Orain aldaketa TRASLAZIOA da\n");
		break;
	case 'r':
		aldaketa = 'r';
        printf("Orain aldaketa ERROTAZIOA da\n");
		break;
	case 'g':
    case 'G':
    if(ald_kamera==1)
    {
        if (ald_lokala==1)
        {
            ald_lokala==0;
            analisi_kokaketa();
            printf("Orain aldaketa ANALISI MODUAN egingo da\n");
        }
        else
        {
            ald_lokala==1;
            printf("Orain aldaketa HEGALDI MODUAN egingo da\n");
        }
        
    }
    else
    {
		if (ald_lokala == 1)
        {
            ald_lokala = 0;
            printf("Orain aldaketa GLOBALA da\n");
        }
		else
        {
            ald_lokala = 1;
            printf("Orain aldaketa LOKALA da\n");
        }
    }
		break;
    case 'x':
        x_aldaketa(1);
        break;
    case 'y':
        y_aldaketa(1);
        break;
    case 'z':
        z_aldaketa(1);
        break;
    case 'X':
        x_aldaketa(0);
        break;
    case 'Y':
        y_aldaketa(0);
        break;
    case 'Z':
        z_aldaketa(0);
        break;
    case 'u':
        undo();
        break;
    case 's':
        eskalatu(0);
        break;
    case 'S':
        eskalatu(1);
        break;
    case 'c':
        if (ald_kamera == 1)
        {
            ald_kamera = 0;
            printf("Orain aldaketa OBJEKTUARENA da\n");
        } 
	    else
        {
            ald_kamera = 1;
            printf("Orain aldaketa KAMERARENA da\n");
            if(ald_lokala==0){
                analisi_kokaketa();
                printf("Aldaketa ANALISI MODUAN egingo da\n");
            }
            else printf("Aldaketa HEGALDI MODUAN egingo da\n");
        }
        break;
    case 'C':
        kamera_obj_ald();
        break;
    case 'b':
        if(back_culling==1) back_culling=0;
        else back_culling=1;
        break;
    case 'n':
        if(normala==1) normala=0;
        else normala=1;
        break;
    case 'p':
        if(paralelo==1){
            paralelo=0;
            printf("Orain PERSPEKTIBAN dago\n");
        }
        else{
            paralelo=1;
            printf("Orain PARALELOAN dago\n");
        } 
        break;
	case 'f':
	        /*Ask for file*/
	    printf("idatzi fitxategi izena\n"); 
	    scanf("%s", &(fitxiz[0]));
	    read_from_file(fitxiz,0);
	    indexx = 0;
            break;
    /* case 'S':  // save to file
	    printf("idatzi fitxategi izena\n"); 
	    scanf("%s", &(fitxiz[0]));
            if ((obj_file = fopen(fitxiz, "w")) == NULL)
                {
                printf("ezin fitxategia ireki\n");
                }
            else
                {
                for (i =0; i < sel_ptr->num_triangles; i++)
                    {
                     fprintf(obj_file,"t %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n",
                     sel_ptr->triptr[i].p1.x-250, sel_ptr->triptr[i].p1.y-250, sel_ptr->triptr[i].p1.z, 
                     sel_ptr->triptr[i].p1.u, sel_ptr->triptr[i].p1.v,
                     sel_ptr->triptr[i].p2.x-250, sel_ptr->triptr[i].p2.y-250, sel_ptr->triptr[i].p2.z, 
                     sel_ptr->triptr[i].p2.u, sel_ptr->triptr[i].p2.v,
                     sel_ptr->triptr[i].p3.x-250, sel_ptr->triptr[i].p3.y-250, sel_ptr->triptr[i].p3.z, 
                     sel_ptr->triptr[i].p3.u, sel_ptr->triptr[i].p3.v );
                    }
                fclose(obj_file);
                }
        break; */
    case 9: /* <TAB> */
        if (foptr != 0) // objekturik gabe ez du ezer egin behar
                            // si no hay objeto no hace nada
            {
            sel_ptr = sel_ptr->hptr;
            /*The selection is circular, thus if we move out of the list we go back to the first element*/
            if (sel_ptr == 0) sel_ptr = foptr;
            indexx =0; // the selected polygon is the first one
            }
        break;
	case 27:  // <ESC>
		exit( 0 );
		break;
	default:
		printf("%d %c\n", key, key );
	}

// The screen must be drawn to show the new triangle
glutPostRedisplay();
}

void hasieraketa_globalak()
{
    denak = 0;
        lineak =1;
        objektuak = 1;
        foptr = 0;
        sel_ptr = 0;
        aldaketa = 'r';
        printf("Orain aldaketa ERROTAZIOA da\n");
        back_culling = 0;
        normala = 0;
        paralelo = 1;
        printf("Orain PARALELOAN dago\n");
        ald_lokala = 1;
        printf("Orain aldaketa LOKALA da\n");
        ald_kamera=0;
        printf("Orain aldaketa OBJEKTUARENA da\n");
}
int main(int argc, char** argv)
{
int retval;

   
	printf(" Triangeluak: barneko puntuak eta testura\n Triángulos con puntos internos y textura \n");
	printf("Press <ESC> to finish\n");
	glutInit(&argc,argv);
	glutInitDisplayMode ( GLUT_RGB|GLUT_DEPTH );
	glutInitWindowSize ( 500, 500 );
	glutInitWindowPosition ( 100, 100 );
	glutCreateWindow( "KBG/GO praktika" );
	
	glutDisplayFunc( marraztu );
	glutKeyboardFunc( teklatua );
	/* we put the information of the texture in the buffer pointed by bufferra. The dimensions of the texture are loaded into dimx and dimy */ 
        retval = load_ppm("testura.ppm", &bufferra, &dimx, &dimy);
        if (retval !=1) 
            {
            printf("Ez dago texturaren fitxategia (testura.ppm)\n");
            exit(-1);
            }
        
	glClearColor( 0.0f, 0.0f, 0.7f, 1.0f );
	glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
        glEnable(GL_DEPTH_TEST); // activar el test de profundidad (Z-buffer)
        read_from_file("cam.txt",1);
        hasieraketa_globalak();
        if (argc>1) read_from_file(argv[1], 0);
        else 
        {
            read_from_file("k.txt",0);
            read_from_file("k.txt",0);
        }    
            //read_from_file("adibideak.txt",0);
        //print_matrizea("Kameraren mobj matrizea", sel_cam_ptr->mptr->mesa);
	glutMainLoop();

	return 0;   
}
