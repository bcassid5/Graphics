/*
 Ben Cassidy - BCASSID5 - 250808910
 CS 3388 Graphics
 Assignment 2
 */

#include <X11/Xlib.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "matrix.h"

#define Ex 15.0
#define Ey 15.0
#define Ez 15.0

#define Gx 0.0
#define Gy 0.0
#define Gz 0.0

#define UPx 0.0
#define UPy 0.0
#define UPz 1.0

#define NP 5.0
#define FP 50.0

#define THETA 90.0
#define ASPECT 1.0

#define W  512
#define H  512

#define FALSE 0
#define TRUE  1
#define COLS  600
#define ROWS  600

#define POSX  0
#define POSY  0

#define TRIANGLE 3


Display *InitX(Display *d, Window *w, int *s) {
    
    d = XOpenDisplay(NULL) ;
    if(d == NULL) {
        printf("Cannot open display\n") ;
        exit(1) ;
    }
    *s = DefaultScreen(d) ;
    *w = XCreateSimpleWindow(d, RootWindow(d, *s), POSX, POSY, COLS, ROWS, 1, BlackPixel(d, *s), WhitePixel(d, *s)) ;
    Atom delWindow = XInternAtom(d, "WM_DELETE_WINDOW", 0) ;
    XSetWMProtocols(d, *w, &delWindow, 1) ;
    XSelectInput(d, *w, ExposureMask | KeyPressMask) ;
    XMapWindow(d, *w) ;
    return(d) ;
}

void SetCurrentColorX(Display *d, GC *gc, unsigned int r, unsigned int g, unsigned int b) {
    
    XSetForeground(d, *gc, r << 16 | g << 8 | b) ;
}

void SetPixelX(Display *d, Window w, int s, int i, int j) {
    
    XDrawPoint(d, w, DefaultGC(d, s), i, j) ;
}

void QuitX(Display *d, Window w) {
    
    XDestroyWindow(d,w) ;
    XCloseDisplay(d) ;
}


dmatrix_t *build_camera_matrix(dmatrix_t *E, dmatrix_t *G) {
    
    dmatrix_t N ; /* Viewing axis */
    
    N = *dmat_normalize(dmat_sub(E,G)) ;
    N.l = 3 ;

    dmatrix_t UP ;
    dmat_alloc(&UP,4,1) ;
    UP.l = 3 ;
    
    UP.m[1][1] = UPx ;
    UP.m[2][1] = UPy ;
    UP.m[3][1] = UPz ;
    UP.m[4][1] = 1.0 ;
    
    dmatrix_t U ;
    
    U = *dmat_normalize(dcross_product(&UP,&N)) ;
    
    dmatrix_t V ;
    V = *dcross_product(&N,&U) ;
    
    dmatrix_t Mv ; /* Build matrix M_v */
    dmat_alloc(&Mv,4,4) ;
    
    Mv.m[1][1] = U.m[1][1] ; 
    Mv.m[1][2] = U.m[2][1] ; 
    Mv.m[1][3] = U.m[3][1] ; 
    Mv.m[1][4] = -1.0*((*E).m[1][1]*U.m[1][1] + (*E).m[2][1]*U.m[2][1] + (*E).m[3][1]*U.m[3][1]) ;
    
    Mv.m[2][1] = V.m[1][1] ; 
    Mv.m[2][2] = V.m[2][1] ; 
    Mv.m[2][3] = V.m[3][1] ; 
    Mv.m[2][4] = -1.0*((*E).m[1][1]*V.m[1][1] + (*E).m[2][1]*V.m[2][1] + (*E).m[3][1]*V.m[3][1]) ;
    
    Mv.m[3][1] = N.m[1][1] ; 
    Mv.m[3][2] = N.m[2][1] ; 
    Mv.m[3][3] = N.m[3][1] ; 
    Mv.m[3][4] = -1.0*((*E).m[1][1]*N.m[1][1] + (*E).m[2][1]*N.m[2][1] + (*E).m[3][1]*N.m[3][1]) ;
    
    Mv.m[4][1] = 0.0 ; 
    Mv.m[4][2] = 0.0 ; 
    Mv.m[4][3] = 0.0 ; 
    Mv.m[4][4] = 1.0 ;
    
    dmatrix_t Mp ; /* Build matrix Mp */
    dmat_alloc(&Mp,4,4) ;
    Mp = *dmat_identity(&Mp) ;
    
    float a = -1.0*(FP + NP)/(FP - NP) ;
    float b = -2.0*(FP*NP)/(FP - NP) ;
    
    Mp.m[1][1] = NP ;
    Mp.m[2][2] = NP ;
    Mp.m[3][3] = a ;
    Mp.m[3][4] = b ;
    Mp.m[4][3] = -1.0 ;
    Mp.m[4][4] = 0.0 ;
    
    /* Build matrices T_1 and S_1 */
    
    /* Work out coordinates of near plane corners */
    
    float top = NP*tan(M_PI/180.0*THETA/2.0) ;
    float right = ASPECT*top ;
    float bottom = -top ;
    float left = -right ;
   
    dmatrix_t T1 ;
    dmat_alloc(&T1,4,4) ;
    
    T1 = *dmat_identity(&T1) ;
    T1.m[1][4] = -(right + left)/2.0 ;
    T1.m[2][4] = -(top + bottom)/2.0 ;

    dmatrix_t S1 ;
    dmat_alloc(&S1,4,4) ;
    
    S1 = *dmat_identity(&S1) ;
    S1.m[1][1] = 2.0/(right - left) ;
    S1.m[2][2] = 2.0/(top - bottom) ;

    /* Build matrices T2, S2, and W2 */
    
    dmatrix_t T2 ;
    dmatrix_t S2 ;
    dmatrix_t W2 ;
    
    dmat_alloc(&T2,4,4) ;
    dmat_alloc(&S2,4,4) ;
    dmat_alloc(&W2,4,4) ;
    
    T2 = *dmat_identity(&T2) ;
    S2 = *dmat_identity(&S2) ;
    W2 = *dmat_identity(&W2) ;
    
    T2.m[1][4] = 1.0 ;
    T2.m[2][4] = 1.0 ;

    S2.m[1][1] = W/2.0 ;
    S2.m[2][2] = H/2.0 ;
    
    W2.m[2][2] = -1.0 ;
    W2.m[2][4] = (double)H ;
    
    return dmat_mult(&W2,dmat_mult(&S2,dmat_mult(&T2,dmat_mult(&S1,dmat_mult(&T1,dmat_mult(&Mp,&Mv)))))) ;
}

void exchangeInt(int *a, int *b){ 
    int t ;
    t = *a ;
    *a = *b ;
    *b = t ;
}

/*
 Bresenham Line Drawing Algorithm
 - takes two points and draws a line between those points (in any direction - symmetric cases solved)
 - from assignment 1
 */
void Bresenham(Display *d, Window w, int s, int x1, int y1, int x2, int y2){

    int Transx, Transy ;
    int Pi, Dx, Dy, Two_Dx, Two_Dy, i, Inc1stcoord, Inc2ndcoord, Exchange ;
    
    Exchange = FALSE ;
    Inc1stcoord = 1 ;
    Inc2ndcoord = 1 ;
    
    Transx = -x1 ;
    Transy = -y1 ;
    
    x1 = 0 ;
    y1 = 0 ;
    
    x2 += Transx ;
    y2 += Transy ;
    
    Dx = x2 ;
    Dy = y2 ;
    
    Two_Dx = 2*x2 ;
    Two_Dy = 2*y2 ;
    
    if (Dy < 0) {
        Inc2ndcoord = -1 ;
        Dy *= -1 ;
        Two_Dy *= -1 ;
    }
    
    if (Dx < 0) {
        Inc1stcoord = -1 ;
        Dx *= -1 ;
        Two_Dx *= -1 ;
    }
    
    
    if (Dy > Dx) {
        Exchange = TRUE ;
        exchangeInt(&Two_Dx,&Two_Dy) ;
        exchangeInt(&Dx,&Dy) ;
        exchangeInt(&Inc1stcoord,&Inc2ndcoord) ;
    }
    
    Pi = Two_Dy - Dx ;
    if (Exchange) {
        SetPixelX(d, w, s, y1 - Transx, x1 - Transy) ;
    }
    else {
        SetPixelX(d, w, s, x1 - Transx, y1 - Transy) ;
    }
    for (i = 0 ; i < Dx ; i++) {
        if (Pi < 0) {
            Pi += Two_Dy ;
        }
        else {
            Pi += Two_Dy - Two_Dx ;
            y1 += Inc2ndcoord ;
        }
        x1 += Inc1stcoord ;
        if (Exchange) {
            SetPixelX(d, w, s, y1 - Transx, x1 - Transy) ;
        }
        else {
            SetPixelX(d, w, s, x1 - Transx, y1 - Transy) ;
        }
    }
}

dmatrix_t *perspective_projection(dmatrix_t *P) {

    (*P).m[1][1] /= (*P).m[4][1] ;
    (*P).m[2][1] /= (*P).m[4][1] ;
    (*P).m[3][1] /= (*P).m[4][1] ;
    (*P).m[4][1] /= (*P).m[4][1] ;

    return P ;
}

/*
 sphere Matrix Points Calculator
 - creates the matrix points for the sphere
 */
dmatrix_t *sphere_matrix(dmatrix_t *P, double r, double theta, double rho) {
    //parametric equations for a sphere to create the P matrix
    (*P).m[1][1] = r*cos(theta)*sin(rho);
    (*P).m[2][1] = r*sin(theta)*sin(rho);
    (*P).m[3][1] = r*cos(rho);
    (*P).m[4][1] = 1.0;
    
    return P;
}

/*
 Torus Matrix Points Calculator
 - creates the matrix points for the torus
 */
dmatrix_t *torus_matrix(dmatrix_t *P, double a, double c, double theta, double rho) {
    //parametric equations for a sphere to create the P matrix
    (*P).m[1][1] = (c + a*cos(rho))*cos(theta);
    (*P).m[2][1] = (c + a*cos(rho))*sin(theta);
    (*P).m[3][1] = a*sin(rho);
    (*P).m[4][1] = 1.0;
    
    return P;
}

/*
 Sphere Wire Mesh Renderer
 - This function takes all necessary paramters and calculates the points needed to draw the wire mesh frame of a sphere. The lines are drawn using bresenham's line drawing algorithm implemented in assignment 1.
 */

void sphere_wiremesh(Display *d, Window w, int s, dmatrix_t *C, double rad, double theta_low, double theta_up, double rho_low, double rho_high, double dTheta, double dRho){
    
    //initialize variables
    int x;
    double theta, rho;
    dmatrix_t P[TRIANGLE];
    
    //dmat_alloc call loop for new matrix P
    for (x = 0; x < TRIANGLE; x++){
        dmat_alloc(&P[x],4,1);
    }
    //main loop - nested for loop for each plane to rotate about different axes
    for (theta = theta_low; theta < theta_up; theta += dTheta){
        for (rho = rho_low; rho < rho_high; rho += dRho) {
            //call projection 3 times to create the points from 2D to 3D
            P[0] = *perspective_projection(dmat_mult(C, sphere_matrix(&P[0], rad, theta, rho)));
            P[1] = *perspective_projection(dmat_mult(C, sphere_matrix(&P[1], rad, theta+dTheta, rho)));
            P[2] = *perspective_projection(dmat_mult(C, sphere_matrix(&P[2], rad, theta, rho+dRho)));
            //draw the lines
            Bresenham(d,w,s,(int)P[0].m[1][1],(int)P[0].m[2][1],(int)P[1].m[1][1],(int)P[1].m[2][1]);
            Bresenham(d,w,s,(int)P[0].m[1][1],(int)P[0].m[2][1],(int)P[2].m[1][1],(int)P[2].m[2][1]);
        }
    }
    //clear the matrices
    for (x=0;x<TRIANGLE;x++){
        free_dmatrix(P[x].m, 1, P[x].l, 1, P[x].c);
    }
    
}

/*
 Torus Wire Mesh Renderer
 - This function takes all necessary paramters and calculates the points needed to draw the wire mesh frame of a torus. The lines are drawn using bresenham's line drawing algorithm implemented in assignment 1.
 */

void torus_wiremesh(Display *d, Window w, int s, dmatrix_t *C, double a, double c, double theta_low, double theta_up, double rho_low, double rho_high, double dTheta, double dRho){
    
    //initialize variables
    int x;
    double theta, rho;
    dmatrix_t P[TRIANGLE];
    
    //dmat_alloc call loop for new matrix P
    for (x = 0; x < TRIANGLE; x++){
        dmat_alloc(&P[x],4,1);
    }
    
    //main loop - nested for loop for each plane to rotate about different axes
    for (theta = theta_low; theta < theta_up; theta += dTheta){
        for (rho = rho_low; rho < rho_high; rho += dRho) {
            //call projection 3 times to create the points from 2D to 3D
            P[0] = *perspective_projection(dmat_mult(C, torus_matrix(&P[0], a, c, theta, rho)));
            P[1] = *perspective_projection(dmat_mult(C, torus_matrix(&P[1], a, c, theta+dTheta, rho)));
            P[2] = *perspective_projection(dmat_mult(C, torus_matrix(&P[2], a, c, theta, rho+dRho)));
            //draw the lines
            Bresenham(d,w,s,(int)P[0].m[1][1],(int)P[0].m[2][1],(int)P[1].m[1][1],(int)P[1].m[2][1]);
            Bresenham(d,w,s,(int)P[0].m[1][1],(int)P[0].m[2][1],(int)P[2].m[1][1],(int)P[2].m[2][1]);
        }
    }
    //clear the matrices
    for (x=0;x<TRIANGLE;x++){
        free_dmatrix(P[x].m, 1, P[x].l, 1, P[x].c);
    }
    
}

int main() {
    
    //create the window variables
    Display *d ;
    Window w ;
    XEvent e ;
    int s ;
    
    unsigned int r, g, b ;
    double rad, u, v, theta, dTheta, rho, dRho, dM_PI;
    
    r = g = b = 0 ;
    
    //open the window
    d = InitX(d, &w, &s) ;
    SetCurrentColorX(d, &(DefaultGC(d, s)), r, g, b) ;
    

    dmatrix_t E ; /* The centre of projection for the camera */
    
    dmat_alloc(&E,4,1) ;
    
    E.m[1][1] = Ex ;
    E.m[2][1] = Ey ;
    E.m[3][1] = Ez ;
    E.m[4][1] = 1.0 ;
    
    dmatrix_t G ; /* Point gazed at by camera */
    
    dmat_alloc(&G,4,1) ;
    
    G.m[1][1] = Gx ;
    G.m[2][1] = Gy ;
    G.m[3][1] = Gz ;
    G.m[4][1] = 1.0 ;

    dmatrix_t C ; /* The camera matrix */

    dmat_alloc(&C,4,4) ;
    C = *build_camera_matrix(&E,&G) ;

    //declare some variables for use in the renderer functions
    dM_PI = 2.0*M_PI;
    dTheta = 2.0*M_PI/36.0;
    dRho = M_PI/18.0;
    rad = 10.0;

    while (1) {
        XNextEvent(d, &e) ; // fixes double tracing issue
        
        //call the wire mesh renderer functions here to draw on the run of the program
        if(e.type == Expose){
            sphere_wiremesh(d, w, s, &C, 10.0, 0.0, dM_PI, 0.0, M_PI, dTheta, dRho);
            torus_wiremesh(d, w, s, &C, 2.0, 15.0, 0.0, dM_PI, 0.0, dM_PI, dTheta, dTheta);
        }
        if(e.type == KeyPress)
            break ;
        if(e.type == ClientMessage)
            break ;
    }
    QuitX(d,w) ;
    
    //clear the matrices after ending the program
    free_dmatrix(E.m, 1, E.l, 1, E.c);
    free_dmatrix(G.m, 1, G.l, 1, G.c);
    free_dmatrix(C.m, 1, C.l, 1, C.c);
}
