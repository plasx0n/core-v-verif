#include <stdint.h>
#include "perfs.h"
#include "trame.h"


//--------------
    // min arrays n-sys 
    int8_t    LLR[2*codw_N]  = {0}; 
    int8_t    PS[codw_N]     = {0};
    int8_t    decode[codw_N] = {0};

// low level polar funcs
int8_t func_f(int8_t la,int8_t lb){
    int8_t min1 , min2 ; 
    int8_t sign = 0 ; 

    min1 = abs(la) ; 
    min2 = abs(lb) ;

    if(min1>min2) min1=min2 ; 
    sign = (la < 0) ^ (lb < 0) ; 

    return (sign == 0 )? min1 : -min1 ; 
}

int8_t func_r(int8_t la,int8_t  froozen){
    if(froozen){
        return 0 ;
    } 
    else {
        return (la < 0) ; 
    }
}

int16_t func_g(int8_t sa,int16_t la,int16_t lb){
    if ( sa==0 ){
        return la+lb ; 
    }
    else{ 
        return lb-la ; 
    }
}

int8_t func_h(int8_t sa,int8_t sb){
    return sa^sb; 
}

int8_t sat( int16_t val){
    // evaluation uniquement sur G en 16 bits 
    if( val >= 127 )
        return 127 ; 
    else if( val <= -127)
        return -127 ; 
    else 
        return val ;
}

////////////////////////// 
// main func 

void node( int8_t* ptr_sum, int8_t *LLR , int N, int8_t *fz_bits,int8_t *decode)
{
    if( N == 1 ){    
		*ptr_sum = func_r(*LLR, *fz_bits ); 
        if ( *fz_bits == 0 )
			*decode = *ptr_sum ;
		else
			*decode = 0;
        return;
    }
        // ON CALCULE LES F
        for( int x = 0 ; x < N/2; x += 1 ){
            (LLR+N)[ x ]= func_f( LLR[ x ], (LLR+N/2)[ x ]);  
        }
 
        // ON CALCULE LA BRANCHE GAUCHE
        node( ptr_sum, (LLR+N), N/2, fz_bits,decode);

        // ON CALCULE LES G
        for( int x = 0;  x < N/2; x += 1 ){
            __int16_t temp = func_g( ptr_sum[x] , (__int16_t) LLR[ x ], (__int16_t) (LLR+N/2)[ x ]) ; 
            (LLR+N)[ x ] =sat( temp)  ; 
        }

        // ON CALCULE LA BRANCHE DROITE
        node( ptr_sum+ N/2, (LLR+N), N/2,  fz_bits+ N/2, decode+N/2);
    
        // ON FAIT LES CALCUL DES H (XOR DES SP)
        for(int x = 0 ; x < N/2 ; x += 1 ){          
            ptr_sum[x] = func_h( ptr_sum[x], ptr_sum[ x + (N/2) ]);     
        }
}


int main() {

    perfs_init();

    node(PS, LLR, codw_N, froozen_bits,decode  );

    perfs_end();
    display_perf();

	return 0 ; 
}

