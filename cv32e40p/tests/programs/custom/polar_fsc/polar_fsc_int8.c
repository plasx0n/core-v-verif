

#include "perfs.h"
#include <stdint.h>
// which trame 
#include "n3_repspc_k512.h"
#define FSC

int8_t sat( int16_t val){
    // evaluation uniquement sur G en 16 bits 

    if( val >= 127 )
        return 127 ; 
    else if( val <= -127)
        return -127 ; 
    else 
        return val ;
    }   

int8_t func_f(int8_t la,int8_t lb){
    int8_t min1 , min2 ; 
    int8_t sign = 0 ; 

    min1 = abs(la) ; 
    min2 = abs(lb) ;

    if(min1>min2) min1=min2 ; 
    sign = (la < 0) ^ (lb < 0) ; 

    return (sign == 0 )? min1 : -min1 ; 
}

int16_t func_g(int8_t sa,int16_t la,int16_t lb){
    if ( sa==0 )
    {
        return la+lb ; 
    }
    else
    { 
        return lb-la ; 
    }
}

int8_t func_r(int8_t la,int8_t  froozen)
{
    if(froozen)
    {
        return 0 ;  
    }
    else 
    {
        return (la < 0) ; 
    }
}

void node_8(int8_t* ptr_sum, int8_t *LLR , int N, char *fz_bits)
{ 
    if( N == 1 )
    {     
	    *ptr_sum = func_r(*LLR, *fz_bits ); 
        return;
    }

    for( int x = 0 ; x < N/2; x += 1 )
    {
        (LLR+N)[ x ] = func_f( LLR[ x ], (LLR+N/2)[ x ]) ;
    }

         
    // ON CALCULE LA BRANCHE GAUCHE
    node_8(ptr_sum, (LLR+N), N/2, fz_bits);

    // ON CALCULE LES G
    for( int x = 0;  x < N/2; x += 1 )
    {
        int16_t temp = func_g( ptr_sum[x] , (int16_t) LLR[ x ], (int16_t) (LLR+N/2)[ x ]) ; 
        (LLR+N)[ x ] =sat( temp)  ;  
    }

    // ON CALCULE LA BRANCHE DROITE
    node_8(ptr_sum+ N/2, (LLR+N), N/2,  fz_bits+ N/2 );
    
    // ON FAIT LES CALCUL DES H (XOR DES SP)
    for(int x = 0 ; x < N/2 ; x += 1 )
    {          
        ptr_sum[x] ^=ptr_sum[ x + (N/2) ]; 
    }
}

void node(int8_t* ptr_sum, int8_t *LLR , int N, char *fz_bits)
{
    // ON CALCULE LA BRANCHE GAUCHE
    // get the node status 
    char enab1 = *fb_table_tileN++  ; 

    if(enab1==0) // r0 
    {
        // update res with G 
        for( int x = 0;  x < N/2; x += 1 )
        {
            int16_t temp = func_g( ptr_sum[x] , (int16_t) LLR[ x ], (int16_t) (LLR+N/2)[ x ]) ; 
            (LLR+N)[ x ] =sat( temp)  ;  
        }

    }else 
    if (enab1==1) // r1 
    {
        // ON CALCULE LES F
        for( int x = 0 ; x < N/2; x += 1 )
        {
            (LLR+N)[ x ] = func_f( LLR[ x ], (LLR+N/2)[ x ]) ; 
            ptr_sum[x] = func_r((LLR+N)[x], 0 );      
        }

        // update Res with G 
        for( int x = 0;  x < N/2; x += 1 )
        {
            int16_t temp = func_g( ptr_sum[x] , (int16_t) LLR[ x ], (int16_t) (LLR+N/2)[ x ]) ; 
            (LLR+N)[ x ] =sat( temp)  ;  
        }

    }else 
    if (enab1==2) // REP
    {
        // Somme des LLR 
        // x > 0 ? PS => 0 
        // x < 0 ? PS => 1 

        // ON CALCULE LES F
        int tot = 0 ; 
        for( int x = 0 ; x < N/2; x += 1 )
        {
            (LLR+N)[ x ] = func_f( LLR[ x ], (LLR+N/2)[ x ]) ; 
            tot += (LLR+N)[x] ;
        }

        if(tot < 0 )
        {
            for( int x = 0 ; x < N/2; x += 1 )
                ptr_sum[x] = 1 ; 
        } 

        // update Res with G for next node 
        for( int x = 0;  x < N/2; x += 1 )
        {
            int16_t temp = func_g( ptr_sum[x] , (int16_t) LLR[ x ], (int16_t) (LLR+N/2)[ x ]) ; 
            (LLR+N)[ x ] =sat( temp)  ;  
        }


    } else 
    if (enab1==3 ) // SPC 
    {

        int8_t sign=0;
        int mina = 100  ;
        int8_t idx_min= 0 ; 
        for( int x = 0 ; x < N/2; x += 1 )
        {
            (LLR+N)[ x ] = func_f( LLR[ x ], (LLR+N/2)[ x ]) ;
            ptr_sum[x] = func_r((LLR+N)[x], 0 );
            
            sign ^= ptr_sum[x] ; 

            int a = abs((LLR+N)[ x ]) ; 

			if(a < mina)
			{
				mina = a;
                idx_min = x ; 
			}
        }
            
        // printf("sign %d min %f  idx %d " ,sign,  mina,idx_min); 
        ptr_sum[idx_min]^=sign ;

        // update Res with G for next node 
        for( int x = 0;  x < N/2; x += 1 )
        {
            int16_t temp = func_g( ptr_sum[x] , (int16_t) LLR[ x ], (int16_t) (LLR+N/2)[ x ]) ; 
            (LLR+N)[ x ] =sat( temp)  ;  
        }


    } 
    else // opération normale 
    {   
        // if ici permet de mieux dérouler le code afterwards ( reduc insn reduc cycles ) ? 
        if((N/2)==8){                          
                                
            // ON CALCULE LES F
            for( int x = 0 ; x < N/2; x += 1 )
            {
                (LLR+N)[ x ] = func_f( LLR[ x ], (LLR+N/2)[ x ]) ; 
            }

            // Recursif Node_8 
            node_8( ptr_sum, (LLR+N), N/2, fz_bits) ; 

            // ON CALCULE LES G
            for( int x = 0;  x < N/2; x += 1 )
            {
                int16_t temp = func_g( ptr_sum[x] , (int16_t) LLR[ x ], (int16_t) (LLR+N/2)[ x ]) ; 
                (LLR+N)[ x ] =sat( temp)  ;  
            }
        }
        else{
            // ON CALCULE LES F
            for( int x = 0 ; x < N/2; x += 1 )
            {
                (LLR+N)[ x ] = func_f( LLR[ x ], (LLR+N/2)[ x ]) ; 
            }
            
            node( ptr_sum, (LLR+N), N/2, fz_bits);
                
            // ON CALCULE LES G
            for( int x = 0;  x < N/2; x += 1 )
            {
                int16_t temp = func_g( ptr_sum[x] , (int16_t) LLR[ x ], (int16_t) (LLR+N/2)[ x ]) ; 
                (LLR+N)[ x ] =sat( temp)  ;  
            }
        }
        
    }
        

    // ON CALCULE LA BRANCHE DROITE
    // get the node status 
    char enab2 = *fb_table_tileN++  ; 
    
    if(enab2==0) // r0 
    {
        //gestion de ptr ... 

    }else 
    if (enab2==1) // r1 
    {
        // copy leaf calc to PS  
        for( int x = 0 ; x < N/2; x += 1 )
        {
            // +N puisque ne saute pas dans un nouveau node 
            // force r1 
            (ptr_sum+N/2)[x] = func_r( (LLR+N)[x], 0 );
        }

    }else 
    if (enab2==2) // REP
    {
        // Somme LLRS et PS update 
        int tot = 0 ; 
        for (int x = 0; x < N/2; x++)
            tot += (LLR+N)[x] ; 

        if(tot <= 0 )
        {
            for( int x = 0 ; x < N/2; x += 1 )
                (ptr_sum+N/2)[x] = 1 ; 
        } 
    }else 
    if (enab2==3) // SPC
    {
        // printf("\n SPC droit \n "); 
        // parité + sign 
        int8_t sign=0 ; 
        int mina = 100  ; 
        int8_t idx_min= 0 ; 
        for( int x = 0 ; x < N/2; x += 1 )
        {
            // +N puisque ne saute pas dans un nouveau node 
            // force R1 (ne tiens pas compte des bits gelés )
            (ptr_sum+N/2)[x] = func_r( (LLR+N)[x], 0 );

            
            sign ^= (ptr_sum+N/2)[x] ; 
            int a = abs((LLR+N)[ x ]) ; 
            
			if(a < mina)
			{	
				mina = a;
                idx_min = x ; 
			}
        }
            
        // printf("sign %d min %f  idx %d " ,sign,  mina,idx_min);
        ptr_sum[idx_min]^=sign ;
    }       
    else // opération normale 
    {
        if( (N/2)==8 ){
            node_8(ptr_sum+ N/2, (LLR+N), N/2,  fz_bits+ N/2 );
        }else{
            node(ptr_sum+ N/2, (LLR+N), N/2,  fz_bits+ N/2 );
        }
    }  

    // ON FAIT LES CALCUL DES H (XOR DES SP)
    for(int x = 0 ; x < N/2 ; x += 1 )
    {          
        ptr_sum[x] ^=ptr_sum[ x + (N/2) ];    
    }
}

// calcule unqiuement le premier f ,  g  et le dernier Xor(h)
void node_top(int8_t* ptr_sum, int8_t *LLR , int N, char *fz_bits)
{

    // if (*fb_table_tileN != 4) {
    //     printf("(EE) TOP LEVEL TILE IS FROZEN !\n");
    //     exit(0);
    // }

    fb_table_tileN+=1 ; 

    // une fois F pr la branche gauche 
    for( int x = 0 ; x < N/2; x += 1 )
    {
        (LLR+N)[ x ] = func_f( LLR[ x ], (LLR+N/2)[ x ]) ; 
    }

    // tt la branche gauche 
    const int not_frozen_value = 4;
    if (*fb_table_tileN++ != not_frozen_value) {
        printf("(EE) Un truc impossible vient de se produire (1:%d)\n", *(fb_table_tileN-1));
        
    } else {
        node( ptr_sum, (LLR+N), N/2, fz_bits);
    }

    // une fois G pour la branche droite   
    for( int x = 0;  x < N/2; x += 1 )
    {
        int16_t temp = func_g( ptr_sum[x] , (int16_t) LLR[ x ], (int16_t) (LLR+N/2)[ x ]) ; 
        (LLR+N)[ x ] =sat( temp)  ;  
    }

    // tt la branche droite 
    if (*fb_table_tileN++ != not_frozen_value) {
        printf("(EE) Un truc impossible vient de se produire (1:%d)\n", *(fb_table_tileN-1));
        
    } else {
         node( ptr_sum+ N/2, (LLR+N), N/2,  fz_bits+ N/2 );
    }

     // ON FAIT LES CALCUL DES H (XOR DES SP)
    for(int x = 0 ; x < N/2 ; x += 1 )
    {          
        // xor
        ptr_sum[x] ^=ptr_sum[ x + (N/2) ]; 
    }

}

// main 
int main() 
{

    int8_t LLR[2*codw_N]       = {0} ; 
    int8_t PS[codw_N]          = {0};    
    int8_t decode_out_[K ]     = {0} ; 

    for (int i = 0; i < codw_N; i++)    
        LLR[i] = codw[i] ;

    perfs_init();
    // ------------
    // node 8   => standard SC 
    // node top => Fast 
    #ifdef FSC
	    node_top( PS, LLR, codw_N, froozen_bits );
    #else 
        node_8( PS, LLR, codw_N, froozen_bits );
    #endif 
    // ------------
	perfs_end();
    display_perf(); 

    return 0 ; 
}
