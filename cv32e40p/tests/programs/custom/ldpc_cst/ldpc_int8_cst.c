#include <stdlib.h>
#include <stdint.h>

#define CODE   	("LDPC")
#define ALGO	("MPA - MIN-SUM")
#define ordo   	("Horizontal layered")

#define iter 	10
#define _8b 

#include "perfs.h"
#include "insn_2r_32b.h"
#include "trame_34_20_BG2.h"

void process()
{
	// accu à conserver pour les iterations 
	int8_t Resu[32]  ;

		for(int l=0;l<iter;l++)
		{
			// garder les ptrs en data type 
			int8_t* ptr_posVn = posVn ;
			int8_t* ptr_c2v   = c2v ;

			// parcours des CN
			for( int idex_Cn = 0 ; idex_Cn < nb_CN ; idex_Cn++)
			{

				int8_t min1    = INT8_MAX ;
				int8_t min2    = INT8_MAX ;
				int8_t sign    =   0 ;

				// parcours des VN liés au Cn courant
				// on conserve une data en int pour la comparaison et 
				// int8_t idex_vn baisse les perfos 
				// compilateur préfère les int pour la loop
				int degCn = deg_Cns[idex_Cn];
				for( int idex_Vn =0 ; idex_Vn < degCn ; idex_Vn++)
				{
					int8_t vAccu ;
					
					int8_t a ;

					int8_t indice = ptr_posVn[ idex_Vn ];
					int8_t pVn  =	accuVn[ indice];
					int8_t msg  =	ptr_c2v [ idex_Vn ];
					

						callSubSat(vAccu,pVn,msg);				

					Resu[idex_Vn] =  vAccu; 

					// check min & signe ;
					// int8_t testacc =  
					// andi	s8,t3,255
					// sb	t3,0(a6)
					// srli	s8,s8,0x7

					sign  ^=  ( vAccu < 0); 
	
					// min casse la séquence car force slli & srai 

						callAbs(a,vAccu,0); 



						int8_t min_temp ;
						callMax(min_temp,min1,a) ; 
						callMin(min2, min2, min_temp )  ;   
						callMin(min1, a,min1) ; 


				}

				// parcours des VN liés au Cn courant
				for( int idex_Vn =0 ; idex_Vn < degCn ; idex_Vn++)
				{
					int8_t nMessage ;
					int8_t eval ; 
					int8_t Rsign; ; 

					int8_t temp = Resu[idex_Vn] ; 

					// idem avec eval qui slli & srai 

						callEval(eval,min1,temp); 
						callRsign(Rsign,sign,temp) ;


					// generation du mask + min à sortir
					int8_t min_t = min1 & ~eval ; 
					int8_t min_u = min2 & eval ; 
					int8_t min_  = min_t | min_u ; 
 
						callNmess(nMessage,Rsign,min_ ) ;

					// maj c2v
					ptr_c2v[idex_Vn] = nMessage ;

						callAddSat(temp,temp, nMessage) ;
 

					int8_t    indice = ptr_posVn[ idex_Vn ];
					accuVn[ indice ] = temp ;
				}
				
				ptr_posVn   += degCn;
				ptr_c2v     += degCn;
			}
		}
}

int main() {

	perfs_init();

   	process() ; 

	perfs_end();
	display_perf();

	return 0; 
}