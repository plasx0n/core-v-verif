// PERFS FUNCTIONS 
    // dealing with 64 bits csr
long cycles(){
	long cycles;
	asm volatile ("rdcycle %0" : "=r"(cycles));
	// printf("[time() -> %d]", cycles);
	return cycles;
}

long insn(){
	long insns;
	asm volatile ("rdinstret %0" : "=r"(insns));
	// printf("[insn() -> %d]", insns);
	return insns;
}

long insn_start,insn_stop,insn_tot ; 
long cycle_start,cycle_stop , cycle_tot ; 


void perfs_init(){
    insn_start  = insn() - 4; 
    cycle_start = cycles() - 4 ; 
}

void perfs_end(){
    cycle_stop= cycles()-4;
	insn_stop = insn()-4; 

	cycle_tot = cycle_stop - cycle_start ; 
	insn_tot = insn_stop - insn_start ; 
}

void display_perf(){
    printf("|=======================|\n");
	printf("|Bench\n"); 
	printf("cycles: %ld\n", cycle_tot) ; 
	printf("insn: %ld\n", insn_tot); 
	printf("|=======================|\n");
}

void displayVectorV1( int32_t v1)
{
    printf("op V[7:0]   %d \n", (v1   <<24)>>24  );
    printf("op V[15:8]  %d \n", (v1   <<16)>>24  );
    printf("op V[23:16] %d \n", (v1   <<8 )>>24  );
    printf("op V[31:24] %d \n", (v1       )>>24  );
}

void displayVector( int32_t v1, int32_t v2, int32_t v3)
{
    printf("op V[7:0]   (%d,%d) = %d \n", (v1   <<24)>>24 , ( v2   <<24)>>24 , ( v3   <<24)>>24 );
    printf("op V[15:8]  (%d,%d) = %d \n", (v1   <<16)>>24 , ( v2   <<16)>>24 , ( v3   <<16)>>24 );
    printf("op V[23:16] (%d,%d) = %d \n", (v1   <<8 )>>24 , ( v2   <<8 )>>24 , ( v3   <<8 )>>24 );
    printf("op V[31:24] (%d,%d) = %d \n", (v1       )>>24 , ( v2       )>>24 , ( v3       )>>24 );
}

