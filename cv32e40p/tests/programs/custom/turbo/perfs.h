// PERFS FUNCTIONS 
    // dealing with 64 bits csr


long insn_start,insn_stop,insn_tot ; 
long cycle_start,cycle_stop , cycle_tot ; 



#ifdef _X86_

uint64_t rdtsc(){
    unsigned int lo,hi;
    __asm__ __volatile__ ("rdtsc" : "=a" (lo), "=d" (hi));
    return ((uint64_t)hi << 32) | lo;
}

long cycles(){
    unsigned int lo,hi;
    __asm__ __volatile__ ("rdtsc" : "=a" (lo), "=d" (hi));
    return ((uint64_t)hi << 32) | lo;
}

void perfs_init(){
    // insn_start  = insn() - 4; 
    cycle_start = cycles() - 4 ; 
}

void perfs_end(){
    cycle_stop= cycles()-4;
	// insn_stop = insn()-4; 

	cycle_tot = cycle_stop - cycle_start ; 
	// insn_tot = insn_stop - insn_start ; 
}

#else 
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

#endif 

void display_perf(){
    printf("#=======================|\n");
	printf("Bench\n"); 
	printf("cycles: %ld\n", cycle_tot) ; 
	printf("insn: %ld\n", insn_tot); 
	printf("#=======================|\n");
}

