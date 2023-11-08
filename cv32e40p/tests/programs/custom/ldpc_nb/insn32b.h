
#define callAdd32sSat64(rd,rs1,rs2) asm volatile("ldn.add32s_sat64 %0,%1,%2" \
	                            : "=r" (rd) \
	                            : "r" (rs1), "r" (rs2)); 

#define callSub32sSat64(rd,rs1,rs2) asm volatile("ldn.sub32s_sat64 %0,%1,%2" \
	                            : "=r" (rd) \
	                            : "r" (rs1), "r" (rs2));         

#define callMin(rd,rs1,rs2) asm volatile("ldn.min %0,%1,%2" \
	                            : "=r" (rd) \
	                            : "r" (rs1), "r" (rs2)); 

#define callAdd32uSat64(rd,rs1,rs2) asm volatile("ldn.add32u_sat64 %0,%1,%2" \
	                            : "=r" (rd) \
	                            : "r" (rs1), "r" (rs2)); 

#define callSat8u(rd,rs1,rs2) asm volatile("ldn.sat8 %0,%1,%2" \
	                            : "=r" (rd) \
	                            : "r" (rs1), "r" (rs2)); 