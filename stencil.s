	.file	"stencil_template_parallel.c"
	.intel_syntax noprefix
# GNU C17 (Ubuntu 13.3.0-6ubuntu2~24.04) version 13.3.0 (x86_64-linux-gnu)
#	compiled by GNU C version 13.3.0, GMP version 6.3.0, MPFR version 4.2.1, MPC version 1.3.1, isl version isl-0.26-GMP

# GGC heuristics: --param ggc-min-expand=100 --param ggc-min-heapsize=131072
# options passed: -march=alderlake -mmmx -mpopcnt -msse -msse2 -msse3 -mssse3 -msse4.1 -msse4.2 -mavx -mavx2 -mno-sse4a -mno-fma4 -mno-xop -mfma -mno-avx512f -mbmi -mbmi2 -maes -mpclmul -mno-avx512vl -mno-avx512bw -mno-avx512dq -mno-avx512cd -mno-avx512er -mno-avx512pf -mno-avx512vbmi -mno-avx512ifma -mno-avx5124vnniw -mno-avx5124fmaps -mno-avx512vpopcntdq -mno-avx512vbmi2 -mgfni -mvpclmulqdq -mno-avx512vnni -mno-avx512bitalg -mno-avx512bf16 -mno-avx512vp2intersect -mno-3dnow -madx -mabm -mno-cldemote -mclflushopt -mclwb -mno-clzero -mcx16 -mno-enqcmd -mf16c -mfsgsbase -mfxsr -mno-hle -msahf -mno-lwp -mlzcnt -mmovbe -mmovdir64b -mmovdiri -mno-mwaitx -mno-pconfig -mno-pku -mno-prefetchwt1 -mprfchw -mno-ptwrite -mrdpid -mrdrnd -mrdseed -mno-rtm -mserialize -mno-sgx -msha -mshstk -mno-tbm -mno-tsxldtrk -mvaes -mwaitpkg -mno-wbnoinvd -mxsave -mxsavec -mxsaveopt -mxsaves -mno-amx-tile -mno-amx-int8 -mno-amx-bf16 -mno-uintr -mno-hreset -mno-kl -mno-widekl -mavxvnni -mno-avx512fp16 -mno-avxifma -mno-avxvnniint8 -mno-avxneconvert -mno-cmpccxadd -mno-amx-fp16 -mno-prefetchi -mno-raoint -mno-amx-complex --param=l1-cache-size=48 --param=l1-cache-line-size=64 --param=l2-cache-size=12288 -mtune=alderlake -masm=intel -O3 -fopenmp -fasynchronous-unwind-tables -fstack-protector-strong -fstack-clash-protection -fcf-protection
	.text
	.p2align 4
	.type	update_plane._omp_fn.0, @function
update_plane._omp_fn.0:
.LFB61:
	.cfi_startproc
	endbr64	
	push	r15	#
	.cfi_def_cfa_offset 16
	.cfi_offset 15, -16
	push	r14	#
	.cfi_def_cfa_offset 24
	.cfi_offset 14, -24
	push	r13	#
	.cfi_def_cfa_offset 32
	.cfi_offset 13, -32
	push	r12	#
	.cfi_def_cfa_offset 40
	.cfi_offset 12, -40
	push	rbp	#
	.cfi_def_cfa_offset 48
	.cfi_offset 6, -48
	push	rbx	#
	.cfi_def_cfa_offset 56
	.cfi_offset 3, -56
	sub	rsp, 24	#,
	.cfi_def_cfa_offset 80
# include/stencil_template_parallel.h:210:     #pragma omp parallel for schedule(static)
	mov	ebx, DWORD PTR 24[rdi]	# ysize, *.omp_data_i_9(D).ysize
	lea	eax, 1[rbx]	# tmp143,
	cmp	eax, 1	# tmp143,
	jbe	.L14	#,
	mov	r12, rdi	# .omp_data_i, tmp175
	call	omp_get_num_threads@PLT	#
	mov	ebp, eax	# _16, tmp176
	call	omp_get_thread_num@PLT	#
	xor	edx, edx	# tt.18_2
	mov	ecx, eax	# _19, tmp177
	mov	eax, ebx	# ysize, ysize
	div	ebp	# _16
	cmp	ecx, edx	# _19, tt.18_2
	jb	.L3	#,
.L8:
	imul	ecx, eax	# tmp146, q.17_1
	add	ecx, edx	# _24, tt.18_2
	add	eax, ecx	# _25, _24
	cmp	ecx, eax	# _24, _25
	jnb	.L14	#,
	mov	ebx, DWORD PTR 20[r12]	# xsize, *.omp_data_i_9(D).xsize
	add	ecx, 1	# j,
	mov	rdx, QWORD PTR [r12]	# old, *.omp_data_i_9(D).old
	lea	r8d, 1[rax]	# _27,
	mov	rbp, QWORD PTR 8[r12]	# new, *.omp_data_i_9(D).new
	mov	r12d, DWORD PTR 16[r12]	# fxsize, *.omp_data_i_9(D).fxsize
	test	ebx, ebx	# xsize
	je	.L14	#,
	mov	edi, r12d	# ivtmp.158, fxsize
# include/stencil_template_parallel.h:212:         for ( uint i = 1; i <= xsize; i++)
	mov	r13d, 1	# tmp169,
	mov	esi, r8d	# _27, _27
	vmovsd	xmm2, QWORD PTR .LC0[rip]	# tmp171,
	imul	edi, ecx	# ivtmp.158, j
	vmovsd	xmm1, QWORD PTR .LC1[rip]	# tmp172,
	.p2align 4,,10
	.p2align 3
.L7:
# include/stencil_template_parallel.h:225:                                 old[IDX(i, j-1)] + old[IDX(i, j+1)] ) /4.0 / 2.0;
	mov	r8d, edi	# _48, ivtmp.158
	add	ecx, 1	# j,
	mov	eax, edi	# ivtmp.158, ivtmp.158
# include/stencil_template_parallel.h:225:                                 old[IDX(i, j-1)] + old[IDX(i, j+1)] ) /4.0 / 2.0;
	add	edi, r12d	# ivtmp.158, fxsize
# include/stencil_template_parallel.h:225:                                 old[IDX(i, j-1)] + old[IDX(i, j+1)] ) /4.0 / 2.0;
	sub	r8d, r12d	# _48, fxsize
# include/stencil_template_parallel.h:212:         for ( uint i = 1; i <= xsize; i++)
	mov	DWORD PTR 12[rsp], ecx	# %sfp, j
	lea	r11d, 1[rdi]	# tmp174,
	mov	r15d, r13d	# tmp168, tmp169
	lea	r14d, 1[r8]	# tmp173,
# include/stencil_template_parallel.h:225:                                 old[IDX(i, j-1)] + old[IDX(i, j+1)] ) /4.0 / 2.0;
	sub	r11d, eax	# tmp159, ivtmp.146
# include/stencil_template_parallel.h:212:         for ( uint i = 1; i <= xsize; i++)
	sub	r15d, eax	# tmp168, ivtmp.146
# include/stencil_template_parallel.h:225:                                 old[IDX(i, j-1)] + old[IDX(i, j+1)] ) /4.0 / 2.0;
	sub	r14d, eax	# tmp154, ivtmp.146
	.p2align 4,,10
	.p2align 3
.L6:
	mov	r9d, eax	#, ivtmp.146
# include/stencil_template_parallel.h:224:             new[ IDX(i,j) ] = old[ IDX(i,j) ] / 2.0 + ( old[IDX(i-1, j)] + old[IDX(i+1, j)] +
	lea	r8d, 1[rax]	#,
	mov	rcx, r9	#,
# include/stencil_template_parallel.h:224:             new[ IDX(i,j) ] = old[ IDX(i,j) ] / 2.0 + ( old[IDX(i-1, j)] + old[IDX(i+1, j)] +
	lea	r10d, 2[r9]	# tmp149,
# include/stencil_template_parallel.h:224:             new[ IDX(i,j) ] = old[ IDX(i,j) ] / 2.0 + ( old[IDX(i-1, j)] + old[IDX(i+1, j)] +
	mov	rax, r8	#,
# include/stencil_template_parallel.h:225:                                 old[IDX(i, j-1)] + old[IDX(i, j+1)] ) /4.0 / 2.0;
	add	ecx, r11d	# tmp161, tmp159
# include/stencil_template_parallel.h:224:             new[ IDX(i,j) ] = old[ IDX(i,j) ] / 2.0 + ( old[IDX(i-1, j)] + old[IDX(i+1, j)] +
	vmovsd	xmm0, QWORD PTR [rdx+r10*8]	# *_44, *_44
	vaddsd	xmm0, xmm0, QWORD PTR [rdx+r9*8]	# tmp151, *_44, *_39
# include/stencil_template_parallel.h:225:                                 old[IDX(i, j-1)] + old[IDX(i, j+1)] ) /4.0 / 2.0;
	lea	r9d, [r14+r9]	# tmp156,
# include/stencil_template_parallel.h:225:                                 old[IDX(i, j-1)] + old[IDX(i, j+1)] ) /4.0 / 2.0;
	mov	ecx, ecx	# tmp161, tmp161
# include/stencil_template_parallel.h:224:             new[ IDX(i,j) ] = old[ IDX(i,j) ] / 2.0 + ( old[IDX(i-1, j)] + old[IDX(i+1, j)] +
	vaddsd	xmm0, xmm0, QWORD PTR [rdx+r9*8]	# tmp157, tmp151, *_52
# include/stencil_template_parallel.h:225:                                 old[IDX(i, j-1)] + old[IDX(i, j+1)] ) /4.0 / 2.0;
	vaddsd	xmm0, xmm0, QWORD PTR [rdx+rcx*8]	# tmp162, tmp157, *_60
# include/stencil_template_parallel.h:212:         for ( uint i = 1; i <= xsize; i++)
	lea	ecx, [r15+r8]	# i,
	cmp	ebx, ecx	# xsize, i
# include/stencil_template_parallel.h:225:                                 old[IDX(i, j-1)] + old[IDX(i, j+1)] ) /4.0 / 2.0;
	vmulsd	xmm0, xmm0, xmm2	# tmp163, tmp162, tmp171
# include/stencil_template_parallel.h:225:                                 old[IDX(i, j-1)] + old[IDX(i, j+1)] ) /4.0 / 2.0;
	vmulsd	xmm0, xmm0, xmm1	# tmp165, tmp163, tmp172
# include/stencil_template_parallel.h:224:             new[ IDX(i,j) ] = old[ IDX(i,j) ] / 2.0 + ( old[IDX(i-1, j)] + old[IDX(i+1, j)] +
	vfmadd231sd	xmm0, xmm1, QWORD PTR [rdx+r8*8]	# _66, tmp172, *_33
# include/stencil_template_parallel.h:224:             new[ IDX(i,j) ] = old[ IDX(i,j) ] / 2.0 + ( old[IDX(i-1, j)] + old[IDX(i+1, j)] +
	vmovsd	QWORD PTR 0[rbp+r8*8], xmm0	# *_65, _66
# include/stencil_template_parallel.h:212:         for ( uint i = 1; i <= xsize; i++)
	jnb	.L6	#,
	mov	ecx, DWORD PTR 12[rsp]	# j, %sfp
	cmp	ecx, esi	# j, _27
	jb	.L7	#,
.L14:
# include/stencil_template_parallel.h:210:     #pragma omp parallel for schedule(static)
	add	rsp, 24	#,
	.cfi_remember_state
	.cfi_def_cfa_offset 56
	pop	rbx	#
	.cfi_def_cfa_offset 48
	pop	rbp	#
	.cfi_def_cfa_offset 40
	pop	r12	#
	.cfi_def_cfa_offset 32
	pop	r13	#
	.cfi_def_cfa_offset 24
	pop	r14	#
	.cfi_def_cfa_offset 16
	pop	r15	#
	.cfi_def_cfa_offset 8
	ret	
.L3:
	.cfi_restore_state
	add	eax, 1	# q.17_1,
# include/stencil_template_parallel.h:210:     #pragma omp parallel for schedule(static)
	xor	edx, edx	# tt.18_2
	jmp	.L8	#
	.cfi_endproc
.LFE61:
	.size	update_plane._omp_fn.0, .-update_plane._omp_fn.0
	.p2align 4
	.type	get_total_energy._omp_fn.0, @function
get_total_energy._omp_fn.0:
.LFB62:
	.cfi_startproc
	endbr64	
	push	r14	#
	.cfi_def_cfa_offset 16
	.cfi_offset 14, -16
	push	r13	#
	.cfi_def_cfa_offset 24
	.cfi_offset 13, -24
	push	r12	#
	.cfi_def_cfa_offset 32
	.cfi_offset 12, -32
	mov	r12, rdi	# .omp_data_i, tmp181
	push	rbp	#
	.cfi_def_cfa_offset 40
	.cfi_offset 6, -40
	push	rbx	#
	.cfi_def_cfa_offset 48
	.cfi_offset 3, -48
# include/stencil_template_parallel.h:288:     #pragma omp parallel for reduction(+:totenergy) schedule(static)
	mov	rbp, QWORD PTR [rdi]	# data, *.omp_data_i_11(D).data
	mov	r13d, DWORD PTR 24[rdi]	# fsize, *.omp_data_i_11(D).fsize
	mov	ebx, DWORD PTR 16[rdi]	# xsize, *.omp_data_i_11(D).xsize
	call	omp_get_num_threads@PLT	#
	mov	r14d, eax	# _16, tmp182
	call	omp_get_thread_num@PLT	#
	mov	ecx, eax	# _17, tmp183
	mov	eax, DWORD PTR 20[r12]	# *.omp_data_i_11(D).ysize, *.omp_data_i_11(D).ysize
	cdq
	idiv	r14d	# _16
	cmp	ecx, edx	# _17, tt.23_5
	jl	.L18	#,
.L28:
	imul	ecx, eax	# tmp160, q.22_4
	vxorpd	xmm0, xmm0, xmm0	# totenergy
	add	ecx, edx	# _22, tt.23_5
	add	eax, ecx	# _23, _22
	cmp	ecx, eax	# _22, _23
	jge	.L19	#,
	add	ecx, 1	# j,
	mov	r8d, ebx	# bnd.167, xsize
	mov	r10d, ebx	# tmp163, xsize
	mov	r11d, ebx	# tmp179, xsize
	mov	esi, ecx	# ivtmp.187, j
	shr	r8d, 2	#,
	and	r10d, -4	# tmp163,
	lea	edi, 1[rax]	# _104,
	imul	esi, r13d	# ivtmp.187, fsize
	lea	r9d, -1[rbx]	# _53,
	sal	r8, 5	# _66,
	add	r10d, 1	# tmp.169,
	and	r11d, 3	# tmp179,
	.p2align 4,,10
	.p2align 3
.L21:
# include/stencil_template_parallel.h:290:         for ( int i = 1; i <= xsize; i++ ) {
	test	ebx, ebx	# xsize
	jle	.L25	#,
	cmp	r9d, 2	# _53,
	jbe	.L30	#,
	movsx	rax, esi	# _26, ivtmp.187
	lea	rax, 8[rbp+rax*8]	# ivtmp.178,
	lea	rdx, [rax+r8]	# _63,
	.p2align 4,,10
	.p2align 3
.L23:
	vaddsd	xmm0, xmm0, QWORD PTR [rax]	# stmp_totenergy_32.173, totenergy, BIT_FIELD_REF <MEM <vector(4) double> [(double *)_56], 64, 0>
	add	rax, 32	# ivtmp.178,
	vaddsd	xmm0, xmm0, QWORD PTR -24[rax]	# stmp_totenergy_32.173, stmp_totenergy_32.173, BIT_FIELD_REF <MEM <vector(4) double> [(double *)_56], 64, 64>
# include/stencil_template_parallel.h:291:             totenergy += data[ IDX(i, j) ];
	vaddsd	xmm0, xmm0, QWORD PTR -16[rax]	# stmp_totenergy_32.173, stmp_totenergy_32.173, BIT_FIELD_REF <MEM <vector(4) double> [(double *)_56], 64, 128>
	vaddsd	xmm0, xmm0, QWORD PTR -8[rax]	# totenergy, stmp_totenergy_32.173, BIT_FIELD_REF <MEM <vector(4) double> [(double *)_56], 64, 192>
	cmp	rax, rdx	# ivtmp.178, _63
	jne	.L23	#,
	test	r11d, r11d	# tmp179
# include/stencil_template_parallel.h:290:         for ( int i = 1; i <= xsize; i++ ) {
	mov	eax, r10d	# i, tmp.169
	je	.L25	#,
.L22:
# include/stencil_template_parallel.h:291:             totenergy += data[ IDX(i, j) ];
	lea	edx, [rsi+rax]	# tmp169,
	movsx	rdx, edx	# tmp170, tmp169
# include/stencil_template_parallel.h:291:             totenergy += data[ IDX(i, j) ];
	vaddsd	xmm0, xmm0, QWORD PTR 0[rbp+rdx*8]	# totenergy, totenergy, *_31
# include/stencil_template_parallel.h:290:         for ( int i = 1; i <= xsize; i++ ) {
	lea	edx, 1[rax]	# i,
# include/stencil_template_parallel.h:290:         for ( int i = 1; i <= xsize; i++ ) {
	cmp	ebx, edx	# xsize, i
	jl	.L25	#,
# include/stencil_template_parallel.h:291:             totenergy += data[ IDX(i, j) ];
	add	edx, esi	# tmp171, ivtmp.187
# include/stencil_template_parallel.h:290:         for ( int i = 1; i <= xsize; i++ ) {
	add	eax, 2	# i,
# include/stencil_template_parallel.h:291:             totenergy += data[ IDX(i, j) ];
	movsx	rdx, edx	# tmp172, tmp171
# include/stencil_template_parallel.h:290:         for ( int i = 1; i <= xsize; i++ ) {
	cmp	ebx, eax	# xsize, i
# include/stencil_template_parallel.h:291:             totenergy += data[ IDX(i, j) ];
	vaddsd	xmm0, xmm0, QWORD PTR 0[rbp+rdx*8]	# totenergy, totenergy, *_90
# include/stencil_template_parallel.h:290:         for ( int i = 1; i <= xsize; i++ ) {
	jl	.L25	#,
# include/stencil_template_parallel.h:291:             totenergy += data[ IDX(i, j) ];
	add	eax, esi	# tmp173, ivtmp.187
	cdqe
# include/stencil_template_parallel.h:291:             totenergy += data[ IDX(i, j) ];
	vaddsd	xmm0, xmm0, QWORD PTR 0[rbp+rax*8]	# totenergy, totenergy, *_59
.L25:
	add	ecx, 1	# j,
	add	esi, r13d	# ivtmp.187, fsize
	cmp	ecx, edi	# j, _104
	jne	.L21	#,
.L19:
	mov	rdx, QWORD PTR 8[r12]	# _8,
# include/stencil_template_parallel.h:288:     #pragma omp parallel for reduction(+:totenergy) schedule(static)
	lea	rcx, 8[r12]	# _34,
.L27:
# include/stencil_template_parallel.h:288:     #pragma omp parallel for reduction(+:totenergy) schedule(static)
	vmovq	xmm2, rdx	# tmp190, _8
	mov	rax, rdx	# tmp176, _8
	vaddsd	xmm1, xmm0, xmm2	# tmp175, totenergy, tmp190
	vmovq	rsi, xmm1	# tmp175, tmp175
	lock cmpxchg	QWORD PTR [rcx], rsi	#,* _34, tmp175
	jne	.L40	#,
# include/stencil_template_parallel.h:288:     #pragma omp parallel for reduction(+:totenergy) schedule(static)
	pop	rbx	#
	.cfi_remember_state
	.cfi_def_cfa_offset 40
	pop	rbp	#
	.cfi_def_cfa_offset 32
	pop	r12	#
	.cfi_def_cfa_offset 24
	pop	r13	#
	.cfi_def_cfa_offset 16
	pop	r14	#
	.cfi_def_cfa_offset 8
	ret	
.L30:
	.cfi_restore_state
# include/stencil_template_parallel.h:290:         for ( int i = 1; i <= xsize; i++ ) {
	mov	eax, 1	# i,
	jmp	.L22	#
.L18:
	add	eax, 1	# q.22_4,
# include/stencil_template_parallel.h:288:     #pragma omp parallel for reduction(+:totenergy) schedule(static)
	xor	edx, edx	# tt.23_5
	jmp	.L28	#
.L40:
# include/stencil_template_parallel.h:288:     #pragma omp parallel for reduction(+:totenergy) schedule(static)
	mov	rdx, rax	# _8, tmp176
	jmp	.L27	#
	.cfi_endproc
.LFE62:
	.size	get_total_energy._omp_fn.0, .-get_total_energy._omp_fn.0
	.p2align 4
	.globl	inject_energy
	.type	inject_energy, @function
inject_energy:
.LFB48:
	.cfi_startproc
	endbr64	
	movsx	rax, esi	#, tmp187
	push	r12	#
	.cfi_def_cfa_offset 16
	.cfi_offset 12, -16
# include/stencil_template_parallel.h:129:     for (int s = 0; s < Nsources; s++) {
	test	eax, eax	# Nsources
# include/stencil_template_parallel.h:120: {
	push	rbp	#
	.cfi_def_cfa_offset 24
	.cfi_offset 6, -24
	push	rbx	#
	.cfi_def_cfa_offset 32
	.cfi_offset 3, -32
# include/stencil_template_parallel.h:122:     const uint register sizex = plane->size[_x_]+2; // interior + halos
	mov	ebp, DWORD PTR 8[rcx]	# _1, plane_50(D)->size[0]
# include/stencil_template_parallel.h:123:     double * restrict data = plane->data;
	mov	rsi, QWORD PTR [rcx]	# data, plane_50(D)->data
# include/stencil_template_parallel.h:129:     for (int s = 0; s < Nsources; s++) {
	jle	.L57	#,
	test	edi, edi	# periodic
	lea	r9d, 2[rbp]	# sizex,
	je	.L59	#,
# include/stencil_template_parallel.h:140:                 if ( (N[_x_] == 1)  ) {
	mov	r11d, DWORD PTR [r8]	# _14, *N_58(D)
# include/stencil_template_parallel.h:157:                 if ( (N[_y_] == 1) )
	mov	r12d, DWORD PTR 4[r8]	# _27, MEM[(const uint *)N_58(D) + 4B]
	mov	r10, rcx	# plane, tmp190
	lea	rbx, [rdx+rax*8]	# _79,
	jmp	.L53	#
	.p2align 4,,10
	.p2align 3
.L49:
# include/stencil_template_parallel.h:157:                 if ( (N[_y_] == 1) )
	cmp	r12d, 1	# _27,
	je	.L60	#,
.L48:
# include/stencil_template_parallel.h:129:     for (int s = 0; s < Nsources; s++) {
	add	rdx, 8	# ivtmp.202,
	cmp	rdx, rbx	# ivtmp.202, _79
	je	.L57	#,
.L53:
# include/stencil_template_parallel.h:133:             int y = Sources[s][_y_];
	mov	edi, DWORD PTR 4[rdx]	# _6, MEM[(unsigned int *)_69 + 4B]
# include/stencil_template_parallel.h:132:             int x = Sources[s][_x_];
	mov	ecx, DWORD PTR [rdx]	#, MEM[(unsigned int *)_69]
# include/stencil_template_parallel.h:136:             data[ IDX(x,y) ] += energy;
	mov	r8d, edi	# _7, _6
	imul	r8d, r9d	#, sizex
# include/stencil_template_parallel.h:140:                 if ( (N[_x_] == 1)  ) {
	cmp	r11d, 1	# _14,
# include/stencil_template_parallel.h:136:             data[ IDX(x,y) ] += energy;
	lea	eax, [rcx+r8]	# tmp161,
	lea	rax, [rsi+rax*8]	# _11,
# include/stencil_template_parallel.h:136:             data[ IDX(x,y) ] += energy;
	vaddsd	xmm1, xmm0, QWORD PTR [rax]	# tmp163, energy, *_11
	vmovsd	QWORD PTR [rax], xmm1	# *_11, tmp163
# include/stencil_template_parallel.h:140:                 if ( (N[_x_] == 1)  ) {
	jne	.L49	#,
# include/stencil_template_parallel.h:150:                     if ( x == 1 ){ // source on easr edge case
	cmp	ecx, 1	# _5,
	jne	.L50	#,
# include/stencil_template_parallel.h:151:                         data[IDX(plane->size[_x_]+1, y)] += energy;} // propagate on west halo
	lea	eax, 1[rbp+r8]	# tmp167,
	lea	rax, [rsi+rax*8]	# _19,
# include/stencil_template_parallel.h:151:                         data[IDX(plane->size[_x_]+1, y)] += energy;} // propagate on west halo
	vaddsd	xmm1, xmm0, QWORD PTR [rax]	# tmp169, energy, *_19
	vmovsd	QWORD PTR [rax], xmm1	# *_19, tmp169
.L50:
# include/stencil_template_parallel.h:153:                     if ( x == plane->size[_x_] ) { // source on west edge case
	cmp	ebp, ecx	# _1, _5
	jne	.L49	#,
# include/stencil_template_parallel.h:154:                         data[IDX(0, y)] += energy;} // propagate on east halo
	lea	rax, [rsi+r8*8]	# _24,
# include/stencil_template_parallel.h:157:                 if ( (N[_y_] == 1) )
	cmp	r12d, 1	# _27,
# include/stencil_template_parallel.h:154:                         data[IDX(0, y)] += energy;} // propagate on east halo
	vaddsd	xmm1, xmm0, QWORD PTR [rax]	# tmp173, energy, *_24
	vmovsd	QWORD PTR [rax], xmm1	# *_24, tmp173
# include/stencil_template_parallel.h:157:                 if ( (N[_y_] == 1) )
	jne	.L48	#,
	.p2align 4,,10
	.p2align 3
.L60:
# include/stencil_template_parallel.h:167:                     if ( y == 1 ){ // source on easr edge case
	cmp	edi, 1	# _6,
# include/stencil_template_parallel.h:168:                         data[IDX(x, plane->size[_y_]+1)] += energy;
	mov	eax, DWORD PTR 12[r10]	# pretmp_113, plane_50(D)->size[1]
# include/stencil_template_parallel.h:167:                     if ( y == 1 ){ // source on easr edge case
	jne	.L52	#,
# include/stencil_template_parallel.h:168:                         data[IDX(x, plane->size[_y_]+1)] += energy;
	lea	r8d, 1[rax]	# tmp175,
	imul	r8d, r9d	# tmp176, sizex
	add	r8d, ecx	# tmp178, _5
	mov	r8d, r8d	# tmp178, tmp178
	lea	r8, [rsi+r8*8]	# _34,
# include/stencil_template_parallel.h:168:                         data[IDX(x, plane->size[_y_]+1)] += energy;
	vaddsd	xmm1, xmm0, QWORD PTR [r8]	# tmp180, energy, *_34
	vmovsd	QWORD PTR [r8], xmm1	# *_34, tmp180
.L52:
# include/stencil_template_parallel.h:171:                     if ( y == plane->size[_y_] ) { // source on west edge case
	cmp	edi, eax	# _6, pretmp_113
	jne	.L48	#,
# include/stencil_template_parallel.h:172:                         data[IDX(x, 0)] += energy;
	lea	rax, [rsi+rcx*8]	# _40,
# include/stencil_template_parallel.h:129:     for (int s = 0; s < Nsources; s++) {
	add	rdx, 8	# ivtmp.202,
# include/stencil_template_parallel.h:172:                         data[IDX(x, 0)] += energy;
	vaddsd	xmm1, xmm0, QWORD PTR [rax]	# tmp184, energy, *_40
# include/stencil_template_parallel.h:129:     for (int s = 0; s < Nsources; s++) {
	cmp	rdx, rbx	# ivtmp.202, _79
# include/stencil_template_parallel.h:172:                         data[IDX(x, 0)] += energy;
	vmovsd	QWORD PTR [rax], xmm1	# *_40, tmp184
# include/stencil_template_parallel.h:129:     for (int s = 0; s < Nsources; s++) {
	jne	.L53	#,
.L57:
# include/stencil_template_parallel.h:180: }
	pop	rbx	#
	.cfi_remember_state
	.cfi_def_cfa_offset 24
	xor	eax, eax	#
	pop	rbp	#
	.cfi_def_cfa_offset 16
	pop	r12	#
	.cfi_def_cfa_offset 8
	ret	
	.p2align 4,,10
	.p2align 3
.L59:
	.cfi_restore_state
	lea	rcx, [rdx+rax*8]	# _74,
	.p2align 4,,10
	.p2align 3
.L44:
# include/stencil_template_parallel.h:136:             data[ IDX(x,y) ] += energy;
	mov	eax, DWORD PTR 4[rdx]	# tmp151, MEM[(unsigned int *)_28 + 4B]
# include/stencil_template_parallel.h:129:     for (int s = 0; s < Nsources; s++) {
	add	rdx, 8	# ivtmp.195,
# include/stencil_template_parallel.h:136:             data[ IDX(x,y) ] += energy;
	imul	eax, r9d	# tmp151, sizex
	add	eax, DWORD PTR -8[rdx]	# tmp154, MEM[(unsigned int *)_28]
# include/stencil_template_parallel.h:129:     for (int s = 0; s < Nsources; s++) {
	cmp	rdx, rcx	# ivtmp.195, _74
# include/stencil_template_parallel.h:136:             data[ IDX(x,y) ] += energy;
	lea	rax, [rsi+rax*8]	# _99,
# include/stencil_template_parallel.h:136:             data[ IDX(x,y) ] += energy;
	vaddsd	xmm1, xmm0, QWORD PTR [rax]	# tmp156, energy, *_99
	vmovsd	QWORD PTR [rax], xmm1	# *_99, tmp156
# include/stencil_template_parallel.h:129:     for (int s = 0; s < Nsources; s++) {
	jne	.L44	#,
# include/stencil_template_parallel.h:180: }
	pop	rbx	#
	.cfi_def_cfa_offset 24
	xor	eax, eax	#
	pop	rbp	#
	.cfi_def_cfa_offset 16
	pop	r12	#
	.cfi_def_cfa_offset 8
	ret	
	.cfi_endproc
.LFE48:
	.size	inject_energy, .-inject_energy
	.p2align 4
	.globl	update_plane
	.type	update_plane, @function
update_plane:
.LFB49:
	.cfi_startproc
	endbr64	
	push	r15	#
	.cfi_def_cfa_offset 16
	.cfi_offset 15, -16
	mov	r15d, edi	# periodic, tmp160
	lea	rdi, update_plane._omp_fn.0[rip]	# tmp141,
	push	r14	#
	.cfi_def_cfa_offset 24
	.cfi_offset 14, -24
	push	r13	#
	.cfi_def_cfa_offset 32
	.cfi_offset 13, -32
	push	r12	#
	.cfi_def_cfa_offset 40
	.cfi_offset 12, -40
	mov	r12, rsi	# N, tmp161
	push	rbp	#
	.cfi_def_cfa_offset 48
	.cfi_offset 6, -48
	push	rbx	#
	.cfi_def_cfa_offset 56
	.cfi_offset 3, -56
	sub	rsp, 56	#,
	.cfi_def_cfa_offset 112
# include/stencil_template_parallel.h:190:     uint register fxsize = oldplane->size[_x_]+2; // interior + halo
	mov	ebp, DWORD PTR 8[rdx]	# _1, oldplane_50(D)->size[0]
# include/stencil_template_parallel.h:208:     double * restrict new = newplane->data; // and write on this
	mov	rbx, QWORD PTR [rcx]	# new, newplane_53(D)->data
	xor	ecx, ecx	#
# include/stencil_template_parallel.h:189: {
	mov	rax, QWORD PTR fs:40	# tmp164, MEM[(<address-space-1> long unsigned int *)40B]
	mov	QWORD PTR 40[rsp], rax	# D.12225, tmp164
	xor	eax, eax	# tmp164
# include/stencil_template_parallel.h:210:     #pragma omp parallel for schedule(static)
	vmovq	xmm0, QWORD PTR [rdx]	# oldplane_50(D)->data, oldplane_50(D)->data
# include/stencil_template_parallel.h:191:     uint register fysize = oldplane->size[_y_]+2;
	mov	r13d, DWORD PTR 12[rdx]	# _2, oldplane_50(D)->size[1]
	mov	rsi, rsp	# tmp140,
# include/stencil_template_parallel.h:190:     uint register fxsize = oldplane->size[_x_]+2; // interior + halo
	lea	r14d, 2[rbp]	# fxsize,
	xor	edx, edx	#
# include/stencil_template_parallel.h:210:     #pragma omp parallel for schedule(static)
	vpinsrq	xmm0, xmm0, rbx, 1	# tmp137, oldplane_50(D)->data, new
	vmovd	xmm1, r14d	# fxsize, fxsize
	mov	DWORD PTR 24[rsp], r13d	# .omp_data_o.16.ysize, _2
	vmovdqa	XMMWORD PTR [rsp], xmm0	# MEM <vector(2) long unsigned int> [(double * *)&.omp_data_o.16], tmp137
	vpinsrd	xmm0, xmm1, ebp, 1	# tmp139, fxsize, _1
	vmovq	QWORD PTR 16[rsp], xmm0	# MEM <vector(2) unsigned int> [(unsigned int *)&.omp_data_o.16 + 16B], tmp139
	call	GOMP_parallel@PLT	#
# include/stencil_template_parallel.h:230:     if ( periodic ) {
	test	r15d, r15d	# periodic
	je	.L62	#,
# include/stencil_template_parallel.h:234:         if ( N[_x_] == 1 ) {
	mov	edx, DWORD PTR [r12]	# j, *N_63(D)
# include/stencil_template_parallel.h:234:         if ( N[_x_] == 1 ) {
	cmp	edx, 1	# j,
	je	.L79	#,
.L63:
# include/stencil_template_parallel.h:243:         if ( N[_y_] == 1 ) {
	mov	eax, DWORD PTR 4[r12]	# i, MEM[(const uint *)N_63(D) + 4B]
# include/stencil_template_parallel.h:243:         if ( N[_y_] == 1 ) {
	cmp	eax, 1	# i,
	je	.L80	#,
.L62:
# include/stencil_template_parallel.h:257: }
	mov	rax, QWORD PTR 40[rsp]	# tmp165, D.12225
	sub	rax, QWORD PTR fs:40	# tmp165, MEM[(<address-space-1> long unsigned int *)40B]
	jne	.L81	#,
	add	rsp, 56	#,
	.cfi_remember_state
	.cfi_def_cfa_offset 56
	xor	eax, eax	#
	pop	rbx	#
	.cfi_def_cfa_offset 48
	pop	rbp	#
	.cfi_def_cfa_offset 40
	pop	r12	#
	.cfi_def_cfa_offset 32
	pop	r13	#
	.cfi_def_cfa_offset 24
	pop	r14	#
	.cfi_def_cfa_offset 16
	pop	r15	#
	.cfi_def_cfa_offset 8
	ret	
	.p2align 4,,10
	.p2align 3
.L80:
	.cfi_restore_state
# include/stencil_template_parallel.h:246:             for (uint i = 1; i <= xsize; i++) {
	test	ebp, ebp	# _1
	je	.L62	#,
# include/stencil_template_parallel.h:247:                 new[ IDX(i, 0) ] = new[ IDX(i, ysize) ]; // propagate south edge on the north halo
	imul	r13d, r14d	# _23, fxsize
# include/stencil_template_parallel.h:248:                 new[ IDX(i, ysize+1) ] = new[ IDX(i,1) ]; // propagate north edge on the south halo
	lea	ecx, 0[r13+r14]	# _37,
	.p2align 4,,10
	.p2align 3
.L65:
# include/stencil_template_parallel.h:247:                 new[ IDX(i, 0) ] = new[ IDX(i, ysize) ]; // propagate south edge on the north halo
	lea	edx, 0[r13+rax]	# tmp151,
# include/stencil_template_parallel.h:247:                 new[ IDX(i, 0) ] = new[ IDX(i, ysize) ]; // propagate south edge on the north halo
	vmovsd	xmm0, QWORD PTR [rbx+rdx*8]	# _31, *_27
# include/stencil_template_parallel.h:247:                 new[ IDX(i, 0) ] = new[ IDX(i, ysize) ]; // propagate south edge on the north halo
	mov	edx, eax	# i, i
# include/stencil_template_parallel.h:247:                 new[ IDX(i, 0) ] = new[ IDX(i, ysize) ]; // propagate south edge on the north halo
	vmovsd	QWORD PTR [rbx+rdx*8], xmm0	# *_30, _31
# include/stencil_template_parallel.h:248:                 new[ IDX(i, ysize+1) ] = new[ IDX(i,1) ]; // propagate north edge on the south halo
	lea	edx, [r14+rax]	# tmp154,
# include/stencil_template_parallel.h:248:                 new[ IDX(i, ysize+1) ] = new[ IDX(i,1) ]; // propagate north edge on the south halo
	vmovsd	xmm0, QWORD PTR [rbx+rdx*8]	# _42, *_35
# include/stencil_template_parallel.h:248:                 new[ IDX(i, ysize+1) ] = new[ IDX(i,1) ]; // propagate north edge on the south halo
	lea	edx, [rcx+rax]	# tmp156,
# include/stencil_template_parallel.h:246:             for (uint i = 1; i <= xsize; i++) {
	add	eax, 1	# i,
# include/stencil_template_parallel.h:246:             for (uint i = 1; i <= xsize; i++) {
	cmp	ebp, eax	# _1, i
# include/stencil_template_parallel.h:248:                 new[ IDX(i, ysize+1) ] = new[ IDX(i,1) ]; // propagate north edge on the south halo
	vmovsd	QWORD PTR [rbx+rdx*8], xmm0	# *_41, _42
# include/stencil_template_parallel.h:246:             for (uint i = 1; i <= xsize; i++) {
	jnb	.L65	#,
	jmp	.L62	#
	.p2align 4,,10
	.p2align 3
.L79:
# include/stencil_template_parallel.h:237:             for (uint j = 1; j <= ysize; j++) {
	test	r13d, r13d	# _2
	je	.L63	#,
	mov	eax, r14d	# ivtmp.230, fxsize
	lea	esi, 1[rbp]	# tmp159,
	.p2align 4,,10
	.p2align 3
.L64:
# include/stencil_template_parallel.h:238:                 new[ IDX(0,j) ] = new[ IDX(xsize,j) ]; // we propagate east edge on west halo
	lea	ecx, 0[rbp+rax]	# tmp143,
# include/stencil_template_parallel.h:237:             for (uint j = 1; j <= ysize; j++) {
	add	edx, 1	# j,
# include/stencil_template_parallel.h:238:                 new[ IDX(0,j) ] = new[ IDX(xsize,j) ]; // we propagate east edge on west halo
	vmovsd	xmm0, QWORD PTR [rbx+rcx*8]	# _12, *_8
# include/stencil_template_parallel.h:238:                 new[ IDX(0,j) ] = new[ IDX(xsize,j) ]; // we propagate east edge on west halo
	mov	ecx, eax	# ivtmp.230, ivtmp.230
# include/stencil_template_parallel.h:238:                 new[ IDX(0,j) ] = new[ IDX(xsize,j) ]; // we propagate east edge on west halo
	vmovsd	QWORD PTR [rbx+rcx*8], xmm0	# *_11, _12
# include/stencil_template_parallel.h:239:                 new[ IDX(xsize+1,j) ] = new[ IDX(1,j) ]; // we propagate west edge on east halo
	lea	ecx, 1[rax]	# tmp146,
# include/stencil_template_parallel.h:239:                 new[ IDX(xsize+1,j) ] = new[ IDX(1,j) ]; // we propagate west edge on east halo
	vmovsd	xmm0, QWORD PTR [rbx+rcx*8]	# _21, *_16
# include/stencil_template_parallel.h:239:                 new[ IDX(xsize+1,j) ] = new[ IDX(1,j) ]; // we propagate west edge on east halo
	lea	ecx, [rsi+rax]	# tmp149,
# include/stencil_template_parallel.h:237:             for (uint j = 1; j <= ysize; j++) {
	add	eax, r14d	# ivtmp.230, fxsize
	cmp	r13d, edx	# _2, j
# include/stencil_template_parallel.h:239:                 new[ IDX(xsize+1,j) ] = new[ IDX(1,j) ]; // we propagate west edge on east halo
	vmovsd	QWORD PTR [rbx+rcx*8], xmm0	# *_20, _21
# include/stencil_template_parallel.h:237:             for (uint j = 1; j <= ysize; j++) {
	jnb	.L64	#,
# include/stencil_template_parallel.h:243:         if ( N[_y_] == 1 ) {
	mov	eax, DWORD PTR 4[r12]	# i, MEM[(const uint *)N_63(D) + 4B]
# include/stencil_template_parallel.h:243:         if ( N[_y_] == 1 ) {
	cmp	eax, 1	# i,
	jne	.L62	#,
	jmp	.L80	#
.L81:
# include/stencil_template_parallel.h:257: }
	call	__stack_chk_fail@PLT	#
	.cfi_endproc
.LFE49:
	.size	update_plane, .-update_plane
	.p2align 4
	.globl	get_total_energy
	.type	get_total_energy, @function
get_total_energy:
.LFB50:
	.cfi_startproc
	endbr64	
	push	rbx	#
	.cfi_def_cfa_offset 16
	.cfi_offset 3, -16
	xor	ecx, ecx	#
	mov	rbx, rsi	# energy, tmp100
	xor	edx, edx	#
	sub	rsp, 48	#,
	.cfi_def_cfa_offset 64
# include/stencil_template_parallel.h:269:     const int register xsize = plane->size[_x_];
	vmovq	xmm0, QWORD PTR 8[rdi]	# vect__1.239, MEM <vector(2) unsigned int> [(unsigned int *)plane_4(D) + 8B]
# include/stencil_template_parallel.h:267: {
	mov	rax, QWORD PTR fs:40	# tmp101, MEM[(<address-space-1> long unsigned int *)40B]
	mov	QWORD PTR 40[rsp], rax	# D.12247, tmp101
	xor	eax, eax	# tmp101
# include/stencil_template_parallel.h:273:     double * restrict data = plane->data;
	mov	rax, QWORD PTR [rdi]	# data, plane_4(D)->data
	mov	rsi, rsp	# tmp95,
	lea	rdi, get_total_energy._omp_fn.0[rip]	# tmp96,
# include/stencil_template_parallel.h:288:     #pragma omp parallel for reduction(+:totenergy) schedule(static)
	vmovq	QWORD PTR 16[rsp], xmm0	# MEM <const vector(2) int> [(int *)&.omp_data_o.21 + 16B], vect__1.239
	mov	QWORD PTR [rsp], rax	# .omp_data_o.21.data, data
# include/stencil_template_parallel.h:271:     const int register fsize = xsize+2;
	vmovd	eax, xmm0	# tmp93, vect__1.239
# include/stencil_template_parallel.h:288:     #pragma omp parallel for reduction(+:totenergy) schedule(static)
	mov	QWORD PTR 8[rsp], 0x000000000	# .omp_data_o.21.totenergy,
# include/stencil_template_parallel.h:271:     const int register fsize = xsize+2;
	add	eax, 2	# tmp94,
	mov	DWORD PTR 24[rsp], eax	# .omp_data_o.21.fsize, tmp94
	call	GOMP_parallel@PLT	#
# include/stencil_template_parallel.h:288:     #pragma omp parallel for reduction(+:totenergy) schedule(static)
	vmovsd	xmm0, QWORD PTR 8[rsp]	# totenergy, .omp_data_o.21.totenergy
# include/stencil_template_parallel.h:297:     *energy = (double)totenergy;
	vmovsd	QWORD PTR [rbx], xmm0	# *energy_17(D), totenergy
# include/stencil_template_parallel.h:299: }
	mov	rax, QWORD PTR 40[rsp]	# tmp102, D.12247
	sub	rax, QWORD PTR fs:40	# tmp102, MEM[(<address-space-1> long unsigned int *)40B]
	jne	.L85	#,
	add	rsp, 48	#,
	.cfi_remember_state
	.cfi_def_cfa_offset 16
	xor	eax, eax	#
	pop	rbx	#
	.cfi_def_cfa_offset 8
	ret	
.L85:
	.cfi_restore_state
	call	__stack_chk_fail@PLT	#
	.cfi_endproc
.LFE50:
	.size	get_total_energy, .-get_total_energy
	.p2align 4
	.globl	fill_buffers
	.type	fill_buffers, @function
fill_buffers:
.LFB52:
	.cfi_startproc
	endbr64	
	push	r13	#
	.cfi_def_cfa_offset 16
	.cfi_offset 13, -16
	push	r12	#
	.cfi_def_cfa_offset 24
	.cfi_offset 12, -24
	push	rbp	#
	.cfi_def_cfa_offset 32
	.cfi_offset 6, -32
	push	rbx	#
	.cfi_def_cfa_offset 40
	.cfi_offset 3, -40
# src/stencil_template_parallel.c:217:                     const vec2_t N) {
	mov	ebx, edx	# tmp303, periodic
	mov	rdx, rcx	# N, tmp304
# src/stencil_template_parallel.c:221:   const size_t fx = (size_t)nx + 2;
	mov	ecx, DWORD PTR 8[rsi]	# _1, plane_58(D)->size[0]
# src/stencil_template_parallel.c:220:   const unsigned int ny = plane->size[_y_];
	mov	r8d, DWORD PTR 12[rsi]	#, plane_58(D)->size[1]
# src/stencil_template_parallel.c:226:   buffers[SEND][NORTH] = &plane->data[IDX(1,1)]; // north starts from position (1,1) (since we also have the halos in data)
	mov	rax, QWORD PTR [rsi]	# tmp216, plane_58(D)->data
# src/stencil_template_parallel.c:226:   buffers[SEND][NORTH] = &plane->data[IDX(1,1)]; // north starts from position (1,1) (since we also have the halos in data)
	lea	rbp, 3[rcx]	# _3,
# src/stencil_template_parallel.c:221:   const size_t fx = (size_t)nx + 2;
	lea	r9, 2[rcx]	# fx,
# src/stencil_template_parallel.c:226:   buffers[SEND][NORTH] = &plane->data[IDX(1,1)]; // north starts from position (1,1) (since we also have the halos in data)
	lea	r11, 0[0+rbp*8]	# _4,
# src/stencil_template_parallel.c:227:   buffers[SEND][SOUTH] = &plane->data[IDX(1,ny)]; // north starts from position (1,ny)
	mov	r12, r8	# tmp218, _7
# src/stencil_template_parallel.c:220:   const unsigned int ny = plane->size[_y_];
	mov	r10, r8	#,
# src/stencil_template_parallel.c:227:   buffers[SEND][SOUTH] = &plane->data[IDX(1,ny)]; // north starts from position (1,ny)
	imul	r12, r9	# tmp218, fx
# src/stencil_template_parallel.c:226:   buffers[SEND][NORTH] = &plane->data[IDX(1,1)]; // north starts from position (1,1) (since we also have the halos in data)
	add	rax, r11	# tmp216, _4
# src/stencil_template_parallel.c:226:   buffers[SEND][NORTH] = &plane->data[IDX(1,1)]; // north starts from position (1,1) (since we also have the halos in data)
	mov	QWORD PTR [rdi], rax	# (*buffers_62(D))[0], tmp216
# src/stencil_template_parallel.c:227:   buffers[SEND][SOUTH] = &plane->data[IDX(1,ny)]; // north starts from position (1,ny)
	mov	rax, QWORD PTR [rsi]	# plane_58(D)->data, plane_58(D)->data
	lea	rax, 8[rax+r12*8]	# tmp223,
# src/stencil_template_parallel.c:227:   buffers[SEND][SOUTH] = &plane->data[IDX(1,ny)]; // north starts from position (1,ny)
	mov	QWORD PTR 8[rdi], rax	# (*buffers_62(D))[1], tmp223
# src/stencil_template_parallel.c:230:   buffers[RECV][NORTH] = &plane->data[IDX(1,0)]; // north starts from position (1,1) (since we also have the halos in data)
	mov	rax, QWORD PTR [rsi]	# tmp308, plane_58(D)->data
	add	rax, 8	# tmp227,
# src/stencil_template_parallel.c:230:   buffers[RECV][NORTH] = &plane->data[IDX(1,0)]; // north starts from position (1,1) (since we also have the halos in data)
	mov	QWORD PTR 32[rdi], rax	# MEM[(double * restrict[4] *)buffers_62(D) + 32B][0], tmp227
# src/stencil_template_parallel.c:231:   buffers[RECV][SOUTH] = &plane->data[IDX(1,ny+1)]; 
	lea	eax, 1[r8]	# tmp230,
# src/stencil_template_parallel.c:231:   buffers[RECV][SOUTH] = &plane->data[IDX(1,ny+1)]; 
	mov	r12, QWORD PTR [rsi]	# tmp235, plane_58(D)->data
# src/stencil_template_parallel.c:231:   buffers[RECV][SOUTH] = &plane->data[IDX(1,ny+1)]; 
	imul	rax, r9	# tmp231, fx
# src/stencil_template_parallel.c:231:   buffers[RECV][SOUTH] = &plane->data[IDX(1,ny+1)]; 
	lea	rax, 8[0+rax*8]	# _19,
# src/stencil_template_parallel.c:231:   buffers[RECV][SOUTH] = &plane->data[IDX(1,ny+1)]; 
	add	r12, rax	# tmp235, _19
# src/stencil_template_parallel.c:233:   if (periodic && N[_y_]==2) {
	test	ebx, ebx	# tmp303
# src/stencil_template_parallel.c:231:   buffers[RECV][SOUTH] = &plane->data[IDX(1,ny+1)]; 
	mov	QWORD PTR 40[rdi], r12	# MEM[(double * restrict[4] *)buffers_62(D) + 32B][1], tmp235
# src/stencil_template_parallel.c:233:   if (periodic && N[_y_]==2) {
	je	.L87	#,
# src/stencil_template_parallel.c:233:   if (periodic && N[_y_]==2) {
	cmp	DWORD PTR 4[rdx], 2	# MEM[(const uint *)N_68(D) + 4B],
	je	.L138	#,
.L87:
# src/stencil_template_parallel.c:239:   if (buffers[SEND][WEST]) {
	mov	rax, QWORD PTR 24[rdi]	# _26, (*buffers_62(D))[3]
# src/stencil_template_parallel.c:239:   if (buffers[SEND][WEST]) {
	test	rax, rax	# _26
	je	.L88	#,
# src/stencil_template_parallel.c:240:     for (uint j = 0; j < ny; j++)
	test	r10d, r10d	# ny
	je	.L115	#,
# src/stencil_template_parallel.c:241:       buffers[SEND][WEST][j] = plane->data[IDX(1, j+1)];
	mov	rsi, QWORD PTR [rsi]	# _136, plane_58(D)->data
	lea	r12d, -1[r10]	# tmp298,
	cmp	r12d, 14	# tmp298,
	lea	rdx, [rsi+r11]	# ivtmp.282,
	jbe	.L90	#,
	lea	rbx, -1[r8]	# _100,
	mov	r13, r9	# tmp246, fx
	lea	rdx, [rsi+r11]	# ivtmp.282,
	imul	r13, rbx	# tmp246, _100
	add	rbp, r13	# tmp247, tmp246
	lea	rbp, [rsi+rbp*8]	# tmp249,
	cmp	rax, rbp	# _26, tmp249
	ja	.L106	#,
	lea	rbp, [rax+rbx*8]	# tmp253,
	cmp	rbp, rdx	# tmp253, ivtmp.282
	jb	.L106	#,
.L90:
	lea	rbp, -8[r11]	# _30,
	lea	rbx, [rax+r8*8]	# _43,
	.p2align 4,,10
	.p2align 3
.L95:
# src/stencil_template_parallel.c:241:       buffers[SEND][WEST][j] = plane->data[IDX(1, j+1)];
	vmovsd	xmm0, QWORD PTR [rdx]	# _190, MEM[(double *)_37]
# src/stencil_template_parallel.c:240:     for (uint j = 0; j < ny; j++)
	add	rax, 8	# ivtmp.284,
	add	rdx, rbp	# ivtmp.282, _30
# src/stencil_template_parallel.c:241:       buffers[SEND][WEST][j] = plane->data[IDX(1, j+1)];
	vmovsd	QWORD PTR -8[rax], xmm0	# MEM[(double *)_78], _190
# src/stencil_template_parallel.c:240:     for (uint j = 0; j < ny; j++)
	cmp	rax, rbx	# ivtmp.284, _43
	jne	.L95	#,
# src/stencil_template_parallel.c:245:   if (buffers[SEND][EAST]) {
	mov	rax, QWORD PTR 16[rdi]	# _194, (*buffers_62(D))[2]
# src/stencil_template_parallel.c:245:   if (buffers[SEND][EAST]) {
	test	rax, rax	# _194
	jne	.L99	#,
.L115:
# src/stencil_template_parallel.c:252: }
	pop	rbx	#
	.cfi_remember_state
	.cfi_def_cfa_offset 32
	xor	eax, eax	#
	pop	rbp	#
	.cfi_def_cfa_offset 24
	pop	r12	#
	.cfi_def_cfa_offset 16
	pop	r13	#
	.cfi_def_cfa_offset 8
	ret	
	.p2align 4,,10
	.p2align 3
.L106:
	.cfi_restore_state
	mov	r12d, r10d	# bnd.255, ny
	mov	r13, r9	# _216, fx
	mov	rbp, rax	# ivtmp.290, _26
	shr	r12d	#
	sal	r13, 4	# _216,
	sal	r12, 4	# tmp258,
	add	r12, rax	# _154, _26
	.p2align 4,,10
	.p2align 3
.L92:
# src/stencil_template_parallel.c:241:       buffers[SEND][WEST][j] = plane->data[IDX(1, j+1)];
	vmovsd	xmm0, QWORD PTR [rdx]	# MEM[(double *)_49], MEM[(double *)_49]
	add	rbp, 16	# ivtmp.290,
	vmovhpd	xmm0, xmm0, QWORD PTR 16[rdx+rcx*8]	# tmp259, MEM[(double *)_49], MEM[(double *)_49 + 16B + _1 * 8]
	add	rdx, r13	# ivtmp.287, _216
# src/stencil_template_parallel.c:241:       buffers[SEND][WEST][j] = plane->data[IDX(1, j+1)];
	vmovupd	XMMWORD PTR -16[rbp], xmm0	# MEM <vector(2) double> [(double *)_107], tmp259
	cmp	rbp, r12	# ivtmp.290, _154
	jne	.L92	#,
	mov	edx, r10d	# niters_vector_mult_vf.256, ny
	and	edx, -2	#,
	test	r10b, 1	# ny,
	je	.L93	#,
# src/stencil_template_parallel.c:241:       buffers[SEND][WEST][j] = plane->data[IDX(1, j+1)];
	lea	ebp, 1[rdx]	# tmp264,
	imul	rbp, r9	# tmp265, fx
# src/stencil_template_parallel.c:241:       buffers[SEND][WEST][j] = plane->data[IDX(1, j+1)];
	vmovsd	xmm0, QWORD PTR 8[rsi+rbp*8]	# _205, *_201
# src/stencil_template_parallel.c:241:       buffers[SEND][WEST][j] = plane->data[IDX(1, j+1)];
	vmovsd	QWORD PTR [rax+rdx*8], xmm0	# *_204, _205
# src/stencil_template_parallel.c:245:   if (buffers[SEND][EAST]) {
	mov	rax, QWORD PTR 16[rdi]	# _194, (*buffers_62(D))[2]
# src/stencil_template_parallel.c:245:   if (buffers[SEND][EAST]) {
	test	rax, rax	# _194
	je	.L115	#,
.L98:
	mov	r12, r9	# tmp270, fx
	lea	rbp, 1[rcx]	# _73,
	imul	r12, rbx	# tmp270, _100
	mov	rdi, rbp	# _86, _73
	sal	rdi, 4	# _86,
	lea	rdx, [rsi+rdi]	# _85,
	lea	rbp, [r12+rbp*2]	# tmp272,
	lea	rbp, [rsi+rbp*8]	# tmp274,
	cmp	rax, rbp	# _194, tmp274
	ja	.L107	#,
	lea	rbx, [rax+rbx*8]	# tmp278,
	cmp	rdx, rbx	# _85, tmp278
	ja	.L107	#,
.L100:
	add	rcx, r9	# tmp293, fx
	sub	r11, 8	# _174,
	lea	rdx, [rsi+rcx*8]	# ivtmp.266,
	lea	rcx, [rax+r8*8]	# _103,
	.p2align 4,,10
	.p2align 3
.L104:
# src/stencil_template_parallel.c:247:       buffers[SEND][EAST][j] = plane->data[IDX(nx, j+1)];
	vmovsd	xmm0, QWORD PTR [rdx]	# _113, MEM[(double *)_213]
# src/stencil_template_parallel.c:246:     for (uint j = 0; j < ny; j++)
	add	rax, 8	# ivtmp.267,
	add	rdx, r11	# ivtmp.266, _174
# src/stencil_template_parallel.c:247:       buffers[SEND][EAST][j] = plane->data[IDX(nx, j+1)];
	vmovsd	QWORD PTR -8[rax], xmm0	# MEM[(double *)_214], _113
# src/stencil_template_parallel.c:246:     for (uint j = 0; j < ny; j++)
	cmp	rax, rcx	# ivtmp.267, _103
	jne	.L104	#,
# src/stencil_template_parallel.c:252: }
	pop	rbx	#
	.cfi_remember_state
	.cfi_def_cfa_offset 32
	xor	eax, eax	#
	pop	rbp	#
	.cfi_def_cfa_offset 24
	pop	r12	#
	.cfi_def_cfa_offset 16
	pop	r13	#
	.cfi_def_cfa_offset 8
	ret	
	.p2align 4,,10
	.p2align 3
.L88:
	.cfi_restore_state
# src/stencil_template_parallel.c:245:   if (buffers[SEND][EAST]) {
	mov	rax, QWORD PTR 16[rdi]	# _194, (*buffers_62(D))[2]
# src/stencil_template_parallel.c:245:   if (buffers[SEND][EAST]) {
	test	rax, rax	# _194
	je	.L115	#,
# src/stencil_template_parallel.c:246:     for (uint j = 0; j < ny; j++)
	test	r10d, r10d	# ny
	je	.L115	#,
	mov	rsi, QWORD PTR [rsi]	# _136, plane_58(D)->data
	lea	r12d, -1[r10]	# tmp298,
.L99:
	cmp	r12d, 14	# tmp298,
	jbe	.L100	#,
	lea	rbx, -1[r8]	# _100,
	jmp	.L98	#
	.p2align 4,,10
	.p2align 3
.L107:
	mov	r8d, r10d	# bnd.249, ny
	lea	r11, 16[rdi]	# _140,
	mov	rdi, rax	# ivtmp.273, _194
	shr	r8d	#
	sal	r8, 4	# tmp283,
	add	r8, rax	# _117, _194
	.p2align 4,,10
	.p2align 3
.L102:
# src/stencil_template_parallel.c:247:       buffers[SEND][EAST][j] = plane->data[IDX(nx, j+1)];
	vmovsd	xmm0, QWORD PTR [rdx]	# MEM[(double *)_116], MEM[(double *)_116]
	add	rdi, 16	# ivtmp.273,
	vmovhpd	xmm0, xmm0, QWORD PTR 16[rdx+rcx*8]	# tmp284, MEM[(double *)_116], MEM[(double *)_116 + 16B + _1 * 8]
	add	rdx, r11	# ivtmp.270, _140
# src/stencil_template_parallel.c:247:       buffers[SEND][EAST][j] = plane->data[IDX(nx, j+1)];
	vmovupd	XMMWORD PTR -16[rdi], xmm0	# MEM <vector(2) double> [(double *)_137], tmp284
	cmp	rdi, r8	# ivtmp.273, _117
	jne	.L102	#,
	mov	edi, r10d	# niters_vector_mult_vf.250, ny
	and	edi, -2	#,
	and	r10d, 1	# ny,
	je	.L115	#,
# src/stencil_template_parallel.c:247:       buffers[SEND][EAST][j] = plane->data[IDX(nx, j+1)];
	lea	edx, 1[rdi]	# tmp289,
	imul	rdx, r9	# tmp290, fx
# src/stencil_template_parallel.c:247:       buffers[SEND][EAST][j] = plane->data[IDX(nx, j+1)];
	add	rdx, rcx	# tmp291, _1
	vmovsd	xmm0, QWORD PTR [rsi+rdx*8]	# _128, *_124
# src/stencil_template_parallel.c:247:       buffers[SEND][EAST][j] = plane->data[IDX(nx, j+1)];
	vmovsd	QWORD PTR [rax+rdi*8], xmm0	# *_127, _128
# src/stencil_template_parallel.c:252: }
	xor	eax, eax	#
	pop	rbx	#
	.cfi_remember_state
	.cfi_def_cfa_offset 32
	pop	rbp	#
	.cfi_def_cfa_offset 24
	pop	r12	#
	.cfi_def_cfa_offset 16
	pop	r13	#
	.cfi_def_cfa_offset 8
	ret	
	.p2align 4,,10
	.p2align 3
.L138:
	.cfi_restore_state
# src/stencil_template_parallel.c:234:     buffers[RECV][NORTH] = &plane->data[IDX(1,ny+1)];  // north starts from position (1,1) (since we also have the halos in data)
	add	rax, QWORD PTR [rsi]	# tmp239, plane_58(D)->data
# src/stencil_template_parallel.c:234:     buffers[RECV][NORTH] = &plane->data[IDX(1,ny+1)];  // north starts from position (1,1) (since we also have the halos in data)
	mov	QWORD PTR 32[rdi], rax	# MEM[(double * restrict[4] *)buffers_62(D) + 32B][0], tmp239
# src/stencil_template_parallel.c:235:     buffers[RECV][SOUTH] = &plane->data[IDX(1,0)]; 
	mov	rax, QWORD PTR [rsi]	# tmp310, plane_58(D)->data
	add	rax, 8	# tmp243,
# src/stencil_template_parallel.c:235:     buffers[RECV][SOUTH] = &plane->data[IDX(1,0)]; 
	mov	QWORD PTR 40[rdi], rax	# MEM[(double * restrict[4] *)buffers_62(D) + 32B][1], tmp243
	jmp	.L87	#
	.p2align 4,,10
	.p2align 3
.L93:
# src/stencil_template_parallel.c:245:   if (buffers[SEND][EAST]) {
	mov	rax, QWORD PTR 16[rdi]	# _194, (*buffers_62(D))[2]
# src/stencil_template_parallel.c:245:   if (buffers[SEND][EAST]) {
	test	rax, rax	# _194
	jne	.L98	#,
	jmp	.L115	#
	.cfi_endproc
.LFE52:
	.size	fill_buffers, .-fill_buffers
	.p2align 4
	.globl	post_MPI_reqs
	.type	post_MPI_reqs, @function
post_MPI_reqs:
.LFB53:
	.cfi_startproc
	endbr64	
	push	r15	#
	.cfi_def_cfa_offset 16
	.cfi_offset 15, -16
	push	r14	#
	.cfi_def_cfa_offset 24
	.cfi_offset 14, -24
	push	r13	#
	.cfi_def_cfa_offset 32
	.cfi_offset 13, -32
	mov	r13, rsi	# buffers, tmp146
	push	r12	#
	.cfi_def_cfa_offset 40
	.cfi_offset 12, -40
	mov	r12, r8	# comm, tmp149
	push	rbp	#
	.cfi_def_cfa_offset 48
	.cfi_offset 6, -48
	mov	rbp, rdi	# reqs, tmp145
	push	rbx	#
	.cfi_def_cfa_offset 56
	.cfi_offset 3, -56
	mov	rbx, rcx	# my_neigh, tmp148
	sub	rsp, 8	#,
	.cfi_def_cfa_offset 64
# src/stencil_template_parallel.c:267:     if (my_neigh[d] != MPI_PROC_NULL) {
	mov	ecx, DWORD PTR [rcx]	# _83, *my_neigh_39(D)
# src/stencil_template_parallel.c:261:   const unsigned int nx = plane->size[_x_];
	mov	r14d, DWORD PTR 8[rdx]	# nx, plane_36(D)->size[0]
# src/stencil_template_parallel.c:262:   const unsigned int ny = plane->size[_y_];
	mov	r15d, DWORD PTR 12[rdx]	# ny, plane_36(D)->size[1]
# src/stencil_template_parallel.c:267:     if (my_neigh[d] != MPI_PROC_NULL) {
	cmp	ecx, -2	# _83,
	je	.L140	#,
# src/stencil_template_parallel.c:268:       MPI_Irecv(buffers[RECV][d], nx, MPI_DOUBLE, my_neigh[d], 0,  comm, &reqs[d+4]);
	sub	rsp, 8	#,
	.cfi_def_cfa_offset 72
	lea	rax, 32[rbp]	# tmp118,
	mov	rdi, QWORD PTR 32[rsi]	# MEM[(double * restrict[4] *)buffers_41(D) + 32B][0], MEM[(double * restrict[4] *)buffers_41(D) + 32B][0]
	mov	r9, r8	#, comm
	push	rax	# tmp118
	.cfi_def_cfa_offset 80
	xor	r8d, r8d	#
	lea	rdx, ompi_mpi_double[rip]	# tmp119,
	mov	esi, r14d	#, nx
	call	MPI_Irecv@PLT	#
# src/stencil_template_parallel.c:269:       MPI_Isend(buffers[SEND][d], nx, MPI_DOUBLE, my_neigh[d], 0,  comm, &reqs[d]); 
	mov	ecx, DWORD PTR [rbx]	# *my_neigh_39(D), *my_neigh_39(D)
	mov	rdi, QWORD PTR 0[r13]	# (*buffers_41(D))[0], (*buffers_41(D))[0]
	mov	r9, r12	#, comm
	xor	r8d, r8d	#
	lea	rdx, ompi_mpi_double[rip]	# tmp119,
	mov	esi, r14d	#, nx
	mov	QWORD PTR [rsp], rbp	#, reqs
	call	MPI_Isend@PLT	#
	pop	r9	#
	.cfi_def_cfa_offset 72
	pop	r10	#
	.cfi_def_cfa_offset 64
.L140:
# src/stencil_template_parallel.c:267:     if (my_neigh[d] != MPI_PROC_NULL) {
	mov	ecx, DWORD PTR 4[rbx]	# _4, MEM[(int *)my_neigh_39(D) + 4B]
# src/stencil_template_parallel.c:267:     if (my_neigh[d] != MPI_PROC_NULL) {
	cmp	ecx, -2	# _4,
	je	.L141	#,
# src/stencil_template_parallel.c:268:       MPI_Irecv(buffers[RECV][d], nx, MPI_DOUBLE, my_neigh[d], 0,  comm, &reqs[d+4]);
	sub	rsp, 8	#,
	.cfi_def_cfa_offset 72
	lea	rax, 40[rbp]	# tmp124,
	mov	rdi, QWORD PTR 40[r13]	# MEM[(double * restrict[4] *)buffers_41(D) + 32B][1], MEM[(double * restrict[4] *)buffers_41(D) + 32B][1]
	mov	r9, r12	#, comm
	push	rax	# tmp124
	.cfi_def_cfa_offset 80
	xor	r8d, r8d	#
	lea	rdx, ompi_mpi_double[rip]	# tmp125,
	mov	esi, r14d	#, nx
	call	MPI_Irecv@PLT	#
# src/stencil_template_parallel.c:269:       MPI_Isend(buffers[SEND][d], nx, MPI_DOUBLE, my_neigh[d], 0,  comm, &reqs[d]); 
	mov	rdi, QWORD PTR 8[r13]	# (*buffers_41(D))[1], (*buffers_41(D))[1]
	mov	ecx, DWORD PTR 4[rbx]	# MEM[(int *)my_neigh_39(D) + 4B], MEM[(int *)my_neigh_39(D) + 4B]
	lea	rax, 8[rbp]	# tmp128,
	xor	r8d, r8d	#
	mov	r9, r12	#, comm
	lea	rdx, ompi_mpi_double[rip]	# tmp125,
	mov	esi, r14d	#, nx
	mov	QWORD PTR [rsp], rax	#, tmp128
	call	MPI_Isend@PLT	#
	pop	rdi	#
	.cfi_def_cfa_offset 72
	pop	r8	#
	.cfi_def_cfa_offset 64
.L141:
# src/stencil_template_parallel.c:275:     if (my_neigh[d] != MPI_PROC_NULL) {
	mov	ecx, DWORD PTR 8[rbx]	# _32, MEM[(int *)my_neigh_39(D) + 8B]
# src/stencil_template_parallel.c:275:     if (my_neigh[d] != MPI_PROC_NULL) {
	cmp	ecx, -2	# _32,
	je	.L142	#,
# src/stencil_template_parallel.c:276:       MPI_Irecv(buffers[RECV][d], ny, MPI_DOUBLE, my_neigh[d], 0,  comm, &reqs[d+4]); 
	sub	rsp, 8	#,
	.cfi_def_cfa_offset 72
	lea	rax, 48[rbp]	# tmp131,
	mov	rdi, QWORD PTR 48[r13]	# MEM[(double * restrict[4] *)buffers_41(D) + 32B][2], MEM[(double * restrict[4] *)buffers_41(D) + 32B][2]
	mov	r9, r12	#, comm
	push	rax	# tmp131
	.cfi_def_cfa_offset 80
	lea	r14, ompi_mpi_double[rip]	# tmp132,
	xor	r8d, r8d	#
	mov	esi, r15d	#, ny
	mov	rdx, r14	#, tmp132
	call	MPI_Irecv@PLT	#
# src/stencil_template_parallel.c:277:       MPI_Isend(buffers[SEND][d], ny, MPI_DOUBLE, my_neigh[d], 0,  comm, &reqs[d]);  
	mov	ecx, DWORD PTR 8[rbx]	# MEM[(int *)my_neigh_39(D) + 8B], MEM[(int *)my_neigh_39(D) + 8B]
	mov	rdi, QWORD PTR 16[r13]	# (*buffers_41(D))[2], (*buffers_41(D))[2]
	lea	rax, 16[rbp]	# tmp135,
	mov	esi, r15d	#, ny
	mov	r9, r12	#, comm
	xor	r8d, r8d	#
	mov	rdx, r14	#, tmp132
	mov	QWORD PTR [rsp], rax	#, tmp135
	call	MPI_Isend@PLT	#
	pop	rcx	#
	.cfi_def_cfa_offset 72
	pop	rsi	#
	.cfi_def_cfa_offset 64
.L142:
# src/stencil_template_parallel.c:275:     if (my_neigh[d] != MPI_PROC_NULL) {
	mov	ecx, DWORD PTR 12[rbx]	# _18, MEM[(int *)my_neigh_39(D) + 12B]
# src/stencil_template_parallel.c:275:     if (my_neigh[d] != MPI_PROC_NULL) {
	cmp	ecx, -2	# _18,
	jne	.L157	#,
.L143:
# src/stencil_template_parallel.c:282: }
	add	rsp, 8	#,
	.cfi_remember_state
	.cfi_def_cfa_offset 56
	xor	eax, eax	#
	pop	rbx	#
	.cfi_def_cfa_offset 48
	pop	rbp	#
	.cfi_def_cfa_offset 40
	pop	r12	#
	.cfi_def_cfa_offset 32
	pop	r13	#
	.cfi_def_cfa_offset 24
	pop	r14	#
	.cfi_def_cfa_offset 16
	pop	r15	#
	.cfi_def_cfa_offset 8
	ret	
	.p2align 4,,10
	.p2align 3
.L157:
	.cfi_restore_state
# src/stencil_template_parallel.c:276:       MPI_Irecv(buffers[RECV][d], ny, MPI_DOUBLE, my_neigh[d], 0,  comm, &reqs[d+4]); 
	sub	rsp, 8	#,
	.cfi_def_cfa_offset 72
	lea	rax, 56[rbp]	# tmp138,
	mov	rdi, QWORD PTR 56[r13]	# MEM[(double * restrict[4] *)buffers_41(D) + 32B][3], MEM[(double * restrict[4] *)buffers_41(D) + 32B][3]
	mov	r9, r12	#, comm
	push	rax	# tmp138
	.cfi_def_cfa_offset 80
	lea	r14, ompi_mpi_double[rip]	# tmp139,
	xor	r8d, r8d	#
	mov	esi, r15d	#, ny
	mov	rdx, r14	#, tmp139
# src/stencil_template_parallel.c:277:       MPI_Isend(buffers[SEND][d], ny, MPI_DOUBLE, my_neigh[d], 0,  comm, &reqs[d]);  
	add	rbp, 24	# tmp142,
# src/stencil_template_parallel.c:276:       MPI_Irecv(buffers[RECV][d], ny, MPI_DOUBLE, my_neigh[d], 0,  comm, &reqs[d+4]); 
	call	MPI_Irecv@PLT	#
# src/stencil_template_parallel.c:277:       MPI_Isend(buffers[SEND][d], ny, MPI_DOUBLE, my_neigh[d], 0,  comm, &reqs[d]);  
	mov	ecx, DWORD PTR 12[rbx]	# MEM[(int *)my_neigh_39(D) + 12B], MEM[(int *)my_neigh_39(D) + 12B]
	mov	rdi, QWORD PTR 24[r13]	# (*buffers_41(D))[3], (*buffers_41(D))[3]
	mov	rdx, r14	#, tmp139
	mov	r9, r12	#, comm
	xor	r8d, r8d	#
	mov	esi, r15d	#, ny
	mov	QWORD PTR [rsp], rbp	#, tmp142
	call	MPI_Isend@PLT	#
	pop	rax	#
	.cfi_def_cfa_offset 72
	pop	rdx	#
	.cfi_def_cfa_offset 64
	jmp	.L143	#
	.cfi_endproc
.LFE53:
	.size	post_MPI_reqs, .-post_MPI_reqs
	.p2align 4
	.globl	copy_halos
	.type	copy_halos, @function
copy_halos:
.LFB54:
	.cfi_startproc
	endbr64	
	push	r15	#
	.cfi_def_cfa_offset 16
	.cfi_offset 15, -16
	mov	r11, rdi	# buffers, tmp179
	mov	rdi, rsi	# plane, tmp180
	push	r14	#
	.cfi_def_cfa_offset 24
	.cfi_offset 14, -24
	push	r13	#
	.cfi_def_cfa_offset 32
	.cfi_offset 13, -32
	push	r12	#
	.cfi_def_cfa_offset 40
	.cfi_offset 12, -40
	push	rbp	#
	.cfi_def_cfa_offset 48
	.cfi_offset 6, -48
	push	rbx	#
	.cfi_def_cfa_offset 56
	.cfi_offset 3, -56
# src/stencil_template_parallel.c:291:   const unsigned int nx = plane->size[_x_];
	mov	r10d, DWORD PTR 8[rsi]	#, plane_68(D)->size[0]
# src/stencil_template_parallel.c:289:               const vec2_t N) {
	mov	ebx, ecx	# periodic, tmp182
# src/stencil_template_parallel.c:297:   if (neigh[WEST] != MPI_PROC_NULL) {
	mov	r15d, DWORD PTR 12[rdx]	# _2, MEM[(int *)neigh_72(D) + 12B]
# src/stencil_template_parallel.c:292:   const unsigned int ny = plane->size[_y_];
	mov	ebp, DWORD PTR 12[rsi]	# ny, plane_68(D)->size[1]
# src/stencil_template_parallel.c:291:   const unsigned int nx = plane->size[_x_];
	mov	r12, r10	#,
# src/stencil_template_parallel.c:293:   const size_t fx = (size_t)nx + 2;
	add	r10, 2	# fx,
# src/stencil_template_parallel.c:303:   if (neigh[EAST] != MPI_PROC_NULL) {
	mov	r14d, DWORD PTR 8[rdx]	# _61, MEM[(int *)neigh_72(D) + 8B]
# src/stencil_template_parallel.c:297:   if (neigh[WEST] != MPI_PROC_NULL) {
	cmp	r15d, -2	# _2,
# src/stencil_template_parallel.c:289:               const vec2_t N) {
	mov	QWORD PTR -8[rsp], r8	# %sfp, tmp183
# src/stencil_template_parallel.c:297:   if (neigh[WEST] != MPI_PROC_NULL) {
	je	.L159	#,
# src/stencil_template_parallel.c:298:   for (uint j=0; j<ny; j++){
	test	ebp, ebp	# ny
	je	.L160	#,
	mov	r9, QWORD PTR [rsi]	# plane_68(D)->data, plane_68(D)->data
	mov	rax, QWORD PTR 56[r11]	# ivtmp.327, MEM[(double * restrict[4] *)buffers_73(D) + 32B][3]
	lea	rcx, 0[0+r10*8]	# _146,
	mov	r13d, ebp	# ny, ny
	lea	rdx, [r9+rcx]	# ivtmp.328,
	lea	rsi, [rax+r13*8]	# _156,
	.p2align 4,,10
	.p2align 3
.L161:
# src/stencil_template_parallel.c:299:     plane->data[IDX(0, j+1)] = buffers[RECV][WEST][j];
	vmovsd	xmm0, QWORD PTR [rax]	# _13, MEM[(double *)_23]
# src/stencil_template_parallel.c:298:   for (uint j=0; j<ny; j++){
	add	rax, 8	# ivtmp.327,
# src/stencil_template_parallel.c:299:     plane->data[IDX(0, j+1)] = buffers[RECV][WEST][j];
	vmovsd	QWORD PTR [rdx], xmm0	# MEM[(double *)_24], _13
# src/stencil_template_parallel.c:298:   for (uint j=0; j<ny; j++){
	add	rdx, rcx	# ivtmp.328, _146
	cmp	rax, rsi	# ivtmp.327, _156
	jne	.L161	#,
# src/stencil_template_parallel.c:303:   if (neigh[EAST] != MPI_PROC_NULL) {
	cmp	r14d, -2	# _61,
	je	.L163	#,
.L166:
	mov	r8, QWORD PTR 48[r11]	# ivtmp.312, MEM[(double * restrict[4] *)buffers_73(D) + 32B][2]
# src/stencil_template_parallel.c:305:     plane->data[IDX(nx+1, j+1)] = buffers[RECV][EAST][j];
	lea	edx, 1[r12]	# tmp159,
	add	rdx, r10	# tmp160, fx
	mov	rax, r8	# ivtmp.322, ivtmp.312
	lea	rdx, [r9+rdx*8]	# ivtmp.323,
	lea	rsi, [r8+r13*8]	# _86,
	.p2align 4,,10
	.p2align 3
.L169:
# src/stencil_template_parallel.c:305:     plane->data[IDX(nx+1, j+1)] = buffers[RECV][EAST][j];
	vmovsd	xmm0, QWORD PTR [rax]	# _28, MEM[(double *)_102]
# src/stencil_template_parallel.c:304:     for (uint j=0; j<ny; j++){
	add	rax, 8	# ivtmp.322,
# src/stencil_template_parallel.c:305:     plane->data[IDX(nx+1, j+1)] = buffers[RECV][EAST][j];
	vmovsd	QWORD PTR [rdx], xmm0	# MEM[(double *)_92], _28
# src/stencil_template_parallel.c:304:     for (uint j=0; j<ny; j++){
	add	rdx, rcx	# ivtmp.323, _146
	cmp	rsi, rax	# _86, ivtmp.322
	jne	.L169	#,
# src/stencil_template_parallel.c:310:   if (periodic && N[_x_] == 2) {
	test	ebx, ebx	# periodic
	je	.L179	#,
# src/stencil_template_parallel.c:310:   if (periodic && N[_x_] == 2) {
	mov	rax, QWORD PTR -8[rsp]	# N, %sfp
	cmp	DWORD PTR [rax], 2	# *N_77(D),
	je	.L196	#,
.L179:
# src/stencil_template_parallel.c:326: }
	pop	rbx	#
	.cfi_remember_state
	.cfi_def_cfa_offset 48
	xor	eax, eax	#
	pop	rbp	#
	.cfi_def_cfa_offset 40
	pop	r12	#
	.cfi_def_cfa_offset 32
	pop	r13	#
	.cfi_def_cfa_offset 24
	pop	r14	#
	.cfi_def_cfa_offset 16
	pop	r15	#
	.cfi_def_cfa_offset 8
	ret	
	.p2align 4,,10
	.p2align 3
.L160:
	.cfi_restore_state
# src/stencil_template_parallel.c:303:   if (neigh[EAST] != MPI_PROC_NULL) {
	cmp	r14d, -2	# _61,
	jne	.L179	#,
	.p2align 4,,10
	.p2align 3
.L163:
# src/stencil_template_parallel.c:310:   if (periodic && N[_x_] == 2) {
	test	ebx, ebx	# periodic
	je	.L179	#,
# src/stencil_template_parallel.c:311:     if (neigh[WEST] != MPI_PROC_NULL) {
	mov	rax, QWORD PTR -8[rsp]	# N, %sfp
	cmp	DWORD PTR [rax], 2	# *N_77(D),
	jne	.L179	#,
	mov	r14d, -2	# _61,
.L165:
# src/stencil_template_parallel.c:312:       for (uint j=0; j<ny; j++){
	test	ebp, ebp	# ny
	je	.L179	#,
	mov	rax, QWORD PTR 56[r11]	# ivtmp.317, MEM[(double * restrict[4] *)buffers_73(D) + 32B][3]
	mov	r9, QWORD PTR [rdi]	# plane_68(D)->data, plane_68(D)->data
# src/stencil_template_parallel.c:313:       plane->data[IDX(nx+1, j+1)] = buffers[RECV][WEST][j];
	lea	edx, 1[r12]	# tmp166,
	mov	r13d, ebp	# ny, ny
	add	rdx, r10	# tmp167, fx
	lea	rcx, 0[0+r10*8]	# _146,
	lea	rdx, [r9+rdx*8]	# ivtmp.318,
	lea	rsi, [rax+r13*8]	# _115,
	.p2align 4,,10
	.p2align 3
.L172:
# src/stencil_template_parallel.c:313:       plane->data[IDX(nx+1, j+1)] = buffers[RECV][WEST][j];
	vmovsd	xmm0, QWORD PTR [rax]	# _43, MEM[(double *)_120]
# src/stencil_template_parallel.c:312:       for (uint j=0; j<ny; j++){
	add	rax, 8	# ivtmp.317,
# src/stencil_template_parallel.c:313:       plane->data[IDX(nx+1, j+1)] = buffers[RECV][WEST][j];
	vmovsd	QWORD PTR [rdx], xmm0	# MEM[(double *)_119], _43
# src/stencil_template_parallel.c:312:       for (uint j=0; j<ny; j++){
	add	rdx, rcx	# ivtmp.318, _146
	cmp	rsi, rax	# _115, ivtmp.317
	jne	.L172	#,
# src/stencil_template_parallel.c:317:     if (neigh[EAST] != MPI_PROC_NULL) {
	cmp	r14d, -2	# _61,
	je	.L179	#,
	mov	r8, QWORD PTR 48[r11]	# ivtmp.312, MEM[(double * restrict[4] *)buffers_73(D) + 32B][2]
.L171:
	add	r9, rcx	# ivtmp.313, _146
	lea	rax, [r8+r13*8]	# _135,
	.p2align 4,,10
	.p2align 3
.L173:
# src/stencil_template_parallel.c:319:       plane->data[IDX(0, j+1)] = buffers[RECV][EAST][j];
	vmovsd	xmm0, QWORD PTR [r8]	# _54, MEM[(double *)_142]
# src/stencil_template_parallel.c:318:       for (uint j=0; j<ny; j++){
	add	r8, 8	# ivtmp.312,
# src/stencil_template_parallel.c:319:       plane->data[IDX(0, j+1)] = buffers[RECV][EAST][j];
	vmovsd	QWORD PTR [r9], xmm0	# MEM[(double *)_141], _54
# src/stencil_template_parallel.c:318:       for (uint j=0; j<ny; j++){
	add	r9, rcx	# ivtmp.313, _146
	cmp	r8, rax	# ivtmp.312, _135
	jne	.L173	#,
	jmp	.L179	#
	.p2align 4,,10
	.p2align 3
.L159:
# src/stencil_template_parallel.c:303:   if (neigh[EAST] != MPI_PROC_NULL) {
	cmp	r14d, -2	# _61,
	je	.L179	#,
# src/stencil_template_parallel.c:304:     for (uint j=0; j<ny; j++){
	test	ebp, ebp	# ny
	je	.L179	#,
	mov	r9, QWORD PTR [rsi]	# plane_68(D)->data, plane_68(D)->data
	lea	rcx, 0[0+r10*8]	# _146,
	mov	r13d, ebp	# ny, ny
	jmp	.L166	#
	.p2align 4,,10
	.p2align 3
.L196:
# src/stencil_template_parallel.c:311:     if (neigh[WEST] != MPI_PROC_NULL) {
	cmp	r15d, -2	# _2,
	jne	.L165	#,
	jmp	.L171	#
	.cfi_endproc
.LFE54:
	.size	copy_halos, .-copy_halos
	.p2align 4
	.globl	simple_factorization
	.type	simple_factorization, @function
simple_factorization:
.LFB56:
	.cfi_startproc
	endbr64	
# src/stencil_template_parallel.c:680:   while ( f < A ) { // we try every f from 2 to A-1
	cmp	edi, 2	# A,
# src/stencil_template_parallel.c:668: uint simple_factorization( uint A, int *Nfactors, uint **factors ) {
	push	r12	#
	.cfi_def_cfa_offset 16
	.cfi_offset 12, -16
	mov	r12, rdx	# factors, tmp136
	push	rbp	#
	.cfi_def_cfa_offset 24
	.cfi_offset 6, -24
	push	rbx	#
	.cfi_def_cfa_offset 32
	.cfi_offset 3, -32
# src/stencil_template_parallel.c:680:   while ( f < A ) { // we try every f from 2 to A-1
	jbe	.L198	#,
	mov	ebx, edi	# A, tmp134
# src/stencil_template_parallel.c:677:   uint _A_ = A; //working copy
	mov	ecx, edi	# _A_, A
# src/stencil_template_parallel.c:676:   int f = 2; // starting from the smallest prime
	mov	ebp, 2	# f,
# src/stencil_template_parallel.c:675:   int N = 0; // counter for prime factors
	xor	edi, edi	# N
	.p2align 4,,10
	.p2align 3
.L199:
# src/stencil_template_parallel.c:681:       while( _A_ % f == 0 ) { // while f divides what is left
	xor	edx, edx	# tmp115
	mov	eax, ecx	# tmp116, _A_
	div	ebp	# f
# src/stencil_template_parallel.c:681:       while( _A_ % f == 0 ) { // while f divides what is left
	test	edx, edx	# tmp115
	jne	.L202	#,
	.p2align 4,,10
	.p2align 3
.L200:
# src/stencil_template_parallel.c:683: 	      _A_ /= f; 
	mov	eax, ecx	# _A_, _A_
	xor	edx, edx	# tmp107
# src/stencil_template_parallel.c:682: 	      N++; // count this factor
	add	edi, 1	# N,
# src/stencil_template_parallel.c:683: 	      _A_ /= f; 
	div	ebp	# f
# src/stencil_template_parallel.c:681:       while( _A_ % f == 0 ) { // while f divides what is left
	xor	edx, edx	# tmp109
# src/stencil_template_parallel.c:683: 	      _A_ /= f; 
	mov	ecx, eax	# _A_, _A_
# src/stencil_template_parallel.c:681:       while( _A_ % f == 0 ) { // while f divides what is left
	div	ebp	# f
# src/stencil_template_parallel.c:681:       while( _A_ % f == 0 ) { // while f divides what is left
	test	edx, edx	# tmp109
	je	.L200	#,
.L202:
# src/stencil_template_parallel.c:687:       f++; 
	lea	eax, 1[rbp]	# f,
# src/stencil_template_parallel.c:680:   while ( f < A ) { // we try every f from 2 to A-1
	cmp	ebx, eax	# A, f
	je	.L201	#,
	mov	ebp, eax	# f, f
	jmp	.L199	#
	.p2align 4,,10
	.p2align 3
.L201:
# src/stencil_template_parallel.c:690:   *Nfactors = N;
	mov	DWORD PTR [rsi], edi	# *Nfactors_29(D), N
# src/stencil_template_parallel.c:691:   uint *_factors_ = (uint*)malloc( N * sizeof(uint) );
	movsx	rdi, edi	# N, N
	sal	rdi, 2	# tmp119,
	call	malloc@PLT	#
# src/stencil_template_parallel.c:695:   f   = 2;
	mov	ecx, 2	# f,
# src/stencil_template_parallel.c:694:   N   = 0;
	xor	edi, edi	# N
# src/stencil_template_parallel.c:691:   uint *_factors_ = (uint*)malloc( N * sizeof(uint) );
	mov	r8, rax	# _factors_, tmp138
	.p2align 4,,10
	.p2align 3
.L204:
# src/stencil_template_parallel.c:700:       while( _A_ % f == 0 ) {
	xor	edx, edx	# tmp131
	mov	eax, ebx	# tmp132, A
	lea	esi, 1[rdi]	# tmp121,
	div	ecx	# f
	movsx	rsi, esi	# ivtmp.338, tmp121
# src/stencil_template_parallel.c:700:       while( _A_ % f == 0 ) {
	test	edx, edx	# tmp131
	jne	.L207	#,
	.p2align 4,,10
	.p2align 3
.L205:
# src/stencil_template_parallel.c:702: 	      _A_ /= f; }
	mov	eax, ebx	# A, A
	xor	edx, edx	# tmp123
# src/stencil_template_parallel.c:701: 	      _factors_[N++] = f;
	mov	DWORD PTR -4[r8+rsi*4], ecx	# MEM[(uint *)_factors__32 + -4B + ivtmp.338_49 * 4], f
	mov	rdi, rsi	# ivtmp.338, ivtmp.338
# src/stencil_template_parallel.c:702: 	      _A_ /= f; }
	div	ecx	# f
# src/stencil_template_parallel.c:700:       while( _A_ % f == 0 ) {
	xor	edx, edx	# tmp125
# src/stencil_template_parallel.c:700:       while( _A_ % f == 0 ) {
	add	rsi, 1	# ivtmp.338,
# src/stencil_template_parallel.c:702: 	      _A_ /= f; }
	mov	ebx, eax	# A, A
# src/stencil_template_parallel.c:700:       while( _A_ % f == 0 ) {
	div	ecx	# f
# src/stencil_template_parallel.c:700:       while( _A_ % f == 0 ) {
	test	edx, edx	# tmp125
	je	.L205	#,
.L207:
# src/stencil_template_parallel.c:699:   while ( f < A ) {
	cmp	ebp, ecx	# f, f
# src/stencil_template_parallel.c:703:       f++; }
	lea	eax, 1[rcx]	# f,
# src/stencil_template_parallel.c:699:   while ( f < A ) {
	je	.L203	#,
	mov	ecx, eax	# f, f
	jmp	.L204	#
.L198:
# src/stencil_template_parallel.c:690:   *Nfactors = N;
	mov	DWORD PTR [rsi], 0	# *Nfactors_29(D),
# src/stencil_template_parallel.c:691:   uint *_factors_ = (uint*)malloc( N * sizeof(uint) );
	xor	edi, edi	#
	call	malloc@PLT	#
	mov	r8, rax	# _factors_, tmp137
.L203:
# src/stencil_template_parallel.c:705:   *factors = _factors_;
	mov	QWORD PTR [r12], r8	# *factors_33(D), _factors_
# src/stencil_template_parallel.c:707: }
	xor	eax, eax	#
	pop	rbx	#
	.cfi_def_cfa_offset 24
	pop	rbp	#
	.cfi_def_cfa_offset 16
	pop	r12	#
	.cfi_def_cfa_offset 8
	ret	
	.cfi_endproc
.LFE56:
	.size	simple_factorization, .-simple_factorization
	.p2align 4
	.globl	initialize_sources
	.type	initialize_sources, @function
initialize_sources:
.LFB57:
	.cfi_startproc
	endbr64	
	push	rbp	#
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	mov	rbp, rsp	#,
	.cfi_def_cfa_register 6
	push	r15	#
	.cfi_offset 15, -24
	mov	r15d, edi	# Me, tmp244
# src/stencil_template_parallel.c:720:   srand48(time(NULL) ^ Me); // get a different seed per rank
	xor	edi, edi	#
# src/stencil_template_parallel.c:718: 			vec2_t  **Sources ) {
	push	r14	#
	.cfi_offset 14, -32
	mov	r14, rcx	# mysize, tmp247
	push	r13	#
	push	r12	#
	push	rbx	#
	.cfi_offset 13, -40
	.cfi_offset 12, -48
	.cfi_offset 3, -56
	mov	ebx, r8d	# Nsources, tmp248
	and	rsp, -32	#,
	sub	rsp, 32	#,
# src/stencil_template_parallel.c:718: 			vec2_t  **Sources ) {
	mov	DWORD PTR 28[rsp], esi	# %sfp, tmp245
	mov	QWORD PTR 8[rsp], rdx	# %sfp, tmp246
	mov	QWORD PTR [rsp], r9	# %sfp, tmp249
# src/stencil_template_parallel.c:720:   srand48(time(NULL) ^ Me); // get a different seed per rank
	call	time@PLT	#
# src/stencil_template_parallel.c:720:   srand48(time(NULL) ^ Me); // get a different seed per rank
	movsx	rdi, r15d	# Me, Me
	xor	rdi, rax	# tmp182, tmp250
	call	srand48@PLT	#
# src/stencil_template_parallel.c:723:   int *tasks_with_sources = (int*)malloc( Nsources * sizeof(int) );
	movsx	rax, ebx	# Nsources, Nsources
	sal	rax, 2	# _5,
	mov	rdi, rax	#, _5
	mov	QWORD PTR 16[rsp], rax	# %sfp, _5
	call	malloc@PLT	#
# src/stencil_template_parallel.c:726:   if ( Me == 0 ) {
	test	r15d, r15d	# Me
# src/stencil_template_parallel.c:723:   int *tasks_with_sources = (int*)malloc( Nsources * sizeof(int) );
	mov	r13, rax	# tasks_with_sources, tmp251
# src/stencil_template_parallel.c:726:   if ( Me == 0 ) {
	jne	.L219	#,
# src/stencil_template_parallel.c:727:       for ( int i = 0; i < Nsources; i++ ) {
	test	ebx, ebx	# Nsources
	jle	.L220	#,
	mov	r12, rax	# ivtmp.404, tasks_with_sources
	mov	rax, QWORD PTR 16[rsp]	# _5, %sfp
	add	rax, r13	# _5, tasks_with_sources
	mov	QWORD PTR 16[rsp], rax	# %sfp, _5
	.p2align 4,,10
	.p2align 3
.L221:
# src/stencil_template_parallel.c:728: 	      tasks_with_sources[i] = (int)lrand48() % Ntasks; }
	call	lrand48@PLT	#
# src/stencil_template_parallel.c:727:       for ( int i = 0; i < Nsources; i++ ) {
	add	r12, 4	# ivtmp.404,
# src/stencil_template_parallel.c:728: 	      tasks_with_sources[i] = (int)lrand48() % Ntasks; }
	cdq
	idiv	DWORD PTR 28[rsp]	# %sfp
# src/stencil_template_parallel.c:728: 	      tasks_with_sources[i] = (int)lrand48() % Ntasks; }
	mov	DWORD PTR -4[r12], edx	# MEM[(int *)_159], tmp186
# src/stencil_template_parallel.c:727:       for ( int i = 0; i < Nsources; i++ ) {
	cmp	QWORD PTR 16[rsp], r12	# %sfp, ivtmp.404
	jne	.L221	#,
# src/stencil_template_parallel.c:733:   MPI_Bcast( tasks_with_sources, Nsources, MPI_INT, 0, *Comm );
	mov	rax, QWORD PTR 8[rsp]	# Comm, %sfp
	xor	ecx, ecx	#
	lea	rdx, ompi_mpi_int[rip]	# tmp243,
	mov	esi, ebx	#, Nsources
	mov	rdi, r13	#, tasks_with_sources
	mov	r8, QWORD PTR [rax]	#, *Comm_56(D)
	call	MPI_Bcast@PLT	#
.L231:
	lea	eax, -1[rbx]	# tmp192,
	cmp	eax, 6	# tmp192,
	jbe	.L232	#,
	mov	edx, ebx	# bnd.362, Nsources
	vmovd	xmm2, r15d	# vect_cst__114, Me
# src/stencil_template_parallel.c:718: 			vec2_t  **Sources ) {
	vpxor	xmm1, xmm1, xmm1	# vect_nlocal_71.371
	mov	rax, r13	# ivtmp.397, tasks_with_sources
	shr	edx, 3	#,
	vpbroadcastd	ymm2, xmm2	# vect_cst__114, vect_cst__114
	sal	rdx, 5	# tmp196,
	add	rdx, r13	# _135, tasks_with_sources
	.p2align 4,,10
	.p2align 3
.L226:
# src/stencil_template_parallel.c:738:     nlocal += (tasks_with_sources[i] == Me); 
	vpcmpeqd	ymm0, ymm2, YMMWORD PTR [rax]	# tmp199, vect_cst__114, MEM <vector(8) int> [(int *)_179]
	add	rax, 32	# ivtmp.397,
	cmp	rdx, rax	# _135, ivtmp.397
# src/stencil_template_parallel.c:738:     nlocal += (tasks_with_sources[i] == Me); 
	vpsubd	ymm1, ymm1, ymm0	# vect_nlocal_71.371, vect_nlocal_71.371, tmp199
	jne	.L226	#,
	vmovdqa	xmm0, xmm1	# tmp200, vect_nlocal_71.371
	vextracti128	xmm1, ymm1, 0x1	# tmp201, vect_nlocal_71.371
	mov	edx, ebx	# tmp.376, Nsources
	vpaddd	xmm0, xmm0, xmm1	# _123, tmp200, tmp201
	and	edx, -8	# tmp.376,
	vpsrldq	xmm1, xmm0, 8	# tmp203, _123,
	cmp	ebx, edx	# Nsources, tmp.376
	mov	ecx, edx	#, tmp.376
	vpaddd	xmm1, xmm0, xmm1	# _125, _123, tmp203
	vpsrldq	xmm2, xmm1, 4	# tmp205, _125,
	vpaddd	xmm1, xmm1, xmm2	# tmp206, _125, tmp205
	vmovd	eax, xmm1	# stmp_nlocal_77.384, tmp206
	je	.L243	#,
	vzeroupper
.L225:
	mov	esi, ebx	# niters.373, Nsources
	sub	esi, ecx	# niters.373, niters_vector_mult_vf.363
	lea	edi, -1[rsi]	# tmp207,
	cmp	edi, 2	# tmp207,
	jbe	.L228	#,
# src/stencil_template_parallel.c:738:     nlocal += (tasks_with_sources[i] == Me); 
	vmovd	xmm3, r15d	# Me, Me
	vpshufd	xmm1, xmm3, 0	# tmp210, Me
	vpcmpeqd	xmm1, xmm1, XMMWORD PTR 0[r13+rcx*4]	# tmp212, tmp210, MEM <vector(4) int> [(int *)vectp_tasks_with_sources.379_160]
	mov	ecx, esi	# niters_vector_mult_vf.375, niters.373
	and	ecx, -4	# niters_vector_mult_vf.375,
	add	edx, ecx	# tmp.376, niters_vector_mult_vf.375
	and	esi, 3	# niters.373,
# src/stencil_template_parallel.c:738:     nlocal += (tasks_with_sources[i] == Me); 
	vpsubd	xmm0, xmm0, xmm1	# vect_nlocal_77.383, _123, tmp212
	vpsrldq	xmm1, xmm0, 8	# tmp214, vect_nlocal_77.383,
	vpaddd	xmm0, xmm0, xmm1	# _174, vect_nlocal_77.383, tmp214
	vpsrldq	xmm1, xmm0, 4	# tmp216, _174,
	vpaddd	xmm0, xmm0, xmm1	# tmp217, _174, tmp216
	vmovd	eax, xmm0	# stmp_nlocal_77.384, tmp217
	je	.L227	#,
.L228:
# src/stencil_template_parallel.c:738:     nlocal += (tasks_with_sources[i] == Me); 
	movsx	rcx, edx	# tmp.376, tmp.376
# src/stencil_template_parallel.c:738:     nlocal += (tasks_with_sources[i] == Me); 
	cmp	DWORD PTR 0[r13+rcx*4], r15d	# *_15, Me
# src/stencil_template_parallel.c:738:     nlocal += (tasks_with_sources[i] == Me); 
	lea	rsi, 0[0+rcx*4]	# _14,
# src/stencil_template_parallel.c:738:     nlocal += (tasks_with_sources[i] == Me); 
	sete	cl	#, tmp221
	movzx	ecx, cl	# tmp221, tmp221
# src/stencil_template_parallel.c:738:     nlocal += (tasks_with_sources[i] == Me); 
	add	eax, ecx	# stmp_nlocal_77.384, tmp221
# src/stencil_template_parallel.c:737:   for ( int i = 0; i < Nsources; i++ ) {
	lea	ecx, 1[rdx]	# i,
# src/stencil_template_parallel.c:737:   for ( int i = 0; i < Nsources; i++ ) {
	cmp	ebx, ecx	# Nsources, i
	jle	.L227	#,
# src/stencil_template_parallel.c:738:     nlocal += (tasks_with_sources[i] == Me); 
	xor	ecx, ecx	# tmp224
	cmp	DWORD PTR 4[r13+rsi], r15d	# *_87, Me
	sete	cl	#, tmp224
# src/stencil_template_parallel.c:737:   for ( int i = 0; i < Nsources; i++ ) {
	add	edx, 2	# i,
# src/stencil_template_parallel.c:738:     nlocal += (tasks_with_sources[i] == Me); 
	add	eax, ecx	# stmp_nlocal_77.384, tmp224
# src/stencil_template_parallel.c:737:   for ( int i = 0; i < Nsources; i++ ) {
	cmp	ebx, edx	# Nsources, i
	jle	.L227	#,
# src/stencil_template_parallel.c:738:     nlocal += (tasks_with_sources[i] == Me); 
	xor	edx, edx	# tmp227
	cmp	DWORD PTR 8[r13+rsi], r15d	# *_138, Me
	sete	dl	#, tmp227
# src/stencil_template_parallel.c:738:     nlocal += (tasks_with_sources[i] == Me); 
	add	eax, edx	# stmp_nlocal_77.384, tmp227
.L227:
# src/stencil_template_parallel.c:740:   *Nsources_local = nlocal;
	mov	rdi, QWORD PTR [rsp]	# Nsources_local, %sfp
# src/stencil_template_parallel.c:743:   if ( nlocal > 0 ) {
	test	eax, eax	# stmp_nlocal_77.384
# src/stencil_template_parallel.c:740:   *Nsources_local = nlocal;
	mov	DWORD PTR [rdi], eax	# *Nsources_local_58(D), stmp_nlocal_77.384
# src/stencil_template_parallel.c:743:   if ( nlocal > 0 ) {
	jg	.L244	#,
.L223:
# src/stencil_template_parallel.c:753:   free( tasks_with_sources );
	mov	rdi, r13	#, tasks_with_sources
	call	free@PLT	#
# src/stencil_template_parallel.c:756: }
	lea	rsp, -40[rbp]	#,
	xor	eax, eax	#
	pop	rbx	#
	pop	r12	#
	pop	r13	#
	pop	r14	#
	pop	r15	#
	pop	rbp	#
	.cfi_remember_state
	.cfi_def_cfa 7, 8
	ret	
	.p2align 4,,10
	.p2align 3
.L219:
	.cfi_restore_state
# src/stencil_template_parallel.c:733:   MPI_Bcast( tasks_with_sources, Nsources, MPI_INT, 0, *Comm );
	mov	rax, QWORD PTR 8[rsp]	# Comm, %sfp
	xor	ecx, ecx	#
	lea	rdx, ompi_mpi_int[rip]	# tmp191,
	mov	esi, ebx	#, Nsources
	mov	rdi, r13	#, tasks_with_sources
	mov	r8, QWORD PTR [rax]	#, *Comm_56(D)
	call	MPI_Bcast@PLT	#
# src/stencil_template_parallel.c:737:   for ( int i = 0; i < Nsources; i++ ) {
	test	ebx, ebx	# Nsources
	jg	.L231	#,
# src/stencil_template_parallel.c:740:   *Nsources_local = nlocal;
	mov	rax, QWORD PTR [rsp]	# Nsources_local, %sfp
	mov	DWORD PTR [rax], 0	# *Nsources_local_58(D),
	jmp	.L223	#
	.p2align 4,,10
	.p2align 3
.L244:
# src/stencil_template_parallel.c:744:       vec2_t * restrict helper = (vec2_t*)malloc( nlocal * sizeof(vec2_t) );      
	movsx	r12, eax	# stmp_nlocal_77.384, stmp_nlocal_77.384
	sal	r12, 3	# _38,
	mov	rdi, r12	#, _38
	call	malloc@PLT	#
	mov	r15, rax	# helper, tmp253
	mov	rbx, rax	# ivtmp.391, helper
	add	r12, rax	# _156, helper
	.p2align 4,,10
	.p2align 3
.L230:
# src/stencil_template_parallel.c:747:         helper[s][_x_] = 1 + lrand48() % mysize[_x_];
	call	lrand48@PLT	#
# src/stencil_template_parallel.c:747:         helper[s][_x_] = 1 + lrand48() % mysize[_x_];
	mov	ecx, DWORD PTR [r14]	# *mysize_65(D), *mysize_65(D)
# src/stencil_template_parallel.c:745:       for ( int s = 0; s < nlocal; s++ ) {
	add	rbx, 8	# ivtmp.391,
# src/stencil_template_parallel.c:747:         helper[s][_x_] = 1 + lrand48() % mysize[_x_];
	cqo
	idiv	rcx	# *mysize_65(D)
# src/stencil_template_parallel.c:747:         helper[s][_x_] = 1 + lrand48() % mysize[_x_];
	add	edx, 1	# tmp235,
	mov	DWORD PTR -8[rbx], edx	# MEM[(unsigned int *)_79], tmp235
# src/stencil_template_parallel.c:748:         helper[s][_y_] = 1 + lrand48() % mysize[_y_]; }
	call	lrand48@PLT	#
# src/stencil_template_parallel.c:748:         helper[s][_y_] = 1 + lrand48() % mysize[_y_]; }
	mov	ecx, DWORD PTR 4[r14]	# MEM[(uint *)mysize_65(D) + 4B], MEM[(uint *)mysize_65(D) + 4B]
# src/stencil_template_parallel.c:748:         helper[s][_y_] = 1 + lrand48() % mysize[_y_]; }
	cqo
	idiv	rcx	# MEM[(uint *)mysize_65(D) + 4B]
# src/stencil_template_parallel.c:748:         helper[s][_y_] = 1 + lrand48() % mysize[_y_]; }
	add	edx, 1	# tmp240,
	mov	DWORD PTR -4[rbx], edx	# MEM[(unsigned int *)_79 + 4B], tmp240
# src/stencil_template_parallel.c:745:       for ( int s = 0; s < nlocal; s++ ) {
	cmp	rbx, r12	# ivtmp.391, _156
	jne	.L230	#,
# src/stencil_template_parallel.c:750:       *Sources = helper;
	mov	rax, QWORD PTR 16[rbp]	# Sources, Sources
	mov	QWORD PTR [rax], r15	# *Sources_62(D), helper
	jmp	.L223	#
.L220:
# src/stencil_template_parallel.c:733:   MPI_Bcast( tasks_with_sources, Nsources, MPI_INT, 0, *Comm );
	mov	rax, QWORD PTR 8[rsp]	# Comm, %sfp
	xor	ecx, ecx	#
	lea	rdx, ompi_mpi_int[rip]	# tmp189,
	mov	esi, ebx	#, Nsources
	mov	rdi, r13	#, tasks_with_sources
	mov	r8, QWORD PTR [rax]	#, *Comm_56(D)
	call	MPI_Bcast@PLT	#
# src/stencil_template_parallel.c:740:   *Nsources_local = nlocal;
	mov	rax, QWORD PTR [rsp]	# Nsources_local, %sfp
	mov	DWORD PTR [rax], 0	# *Nsources_local_58(D),
	jmp	.L223	#
	.p2align 4,,10
	.p2align 3
.L243:
	vzeroupper
	mov	rdi, QWORD PTR [rsp]	# Nsources_local, %sfp
# src/stencil_template_parallel.c:743:   if ( nlocal > 0 ) {
	test	eax, eax	# stmp_nlocal_77.384
# src/stencil_template_parallel.c:740:   *Nsources_local = nlocal;
	mov	DWORD PTR [rdi], eax	# *Nsources_local_58(D), stmp_nlocal_77.384
# src/stencil_template_parallel.c:743:   if ( nlocal > 0 ) {
	jle	.L223	#,
	jmp	.L244	#
.L232:
# src/stencil_template_parallel.c:718: 			vec2_t  **Sources ) {
	vpxor	xmm0, xmm0, xmm0	# _123
	xor	ecx, ecx	#
# src/stencil_template_parallel.c:737:   for ( int i = 0; i < Nsources; i++ ) {
	xor	edx, edx	# tmp.376
# src/stencil_template_parallel.c:736:   int nlocal = 0;
	xor	eax, eax	# stmp_nlocal_77.384
	jmp	.L225	#
	.cfi_endproc
.LFE57:
	.size	initialize_sources, .-initialize_sources
	.section	.rodata.str1.1,"aMS",@progbits,1
.LC3:
	.string	"invalid planes"
.LC4:
	.string	"invalid buffers"
	.section	.rodata.str1.8,"aMS",@progbits,1
	.align 8
.LC5:
	.string	"memory allocation of the planes (OLD) failes"
	.align 8
.LC6:
	.string	"memory allocation of the planes (NEW) failes"
	.align 8
.LC7:
	.string	"Memory allocation for west buffers failed"
	.align 8
.LC8:
	.string	"Memory allocation for east buffers failed"
	.text
	.p2align 4
	.globl	memory_allocate
	.type	memory_allocate, @function
memory_allocate:
.LFB58:
	.cfi_startproc
	endbr64	
# src/stencil_template_parallel.c:802:   if (planes_ptr == NULL ) {
	test	rcx, rcx	# planes_ptr
# src/stencil_template_parallel.c:763:                     plane_t *planes_ptr) {
	push	r14	#
	.cfi_def_cfa_offset 16
	.cfi_offset 14, -16
	push	r13	#
	.cfi_def_cfa_offset 24
	.cfi_offset 13, -24
	push	r12	#
	.cfi_def_cfa_offset 32
	.cfi_offset 12, -32
	push	rbp	#
	.cfi_def_cfa_offset 40
	.cfi_offset 6, -40
	push	rbx	#
	.cfi_def_cfa_offset 48
	.cfi_offset 3, -48
# src/stencil_template_parallel.c:802:   if (planes_ptr == NULL ) {
	je	.L274	#,
# src/stencil_template_parallel.c:808:   if (buffers_ptr == NULL ) {
	test	rdx, rdx	# buffers_ptr
	mov	r12, rdx	# buffers_ptr, tmp158
	je	.L275	#,
# src/stencil_template_parallel.c:818:   unsigned int frame_size = (planes_ptr[OLD].size[_x_]+2) * (planes_ptr[OLD].size[_y_]+2); //(nx+2) x (ny+2)
	mov	eax, DWORD PTR 8[rcx]	# tmp167, planes_ptr_36(D)->size[0]
# src/stencil_template_parallel.c:818:   unsigned int frame_size = (planes_ptr[OLD].size[_x_]+2) * (planes_ptr[OLD].size[_y_]+2); //(nx+2) x (ny+2)
	mov	r14d, DWORD PTR 12[rcx]	#, planes_ptr_36(D)->size[1]
	mov	r13, rdi	# neighbours, tmp157
	mov	rbp, rcx	# planes_ptr, tmp159
# src/stencil_template_parallel.c:818:   unsigned int frame_size = (planes_ptr[OLD].size[_x_]+2) * (planes_ptr[OLD].size[_y_]+2); //(nx+2) x (ny+2)
	lea	ebx, 2[rax]	# tmp127,
# src/stencil_template_parallel.c:818:   unsigned int frame_size = (planes_ptr[OLD].size[_x_]+2) * (planes_ptr[OLD].size[_y_]+2); //(nx+2) x (ny+2)
	lea	eax, 2[r14]	# tmp129,
# src/stencil_template_parallel.c:820:   planes_ptr[OLD].data = (double*)malloc( frame_size * sizeof(double) );
	imul	ebx, eax	# frame_size, tmp129
	sal	rbx, 3	# _8,
	mov	rdi, rbx	#, _8
	call	malloc@PLT	#
# src/stencil_template_parallel.c:822:   if ( planes_ptr[OLD].data == NULL ) {
	test	rax, rax	# tmp132
# src/stencil_template_parallel.c:820:   planes_ptr[OLD].data = (double*)malloc( frame_size * sizeof(double) );
	mov	QWORD PTR 0[rbp], rax	# planes_ptr_36(D)->data, tmp132
# src/stencil_template_parallel.c:820:   planes_ptr[OLD].data = (double*)malloc( frame_size * sizeof(double) );
	mov	rdi, rax	# tmp132, tmp160
# src/stencil_template_parallel.c:822:   if ( planes_ptr[OLD].data == NULL ) {
	je	.L276	#,
# /usr/include/x86_64-linux-gnu/bits/string_fortified.h:59:   return __builtin___memset_chk (__dest, __ch, __len,
	xor	esi, esi	#
	mov	rcx, rbx	#, _8
	mov	rdx, rbx	#, _8
	call	__memset_chk@PLT	#
# src/stencil_template_parallel.c:831:   planes_ptr[NEW].data = (double*)malloc( frame_size * sizeof(double) );
	mov	rdi, rbx	#, _8
	call	malloc@PLT	#
# src/stencil_template_parallel.c:833:   if ( planes_ptr[NEW].data == NULL ) {
	test	rax, rax	# tmp135
# src/stencil_template_parallel.c:831:   planes_ptr[NEW].data = (double*)malloc( frame_size * sizeof(double) );
	mov	QWORD PTR 16[rbp], rax	# MEM[(struct plane_t *)planes_ptr_36(D) + 16B].data, tmp135
# src/stencil_template_parallel.c:831:   planes_ptr[NEW].data = (double*)malloc( frame_size * sizeof(double) );
	mov	rdi, rax	# tmp135, tmp161
# src/stencil_template_parallel.c:833:   if ( planes_ptr[NEW].data == NULL ) {
	je	.L277	#,
# /usr/include/x86_64-linux-gnu/bits/string_fortified.h:59:   return __builtin___memset_chk (__dest, __ch, __len,
	xor	esi, esi	#
	mov	rcx, rbx	#, _8
	mov	rdx, rbx	#, _8
	call	__memset_chk@PLT	#
# src/stencil_template_parallel.c:875:   if (neighbours[WEST] != MPI_PROC_NULL) {
	cmp	DWORD PTR 12[r13], -2	# MEM[(const int *)neighbours_44(D) + 12B],
	jne	.L278	#,
.L252:
# src/stencil_template_parallel.c:890:   if (neighbours[EAST] != MPI_PROC_NULL) {
	cmp	DWORD PTR 8[r13], -2	# MEM[(const int *)neighbours_44(D) + 8B],
	jne	.L255	#,
.L256:
# src/stencil_template_parallel.c:907: }
	pop	rbx	#
	.cfi_remember_state
	.cfi_def_cfa_offset 40
# src/stencil_template_parallel.c:906:   return 0;
	xor	eax, eax	# <retval>
# src/stencil_template_parallel.c:907: }
	pop	rbp	#
	.cfi_def_cfa_offset 32
	pop	r12	#
	.cfi_def_cfa_offset 24
	pop	r13	#
	.cfi_def_cfa_offset 16
	pop	r14	#
	.cfi_def_cfa_offset 8
	ret	
	.p2align 4,,10
	.p2align 3
.L278:
	.cfi_restore_state
# src/stencil_template_parallel.c:877:     buffers_ptr[SEND][WEST] = (double*)malloc( sizey * sizeof(double) );
	mov	ebp, r14d	# _5, _5
	sal	rbp, 3	# _17,
	mov	rdi, rbp	#, _17
	call	malloc@PLT	#
# src/stencil_template_parallel.c:878:     buffers_ptr[RECV][WEST] = (double*)malloc( sizey * sizeof(double) );
	mov	rdi, rbp	#, _17
# src/stencil_template_parallel.c:877:     buffers_ptr[SEND][WEST] = (double*)malloc( sizey * sizeof(double) );
	mov	QWORD PTR 24[r12], rax	# (*buffers_ptr_37(D))[3], tmp139
# src/stencil_template_parallel.c:877:     buffers_ptr[SEND][WEST] = (double*)malloc( sizey * sizeof(double) );
	mov	rbx, rax	# tmp139, tmp162
# src/stencil_template_parallel.c:878:     buffers_ptr[RECV][WEST] = (double*)malloc( sizey * sizeof(double) );
	call	malloc@PLT	#
# src/stencil_template_parallel.c:880:     if (!buffers_ptr[SEND][WEST] || !buffers_ptr[RECV][WEST]) {
	test	rax, rax	# tmp140
# src/stencil_template_parallel.c:878:     buffers_ptr[RECV][WEST] = (double*)malloc( sizey * sizeof(double) );
	mov	QWORD PTR 56[r12], rax	# MEM[(double * restrict[4] *)buffers_ptr_37(D) + 32B][3], tmp140
# src/stencil_template_parallel.c:880:     if (!buffers_ptr[SEND][WEST] || !buffers_ptr[RECV][WEST]) {
	je	.L258	#,
	test	rbx, rbx	# tmp139
	jne	.L252	#,
.L258:
# /usr/include/x86_64-linux-gnu/bits/stdio2.h:79:   return __fprintf_chk (__stream, __USE_FORTIFY_LEVEL - 1, __fmt,
	mov	rcx, QWORD PTR stderr[rip]	#, stderr
	mov	edx, 41	#,
	mov	esi, 1	#,
	lea	rdi, .LC7[rip]	# tmp146,
	call	fwrite@PLT	#
# src/stencil_template_parallel.c:883:       if ( !buffers_ptr[SEND][WEST] && buffers_ptr[SEND][EAST]) { free(buffers_ptr[SEND][EAST]); }
	cmp	QWORD PTR 24[r12], 0	# (*buffers_ptr_37(D))[3],
	je	.L279	#,
.L254:
# src/stencil_template_parallel.c:884:       else if (!buffers_ptr[RECV][EAST] && buffers_ptr[RECV][WEST]) { free(buffers_ptr[RECV][WEST]); }
	cmp	QWORD PTR 48[r12], 0	# MEM[(double * restrict[4] *)buffers_ptr_37(D) + 32B][2],
	jne	.L247	#,
# src/stencil_template_parallel.c:884:       else if (!buffers_ptr[RECV][EAST] && buffers_ptr[RECV][WEST]) { free(buffers_ptr[RECV][WEST]); }
	mov	rdi, QWORD PTR 56[r12]	# _24, MEM[(double * restrict[4] *)buffers_ptr_37(D) + 32B][3]
# src/stencil_template_parallel.c:884:       else if (!buffers_ptr[RECV][EAST] && buffers_ptr[RECV][WEST]) { free(buffers_ptr[RECV][WEST]); }
	test	rdi, rdi	# _24
	je	.L247	#,
.L273:
# src/stencil_template_parallel.c:884:       else if (!buffers_ptr[RECV][EAST] && buffers_ptr[RECV][WEST]) { free(buffers_ptr[RECV][WEST]); }
	call	free@PLT	#
	jmp	.L247	#
	.p2align 4,,10
	.p2align 3
.L255:
# src/stencil_template_parallel.c:892:     buffers_ptr[SEND][EAST] = (double*)malloc( sizey * sizeof(double) );
	sal	r14, 3	# _27,
	mov	rdi, r14	#, _27
	call	malloc@PLT	#
# src/stencil_template_parallel.c:893:     buffers_ptr[RECV][EAST] = (double*)malloc( sizey * sizeof(double) );
	mov	rdi, r14	#, _27
# src/stencil_template_parallel.c:892:     buffers_ptr[SEND][EAST] = (double*)malloc( sizey * sizeof(double) );
	mov	QWORD PTR 16[r12], rax	# (*buffers_ptr_37(D))[2], tmp148
# src/stencil_template_parallel.c:892:     buffers_ptr[SEND][EAST] = (double*)malloc( sizey * sizeof(double) );
	mov	rbx, rax	# tmp148, tmp164
# src/stencil_template_parallel.c:893:     buffers_ptr[RECV][EAST] = (double*)malloc( sizey * sizeof(double) );
	call	malloc@PLT	#
# src/stencil_template_parallel.c:895:     if (!buffers_ptr[SEND][EAST] || !buffers_ptr[RECV][EAST]) {
	test	rbx, rbx	# tmp148
# src/stencil_template_parallel.c:893:     buffers_ptr[RECV][EAST] = (double*)malloc( sizey * sizeof(double) );
	mov	QWORD PTR 48[r12], rax	# MEM[(double * restrict[4] *)buffers_ptr_37(D) + 32B][2], tmp149
# src/stencil_template_parallel.c:895:     if (!buffers_ptr[SEND][EAST] || !buffers_ptr[RECV][EAST]) {
	je	.L259	#,
	test	rax, rax	# tmp149
	jne	.L256	#,
.L259:
# /usr/include/x86_64-linux-gnu/bits/stdio2.h:79:   return __fprintf_chk (__stream, __USE_FORTIFY_LEVEL - 1, __fmt,
	mov	rcx, QWORD PTR stderr[rip]	#, stderr
	mov	edx, 41	#,
	mov	esi, 1	#,
	lea	rdi, .LC8[rip]	# tmp155,
	call	fwrite@PLT	#
# src/stencil_template_parallel.c:901:       return 1;
	jmp	.L247	#
	.p2align 4,,10
	.p2align 3
.L275:
# /usr/include/x86_64-linux-gnu/bits/stdio2.h:79:   return __fprintf_chk (__stream, __USE_FORTIFY_LEVEL - 1, __fmt,
	mov	rcx, QWORD PTR stderr[rip]	#, stderr
	mov	edx, 15	#,
	mov	esi, 1	#,
	lea	rdi, .LC4[rip]	# tmp126,
	call	fwrite@PLT	#
.L247:
# src/stencil_template_parallel.c:907: }
	pop	rbx	#
	.cfi_remember_state
	.cfi_def_cfa_offset 40
# src/stencil_template_parallel.c:806:     return 1;}
	mov	eax, 1	# <retval>,
# src/stencil_template_parallel.c:907: }
	pop	rbp	#
	.cfi_def_cfa_offset 32
	pop	r12	#
	.cfi_def_cfa_offset 24
	pop	r13	#
	.cfi_def_cfa_offset 16
	pop	r14	#
	.cfi_def_cfa_offset 8
	ret	
	.p2align 4,,10
	.p2align 3
.L274:
	.cfi_restore_state
# /usr/include/x86_64-linux-gnu/bits/stdio2.h:79:   return __fprintf_chk (__stream, __USE_FORTIFY_LEVEL - 1, __fmt,
	mov	rcx, QWORD PTR stderr[rip]	#, stderr
	mov	edx, 14	#,
	mov	esi, 1	#,
	lea	rdi, .LC3[rip]	# tmp124,
	call	fwrite@PLT	#
# src/stencil_template_parallel.c:806:     return 1;}
	jmp	.L247	#
.L276:
# /usr/include/x86_64-linux-gnu/bits/stdio2.h:79:   return __fprintf_chk (__stream, __USE_FORTIFY_LEVEL - 1, __fmt,
	mov	rcx, QWORD PTR stderr[rip]	#, stderr
	mov	edx, 44	#,
	mov	esi, 1	#,
	lea	rdi, .LC5[rip]	# tmp134,
	call	fwrite@PLT	#
# src/stencil_template_parallel.c:825:     return 1;
	jmp	.L247	#
.L277:
# /usr/include/x86_64-linux-gnu/bits/stdio2.h:79:   return __fprintf_chk (__stream, __USE_FORTIFY_LEVEL - 1, __fmt,
	mov	rcx, QWORD PTR stderr[rip]	#, stderr
	mov	edx, 44	#,
	mov	esi, 1	#,
	lea	rdi, .LC6[rip]	# tmp137,
	call	fwrite@PLT	#
# src/stencil_template_parallel.c:839:     free(planes_ptr[OLD].data);
	mov	rdi, QWORD PTR 0[rbp]	#, planes_ptr_36(D)->data
	call	free@PLT	#
# src/stencil_template_parallel.c:841:     return 1;
	jmp	.L247	#
.L279:
# src/stencil_template_parallel.c:883:       if ( !buffers_ptr[SEND][WEST] && buffers_ptr[SEND][EAST]) { free(buffers_ptr[SEND][EAST]); }
	mov	rdi, QWORD PTR 16[r12]	# _22, (*buffers_ptr_37(D))[2]
# src/stencil_template_parallel.c:883:       if ( !buffers_ptr[SEND][WEST] && buffers_ptr[SEND][EAST]) { free(buffers_ptr[SEND][EAST]); }
	test	rdi, rdi	# _22
	jne	.L273	#,
	jmp	.L254	#
	.cfi_endproc
.LFE58:
	.size	memory_allocate, .-memory_allocate
	.section	.rodata.str1.8
	.align 8
.LC11:
	.string	"initialisation error: planes pointer is NULL\n"
	.align 8
.LC13:
	.ascii	"\nvalid options are ( values btw [] are the default values"
	.string	" ):\n-x    x size of the plate [10000]\n-y    y size of the plate [10000]\n-e    how many energy sources on the plate [4]\n-E    how many energy sources on the plate [1.0]\n-n    how many iterations [1000]\n-p    whether periodic boundaries applies  [0 = false]\n"
	.align 8
.LC14:
	.string	"option -%c requires an argument\n"
	.align 8
.LC15:
	.string	" -------- help unavailable ----------"
	.section	.rodata.str1.1
.LC16:
	.string	":hx:y:e:E:n:o:p:v:"
	.section	.rodata.str1.8
	.align 8
.LC17:
	.string	"invalid grid input, please select a positive non-zero integer"
	.align 8
.LC18:
	.string	"invalid number of iterations, please select a positive non-zero integer"
	.align 8
.LC19:
	.string	"invalid number of heat sources, please select a positive integer"
	.align 8
.LC20:
	.string	"invalid verbose option, please select a positive integer"
	.align 8
.LC21:
	.string	"Tasks are decomposed in a grid %d x %d\n\n"
	.align 8
.LC22:
	.string	"Task %4d: \tgrid coordinates : %3d, %3d\n\tneighbours: N %4d    E %4d    S %4d    W %4d\n"
	.align 8
.LC23:
	.string	"My (rank %d) patch size is %d x %d\n"
	.text
	.p2align 4
	.globl	initialize
	.type	initialize, @function
initialize:
.LFB55:
	.cfi_startproc
	endbr64	
	push	rbp	#
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	mov	rbp, rsp	#,
	.cfi_def_cfa_register 6
	push	r15	#
	push	r14	#
	push	r13	#
	.cfi_offset 15, -24
	.cfi_offset 14, -32
	.cfi_offset 13, -40
	mov	r13, r8	# argv, tmp407
	push	r12	#
	.cfi_offset 12, -48
	mov	r12d, ecx	# argc, tmp406
	push	r10	#
	.cfi_offset 10, -56
	lea	r10, 16[rbp]	#,
	push	rbx	#
	add	rsp, -128	#,
	.cfi_offset 3, -64
# src/stencil_template_parallel.c:369:                   ) {
	mov	rax, QWORD PTR [r10]	# N, N
	mov	rcx, QWORD PTR 8[r10]	# periodic, periodic
	mov	DWORD PTR -68[rbp], esi	# %sfp, tmp404
	mov	rbx, QWORD PTR 32[r10]	# Niterations, Niterations
	mov	rsi, QWORD PTR 48[r10]	# Nsources_local, Nsources_local
	mov	QWORD PTR -152[rbp], rdi	# %sfp, tmp403
	mov	QWORD PTR -160[rbp], rax	# %sfp, N
	mov	rax, QWORD PTR 16[r10]	# output_energy_stat, output_energy_stat
	mov	rdi, QWORD PTR 40[r10]	# Nsources, Nsources
	mov	DWORD PTR -88[rbp], edx	# %sfp, tmp405
	mov	QWORD PTR -112[rbp], rax	# %sfp, output_energy_stat
	mov	rax, QWORD PTR 24[r10]	# neighbours, neighbours
	mov	QWORD PTR -104[rbp], rcx	# %sfp, periodic
	mov	QWORD PTR -128[rbp], rax	# %sfp, neighbours
	mov	QWORD PTR -96[rbp], rbx	# %sfp, Niterations
	mov	QWORD PTR -80[rbp], rdi	# %sfp, Nsources
	mov	QWORD PTR -168[rbp], rsi	# %sfp, Nsources_local
	mov	rdx, QWORD PTR 56[r10]	# Sources_local, Sources_local
	mov	r11, QWORD PTR 64[r10]	# energy_per_source, energy_per_source
	mov	r8, QWORD PTR 72[r10]	# planes, planes
	mov	rax, QWORD PTR 80[r10]	# buffers, buffers
	mov	QWORD PTR -176[rbp], rdx	# %sfp, Sources_local
	mov	QWORD PTR -120[rbp], r11	# %sfp, energy_per_source
	mov	QWORD PTR -136[rbp], r8	# %sfp, planes
	mov	QWORD PTR -144[rbp], rax	# %sfp, buffers
	mov	rax, QWORD PTR fs:40	# tmp426, MEM[(<address-space-1> long unsigned int *)40B]
	mov	QWORD PTR -56[rbp], rax	# D.12529, tmp426
	xor	eax, eax	# tmp426
# src/stencil_template_parallel.c:378:   (*S)[_x_]         = 10000;
	mov	rax, QWORD PTR .LC9[rip]	# tmp243,
# src/stencil_template_parallel.c:390:   if ( planes == NULL ) {
	test	r8, r8	# planes
# src/stencil_template_parallel.c:378:   (*S)[_x_]         = 10000;
	mov	QWORD PTR [r9], rax	# MEM <vector(2) unsigned int> [(unsigned int *)S_152(D)], tmp243
# src/stencil_template_parallel.c:387:   *energy_per_source = 1.0; // amount of energy per source (injected at every step)
	mov	rax, QWORD PTR .LC10[rip]	# tmp499,
# src/stencil_template_parallel.c:382:   *periodic         = 0; // non-periodic boundary
	mov	DWORD PTR [rcx], 0	# *periodic_155(D),
# src/stencil_template_parallel.c:383:   *Nsources         = 4; // number of global heat sources
	mov	DWORD PTR [rdi], 4	# *Nsources_157(D),
# src/stencil_template_parallel.c:384:   *Nsources_local   = 0; // number of local (rank) heat sources (will be set with initialize_sources)
	mov	DWORD PTR [rsi], 0	# *Nsources_local_159(D),
# src/stencil_template_parallel.c:385:   *Sources_local    = NULL; // positions of local heat sources
	mov	QWORD PTR [rdx], 0	# *Sources_local_161(D),
# src/stencil_template_parallel.c:386:   *Niterations      = 1000; // number of iterations
	mov	DWORD PTR [rbx], 1000	# *Niterations_163(D),
# src/stencil_template_parallel.c:387:   *energy_per_source = 1.0; // amount of energy per source (injected at every step)
	mov	QWORD PTR [r11], rax	# *energy_per_source_165(D), tmp499
# src/stencil_template_parallel.c:390:   if ( planes == NULL ) {
	je	.L382	#,
# src/stencil_template_parallel.c:402:     neighbours[i] = MPI_PROC_NULL;
	mov	eax, -2	# tmp251,
# src/stencil_template_parallel.c:397:   planes[OLD].size[0] = planes[OLD].size[1] = 0; // both 0 and 1??
	mov	QWORD PTR 8[r8], 0	# MEM <vector(2) unsigned int> [(unsigned int *)planes_167(D) + 8B],
	mov	r15, r9	# S, tmp408
	lea	r14, .LC16[rip]	# tmp383,
# src/stencil_template_parallel.c:402:     neighbours[i] = MPI_PROC_NULL;
	vmovd	xmm0, eax	# tmp251, tmp251
	mov	rax, QWORD PTR -128[rbp]	# neighbours, %sfp
# src/stencil_template_parallel.c:398:   planes[NEW].size[0] = planes[NEW].size[1] = 0;
	mov	QWORD PTR 24[r8], 0	# MEM <vector(2) unsigned int> [(unsigned int *)planes_167(D) + 24B],
# src/stencil_template_parallel.c:419:       switch( opt )
	lea	rbx, .L286[rip]	# tmp398,
# src/stencil_template_parallel.c:402:     neighbours[i] = MPI_PROC_NULL;
	vpshufd	xmm0, xmm0, 0	# tmp250, tmp251
# src/stencil_template_parallel.c:373:   int verbose = 0;
	mov	DWORD PTR -84[rbp], 0	# %sfp,
# src/stencil_template_parallel.c:402:     neighbours[i] = MPI_PROC_NULL;
	vmovdqu	XMMWORD PTR [rax], xmm0	# MEM <vector(4) int> [(int *)neighbours_193(D)], tmp250
# src/stencil_template_parallel.c:408:       buffers[b][d] = NULL; 
	mov	rax, QWORD PTR -144[rbp]	# buffers, %sfp
	vpxor	xmm0, xmm0, xmm0	# tmp252
# src/stencil_template_parallel.c:371:   int halt = 0;
	mov	DWORD PTR -72[rbp], 0	# %sfp,
# src/stencil_template_parallel.c:408:       buffers[b][d] = NULL; 
	vmovdqu	YMMWORD PTR [rax], ymm0	# MEM <vector(4) long unsigned int> [(double * *)buffers_223(D)], tmp252
	vmovdqu	YMMWORD PTR 32[rax], ymm0	# MEM <vector(4) long unsigned int> [(double * *)buffers_223(D) + 32B], tmp252
	vzeroupper
	.p2align 4,,10
	.p2align 3
.L370:
# src/stencil_template_parallel.c:418:     while((opt = getopt(argc, argv, ":hx:y:e:E:n:o:p:v:")) != -1) {
	mov	rdx, r14	#, tmp383
	mov	rsi, r13	#, argv
	mov	edi, r12d	#, argc
	call	getopt@PLT	#
# src/stencil_template_parallel.c:418:     while((opt = getopt(argc, argv, ":hx:y:e:E:n:o:p:v:")) != -1) {
	cmp	eax, -1	# opt,
	je	.L383	#,
# src/stencil_template_parallel.c:419:       switch( opt )
	sub	eax, 58	#,
	cmp	eax, 63	# tmp254,
	ja	.L370	#,
	movsx	rax, DWORD PTR [rbx+rax*4]	# tmp258,
	add	rax, rbx	# tmp259, tmp398
	notrack jmp	rax	# tmp259
	.section	.rodata
	.align 4
	.align 4
.L286:
	.long	.L296-.L286
	.long	.L370-.L286
	.long	.L370-.L286
	.long	.L370-.L286
	.long	.L370-.L286
	.long	.L295-.L286
	.long	.L370-.L286
	.long	.L370-.L286
	.long	.L370-.L286
	.long	.L370-.L286
	.long	.L370-.L286
	.long	.L294-.L286
	.long	.L370-.L286
	.long	.L370-.L286
	.long	.L370-.L286
	.long	.L370-.L286
	.long	.L370-.L286
	.long	.L370-.L286
	.long	.L370-.L286
	.long	.L370-.L286
	.long	.L370-.L286
	.long	.L370-.L286
	.long	.L370-.L286
	.long	.L370-.L286
	.long	.L370-.L286
	.long	.L370-.L286
	.long	.L370-.L286
	.long	.L370-.L286
	.long	.L370-.L286
	.long	.L370-.L286
	.long	.L370-.L286
	.long	.L370-.L286
	.long	.L370-.L286
	.long	.L370-.L286
	.long	.L370-.L286
	.long	.L370-.L286
	.long	.L370-.L286
	.long	.L370-.L286
	.long	.L370-.L286
	.long	.L370-.L286
	.long	.L370-.L286
	.long	.L370-.L286
	.long	.L370-.L286
	.long	.L293-.L286
	.long	.L370-.L286
	.long	.L370-.L286
	.long	.L292-.L286
	.long	.L370-.L286
	.long	.L370-.L286
	.long	.L370-.L286
	.long	.L370-.L286
	.long	.L370-.L286
	.long	.L291-.L286
	.long	.L290-.L286
	.long	.L289-.L286
	.long	.L370-.L286
	.long	.L370-.L286
	.long	.L370-.L286
	.long	.L370-.L286
	.long	.L370-.L286
	.long	.L288-.L286
	.long	.L370-.L286
	.long	.L287-.L286
	.long	.L285-.L286
	.text
	.p2align 4,,10
	.p2align 3
.L285:
# /usr/include/stdlib.h:483:   return (int) strtol (__nptr, (char **) NULL, 10);
	mov	rdi, QWORD PTR optarg[rip]	#, optarg
	mov	edx, 10	#,
	xor	esi, esi	#
	call	strtol@PLT	#
# src/stencil_template_parallel.c:424:         case 'y': (*S)[_y_] = (uint)atoi(optarg);
	mov	DWORD PTR 4[r15], eax	# (*S_152(D))[1], tmp410
# src/stencil_template_parallel.c:425:           break;
	jmp	.L370	#
	.p2align 4,,10
	.p2align 3
.L287:
# /usr/include/stdlib.h:483:   return (int) strtol (__nptr, (char **) NULL, 10);
	mov	rdi, QWORD PTR optarg[rip]	#, optarg
	mov	edx, 10	#,
	xor	esi, esi	#
	call	strtol@PLT	#
# src/stencil_template_parallel.c:421:         case 'x': (*S)[_x_] = (uint)atoi(optarg);
	mov	DWORD PTR [r15], eax	# (*S_152(D))[0], tmp409
# src/stencil_template_parallel.c:422:           break;
	jmp	.L370	#
	.p2align 4,,10
	.p2align 3
.L288:
# /usr/include/stdlib.h:483:   return (int) strtol (__nptr, (char **) NULL, 10);
	mov	rdi, QWORD PTR optarg[rip]	#, optarg
	mov	edx, 10	#,
	xor	esi, esi	#
	call	strtol@PLT	#
# /usr/include/stdlib.h:483:   return (int) strtol (__nptr, (char **) NULL, 10);
	mov	DWORD PTR -84[rbp], eax	# %sfp, tmp416
# src/stencil_template_parallel.c:443:           break;
	jmp	.L370	#
	.p2align 4,,10
	.p2align 3
.L289:
# /usr/include/stdlib.h:483:   return (int) strtol (__nptr, (char **) NULL, 10);
	mov	rdi, QWORD PTR optarg[rip]	#, optarg
	xor	esi, esi	#
	mov	edx, 10	#,
	call	strtol@PLT	#
# src/stencil_template_parallel.c:439:         case 'p': *periodic = (atoi(optarg) > 0); // 0 or 1
	mov	rcx, QWORD PTR -104[rbp]	# periodic, %sfp
# src/stencil_template_parallel.c:439:         case 'p': *periodic = (atoi(optarg) > 0); // 0 or 1
	test	eax, eax	# tmp415
	setg	al	#, tmp270
	movzx	eax, al	# tmp270, tmp270
# src/stencil_template_parallel.c:439:         case 'p': *periodic = (atoi(optarg) > 0); // 0 or 1
	mov	DWORD PTR [rcx], eax	#* periodic, tmp270
# src/stencil_template_parallel.c:440:           break;
	jmp	.L370	#
	.p2align 4,,10
	.p2align 3
.L290:
# /usr/include/stdlib.h:483:   return (int) strtol (__nptr, (char **) NULL, 10);
	mov	rdi, QWORD PTR optarg[rip]	#, optarg
	xor	esi, esi	#
	mov	edx, 10	#,
	call	strtol@PLT	#
# src/stencil_template_parallel.c:436:         case 'o': *output_energy_stat = (atoi(optarg) > 0); // 0 or 1
	mov	rcx, QWORD PTR -112[rbp]	# output_energy_stat, %sfp
# src/stencil_template_parallel.c:436:         case 'o': *output_energy_stat = (atoi(optarg) > 0); // 0 or 1
	test	eax, eax	# tmp414
	setg	al	#, tmp267
	movzx	eax, al	# tmp267, tmp267
# src/stencil_template_parallel.c:436:         case 'o': *output_energy_stat = (atoi(optarg) > 0); // 0 or 1
	mov	DWORD PTR [rcx], eax	#* output_energy_stat, tmp267
# src/stencil_template_parallel.c:437:           break;
	jmp	.L370	#
	.p2align 4,,10
	.p2align 3
.L291:
# /usr/include/stdlib.h:483:   return (int) strtol (__nptr, (char **) NULL, 10);
	mov	rdi, QWORD PTR optarg[rip]	#, optarg
	mov	edx, 10	#,
	xor	esi, esi	#
	call	strtol@PLT	#
# /usr/include/stdlib.h:483:   return (int) strtol (__nptr, (char **) NULL, 10);
	mov	rcx, QWORD PTR -96[rbp]	# Niterations, %sfp
	mov	DWORD PTR [rcx], eax	#* Niterations, tmp413
# src/stencil_template_parallel.c:434:           break;
	jmp	.L370	#
	.p2align 4,,10
	.p2align 3
.L292:
# src/stencil_template_parallel.c:446:           if ( Me == 0 )
	mov	eax, DWORD PTR -68[rbp]	#, %sfp
	test	eax, eax	#
	je	.L301	#,
.L302:
# src/stencil_template_parallel.c:454:           halt = 1; }
	mov	DWORD PTR -72[rbp], 1	# %sfp,
	jmp	.L370	#
	.p2align 4,,10
	.p2align 3
.L293:
# /usr/include/stdlib.h:483:   return (int) strtol (__nptr, (char **) NULL, 10);
	mov	rdi, QWORD PTR optarg[rip]	#, optarg
	mov	edx, 10	#,
	xor	esi, esi	#
	call	strtol@PLT	#
# /usr/include/stdlib.h:483:   return (int) strtol (__nptr, (char **) NULL, 10);
	mov	rcx, QWORD PTR -80[rbp]	# Nsources, %sfp
	mov	DWORD PTR [rcx], eax	#* Nsources, tmp411
# src/stencil_template_parallel.c:428:           break;
	jmp	.L370	#
	.p2align 4,,10
	.p2align 3
.L294:
# /usr/include/x86_64-linux-gnu/bits/stdlib-float.h:27:   return strtod (__nptr, (char **) NULL);
	mov	rdi, QWORD PTR optarg[rip]	#, optarg
	xor	esi, esi	#
	call	strtod@PLT	#
# src/stencil_template_parallel.c:430:         case 'E': *energy_per_source = atof(optarg);
	mov	rax, QWORD PTR -120[rbp]	# energy_per_source, %sfp
	vmovsd	QWORD PTR [rax], xmm0	# *energy_per_source_165(D), tmp412
# src/stencil_template_parallel.c:431:           break;
	jmp	.L370	#
	.p2align 4,,10
	.p2align 3
.L295:
# /usr/include/x86_64-linux-gnu/bits/stdio2.h:86:   return __printf_chk (__USE_FORTIFY_LEVEL - 1, __fmt, __va_arg_pack ());
	lea	rdi, .LC15[rip]	#,
	call	puts@PLT	#
	jmp	.L370	#
	.p2align 4,,10
	.p2align 3
.L296:
	mov	edx, DWORD PTR optopt[rip]	#, optopt
	lea	rsi, .LC14[rip]	#,
	mov	edi, 2	#,
	xor	eax, eax	#
	call	__printf_chk@PLT	#
	jmp	.L370	#
	.p2align 4,,10
	.p2align 3
.L301:
	lea	rdi, .LC13[rip]	#,
	call	puts@PLT	#
	jmp	.L302	#
	.p2align 4,,10
	.p2align 3
.L383:
# src/stencil_template_parallel.c:469:   if ( halt )
	mov	r13d, DWORD PTR -72[rbp]	#, %sfp
	test	r13d, r13d	#
	jne	.L282	#,
# src/stencil_template_parallel.c:476:   if ( (*S)[_x_] <= 0 || (*S)[_y_] <=0) {
	mov	r14d, DWORD PTR [r15]	# _27, (*S_152(D))[0]
# src/stencil_template_parallel.c:476:   if ( (*S)[_x_] <= 0 || (*S)[_y_] <=0) {
	test	r14d, r14d	# _27
	je	.L304	#,
# src/stencil_template_parallel.c:476:   if ( (*S)[_x_] <= 0 || (*S)[_y_] <=0) {
	mov	r10d, DWORD PTR 4[r15]	# _28, (*S_152(D))[1]
# src/stencil_template_parallel.c:476:   if ( (*S)[_x_] <= 0 || (*S)[_y_] <=0) {
	test	r10d, r10d	# _28
	je	.L304	#,
# src/stencil_template_parallel.c:481:   if (*Niterations <= 0) {
	mov	rax, QWORD PTR -96[rbp]	# Niterations, %sfp
	mov	r12d, DWORD PTR [rax]	#, *Niterations_163(D)
	test	r12d, r12d	#
	jle	.L384	#,
# src/stencil_template_parallel.c:486:   if (*Nsources < 0) {
	mov	rax, QWORD PTR -80[rbp]	# Nsources, %sfp
	mov	ebx, DWORD PTR [rax]	#, *Nsources_157(D)
	test	ebx, ebx	#
	js	.L385	#,
# src/stencil_template_parallel.c:491:   if (verbose < 0) {
	mov	r11d, DWORD PTR -84[rbp]	#, %sfp
	test	r11d, r11d	#
	js	.L386	#,
	vxorps	xmm2, xmm2, xmm2	# tmp419
# src/stencil_template_parallel.c:517:                     ? (double)(*S)[_x_]/(*S)[_y_] 
	mov	eax, r10d	# _28, _28
# src/stencil_template_parallel.c:518:                     : (double)(*S)[_y_]/(*S)[_x_] );
	cmp	r14d, r10d	# _27, _28
# src/stencil_template_parallel.c:517:                     ? (double)(*S)[_x_]/(*S)[_y_] 
	vcvtsi2sd	xmm1, xmm2, rax	# tmp420, tmp419, _28
# src/stencil_template_parallel.c:517:                     ? (double)(*S)[_x_]/(*S)[_y_] 
	mov	eax, r14d	# _27, _27
	vcvtsi2sd	xmm0, xmm2, rax	# tmp421, tmp419, _27
# src/stencil_template_parallel.c:518:                     : (double)(*S)[_y_]/(*S)[_x_] );
	jb	.L313	#,
# src/stencil_template_parallel.c:518:                     : (double)(*S)[_y_]/(*S)[_x_] );
	vdivsd	xmm1, xmm0, xmm1	# iftmp.70_176, _429, _428
# src/stencil_template_parallel.c:521: int dimensions = 2 - (Ntasks <= ((int)formfactor+1) );
	vcvttsd2si	eax, xmm1	# tmp290, iftmp.70_176
# src/stencil_template_parallel.c:521: int dimensions = 2 - (Ntasks <= ((int)formfactor+1) );
	add	eax, 1	# tmp291,
# src/stencil_template_parallel.c:521: int dimensions = 2 - (Ntasks <= ((int)formfactor+1) );
	cmp	DWORD PTR -88[rbp], eax	# %sfp, tmp291
	jg	.L315	#,
# src/stencil_template_parallel.c:547:   (*N)[_x_] = Grid[_x_];
	mov	ecx, DWORD PTR -88[rbp]	# Ntasks, %sfp
	mov	rax, QWORD PTR -160[rbp]	# N, %sfp
# src/stencil_template_parallel.c:554:   int X = Me % Grid[_x_];
	xor	edx, edx	# _63
# src/stencil_template_parallel.c:547:   (*N)[_x_] = Grid[_x_];
	mov	DWORD PTR [rax], ecx	# (*N_186(D))[0], Ntasks
# src/stencil_template_parallel.c:526:       Grid[_x_] = Ntasks, Grid[_y_] = 1;} // row
	mov	esi, ecx	# Grid$0, Ntasks
# src/stencil_template_parallel.c:548:   (*N)[_y_] = Grid[_y_];
	mov	DWORD PTR 4[rax], 1	# (*N_186(D))[1],
# src/stencil_template_parallel.c:554:   int X = Me % Grid[_x_];
	mov	eax, DWORD PTR -68[rbp]	# tmp321, %sfp
	div	ecx	# Ntasks
# src/stencil_template_parallel.c:555:   int Y = Me / Grid[_x_];
	mov	edi, eax	# _66, tmp321
# src/stencil_template_parallel.c:555:   int Y = Me / Grid[_x_];
	mov	r14d, eax	# Y, _66
# src/stencil_template_parallel.c:526:       Grid[_x_] = Ntasks, Grid[_y_] = 1;} // row
	mov	eax, ecx	# Ntasks, Ntasks
# src/stencil_template_parallel.c:554:   int X = Me % Grid[_x_];
	mov	r8d, edx	# _63, _63
# src/stencil_template_parallel.c:560:   if ( Grid[_x_] > 1 )
	cmp	eax, 1	# Ntasks,
# src/stencil_template_parallel.c:554:   int X = Me % Grid[_x_];
	mov	r13d, edx	# X, _63
# src/stencil_template_parallel.c:526:       Grid[_x_] = Ntasks, Grid[_y_] = 1;} // row
	mov	ecx, 1	# Grid$1,
# src/stencil_template_parallel.c:560:   if ( Grid[_x_] > 1 )
	jbe	.L387	#,
.L330:
# src/stencil_template_parallel.c:562:       if ( *periodic ) {       
	mov	rax, QWORD PTR -104[rbp]	# periodic, %sfp
	mov	r10d, DWORD PTR [rax]	#, *periodic_155(D)
	test	r10d, r10d	#
	je	.L335	#,
# src/stencil_template_parallel.c:563:         neighbours[EAST]  = Y*Grid[_x_] + (Me + 1) % Grid[_x_];
	mov	ebx, DWORD PTR -68[rbp]	# Me, %sfp
# src/stencil_template_parallel.c:563:         neighbours[EAST]  = Y*Grid[_x_] + (Me + 1) % Grid[_x_];
	xor	edx, edx	# tmp337
# src/stencil_template_parallel.c:563:         neighbours[EAST]  = Y*Grid[_x_] + (Me + 1) % Grid[_x_];
	lea	eax, 1[rbx]	# tmp335,
# src/stencil_template_parallel.c:563:         neighbours[EAST]  = Y*Grid[_x_] + (Me + 1) % Grid[_x_];
	div	esi	# Grid$0
# src/stencil_template_parallel.c:563:         neighbours[EAST]  = Y*Grid[_x_] + (Me + 1) % Grid[_x_];
	mov	eax, esi	# tmp339, Grid$0
	imul	eax, edi	# tmp339, _66
# src/stencil_template_parallel.c:564:         neighbours[WEST]  = (X%Grid[_x_] > 0 ? Me-1 : (Y+1)*Grid[_x_]-1); }
	test	r8d, r8d	# _63
# src/stencil_template_parallel.c:563:         neighbours[EAST]  = Y*Grid[_x_] + (Me + 1) % Grid[_x_];
	lea	r9d, [rdx+rax]	# _75,
# src/stencil_template_parallel.c:564:         neighbours[WEST]  = (X%Grid[_x_] > 0 ? Me-1 : (Y+1)*Grid[_x_]-1); }
	lea	edx, -1[rbx]	# iftmp.83_134,
# src/stencil_template_parallel.c:564:         neighbours[WEST]  = (X%Grid[_x_] > 0 ? Me-1 : (Y+1)*Grid[_x_]-1); }
	je	.L388	#,
.L340:
# src/stencil_template_parallel.c:563:         neighbours[EAST]  = Y*Grid[_x_] + (Me + 1) % Grid[_x_];
	mov	rax, QWORD PTR -128[rbp]	# neighbours, %sfp
	mov	DWORD PTR 8[rax], r9d	# MEM[(int *)neighbours_193(D) + 8B], _75
# src/stencil_template_parallel.c:564:         neighbours[WEST]  = (X%Grid[_x_] > 0 ? Me-1 : (Y+1)*Grid[_x_]-1); }
	mov	DWORD PTR 12[rax], edx	# MEM[(int *)neighbours_193(D) + 12B], iftmp.83_134
.L332:
# src/stencil_template_parallel.c:571:   if ( Grid[_y_] > 1 )
	cmp	ecx, 1	# Grid$1,
	jbe	.L331	#,
# src/stencil_template_parallel.c:573:       if ( *periodic ) {      
	mov	rax, QWORD PTR -104[rbp]	# periodic, %sfp
	mov	r9d, DWORD PTR [rax]	#, *periodic_155(D)
	test	r9d, r9d	#
	je	.L342	#,
# src/stencil_template_parallel.c:574:         neighbours[NORTH] = (Ntasks + Me - Grid[_x_]) % Ntasks;
	mov	ebx, DWORD PTR -88[rbp]	# Ntasks, %sfp
	mov	eax, DWORD PTR -68[rbp]	# Me, %sfp
# src/stencil_template_parallel.c:574:         neighbours[NORTH] = (Ntasks + Me - Grid[_x_]) % Ntasks;
	xor	edx, edx	# tmp349
# src/stencil_template_parallel.c:574:         neighbours[NORTH] = (Ntasks + Me - Grid[_x_]) % Ntasks;
	lea	r9d, [rbx+rax]	# _88,
# src/stencil_template_parallel.c:574:         neighbours[NORTH] = (Ntasks + Me - Grid[_x_]) % Ntasks;
	mov	eax, r9d	# tmp347, _88
	sub	eax, esi	# tmp347, Grid$0
# src/stencil_template_parallel.c:574:         neighbours[NORTH] = (Ntasks + Me - Grid[_x_]) % Ntasks;
	div	ebx	# Ntasks
# src/stencil_template_parallel.c:575:         neighbours[SOUTH] = (Ntasks + Me + Grid[_x_]) % Ntasks; }
	lea	eax, [rsi+r9]	# tmp351,
# src/stencil_template_parallel.c:574:         neighbours[NORTH] = (Ntasks + Me - Grid[_x_]) % Ntasks;
	mov	r10d, edx	# tmp349, tmp349
# src/stencil_template_parallel.c:575:         neighbours[SOUTH] = (Ntasks + Me + Grid[_x_]) % Ntasks; }
	xor	edx, edx	# tmp353
	div	ebx	# Ntasks
# src/stencil_template_parallel.c:574:         neighbours[NORTH] = (Ntasks + Me - Grid[_x_]) % Ntasks;
	mov	rax, QWORD PTR -128[rbp]	# neighbours, %sfp
# src/stencil_template_parallel.c:574:         neighbours[NORTH] = (Ntasks + Me - Grid[_x_]) % Ntasks;
	vmovd	xmm4, r10d	# tmp349, tmp349
	vpinsrd	xmm0, xmm4, edx, 1	# vect__92.431, tmp349, tmp353
# src/stencil_template_parallel.c:574:         neighbours[NORTH] = (Ntasks + Me - Grid[_x_]) % Ntasks;
	vmovq	QWORD PTR [rax], xmm0	# MEM <vector(2) int> [(int *)neighbours_193(D)], vect__92.431
.L331:
	mov	eax, DWORD PTR [r15]	# (*S_152(D))[0], (*S_152(D))[0]
	xor	edx, edx	# tmp357
	div	esi	# Grid$0
# src/stencil_template_parallel.c:599:   mysize[_x_] = s + (X < r); // the first r columns get one extra
	cmp	r8d, edx	# _63, tmp357
# src/stencil_template_parallel.c:599:   mysize[_x_] = s + (X < r); // the first r columns get one extra
	adc	eax, 0	# tmp356
	xor	edx, edx	# tmp362
	mov	ebx, eax	# _104, tmp356
# src/stencil_template_parallel.c:599:   mysize[_x_] = s + (X < r); // the first r columns get one extra
	mov	DWORD PTR -64[rbp], eax	# mysize[0], _104
	mov	eax, DWORD PTR 4[r15]	# (*S_152(D))[1], (*S_152(D))[1]
	div	ecx	# Grid$1
# src/stencil_template_parallel.c:606:   mysize[_y_] = s + (Y < r);
	cmp	edi, edx	# _66, tmp362
# src/stencil_template_parallel.c:617:   if ( verbose > 0 ) {// if option was added from commandline
	mov	edi, DWORD PTR -84[rbp]	#, %sfp
# src/stencil_template_parallel.c:606:   mysize[_y_] = s + (Y < r);
	adc	eax, 0	# tmp361
	mov	r12d, eax	# _108, tmp361
# src/stencil_template_parallel.c:606:   mysize[_y_] = s + (Y < r);
	mov	DWORD PTR -60[rbp], eax	# mysize[1], _108
# src/stencil_template_parallel.c:611:   planes[OLD].size[0] = mysize[0];
	mov	rax, QWORD PTR -136[rbp]	# planes, %sfp
# src/stencil_template_parallel.c:617:   if ( verbose > 0 ) {// if option was added from commandline
	test	edi, edi	#
# src/stencil_template_parallel.c:611:   planes[OLD].size[0] = mysize[0];
	mov	DWORD PTR 8[rax], ebx	# planes_167(D)->size[0], _104
# src/stencil_template_parallel.c:612:   planes[OLD].size[1] = mysize[1];
	mov	DWORD PTR 12[rax], r12d	# planes_167(D)->size[1], _108
# src/stencil_template_parallel.c:613:   planes[NEW].size[0] = mysize[0];
	mov	DWORD PTR 24[rax], ebx	# MEM[(struct plane_t *)planes_167(D) + 16B].size[0], _104
# src/stencil_template_parallel.c:614:   planes[NEW].size[1] = mysize[1];
	mov	DWORD PTR 28[rax], r12d	# MEM[(struct plane_t *)planes_167(D) + 16B].size[1], _108
# src/stencil_template_parallel.c:617:   if ( verbose > 0 ) {// if option was added from commandline
	je	.L349	#,
# src/stencil_template_parallel.c:618:     if ( Me == 0 ) {
	mov	edx, DWORD PTR -68[rbp]	#, %sfp
	test	edx, edx	#
	je	.L389	#,
.L348:
# src/stencil_template_parallel.c:623:     MPI_Barrier(*Comm); // so that one rank prints a time
	mov	rax, QWORD PTR -152[rbp]	# Comm, %sfp
	mov	rdi, QWORD PTR [rax]	# *Comm_218(D), *Comm_218(D)
	call	MPI_Barrier@PLT	#
# src/stencil_template_parallel.c:628:     for ( int t = 0; t < Ntasks; t++ ) {
	mov	eax, DWORD PTR -88[rbp]	#, %sfp
	test	eax, eax	#
	jle	.L349	#,
# /usr/include/x86_64-linux-gnu/bits/stdio2.h:86:   return __printf_chk (__USE_FORTIFY_LEVEL - 1, __fmt, __va_arg_pack ());
	xor	eax, eax	# t
	mov	r15d, r12d	# _108, _108
	mov	r12d, ebx	# _104, _104
	mov	ebx, eax	# t, t
	jmp	.L352	#
	.p2align 4,,10
	.p2align 3
.L350:
# src/stencil_template_parallel.c:642:     MPI_Barrier(*Comm); 
	mov	rax, QWORD PTR -152[rbp]	# Comm, %sfp
# src/stencil_template_parallel.c:628:     for ( int t = 0; t < Ntasks; t++ ) {
	add	ebx, 1	# t,
# src/stencil_template_parallel.c:642:     MPI_Barrier(*Comm); 
	mov	rdi, QWORD PTR [rax]	# *Comm_218(D), *Comm_218(D)
	call	MPI_Barrier@PLT	#
# src/stencil_template_parallel.c:628:     for ( int t = 0; t < Ntasks; t++ ) {
	cmp	DWORD PTR -88[rbp], ebx	# %sfp, t
	je	.L349	#,
.L352:
# src/stencil_template_parallel.c:629:       if ( t == Me ) {
	cmp	DWORD PTR -68[rbp], ebx	# %sfp, t
	jne	.L350	#,
# /usr/include/x86_64-linux-gnu/bits/stdio2.h:86:   return __printf_chk (__USE_FORTIFY_LEVEL - 1, __fmt, __va_arg_pack ());
	mov	rcx, QWORD PTR -128[rbp]	# neighbours, %sfp
	sub	rsp, 8	#,
	mov	ebx, DWORD PTR -68[rbp]	# Me, %sfp
	mov	r8d, r14d	#, Y
	lea	rsi, .LC22[rip]	#,
	mov	edi, 2	#,
	mov	eax, DWORD PTR 12[rcx]	# MEM[(int *)neighbours_193(D) + 12B], MEM[(int *)neighbours_193(D) + 12B]
	mov	r9d, DWORD PTR [rcx]	# *neighbours_193(D), *neighbours_193(D)
	mov	edx, ebx	#, Me
	push	rax	# MEM[(int *)neighbours_193(D) + 12B]
	mov	eax, DWORD PTR 4[rcx]	# MEM[(int *)neighbours_193(D) + 4B], MEM[(int *)neighbours_193(D) + 4B]
	push	rax	# MEM[(int *)neighbours_193(D) + 4B]
	mov	eax, DWORD PTR 8[rcx]	# MEM[(int *)neighbours_193(D) + 8B], MEM[(int *)neighbours_193(D) + 8B]
	mov	ecx, r13d	#, X
	push	rax	# MEM[(int *)neighbours_193(D) + 8B]
	xor	eax, eax	#
	call	__printf_chk@PLT	#
	add	rsp, 32	#,
	mov	r8d, r15d	#, _108
	mov	ecx, r12d	#, _104
	mov	edx, ebx	#, Me
	lea	rsi, .LC23[rip]	#,
	mov	edi, 2	#,
	xor	eax, eax	#
	call	__printf_chk@PLT	#
# src/stencil_template_parallel.c:637:         fflush(stdout);
	mov	rdi, QWORD PTR stdout[rip]	#, stdout
	call	fflush@PLT	#
	jmp	.L350	#
.L304:
# /usr/include/x86_64-linux-gnu/bits/stdio2.h:86:   return __printf_chk (__USE_FORTIFY_LEVEL - 1, __fmt, __va_arg_pack ());
	lea	rsi, .LC17[rip]	# tmp278,
	mov	edi, 2	#,
	xor	eax, eax	#
	call	__printf_chk@PLT	#
.L282:
# src/stencil_template_parallel.c:393:     return 1;
	mov	DWORD PTR -72[rbp], 1	# %sfp,
.L280:
# src/stencil_template_parallel.c:662: }
	mov	rax, QWORD PTR -56[rbp]	# tmp427, D.12529
	sub	rax, QWORD PTR fs:40	# tmp427, MEM[(<address-space-1> long unsigned int *)40B]
	jne	.L390	#,
	mov	eax, DWORD PTR -72[rbp]	#, %sfp
	lea	rsp, -48[rbp]	#,
	pop	rbx	#
	pop	r10	#
	pop	r12	#
	pop	r13	#
	pop	r14	#
	pop	r15	#
	pop	rbp	#
	.cfi_remember_state
	.cfi_def_cfa 7, 8
	ret	
.L349:
	.cfi_restore_state
# src/stencil_template_parallel.c:649:   ret = memory_allocate(neighbours, // who the N/S/E/W neighbours are
	mov	rcx, QWORD PTR -136[rbp]	#, %sfp
	mov	rdx, QWORD PTR -144[rbp]	#, %sfp
	mov	rsi, QWORD PTR -160[rbp]	#, %sfp
	mov	rdi, QWORD PTR -128[rbp]	#, %sfp
	call	memory_allocate	#
# src/stencil_template_parallel.c:658:   ret = initialize_sources( Me, Ntasks, Comm, mysize, 
	mov	rax, QWORD PTR -80[rbp]	# Nsources, %sfp
	sub	rsp, 8	#,
	mov	esi, DWORD PTR -88[rbp]	#, %sfp
	mov	r9, QWORD PTR -168[rbp]	#, %sfp
	mov	rdx, QWORD PTR -152[rbp]	#, %sfp
	lea	rcx, -64[rbp]	# tmp367,
	mov	r8d, DWORD PTR [rax]	# *Nsources_157(D), *Nsources_157(D)
	mov	edi, DWORD PTR -68[rbp]	#, %sfp
	push	QWORD PTR -176[rbp]	# %sfp
	call	initialize_sources	#
# src/stencil_template_parallel.c:661:   return 0;  
	pop	rcx	#
	pop	rsi	#
	jmp	.L280	#
.L382:
# /usr/include/x86_64-linux-gnu/bits/stdio2.h:79:   return __fprintf_chk (__stream, __USE_FORTIFY_LEVEL - 1, __fmt,
	mov	rcx, QWORD PTR stderr[rip]	#, stderr
	mov	edx, 45	#,
	mov	esi, 1	#,
	lea	rdi, .LC11[rip]	# tmp246,
	call	fwrite@PLT	#
# src/stencil_template_parallel.c:393:     return 1;
	jmp	.L282	#
.L313:
# src/stencil_template_parallel.c:518:                     : (double)(*S)[_y_]/(*S)[_x_] );
	vdivsd	xmm1, xmm1, xmm0	# iftmp.70_176, _428, _429
# src/stencil_template_parallel.c:521: int dimensions = 2 - (Ntasks <= ((int)formfactor+1) );
	vcvttsd2si	eax, xmm1	# tmp292, iftmp.70_176
# src/stencil_template_parallel.c:521: int dimensions = 2 - (Ntasks <= ((int)formfactor+1) );
	add	eax, 1	# tmp293,
# src/stencil_template_parallel.c:521: int dimensions = 2 - (Ntasks <= ((int)formfactor+1) );
	cmp	DWORD PTR -88[rbp], eax	# %sfp, tmp293
	jg	.L315	#,
# src/stencil_template_parallel.c:547:   (*N)[_x_] = Grid[_x_];
	mov	rax, QWORD PTR -160[rbp]	# N, %sfp
# src/stencil_template_parallel.c:548:   (*N)[_y_] = Grid[_y_];
	mov	ecx, DWORD PTR -88[rbp]	# Ntasks, %sfp
# src/stencil_template_parallel.c:554:   int X = Me % Grid[_x_];
	xor	r13d, r13d	# X
# src/stencil_template_parallel.c:554:   int X = Me % Grid[_x_];
	xor	r8d, r8d	# _63
# src/stencil_template_parallel.c:555:   int Y = Me / Grid[_x_];
	mov	edi, DWORD PTR -68[rbp]	# Me, %sfp
# src/stencil_template_parallel.c:528:       Grid[_x_] = 1, Grid[_y_] = Ntasks;} // column
	mov	esi, 1	# Grid$0,
# src/stencil_template_parallel.c:547:   (*N)[_x_] = Grid[_x_];
	mov	DWORD PTR [rax], 1	# (*N_186(D))[0],
# src/stencil_template_parallel.c:548:   (*N)[_y_] = Grid[_y_];
	mov	DWORD PTR 4[rax], ecx	# (*N_186(D))[1], Ntasks
# src/stencil_template_parallel.c:555:   int Y = Me / Grid[_x_];
	mov	r14d, edi	# Y, Me
	jmp	.L332	#
.L315:
# src/stencil_template_parallel.c:677:   uint _A_ = A; //working copy
	mov	eax, DWORD PTR -88[rbp]	# Ntasks, %sfp
# src/stencil_template_parallel.c:675:   int N = 0; // counter for prime factors
	xor	r12d, r12d	# N
# src/stencil_template_parallel.c:676:   int f = 2; // starting from the smallest prime
	mov	ebx, 2	# f,
# src/stencil_template_parallel.c:680:   while ( f < A ) { // we try every f from 2 to A-1
	cmp	eax, 2	# Ntasks,
# src/stencil_template_parallel.c:677:   uint _A_ = A; //working copy
	mov	ecx, eax	# _A_, Ntasks
	mov	esi, eax	# Ntasks, Ntasks
# src/stencil_template_parallel.c:680:   while ( f < A ) { // we try every f from 2 to A-1
	jbe	.L318	#,
	.p2align 4,,10
	.p2align 3
.L317:
# src/stencil_template_parallel.c:681:       while( _A_ % f == 0 ) { // while f divides what is left
	xor	edx, edx	# tmp303
	mov	eax, ecx	# tmp304, _A_
	div	ebx	# f
# src/stencil_template_parallel.c:681:       while( _A_ % f == 0 ) { // while f divides what is left
	test	edx, edx	# tmp303
	jne	.L321	#,
	.p2align 4,,10
	.p2align 3
.L319:
# src/stencil_template_parallel.c:683: 	      _A_ /= f; 
	mov	eax, ecx	# _A_, _A_
	xor	edx, edx	# tmp295
# src/stencil_template_parallel.c:682: 	      N++; // count this factor
	add	r12d, 1	# N,
# src/stencil_template_parallel.c:683: 	      _A_ /= f; 
	div	ebx	# f
# src/stencil_template_parallel.c:681:       while( _A_ % f == 0 ) { // while f divides what is left
	xor	edx, edx	# tmp297
# src/stencil_template_parallel.c:683: 	      _A_ /= f; 
	mov	ecx, eax	# _A_, _A_
# src/stencil_template_parallel.c:681:       while( _A_ % f == 0 ) { // while f divides what is left
	div	ebx	# f
# src/stencil_template_parallel.c:681:       while( _A_ % f == 0 ) { // while f divides what is left
	test	edx, edx	# tmp297
	je	.L319	#,
.L321:
# src/stencil_template_parallel.c:687:       f++; 
	lea	r13d, 1[rbx]	# f,
# src/stencil_template_parallel.c:680:   while ( f < A ) { // we try every f from 2 to A-1
	cmp	esi, r13d	# Ntasks, f
	je	.L320	#,
	mov	ebx, r13d	# f, f
	jmp	.L317	#
.L335:
# src/stencil_template_parallel.c:567:         neighbours[EAST]  = ( X < Grid[_x_]-1 ? Me+1 : MPI_PROC_NULL );
	mov	ebx, DWORD PTR -68[rbp]	# Me, %sfp
# src/stencil_template_parallel.c:567:         neighbours[EAST]  = ( X < Grid[_x_]-1 ? Me+1 : MPI_PROC_NULL );
	lea	eax, -1[rsi]	# tmp342,
# src/stencil_template_parallel.c:568:         neighbours[WEST]  = ( X > 0 ? (Me-1)%Ntasks : MPI_PROC_NULL ); }  
	mov	edx, -2	# iftmp.83_134,
# src/stencil_template_parallel.c:567:         neighbours[EAST]  = ( X < Grid[_x_]-1 ? Me+1 : MPI_PROC_NULL );
	cmp	r8d, eax	# _63, tmp342
	mov	eax, -2	# tmp387,
	lea	r9d, 1[rbx]	# tmp386,
	cmovnb	r9d, eax	# tmp386,, _75, tmp387
# src/stencil_template_parallel.c:568:         neighbours[WEST]  = ( X > 0 ? (Me-1)%Ntasks : MPI_PROC_NULL ); }  
	test	r13d, r13d	# X
	jle	.L340	#,
# src/stencil_template_parallel.c:568:         neighbours[WEST]  = ( X > 0 ? (Me-1)%Ntasks : MPI_PROC_NULL ); }  
	lea	eax, -1[rbx]	# tmp343,
# src/stencil_template_parallel.c:568:         neighbours[WEST]  = ( X > 0 ? (Me-1)%Ntasks : MPI_PROC_NULL ); }  
	cdq
	idiv	DWORD PTR -88[rbp]	# %sfp
	jmp	.L340	#
.L342:
# src/stencil_template_parallel.c:578:         neighbours[NORTH] = ( Y > 0 ? Me - Grid[_x_]: MPI_PROC_NULL );
	mov	ebx, DWORD PTR -68[rbp]	# Me, %sfp
	mov	edx, -2	# tmp389,
# src/stencil_template_parallel.c:578:         neighbours[NORTH] = ( Y > 0 ? Me - Grid[_x_]: MPI_PROC_NULL );
	mov	r11, QWORD PTR -128[rbp]	# neighbours, %sfp
# src/stencil_template_parallel.c:579:         neighbours[SOUTH] = ( Y < Grid[_y_]-1 ? Me + Grid[_x_] : MPI_PROC_NULL ); }
	lea	r9d, -1[rcx]	# tmp355,
# src/stencil_template_parallel.c:578:         neighbours[NORTH] = ( Y > 0 ? Me - Grid[_x_]: MPI_PROC_NULL );
	mov	eax, ebx	# tmp388, Me
	sub	eax, esi	# tmp388, Grid$0
	test	r14d, r14d	# Y
	cmovle	eax, edx	# tmp388,, iftmp.90_137, tmp389
# src/stencil_template_parallel.c:579:         neighbours[SOUTH] = ( Y < Grid[_y_]-1 ? Me + Grid[_x_] : MPI_PROC_NULL ); }
	cmp	edi, r9d	# _66, tmp355
# src/stencil_template_parallel.c:578:         neighbours[NORTH] = ( Y > 0 ? Me - Grid[_x_]: MPI_PROC_NULL );
	mov	DWORD PTR [r11], eax	# *neighbours_193(D), iftmp.90_137
# src/stencil_template_parallel.c:579:         neighbours[SOUTH] = ( Y < Grid[_y_]-1 ? Me + Grid[_x_] : MPI_PROC_NULL ); }
	lea	eax, [rsi+rbx]	# tmp392,
	cmovnb	eax, edx	# tmp392,, iftmp.92_138, tmp389
# src/stencil_template_parallel.c:579:         neighbours[SOUTH] = ( Y < Grid[_y_]-1 ? Me + Grid[_x_] : MPI_PROC_NULL ); }
	mov	DWORD PTR 4[r11], eax	# MEM[(int *)neighbours_193(D) + 4B], iftmp.92_138
	jmp	.L331	#
.L320:
# src/stencil_template_parallel.c:691:   uint *_factors_ = (uint*)malloc( N * sizeof(uint) );
	movsx	rdi, r12d	# N, N
	mov	DWORD PTR -96[rbp], r10d	# %sfp, _28
	sal	rdi, 2	# tmp306,
	vmovsd	QWORD PTR -112[rbp], xmm1	# %sfp, iftmp.70_176
	call	malloc@PLT	#
# src/stencil_template_parallel.c:696:   _A_ = A;
	mov	esi, DWORD PTR -88[rbp]	# _A_, %sfp
	mov	r10d, DWORD PTR -96[rbp]	# _28, %sfp
# src/stencil_template_parallel.c:694:   N   = 0;
	xor	r8d, r8d	# N
	vmovsd	xmm1, QWORD PTR -112[rbp]	# iftmp.70_176, %sfp
# src/stencil_template_parallel.c:691:   uint *_factors_ = (uint*)malloc( N * sizeof(uint) );
	mov	rdi, rax	# _factors_, tmp418
# src/stencil_template_parallel.c:695:   f   = 2;
	mov	ecx, 2	# f,
	vxorps	xmm2, xmm2, xmm2	# tmp419
	.p2align 4,,10
	.p2align 3
.L322:
# src/stencil_template_parallel.c:700:       while( _A_ % f == 0 ) {
	xor	edx, edx	# tmp318
	mov	eax, esi	# tmp319, _A_
	lea	r9d, 1[r8]	# tmp308,
	div	ecx	# f
	movsx	r9, r9d	# ivtmp.445, tmp308
# src/stencil_template_parallel.c:700:       while( _A_ % f == 0 ) {
	test	edx, edx	# tmp318
	jne	.L326	#,
	.p2align 4,,10
	.p2align 3
.L323:
# src/stencil_template_parallel.c:702: 	      _A_ /= f; }
	mov	eax, esi	# _A_, _A_
	xor	edx, edx	# tmp310
# src/stencil_template_parallel.c:701: 	      _factors_[N++] = f;
	mov	DWORD PTR -4[rdi+r9*4], ecx	# MEM[(uint *)_factors__285 + -4B + ivtmp.445_206 * 4], f
	mov	r8, r9	# ivtmp.445, ivtmp.445
# src/stencil_template_parallel.c:702: 	      _A_ /= f; }
	div	ecx	# f
# src/stencil_template_parallel.c:700:       while( _A_ % f == 0 ) {
	xor	edx, edx	# tmp312
# src/stencil_template_parallel.c:700:       while( _A_ % f == 0 ) {
	add	r9, 1	# ivtmp.445,
# src/stencil_template_parallel.c:702: 	      _A_ /= f; }
	mov	esi, eax	# _A_, _A_
# src/stencil_template_parallel.c:700:       while( _A_ % f == 0 ) {
	div	ecx	# f
# src/stencil_template_parallel.c:700:       while( _A_ % f == 0 ) {
	test	edx, edx	# tmp312
	je	.L323	#,
.L326:
# src/stencil_template_parallel.c:699:   while ( f < A ) {
	cmp	ebx, ecx	# f, f
# src/stencil_template_parallel.c:703:       f++; }
	lea	eax, 1[rcx]	# f,
# src/stencil_template_parallel.c:699:   while ( f < A ) {
	je	.L324	#,
	mov	ecx, eax	# f, f
	jmp	.L322	#
.L324:
# src/stencil_template_parallel.c:538:   for ( int i = 0; (i < Nf) && ((Ntasks/first)/first > formfactor); i++ ) 
	test	r12d, r12d	# N
	jle	.L318	#,
# src/stencil_template_parallel.c:538:   for ( int i = 0; (i < Nf) && ((Ntasks/first)/first > formfactor); i++ ) 
	vcvtsi2sd	xmm0, xmm2, r13d	# tmp422, tmp419, f
# src/stencil_template_parallel.c:538:   for ( int i = 0; (i < Nf) && ((Ntasks/first)/first > formfactor); i++ ) 
	vcomisd	xmm0, xmm1	# tmp320, iftmp.70_176
	jbe	.L318	#,
	lea	eax, -1[r12]	# tmp324,
# src/stencil_template_parallel.c:535:   uint  first = 1;
	mov	r9d, DWORD PTR -88[rbp]	# Ntasks, %sfp
	mov	ecx, 1	# Grid$1,
	lea	r8, [rdi+rax*4]	# _387,
	jmp	.L333	#
	.p2align 4,,10
	.p2align 3
.L391:
# src/stencil_template_parallel.c:538:   for ( int i = 0; (i < Nf) && ((Ntasks/first)/first > formfactor); i++ ) 
	xor	edx, edx	# tmp330
# src/stencil_template_parallel.c:538:   for ( int i = 0; (i < Nf) && ((Ntasks/first)/first > formfactor); i++ ) 
	add	rdi, 4	# ivtmp.440,
# src/stencil_template_parallel.c:538:   for ( int i = 0; (i < Nf) && ((Ntasks/first)/first > formfactor); i++ ) 
	div	ecx	# Grid$1
# src/stencil_template_parallel.c:538:   for ( int i = 0; (i < Nf) && ((Ntasks/first)/first > formfactor); i++ ) 
	vcvtsi2sd	xmm0, xmm2, eax	# tmp423, tmp419, tmp329
# src/stencil_template_parallel.c:538:   for ( int i = 0; (i < Nf) && ((Ntasks/first)/first > formfactor); i++ ) 
	vcomisd	xmm0, xmm1	# tmp331, iftmp.70_176
	jbe	.L328	#,
.L333:
# src/stencil_template_parallel.c:539:     first *= factors[i];
	imul	ecx, DWORD PTR [rdi]	# Grid$1, MEM[(uint *)_381]
# src/stencil_template_parallel.c:538:   for ( int i = 0; (i < Nf) && ((Ntasks/first)/first > formfactor); i++ ) 
	xor	edx, edx	# tmp327
	mov	eax, r9d	# Grid$0, Ntasks
	div	ecx	# Grid$1
# src/stencil_template_parallel.c:538:   for ( int i = 0; (i < Nf) && ((Ntasks/first)/first > formfactor); i++ ) 
	cmp	rdi, r8	# ivtmp.440, _387
# src/stencil_template_parallel.c:538:   for ( int i = 0; (i < Nf) && ((Ntasks/first)/first > formfactor); i++ ) 
	mov	esi, eax	# Grid$0, Grid$0
# src/stencil_template_parallel.c:538:   for ( int i = 0; (i < Nf) && ((Ntasks/first)/first > formfactor); i++ ) 
	jne	.L391	#,
.L328:
# src/stencil_template_parallel.c:541:     if ( (*S)[_x_] > (*S)[_y_] ) 
	cmp	r10d, r14d	# _28, _27
	jnb	.L392	#,
.L334:
# src/stencil_template_parallel.c:547:   (*N)[_x_] = Grid[_x_];
	mov	rax, QWORD PTR -160[rbp]	# N, %sfp
	vmovd	xmm3, esi	# Grid$0, Grid$0
# src/stencil_template_parallel.c:554:   int X = Me % Grid[_x_];
	xor	edx, edx	# _63
# src/stencil_template_parallel.c:547:   (*N)[_x_] = Grid[_x_];
	vpinsrd	xmm0, xmm3, ecx, 1	# tmp332, Grid$0, Grid$1
	vmovq	QWORD PTR [rax], xmm0	# MEM <vector(2) unsigned int> [(unsigned int *)N_186(D)], tmp332
# src/stencil_template_parallel.c:554:   int X = Me % Grid[_x_];
	mov	eax, DWORD PTR -68[rbp]	# tmp333, %sfp
	div	esi	# Grid$0
# src/stencil_template_parallel.c:560:   if ( Grid[_x_] > 1 )
	cmp	esi, 1	# Grid$0,
# src/stencil_template_parallel.c:554:   int X = Me % Grid[_x_];
	mov	r8d, edx	# _63, _63
# src/stencil_template_parallel.c:554:   int X = Me % Grid[_x_];
	mov	r13d, edx	# X, _63
# src/stencil_template_parallel.c:555:   int Y = Me / Grid[_x_];
	mov	edi, eax	# _66, tmp333
# src/stencil_template_parallel.c:555:   int Y = Me / Grid[_x_];
	mov	r14d, eax	# Y, _66
# src/stencil_template_parallel.c:560:   if ( Grid[_x_] > 1 )
	ja	.L330	#,
	xor	r13d, r13d	# X
	xor	r8d, r8d	# _63
	jmp	.L332	#
.L386:
# /usr/include/x86_64-linux-gnu/bits/stdio2.h:86:   return __printf_chk (__USE_FORTIFY_LEVEL - 1, __fmt, __va_arg_pack ());
	lea	rsi, .LC20[rip]	# tmp281,
	mov	edi, 2	#,
	xor	eax, eax	#
	call	__printf_chk@PLT	#
# src/stencil_template_parallel.c:493:     return 1;
	jmp	.L282	#
.L392:
# src/stencil_template_parallel.c:541:     if ( (*S)[_x_] > (*S)[_y_] ) 
	mov	eax, esi	# Grid$0, Grid$0
# src/stencil_template_parallel.c:544:       Grid[_x_] = first, Grid[_y_] = Ntasks/first;
	mov	esi, ecx	# Grid$0, Grid$1
# src/stencil_template_parallel.c:544:       Grid[_x_] = first, Grid[_y_] = Ntasks/first;
	mov	ecx, eax	# Grid$1, Grid$0
	jmp	.L334	#
.L388:
# src/stencil_template_parallel.c:564:         neighbours[WEST]  = (X%Grid[_x_] > 0 ? Me-1 : (Y+1)*Grid[_x_]-1); }
	lea	edx, 1[r14]	# tmp340,
# src/stencil_template_parallel.c:564:         neighbours[WEST]  = (X%Grid[_x_] > 0 ? Me-1 : (Y+1)*Grid[_x_]-1); }
	imul	edx, esi	# tmp341, Grid$0
# src/stencil_template_parallel.c:564:         neighbours[WEST]  = (X%Grid[_x_] > 0 ? Me-1 : (Y+1)*Grid[_x_]-1); }
	sub	edx, 1	# iftmp.83_134,
	jmp	.L340	#
.L384:
# /usr/include/x86_64-linux-gnu/bits/stdio2.h:86:   return __printf_chk (__USE_FORTIFY_LEVEL - 1, __fmt, __va_arg_pack ());
	lea	rsi, .LC18[rip]	# tmp279,
	mov	edi, 2	#,
	xor	eax, eax	#
	call	__printf_chk@PLT	#
# src/stencil_template_parallel.c:483:     return 1;
	jmp	.L282	#
.L385:
# /usr/include/x86_64-linux-gnu/bits/stdio2.h:86:   return __printf_chk (__USE_FORTIFY_LEVEL - 1, __fmt, __va_arg_pack ());
	lea	rsi, .LC19[rip]	# tmp280,
	mov	edi, 2	#,
	xor	eax, eax	#
	call	__printf_chk@PLT	#
# src/stencil_template_parallel.c:488:     return 1;
	jmp	.L282	#
.L389:
# /usr/include/x86_64-linux-gnu/bits/stdio2.h:86:   return __printf_chk (__USE_FORTIFY_LEVEL - 1, __fmt, __va_arg_pack ());
	mov	edx, esi	#, Grid$0
	mov	edi, 2	#,
	lea	rsi, .LC21[rip]	# tmp368,
	xor	eax, eax	#
	call	__printf_chk@PLT	#
# src/stencil_template_parallel.c:620:       fflush(stdout);
	mov	rdi, QWORD PTR stdout[rip]	#, stdout
	call	fflush@PLT	#
	jmp	.L348	#
.L387:
# src/stencil_template_parallel.c:554:   int X = Me % Grid[_x_];
	xor	r13d, r13d	# X
# src/stencil_template_parallel.c:554:   int X = Me % Grid[_x_];
	xor	r8d, r8d	# _63
	jmp	.L331	#
.L318:
# src/stencil_template_parallel.c:535:   uint  first = 1;
	mov	esi, DWORD PTR -88[rbp]	# Grid$0, %sfp
	mov	ecx, 1	# Grid$1,
	jmp	.L328	#
.L390:
# src/stencil_template_parallel.c:662: }
	call	__stack_chk_fail@PLT	#
	.cfi_endproc
.LFE55:
	.size	initialize, .-initialize
	.p2align 4
	.globl	memory_release
	.type	memory_release, @function
memory_release:
.LFB59:
	.cfi_startproc
	endbr64	
	push	rbp	#
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	push	rbx	#
	.cfi_def_cfa_offset 24
	.cfi_offset 3, -24
	mov	rbx, rdi	# buffer_ptr, tmp92
	sub	rsp, 8	#,
	.cfi_def_cfa_offset 32
# src/stencil_template_parallel.c:914:   if ( planes != NULL ) {
	test	rsi, rsi	# planes
	je	.L395	#,
# src/stencil_template_parallel.c:915:       if ( planes[OLD].data != NULL ){
	mov	rdi, QWORD PTR [rsi]	# _1, planes_12(D)->data
	mov	rbp, rsi	# planes, tmp93
# src/stencil_template_parallel.c:915:       if ( planes[OLD].data != NULL ){
	test	rdi, rdi	# _1
	je	.L396	#,
# src/stencil_template_parallel.c:916: 	      free(planes[OLD].data);}
	call	free@PLT	#
.L396:
# src/stencil_template_parallel.c:918:       if ( planes[NEW].data != NULL ){
	mov	rdi, QWORD PTR 16[rbp]	# _2, MEM[(struct plane_t *)planes_12(D) + 16B].data
# src/stencil_template_parallel.c:918:       if ( planes[NEW].data != NULL ){
	test	rdi, rdi	# _2
	je	.L395	#,
# src/stencil_template_parallel.c:919: 	      free(planes[NEW].data);}
	call	free@PLT	#
.L395:
# src/stencil_template_parallel.c:926:     if ( buffer_ptr[SEND][d] != NULL ){
	mov	rdi, QWORD PTR 16[rbx]	# _21, (*buffer_ptr_16(D))[2]
# src/stencil_template_parallel.c:926:     if ( buffer_ptr[SEND][d] != NULL ){
	test	rdi, rdi	# _21
	je	.L398	#,
# src/stencil_template_parallel.c:927:       free(buffer_ptr[SEND][d]);}
	call	free@PLT	#
.L398:
# src/stencil_template_parallel.c:929:     if ( buffer_ptr[RECV][d] != NULL ){
	mov	rdi, QWORD PTR 48[rbx]	# _24, MEM[(double * restrict[4] *)buffer_ptr_16(D) + 32B][2]
# src/stencil_template_parallel.c:929:     if ( buffer_ptr[RECV][d] != NULL ){
	test	rdi, rdi	# _24
	je	.L399	#,
# src/stencil_template_parallel.c:930:       free(buffer_ptr[RECV][d]);}
	call	free@PLT	#
.L399:
# src/stencil_template_parallel.c:926:     if ( buffer_ptr[SEND][d] != NULL ){
	mov	rdi, QWORD PTR 24[rbx]	# _30, (*buffer_ptr_16(D))[3]
# src/stencil_template_parallel.c:926:     if ( buffer_ptr[SEND][d] != NULL ){
	test	rdi, rdi	# _30
	je	.L400	#,
# src/stencil_template_parallel.c:927:       free(buffer_ptr[SEND][d]);}
	call	free@PLT	#
.L400:
# src/stencil_template_parallel.c:929:     if ( buffer_ptr[RECV][d] != NULL ){
	mov	rdi, QWORD PTR 56[rbx]	# _33, MEM[(double * restrict[4] *)buffer_ptr_16(D) + 32B][3]
# src/stencil_template_parallel.c:929:     if ( buffer_ptr[RECV][d] != NULL ){
	test	rdi, rdi	# _33
	je	.L401	#,
# src/stencil_template_parallel.c:930:       free(buffer_ptr[RECV][d]);}
	call	free@PLT	#
.L401:
# src/stencil_template_parallel.c:934: }
	add	rsp, 8	#,
	.cfi_def_cfa_offset 24
	xor	eax, eax	#
	pop	rbx	#
	.cfi_def_cfa_offset 16
	pop	rbp	#
	.cfi_def_cfa_offset 8
	ret	
	.cfi_endproc
.LFE59:
	.size	memory_release, .-memory_release
	.section	.rodata.str1.1
.LC25:
	.string	" [ step %4d ] "
	.section	.rodata.str1.8
	.align 8
.LC26:
	.string	"total injected energy is %g, system energy is %g ( in avg %g per grid point)\n"
	.text
	.p2align 4
	.globl	output_energy_stat
	.type	output_energy_stat, @function
output_energy_stat:
.LFB60:
	.cfi_startproc
	endbr64	
	push	r14	#
	.cfi_def_cfa_offset 16
	.cfi_offset 14, -16
	vmovq	r14, xmm0	# budget, tmp128
	push	r13	#
	.cfi_def_cfa_offset 24
	.cfi_offset 13, -24
	mov	r13d, edi	# step, tmp126
	lea	rdi, get_total_energy._omp_fn.0[rip]	# tmp109,
	push	r12	#
	.cfi_def_cfa_offset 32
	.cfi_offset 12, -32
	mov	r12, rcx	# Comm, tmp130
	xor	ecx, ecx	#
	push	rbp	#
	.cfi_def_cfa_offset 40
	.cfi_offset 6, -40
	mov	ebp, edx	# Me, tmp129
	xor	edx, edx	#
	push	rbx	#
	.cfi_def_cfa_offset 48
	.cfi_offset 3, -48
	mov	rbx, rsi	# plane, tmp127
	sub	rsp, 64	#,
	.cfi_def_cfa_offset 112
# include/stencil_template_parallel.h:269:     const int register xsize = plane->size[_x_];
	vmovq	xmm0, QWORD PTR 8[rsi]	# vect__26.469, MEM <vector(2) unsigned int> [(unsigned int *)plane_14(D) + 8B]
# src/stencil_template_parallel.c:939: {
	mov	rax, QWORD PTR fs:40	# tmp133, MEM[(<address-space-1> long unsigned int *)40B]
	mov	QWORD PTR 56[rsp], rax	# D.12571, tmp133
	xor	eax, eax	# tmp133
# include/stencil_template_parallel.h:273:     double * restrict data = plane->data;
	mov	rax, QWORD PTR [rsi]	# data, plane_14(D)->data
	lea	rsi, 16[rsp]	# tmp108,
# src/stencil_template_parallel.c:940:   double system_energy = 0;
	mov	QWORD PTR [rsp], 0x000000000	# system_energy,
# include/stencil_template_parallel.h:288:     #pragma omp parallel for reduction(+:totenergy) schedule(static)
	vmovq	QWORD PTR 32[rsp], xmm0	# MEM <const vector(2) int> [(int *)&.omp_data_o.21 + 16B], vect__26.469
	mov	QWORD PTR 16[rsp], rax	# .omp_data_o.21.data, data
# include/stencil_template_parallel.h:271:     const int register fsize = xsize+2;
	vmovd	eax, xmm0	# tmp106, vect__26.469
# src/stencil_template_parallel.c:941:   double tot_system_energy = 0;
	mov	QWORD PTR 8[rsp], 0x000000000	# tot_system_energy,
# include/stencil_template_parallel.h:271:     const int register fsize = xsize+2;
	add	eax, 2	# tmp107,
# include/stencil_template_parallel.h:288:     #pragma omp parallel for reduction(+:totenergy) schedule(static)
	mov	QWORD PTR 24[rsp], 0x000000000	# .omp_data_o.21.totenergy,
# include/stencil_template_parallel.h:271:     const int register fsize = xsize+2;
	mov	DWORD PTR 40[rsp], eax	# .omp_data_o.21.fsize, tmp107
	call	GOMP_parallel@PLT	#
# src/stencil_template_parallel.c:944:   MPI_Reduce ( &system_energy, &tot_system_energy, 1, MPI_DOUBLE, MPI_SUM, 0, *Comm );
	lea	rsi, 8[rsp]	# tmp110,
	mov	rdi, rsp	# tmp111,
	xor	r9d, r9d	#
# include/stencil_template_parallel.h:288:     #pragma omp parallel for reduction(+:totenergy) schedule(static)
	vmovsd	xmm0, QWORD PTR 24[rsp]	# totenergy, .omp_data_o.21.totenergy
# src/stencil_template_parallel.c:944:   MPI_Reduce ( &system_energy, &tot_system_energy, 1, MPI_DOUBLE, MPI_SUM, 0, *Comm );
	sub	rsp, 8	#,
	.cfi_def_cfa_offset 120
	mov	edx, 1	#,
	lea	r8, ompi_mpi_op_sum[rip]	#,
	lea	rcx, ompi_mpi_double[rip]	# tmp112,
# include/stencil_template_parallel.h:297:     *energy = (double)totenergy;
	vmovsd	QWORD PTR 8[rsp], xmm0	# system_energy, totenergy
# src/stencil_template_parallel.c:944:   MPI_Reduce ( &system_energy, &tot_system_energy, 1, MPI_DOUBLE, MPI_SUM, 0, *Comm );
	push	QWORD PTR [r12]	# *Comm_16(D)
	.cfi_def_cfa_offset 128
	call	MPI_Reduce@PLT	#
# src/stencil_template_parallel.c:950:   if ( Me == 0 ) {
	test	ebp, ebp	# Me
	pop	rax	#
	.cfi_def_cfa_offset 120
	pop	rdx	#
	.cfi_def_cfa_offset 112
	jne	.L422	#,
# src/stencil_template_parallel.c:951:       if ( step >= 0 ){
	test	r13d, r13d	# step
	jns	.L428	#,
.L423:
# src/stencil_template_parallel.c:959: 	            tot_system_energy / (plane->size[_x_]*plane->size[_y_]) );}
	mov	eax, DWORD PTR 8[rbx]	# plane_14(D)->size[0], plane_14(D)->size[0]
# src/stencil_template_parallel.c:954:       printf( "total injected energy is %g, "
	imul	eax, DWORD PTR 12[rbx]	# tmp118, plane_14(D)->size[1]
	vxorps	xmm2, xmm2, xmm2	# tmp131
# /usr/include/x86_64-linux-gnu/bits/stdio2.h:86:   return __printf_chk (__USE_FORTIFY_LEVEL - 1, __fmt, __va_arg_pack ());
	vmovq	xmm0, r14	#, budget
# src/stencil_template_parallel.c:954:       printf( "total injected energy is %g, "
	vmovsd	xmm1, QWORD PTR 8[rsp]	# tot_system_energy.114_3, tot_system_energy
# /usr/include/x86_64-linux-gnu/bits/stdio2.h:86:   return __printf_chk (__USE_FORTIFY_LEVEL - 1, __fmt, __va_arg_pack ());
	lea	rsi, .LC26[rip]	# tmp123,
	mov	edi, 2	#,
# src/stencil_template_parallel.c:954:       printf( "total injected energy is %g, "
	vcvtsi2sd	xmm2, xmm2, rax	# tmp132, tmp131, tmp118
# /usr/include/x86_64-linux-gnu/bits/stdio2.h:86:   return __printf_chk (__USE_FORTIFY_LEVEL - 1, __fmt, __va_arg_pack ());
	mov	eax, 3	#,
	vdivsd	xmm2, xmm1, xmm2	#, tot_system_energy.114_3, tmp117
	call	__printf_chk@PLT	#
.L422:
# src/stencil_template_parallel.c:962: }
	mov	rax, QWORD PTR 56[rsp]	# tmp134, D.12571
	sub	rax, QWORD PTR fs:40	# tmp134, MEM[(<address-space-1> long unsigned int *)40B]
	jne	.L429	#,
	add	rsp, 64	#,
	.cfi_remember_state
	.cfi_def_cfa_offset 48
	xor	eax, eax	#
	pop	rbx	#
	.cfi_def_cfa_offset 40
	pop	rbp	#
	.cfi_def_cfa_offset 32
	pop	r12	#
	.cfi_def_cfa_offset 24
	pop	r13	#
	.cfi_def_cfa_offset 16
	pop	r14	#
	.cfi_def_cfa_offset 8
	ret	
	.p2align 4,,10
	.p2align 3
.L428:
	.cfi_restore_state
# /usr/include/x86_64-linux-gnu/bits/stdio2.h:86:   return __printf_chk (__USE_FORTIFY_LEVEL - 1, __fmt, __va_arg_pack ());
	mov	edx, r13d	#, step
	lea	rsi, .LC25[rip]	# tmp113,
	mov	edi, 2	#,
	xor	eax, eax	#
	call	__printf_chk@PLT	#
# src/stencil_template_parallel.c:952: 	      printf(" [ step %4d ] ", step ); fflush(stdout);}
	mov	rdi, QWORD PTR stdout[rip]	#, stdout
	call	fflush@PLT	#
	jmp	.L423	#
.L429:
# src/stencil_template_parallel.c:962: }
	call	__stack_chk_fail@PLT	#
	.cfi_endproc
.LFE60:
	.size	output_energy_stat, .-output_energy_stat
	.section	.rodata.str1.8
	.align 8
.LC27:
	.string	"MPI_thread level obtained is %d instead of %d\n"
	.align 8
.LC28:
	.string	"Rank %d signals initialization completed. Beginning computation part.\n"
	.align 8
.LC29:
	.string	"task %d is opting out with termination code %d\n"
	.align 8
.LC30:
	.string	"WEAK rank=0 ranks=%d total =%.6fperstep=%.6f comm=%.6f (pack=%.6f wait=%.6f unpack=%.6f) compute=%.6f inject=%.6f other=%.6f\n"
	.align 8
.LC32:
	.string	"WEAK_PCT comm=%.1f%% compute=%.1f%% inject=%.1f%% other=%.1f%%\n"
	.section	.text.startup,"ax",@progbits
	.p2align 4
	.globl	main
	.type	main, @function
main:
.LFB51:
	.cfi_startproc
	endbr64	
	lea	r10, 8[rsp]	#,
	.cfi_def_cfa 10, 0
	and	rsp, -32	#,
# src/stencil_template_parallel.c:49:   MPI_Init_thread( &argc, &argv, MPI_THREAD_FUNNELED, &level_obtained );
	mov	edx, 1	#,
# src/stencil_template_parallel.c:19: int main(int argc, char **argv) {
	push	QWORD PTR -8[r10]	#
	push	rbp	#
	mov	rbp, rsp	#,
	.cfi_escape 0x10,0x6,0x2,0x76,0
	push	r15	#
	push	r14	#
# src/stencil_template_parallel.c:49:   MPI_Init_thread( &argc, &argv, MPI_THREAD_FUNNELED, &level_obtained );
	lea	rcx, -332[rbp]	# tmp301,
# src/stencil_template_parallel.c:19: int main(int argc, char **argv) {
	push	r13	#
	push	r12	#
	push	r10	#
	.cfi_escape 0xf,0x3,0x76,0x58,0x6
	.cfi_escape 0x10,0xf,0x2,0x76,0x78
	.cfi_escape 0x10,0xe,0x2,0x76,0x70
	.cfi_escape 0x10,0xd,0x2,0x76,0x68
	.cfi_escape 0x10,0xc,0x2,0x76,0x60
	push	rbx	#
	sub	rsp, 608	#,
	.cfi_escape 0x10,0x3,0x2,0x76,0x50
# src/stencil_template_parallel.c:19: int main(int argc, char **argv) {
	mov	DWORD PTR -372[rbp], edi	# argc, argc
# src/stencil_template_parallel.c:49:   MPI_Init_thread( &argc, &argv, MPI_THREAD_FUNNELED, &level_obtained );
	lea	rdi, -372[rbp]	# tmp303,
# src/stencil_template_parallel.c:19: int main(int argc, char **argv) {
	mov	QWORD PTR -384[rbp], rsi	# argv, argv
# src/stencil_template_parallel.c:49:   MPI_Init_thread( &argc, &argv, MPI_THREAD_FUNNELED, &level_obtained );
	lea	rsi, -384[rbp]	# tmp302,
# src/stencil_template_parallel.c:19: int main(int argc, char **argv) {
	mov	rax, QWORD PTR fs:40	# tmp523, MEM[(<address-space-1> long unsigned int *)40B]
	mov	QWORD PTR -56[rbp], rax	# D.12654, tmp523
	xor	eax, eax	# tmp523
# src/stencil_template_parallel.c:49:   MPI_Init_thread( &argc, &argv, MPI_THREAD_FUNNELED, &level_obtained );
	call	MPI_Init_thread@PLT	#
# src/stencil_template_parallel.c:52:   if ( level_obtained < MPI_THREAD_FUNNELED ) {
	mov	edx, DWORD PTR -332[rbp]	# level_obtained.24_1, level_obtained
# src/stencil_template_parallel.c:52:   if ( level_obtained < MPI_THREAD_FUNNELED ) {
	test	edx, edx	# level_obtained.24_1
	jle	.L486	#,
# src/stencil_template_parallel.c:59:   MPI_Comm_rank(MPI_COMM_WORLD, &Rank); // get my ID (of the rank)
	lea	rbx, ompi_mpi_comm_world[rip]	# tmp306,
	lea	rsi, -360[rbp]	# tmp305,
	mov	rdi, rbx	#, tmp306
# src/stencil_template_parallel.c:61:   MPI_Comm_dup(MPI_COMM_WORLD, &myCOMM_WORLD); // make a private communicator
	lea	r15, -328[rbp]	# tmp492,
# src/stencil_template_parallel.c:59:   MPI_Comm_rank(MPI_COMM_WORLD, &Rank); // get my ID (of the rank)
	call	MPI_Comm_rank@PLT	#
# src/stencil_template_parallel.c:60:   MPI_Comm_size(MPI_COMM_WORLD, &Ntasks); // get total number of ranks
	lea	rsi, -356[rbp]	# tmp307,
	mov	rdi, rbx	#, tmp306
	call	MPI_Comm_size@PLT	#
# src/stencil_template_parallel.c:61:   MPI_Comm_dup(MPI_COMM_WORLD, &myCOMM_WORLD); // make a private communicator
	mov	rsi, r15	#, tmp492
	mov	rdi, rbx	#, tmp306
	mov	QWORD PTR -648[rbp], r15	# %sfp, tmp492
	call	MPI_Comm_dup@PLT	#
# src/stencil_template_parallel.c:70:   int ret = initialize( &myCOMM_WORLD, Rank, Ntasks, argc, argv, &S, &N, &periodic, 
	sub	rsp, 8	#,
	lea	rax, -128[rbp]	# tmp494,
	mov	rdi, r15	#, tmp492
	mov	QWORD PTR -400[rbp], rax	# %sfp, tmp494
	mov	r8, QWORD PTR -384[rbp]	#, argv
	lea	r9, -272[rbp]	#,
	mov	ecx, DWORD PTR -372[rbp]	#, argc
	mov	edx, DWORD PTR -356[rbp]	#, Ntasks
	push	rax	# tmp494
	lea	rax, -240[rbp]	# tmp497,
	mov	QWORD PTR -560[rbp], rax	# %sfp, tmp497
	mov	esi, DWORD PTR -360[rbp]	#, Rank
	push	rax	# tmp497
	lea	rax, -312[rbp]	# tmp319,
	push	rax	# tmp319
	lea	rax, -320[rbp]	# tmp320,
	push	rax	# tmp320
	lea	rax, -340[rbp]	# tmp321,
	push	rax	# tmp321
	lea	rax, -344[rbp]	# tmp322,
	push	rax	# tmp322
	lea	rax, -352[rbp]	# tmp323,
	push	rax	# tmp323
	lea	rax, -256[rbp]	# tmp496,
	mov	QWORD PTR -424[rbp], rax	# %sfp, tmp496
	push	rax	# tmp496
	lea	rax, -336[rbp]	# tmp325,
	push	rax	# tmp325
	lea	rax, -348[rbp]	# tmp326,
	push	rax	# tmp326
	lea	rax, -264[rbp]	# tmp495,
	mov	QWORD PTR -416[rbp], rax	# %sfp, tmp495
	push	rax	# tmp495
	call	initialize	#
# /usr/include/x86_64-linux-gnu/bits/stdio2.h:86:   return __printf_chk (__USE_FORTIFY_LEVEL - 1, __fmt, __va_arg_pack ());
	mov	edx, DWORD PTR -360[rbp]	#, Rank
	add	rsp, 96	#,
	lea	rsi, .LC28[rip]	# tmp329,
# src/stencil_template_parallel.c:70:   int ret = initialize( &myCOMM_WORLD, Rank, Ntasks, argc, argv, &S, &N, &periodic, 
	mov	r15d, eax	# current, tmp507
# /usr/include/x86_64-linux-gnu/bits/stdio2.h:86:   return __printf_chk (__USE_FORTIFY_LEVEL - 1, __fmt, __va_arg_pack ());
	mov	edi, 2	#,
	xor	eax, eax	#
	call	__printf_chk@PLT	#
# src/stencil_template_parallel.c:78:   if ( ret ) {
	test	r15d, r15d	# current
	jne	.L487	#,
# src/stencil_template_parallel.c:92:   int M = Niterations - warmup;
	mov	ebx, DWORD PTR -352[rbp]	# Niterations.32_8, Niterations
# src/stencil_template_parallel.c:95:   MPI_Barrier(myCOMM_WORLD);
	mov	rdi, QWORD PTR -328[rbp]	#, myCOMM_WORLD
# src/stencil_template_parallel.c:92:   int M = Niterations - warmup;
	mov	DWORD PTR -568[rbp], ebx	# %sfp, Niterations.32_8
# src/stencil_template_parallel.c:95:   MPI_Barrier(myCOMM_WORLD);
	call	MPI_Barrier@PLT	#
# src/stencil_template_parallel.c:96:   double t_loop0 = MPI_Wtime();
	call	MPI_Wtime@PLT	#
# src/stencil_template_parallel.c:195:   output_energy_stat ( -1, &planes[!current], Niterations * Nsources*energy_per_source, Rank, &myCOMM_WORLD );
	mov	edi, DWORD PTR -344[rbp]	# pretmp_432, Nsources
# src/stencil_template_parallel.c:99:   for (int iter = 0; iter < Niterations; ++iter) { 
	test	ebx, ebx	# Niterations.32_8
# src/stencil_template_parallel.c:195:   output_energy_stat ( -1, &planes[!current], Niterations * Nsources*energy_per_source, Rank, &myCOMM_WORLD );
	vmovsd	xmm8, QWORD PTR -312[rbp]	# pretmp_431, energy_per_source
# src/stencil_template_parallel.c:96:   double t_loop0 = MPI_Wtime();
	vmovsd	QWORD PTR -656[rbp], xmm0	# %sfp, tmp508
# src/stencil_template_parallel.c:195:   output_energy_stat ( -1, &planes[!current], Niterations * Nsources*energy_per_source, Rank, &myCOMM_WORLD );
	mov	DWORD PTR -580[rbp], edi	# %sfp, pretmp_432
# src/stencil_template_parallel.c:99:   for (int iter = 0; iter < Niterations; ++iter) { 
	jle	.L434	#,
# src/stencil_template_parallel.c:110:   inject_energy(periodic, Nsources_local, Sources_local, energy_per_source, &planes[current], N);
	mov	ebx, DWORD PTR -348[rbp]	# pretmp_377, periodic
	lea	rax, ompi_request_null[rip]	# _387,
# src/stencil_template_parallel.c:99:   for (int iter = 0; iter < Niterations; ++iter) { 
	xor	ecx, ecx	#
# src/stencil_template_parallel.c:141:     if ( output_energy_stat_perstep )
	mov	DWORD PTR -564[rbp], edi	# %sfp, pretmp_432
	vmovq	xmm1, rax	# _387, _387
# src/stencil_template_parallel.c:110:   inject_energy(periodic, Nsources_local, Sources_local, energy_per_source, &planes[current], N);
	movsx	rax, DWORD PTR -340[rbp]	#, Nsources_local
# src/stencil_template_parallel.c:99:   for (int iter = 0; iter < Niterations; ++iter) { 
	mov	DWORD PTR -408[rbp], ecx	# %sfp,
# src/stencil_template_parallel.c:110:   inject_energy(periodic, Nsources_local, Sources_local, energy_per_source, &planes[current], N);
	mov	DWORD PTR -392[rbp], ebx	# %sfp, pretmp_377
	vpbroadcastq	ymm4, xmm1	# _383, _387
# src/stencil_template_parallel.c:141:     if ( output_energy_stat_perstep )
	mov	ebx, DWORD PTR -336[rbp]	# pretmp_379, output_energy_stat_perstep
# src/stencil_template_parallel.c:110:   inject_energy(periodic, Nsources_local, Sources_local, energy_per_source, &planes[current], N);
	mov	DWORD PTR -584[rbp], eax	# %sfp, pretmp_375
# src/stencil_template_parallel.c:141:     if ( output_energy_stat_perstep )
	mov	DWORD PTR -588[rbp], ebx	# %sfp, pretmp_379
	mov	rbx, QWORD PTR -320[rbp]	# _308, Sources_local
# src/stencil_template_parallel.c:90:   double t_pack=0.0, t_wait=0.0, t_unpack=0.0, t_inject=0.0, t_compute=0.0;
	mov	QWORD PTR -464[rbp], 0x000000000	# %sfp,
	lea	r13, [rbx+rax*8]	# _269,
	lea	rax, -208[rbp]	# tmp500,
	mov	QWORD PTR -640[rbp], rbx	# %sfp, _308
	mov	QWORD PTR -472[rbp], rax	# %sfp, tmp500
	lea	rax, -304[rbp]	# tmp493,
# src/stencil_template_parallel.c:90:   double t_pack=0.0, t_wait=0.0, t_unpack=0.0, t_inject=0.0, t_compute=0.0;
	mov	QWORD PTR -440[rbp], 0x000000000	# %sfp,
# src/stencil_template_parallel.c:90:   double t_pack=0.0, t_wait=0.0, t_unpack=0.0, t_inject=0.0, t_compute=0.0;
	mov	QWORD PTR -456[rbp], 0x000000000	# %sfp,
# src/stencil_template_parallel.c:90:   double t_pack=0.0, t_wait=0.0, t_unpack=0.0, t_inject=0.0, t_compute=0.0;
	mov	QWORD PTR -432[rbp], 0x000000000	# %sfp,
# src/stencil_template_parallel.c:90:   double t_pack=0.0, t_wait=0.0, t_unpack=0.0, t_inject=0.0, t_compute=0.0;
	mov	QWORD PTR -448[rbp], 0x000000000	# %sfp,
	mov	QWORD PTR -632[rbp], rax	# %sfp, tmp493
	vmovdqa	YMMWORD PTR -624[rbp], ymm4	# %sfp, _383
	vmovsd	QWORD PTR -576[rbp], xmm8	# %sfp, pretmp_431
	.p2align 4,,10
	.p2align 3
.L435:
# src/stencil_template_parallel.c:106:     for (int k = 0; k < 8; ++k) reqs[k] = MPI_REQUEST_NULL;
	vmovdqa	ymm4, YMMWORD PTR -624[rbp]	# _383, %sfp
	vmovdqa	YMMWORD PTR -208[rbp], ymm4	# MEM <vector(4) long unsigned int> [(struct ompi_request_t * *)&reqs], _383
	vmovdqa	YMMWORD PTR -176[rbp], ymm4	# MEM <vector(4) long unsigned int> [(struct ompi_request_t * *)&reqs + 32B], _383
# src/stencil_template_parallel.c:109:     double ti0 = MPI_Wtime();
	vzeroupper
	call	MPI_Wtime@PLT	#
# src/stencil_template_parallel.c:110:   inject_energy(periodic, Nsources_local, Sources_local, energy_per_source, &planes[current], N);
	movsx	r11, r15d	# current, current
	mov	rdi, QWORD PTR -560[rbp]	# tmp497, %sfp
	mov	rax, r11	# tmp444, current
# include/stencil_template_parallel.h:129:     for (int s = 0; s < Nsources; s++) {
	mov	esi, DWORD PTR -584[rbp]	#, %sfp
# src/stencil_template_parallel.c:109:     double ti0 = MPI_Wtime();
	vmovsd	QWORD PTR -480[rbp], xmm0	# %sfp, tmp516
# src/stencil_template_parallel.c:110:   inject_energy(periodic, Nsources_local, Sources_local, energy_per_source, &planes[current], N);
	sal	rax, 4	# tmp444,
# include/stencil_template_parallel.h:122:     const uint register sizex = plane->size[_x_]+2; // interior + halos
	mov	r9d, DWORD PTR -232[rbp+rax]	# _209, MEM <struct plane_t[2]> [(struct plane_t *)&planes][current_106].size[0]
# src/stencil_template_parallel.c:110:   inject_energy(periodic, Nsources_local, Sources_local, energy_per_source, &planes[current], N);
	lea	r14, [rdi+rax]	# _10,
# include/stencil_template_parallel.h:129:     for (int s = 0; s < Nsources; s++) {
	test	esi, esi	#
# include/stencil_template_parallel.h:122:     const uint register sizex = plane->size[_x_]+2; // interior + halos
	lea	rdi, -48[rax]	# tmp649,
# include/stencil_template_parallel.h:123:     double * restrict data = plane->data;
	mov	r12, QWORD PTR -240[rbp+rax]	# data, MEM <struct plane_t[2]> [(struct plane_t *)&planes][current_106].data
# include/stencil_template_parallel.h:122:     const uint register sizex = plane->size[_x_]+2; // interior + halos
	lea	r10, [rdi+rbp]	# tmp447,
# include/stencil_template_parallel.h:122:     const uint register sizex = plane->size[_x_]+2; // interior + halos
	lea	ebx, 2[r9]	# sizex,
	vmovd	xmm7, ebx	# sizex, sizex
	vpinsrd	xmm2, xmm7, r9d, 1	# _397, sizex, _209
# include/stencil_template_parallel.h:129:     for (int s = 0; s < Nsources; s++) {
	jle	.L463	#,
	mov	edx, DWORD PTR -392[rbp]	#, %sfp
# include/stencil_template_parallel.h:140:                 if ( (N[_x_] == 1)  ) {
	mov	eax, DWORD PTR -264[rbp]	# _225, MEM[(const uint *)&N]
# include/stencil_template_parallel.h:157:                 if ( (N[_y_] == 1) )
	mov	edi, DWORD PTR -260[rbp]	# _238, MEM[(const uint *)&N + 4B]
	test	edx, edx	#
	jne	.L436	#,
	mov	rax, QWORD PTR -640[rbp]	# ivtmp.509, %sfp
	vmovsd	xmm8, QWORD PTR -576[rbp]	# pretmp_431, %sfp
	.p2align 4,,10
	.p2align 3
.L437:
# include/stencil_template_parallel.h:136:             data[ IDX(x,y) ] += energy;
	mov	edx, DWORD PTR 4[rax]	# tmp336, MEM[(unsigned int *)_361 + 4B]
# include/stencil_template_parallel.h:129:     for (int s = 0; s < Nsources; s++) {
	add	rax, 8	# ivtmp.509,
# include/stencil_template_parallel.h:136:             data[ IDX(x,y) ] += energy;
	imul	edx, ebx	# tmp336, sizex
	add	edx, DWORD PTR -8[rax]	# tmp339, MEM[(unsigned int *)_361]
# include/stencil_template_parallel.h:129:     for (int s = 0; s < Nsources; s++) {
	cmp	r13, rax	# _269, ivtmp.509
# include/stencil_template_parallel.h:136:             data[ IDX(x,y) ] += energy;
	lea	rdx, [r12+rdx*8]	# _404,
# include/stencil_template_parallel.h:136:             data[ IDX(x,y) ] += energy;
	vaddsd	xmm0, xmm8, QWORD PTR [rdx]	# tmp341, pretmp_431, *_404
	vmovsd	QWORD PTR [rdx], xmm0	# *_404, tmp341
# include/stencil_template_parallel.h:129:     for (int s = 0; s < Nsources; s++) {
	jne	.L437	#,
.L463:
	mov	DWORD PTR -376[rbp], r9d	# %sfp, _209
	mov	QWORD PTR -536[rbp], r11	# %sfp, current
	vmovq	QWORD PTR -552[rbp], xmm2	# %sfp, _397
# src/stencil_template_parallel.c:111:   double ti1 = MPI_Wtime();
	call	MPI_Wtime@PLT	#
	vmovsd	QWORD PTR -488[rbp], xmm0	# %sfp, tmp509
# src/stencil_template_parallel.c:114:   double tA = MPI_Wtime();
	call	MPI_Wtime@PLT	#
# src/stencil_template_parallel.c:115:   fill_buffers(buffers, &planes[current], periodic, N);
	mov	rcx, QWORD PTR -416[rbp]	#, %sfp
	mov	edx, DWORD PTR -392[rbp]	#, %sfp
	mov	rsi, r14	#, _10
	mov	rdi, QWORD PTR -400[rbp]	#, %sfp
# src/stencil_template_parallel.c:114:   double tA = MPI_Wtime();
	vmovsd	QWORD PTR -496[rbp], xmm0	# %sfp, tmp510
# src/stencil_template_parallel.c:115:   fill_buffers(buffers, &planes[current], periodic, N);
	call	fill_buffers	#
# src/stencil_template_parallel.c:116:   post_MPI_reqs(reqs, buffers, &planes[current], neighbours, myCOMM_WORLD);
	mov	r8, QWORD PTR -328[rbp]	#, myCOMM_WORLD
	mov	rdx, r14	#, _10
	mov	rcx, QWORD PTR -424[rbp]	#, %sfp
	mov	rsi, QWORD PTR -400[rbp]	#, %sfp
	mov	rdi, QWORD PTR -472[rbp]	#, %sfp
	vzeroupper
	call	post_MPI_reqs	#
# src/stencil_template_parallel.c:117:   double tB = MPI_Wtime();
	call	MPI_Wtime@PLT	#
# src/stencil_template_parallel.c:120:   MPI_Waitall(8, reqs, MPI_STATUSES_IGNORE);
	mov	rsi, QWORD PTR -472[rbp]	#, %sfp
	xor	edx, edx	#
	mov	edi, 8	#,
# src/stencil_template_parallel.c:117:   double tB = MPI_Wtime();
	vmovsd	QWORD PTR -504[rbp], xmm0	# %sfp, tmp511
# src/stencil_template_parallel.c:120:   MPI_Waitall(8, reqs, MPI_STATUSES_IGNORE);
	call	MPI_Waitall@PLT	#
# src/stencil_template_parallel.c:121:   double tC = MPI_Wtime();
	call	MPI_Wtime@PLT	#
# src/stencil_template_parallel.c:124:   copy_halos(buffers, &planes[current], neighbours, periodic, N);
	mov	r8, QWORD PTR -416[rbp]	#, %sfp
	mov	rsi, r14	#, _10
	mov	ecx, DWORD PTR -392[rbp]	#, %sfp
	mov	rdx, QWORD PTR -424[rbp]	#, %sfp
	mov	rdi, QWORD PTR -400[rbp]	#, %sfp
# src/stencil_template_parallel.c:121:   double tC = MPI_Wtime();
	vmovsd	QWORD PTR -512[rbp], xmm0	# %sfp, tmp512
# src/stencil_template_parallel.c:124:   copy_halos(buffers, &planes[current], neighbours, periodic, N);
	call	copy_halos	#
# src/stencil_template_parallel.c:125:   double tD = MPI_Wtime();
	vzeroupper
	call	MPI_Wtime@PLT	#
	xor	r15d, 1	# _34,
# include/stencil_template_parallel.h:208:     double * restrict new = newplane->data; // and write on this
	movsx	r14, r15d	# _34, _34
# src/stencil_template_parallel.c:125:   double tD = MPI_Wtime();
	vmovsd	QWORD PTR -520[rbp], xmm0	# %sfp, tmp513
# src/stencil_template_parallel.c:128:   double tc0 = MPI_Wtime();
	call	MPI_Wtime@PLT	#
# include/stencil_template_parallel.h:191:     uint register fysize = oldplane->size[_y_]+2;
	mov	r10, QWORD PTR -536[rbp]	# current, %sfp
	xor	ecx, ecx	#
	xor	edx, edx	#
# include/stencil_template_parallel.h:208:     double * restrict new = newplane->data; // and write on this
	mov	rax, r14	# tmp408, _34
# include/stencil_template_parallel.h:210:     #pragma omp parallel for schedule(static)
	vmovq	xmm2, QWORD PTR -552[rbp]	# _397, %sfp
	vmovq	xmm5, r12	# data, data
	mov	rsi, QWORD PTR -632[rbp]	#, %sfp
# include/stencil_template_parallel.h:208:     double * restrict new = newplane->data; // and write on this
	sal	rax, 4	# tmp408,
# include/stencil_template_parallel.h:191:     uint register fysize = oldplane->size[_y_]+2;
	sal	r10, 4	# current,
# src/stencil_template_parallel.c:128:   double tc0 = MPI_Wtime();
	vmovsd	QWORD PTR -528[rbp], xmm0	# %sfp, tmp514
	lea	rdi, update_plane._omp_fn.0[rip]	#,
# include/stencil_template_parallel.h:208:     double * restrict new = newplane->data; // and write on this
	mov	rax, QWORD PTR -240[rbp+rax]	# new, MEM <struct plane_t[2]> [(struct plane_t *)&planes][_34].data
# include/stencil_template_parallel.h:191:     uint register fysize = oldplane->size[_y_]+2;
	mov	r8d, DWORD PTR -228[rbp+r10]	# _162, MEM <struct plane_t[2]> [(const struct plane_t *)&planes][current_106].size[1]
# include/stencil_template_parallel.h:210:     #pragma omp parallel for schedule(static)
	vmovq	QWORD PTR -288[rbp], xmm2	# MEM <vector(2) unsigned int> [(unsigned int *)_201], _397
	vpinsrq	xmm0, xmm5, rax, 1	# tmp411, data, new
	mov	QWORD PTR -544[rbp], rax	# %sfp, new
	mov	DWORD PTR -280[rbp], r8d	# MEM[(struct .omp_data_s.14 *)_201].ysize, _162
	mov	DWORD PTR -536[rbp], r8d	# %sfp, _162
	vmovdqa	XMMWORD PTR -304[rbp], xmm0	# MEM <vector(2) long unsigned int> [(double * *)_201], tmp411
	call	GOMP_parallel@PLT	#
# include/stencil_template_parallel.h:230:     if ( periodic ) {
	mov	r8d, DWORD PTR -392[rbp]	#, %sfp
	mov	rax, QWORD PTR -544[rbp]	# new, %sfp
	mov	r9d, DWORD PTR -376[rbp]	# _209, %sfp
	test	r8d, r8d	#
	mov	r8d, DWORD PTR -536[rbp]	# _162, %sfp
	je	.L456	#,
# include/stencil_template_parallel.h:234:         if ( N[_x_] == 1 ) {
	mov	ecx, DWORD PTR -264[rbp]	# j, MEM[(const uint *)&N]
# include/stencil_template_parallel.h:234:         if ( N[_x_] == 1 ) {
	cmp	ecx, 1	# j,
	je	.L454	#,
.L457:
# include/stencil_template_parallel.h:243:         if ( N[_y_] == 1 ) {
	mov	edx, DWORD PTR -260[rbp]	# i, MEM[(const uint *)&N + 4B]
# include/stencil_template_parallel.h:243:         if ( N[_y_] == 1 ) {
	cmp	edx, 1	# i,
	je	.L488	#,
.L456:
# src/stencil_template_parallel.c:130:   double tc1 = MPI_Wtime();
	call	MPI_Wtime@PLT	#
# src/stencil_template_parallel.c:133:     t_inject += (ti1 - ti0);
	vmovsd	xmm3, QWORD PTR -488[rbp]	# ti1, %sfp
	vsubsd	xmm1, xmm3, QWORD PTR -480[rbp]	# tmp414, ti1, %sfp
# src/stencil_template_parallel.c:133:     t_inject += (ti1 - ti0);
	vaddsd	xmm3, xmm1, QWORD PTR -440[rbp]	# t_inject, tmp414, %sfp
# src/stencil_template_parallel.c:134:     t_pack   += (tB  - tA);
	vmovsd	xmm6, QWORD PTR -504[rbp]	# tB, %sfp
# src/stencil_template_parallel.c:135:     t_wait   += (tC  - tB);
	vmovsd	xmm7, QWORD PTR -512[rbp]	# tC, %sfp
# src/stencil_template_parallel.c:136:     t_unpack += (tD  - tC);
	vmovsd	xmm2, QWORD PTR -520[rbp]	# tD, %sfp
# src/stencil_template_parallel.c:134:     t_pack   += (tB  - tA);
	vsubsd	xmm1, xmm6, QWORD PTR -496[rbp]	# tmp415, tB, %sfp
# src/stencil_template_parallel.c:137:     t_compute+= (tc1 - tc0);
	vsubsd	xmm0, xmm0, QWORD PTR -528[rbp]	# tmp418, tmp515, %sfp
# src/stencil_template_parallel.c:134:     t_pack   += (tB  - tA);
	vaddsd	xmm5, xmm1, QWORD PTR -448[rbp]	# t_pack, tmp415, %sfp
# src/stencil_template_parallel.c:141:     if ( output_energy_stat_perstep )
	mov	edi, DWORD PTR -588[rbp]	#, %sfp
# src/stencil_template_parallel.c:135:     t_wait   += (tC  - tB);
	vsubsd	xmm1, xmm7, xmm6	# tmp416, tC, tB
# src/stencil_template_parallel.c:133:     t_inject += (ti1 - ti0);
	vmovsd	QWORD PTR -440[rbp], xmm3	# %sfp, t_inject
# src/stencil_template_parallel.c:135:     t_wait   += (tC  - tB);
	vaddsd	xmm6, xmm1, QWORD PTR -432[rbp]	# t_wait, tmp416, %sfp
# src/stencil_template_parallel.c:136:     t_unpack += (tD  - tC);
	vsubsd	xmm1, xmm2, xmm7	# tmp417, tD, tC
# src/stencil_template_parallel.c:142:       output_energy_stat( iter, &planes[!current], (iter+1) * Nsources*energy_per_source, Rank, &myCOMM_WORLD);
	mov	eax, DWORD PTR -408[rbp]	# iter, %sfp
# src/stencil_template_parallel.c:136:     t_unpack += (tD  - tC);
	vaddsd	xmm7, xmm1, QWORD PTR -456[rbp]	# t_unpack, tmp417, %sfp
# src/stencil_template_parallel.c:137:     t_compute+= (tc1 - tc0);
	vaddsd	xmm3, xmm0, QWORD PTR -464[rbp]	# t_compute, tmp418, %sfp
# src/stencil_template_parallel.c:141:     if ( output_energy_stat_perstep )
	test	edi, edi	#
# src/stencil_template_parallel.c:134:     t_pack   += (tB  - tA);
	vmovsd	QWORD PTR -448[rbp], xmm5	# %sfp, t_pack
# src/stencil_template_parallel.c:142:       output_energy_stat( iter, &planes[!current], (iter+1) * Nsources*energy_per_source, Rank, &myCOMM_WORLD);
	lea	ebx, 1[rax]	# _429,
# src/stencil_template_parallel.c:135:     t_wait   += (tC  - tB);
	vmovsd	QWORD PTR -432[rbp], xmm6	# %sfp, t_wait
# src/stencil_template_parallel.c:136:     t_unpack += (tD  - tC);
	vmovsd	QWORD PTR -456[rbp], xmm7	# %sfp, t_unpack
# src/stencil_template_parallel.c:137:     t_compute+= (tc1 - tc0);
	vmovsd	QWORD PTR -464[rbp], xmm3	# %sfp, t_compute
# src/stencil_template_parallel.c:141:     if ( output_energy_stat_perstep )
	jne	.L489	#,
# src/stencil_template_parallel.c:99:   for (int iter = 0; iter < Niterations; ++iter) { 
	mov	edi, DWORD PTR -580[rbp]	# pretmp_432, %sfp
	add	DWORD PTR -564[rbp], edi	# %sfp, pretmp_432
	cmp	DWORD PTR -568[rbp], ebx	# %sfp, _429
	je	.L484	#,
.L460:
# include/stencil_template_parallel.h:237:             for (uint j = 1; j <= ysize; j++) {
	mov	DWORD PTR -408[rbp], ebx	# %sfp, _429
	jmp	.L435	#
	.p2align 4,,10
	.p2align 3
.L436:
	cmp	eax, 1	# _225,
# include/stencil_template_parallel.h:157:                 if ( (N[_y_] == 1) )
	mov	rax, QWORD PTR -640[rbp]	# ivtmp.516, %sfp
	je	.L469	#,
	mov	DWORD PTR -488[rbp], r9d	# %sfp, _209
	vmovsd	xmm8, QWORD PTR -576[rbp]	# pretmp_431, %sfp
	jmp	.L444	#
	.p2align 4,,10
	.p2align 3
.L442:
# include/stencil_template_parallel.h:129:     for (int s = 0; s < Nsources; s++) {
	add	rax, 8	# ivtmp.516,
	cmp	r13, rax	# _269, ivtmp.516
	je	.L490	#,
.L444:
# include/stencil_template_parallel.h:133:             int y = Sources[s][_y_];
	mov	edx, DWORD PTR 4[rax]	# _345, MEM[(unsigned int *)_318 + 4B]
# include/stencil_template_parallel.h:136:             data[ IDX(x,y) ] += energy;
	mov	esi, ebx	# _344, sizex
# include/stencil_template_parallel.h:132:             int x = Sources[s][_x_];
	mov	ecx, DWORD PTR [rax]	#, MEM[(unsigned int *)_318]
# include/stencil_template_parallel.h:136:             data[ IDX(x,y) ] += energy;
	imul	esi, edx	# _344, _345
# include/stencil_template_parallel.h:157:                 if ( (N[_y_] == 1) )
	cmp	edi, 1	# _238,
# include/stencil_template_parallel.h:136:             data[ IDX(x,y) ] += energy;
	lea	r8d, [rsi+rcx]	# tmp344,
	lea	r8, [r12+r8*8]	# _340,
# include/stencil_template_parallel.h:136:             data[ IDX(x,y) ] += energy;
	vaddsd	xmm0, xmm8, QWORD PTR [r8]	# tmp346, pretmp_431, *_340
	vmovsd	QWORD PTR [r8], xmm0	# *_340, tmp346
# include/stencil_template_parallel.h:157:                 if ( (N[_y_] == 1) )
	jne	.L442	#,
# include/stencil_template_parallel.h:167:                     if ( y == 1 ){ // source on easr edge case
	cmp	edx, 1	# _345,
# include/stencil_template_parallel.h:168:                         data[IDX(x, plane->size[_y_]+1)] += energy;
	mov	r8d, DWORD PTR -180[r10]	# pretmp_335, MEM <struct plane_t[2]> [(struct plane_t *)&planes][current_106].size[1]
# include/stencil_template_parallel.h:167:                     if ( y == 1 ){ // source on easr edge case
	jne	.L443	#,
# include/stencil_template_parallel.h:168:                         data[IDX(x, plane->size[_y_]+1)] += energy;
	lea	r9d, 1[r8]	# tmp356,
	imul	esi, r9d	# tmp357, tmp356
	add	esi, ecx	# tmp359, _346
	mov	esi, esi	# tmp359, tmp359
	lea	rsi, [r12+rsi*8]	# _322,
# include/stencil_template_parallel.h:168:                         data[IDX(x, plane->size[_y_]+1)] += energy;
	vaddsd	xmm0, xmm8, QWORD PTR [rsi]	# tmp361, pretmp_431, *_322
	vmovsd	QWORD PTR [rsi], xmm0	# *_322, tmp361
.L443:
# include/stencil_template_parallel.h:171:                     if ( y == plane->size[_y_] ) { // source on west edge case
	cmp	r8d, edx	# pretmp_335, _345
	jne	.L442	#,
# include/stencil_template_parallel.h:172:                         data[IDX(x, 0)] += energy;
	lea	rdx, [r12+rcx*8]	# _331,
# include/stencil_template_parallel.h:129:     for (int s = 0; s < Nsources; s++) {
	add	rax, 8	# ivtmp.516,
# include/stencil_template_parallel.h:172:                         data[IDX(x, 0)] += energy;
	vaddsd	xmm0, xmm8, QWORD PTR [rdx]	# tmp354, pretmp_431, *_331
# include/stencil_template_parallel.h:129:     for (int s = 0; s < Nsources; s++) {
	cmp	r13, rax	# _269, ivtmp.516
# include/stencil_template_parallel.h:172:                         data[IDX(x, 0)] += energy;
	vmovsd	QWORD PTR [rdx], xmm0	# *_331, tmp354
# include/stencil_template_parallel.h:129:     for (int s = 0; s < Nsources; s++) {
	jne	.L444	#,
	.p2align 4,,10
	.p2align 3
.L490:
	mov	r9d, DWORD PTR -488[rbp]	# _209, %sfp
	jmp	.L463	#
	.p2align 4,,10
	.p2align 3
.L489:
# src/stencil_template_parallel.c:142:       output_energy_stat( iter, &planes[!current], (iter+1) * Nsources*energy_per_source, Rank, &myCOMM_WORLD);
	vxorpd	xmm4, xmm4, xmm4	# tmp635
	mov	edi, eax	#, iter
# src/stencil_template_parallel.c:129:   update_plane(periodic, N, &planes[current], &planes[!current]);
	sal	r14, 4	# tmp440,
# src/stencil_template_parallel.c:142:       output_energy_stat( iter, &planes[!current], (iter+1) * Nsources*energy_per_source, Rank, &myCOMM_WORLD);
	mov	edx, DWORD PTR -360[rbp]	#, Rank
	mov	r12d, DWORD PTR -564[rbp]	# ivtmp.530, %sfp
# src/stencil_template_parallel.c:129:   update_plane(periodic, N, &planes[current], &planes[!current]);
	mov	rax, QWORD PTR -560[rbp]	# tmp497, %sfp
# src/stencil_template_parallel.c:142:       output_energy_stat( iter, &planes[!current], (iter+1) * Nsources*energy_per_source, Rank, &myCOMM_WORLD);
	mov	rcx, QWORD PTR -648[rbp]	#, %sfp
	vcvtsi2sd	xmm0, xmm4, r12d	# tmp520, tmp635, ivtmp.530
	vmulsd	xmm0, xmm0, QWORD PTR -576[rbp]	# tmp437, tmp436, %sfp
# src/stencil_template_parallel.c:129:   update_plane(periodic, N, &planes[current], &planes[!current]);
	lea	rsi, [rax+r14]	# tmp441,
# src/stencil_template_parallel.c:142:       output_energy_stat( iter, &planes[!current], (iter+1) * Nsources*energy_per_source, Rank, &myCOMM_WORLD);
	call	output_energy_stat	#
# src/stencil_template_parallel.c:99:   for (int iter = 0; iter < Niterations; ++iter) { 
	mov	eax, DWORD PTR -580[rbp]	# pretmp_432, %sfp
	add	r12d, eax	# ivtmp.530, pretmp_432
	cmp	DWORD PTR -568[rbp], ebx	# %sfp, _429
	mov	DWORD PTR -564[rbp], r12d	# %sfp, ivtmp.530
	jne	.L460	#,
.L484:
	vmovsd	xmm8, QWORD PTR -576[rbp]	# pretmp_431, %sfp
# src/stencil_template_parallel.c:162:   MPI_Barrier(myCOMM_WORLD);                      // end of measured region
	mov	rdi, QWORD PTR -328[rbp]	#, myCOMM_WORLD
	vmovsd	QWORD PTR -392[rbp], xmm8	# %sfp, pretmp_431
	call	MPI_Barrier@PLT	#
# src/stencil_template_parallel.c:164:   double my_total = MPI_Wtime() - t_loop0;
	call	MPI_Wtime@PLT	#
# src/stencil_template_parallel.c:165:   double denom    = (M > 0 ? (double)M : 1.0);
	vxorpd	xmm2, xmm2, xmm2	# tmp653
# src/stencil_template_parallel.c:164:   double my_total = MPI_Wtime() - t_loop0;
	vsubsd	xmm9, xmm0, QWORD PTR -656[rbp]	# my_total, tmp517, %sfp
	vmovsd	xmm8, QWORD PTR -392[rbp]	# pretmp_431, %sfp
# src/stencil_template_parallel.c:165:   double denom    = (M > 0 ? (double)M : 1.0);
	vcvtsi2sd	xmm0, xmm2, DWORD PTR -568[rbp]	# tmp521, tmp653, %sfp
.L467:
# src/stencil_template_parallel.c:176:   if (Rank == 0) {
	cmp	DWORD PTR -360[rbp], 0	# Rank,
	je	.L491	#,
.L464:
# src/stencil_template_parallel.c:195:   output_energy_stat ( -1, &planes[!current], Niterations * Nsources*energy_per_source, Rank, &myCOMM_WORLD );
	mov	edi, DWORD PTR -580[rbp]	# pretmp_432, %sfp
	mov	eax, DWORD PTR -568[rbp]	# Niterations.32_8, %sfp
# src/stencil_template_parallel.c:195:   output_energy_stat ( -1, &planes[!current], Niterations * Nsources*energy_per_source, Rank, &myCOMM_WORLD );
	vxorpd	xmm4, xmm4, xmm4	# tmp673
	mov	esi, r15d	# current, current
	xor	esi, 1	# current,
	mov	rbx, QWORD PTR -560[rbp]	# tmp497, %sfp
	mov	rcx, QWORD PTR -648[rbp]	#, %sfp
# src/stencil_template_parallel.c:195:   output_energy_stat ( -1, &planes[!current], Niterations * Nsources*energy_per_source, Rank, &myCOMM_WORLD );
	imul	eax, edi	# Niterations.32_8, pretmp_432
# src/stencil_template_parallel.c:195:   output_energy_stat ( -1, &planes[!current], Niterations * Nsources*energy_per_source, Rank, &myCOMM_WORLD );
	movsx	rsi, esi	# tmp483, tmp482
	mov	edx, DWORD PTR -360[rbp]	#, Rank
	or	edi, -1	#,
	sal	rsi, 4	# tmp484,
	add	rsi, rbx	# tmp485, tmp497
	vcvtsi2sd	xmm0, xmm4, eax	# tmp522, tmp673, tmp478
	vmulsd	xmm0, xmm0, xmm8	# tmp480, tmp479, pretmp_431
	call	output_energy_stat	#
# src/stencil_template_parallel.c:198:   memory_release( buffers, planes );
	mov	rdi, QWORD PTR -400[rbp]	#, %sfp
	mov	rsi, rbx	#, tmp497
	call	memory_release	#
# src/stencil_template_parallel.c:200:   MPI_Finalize();
	call	MPI_Finalize@PLT	#
.L433:
# src/stencil_template_parallel.c:202: }
	mov	rax, QWORD PTR -56[rbp]	# tmp524, D.12654
	sub	rax, QWORD PTR fs:40	# tmp524, MEM[(<address-space-1> long unsigned int *)40B]
	jne	.L492	#,
	lea	rsp, -48[rbp]	#,
	xor	eax, eax	#
	pop	rbx	#
	pop	r10	#
	.cfi_remember_state
	.cfi_def_cfa 10, 0
	pop	r12	#
	pop	r13	#
	pop	r14	#
	pop	r15	#
	pop	rbp	#
	lea	rsp, -8[r10]	#,
	.cfi_def_cfa 7, 8
	ret	
	.p2align 4,,10
	.p2align 3
.L488:
	.cfi_restore_state
# include/stencil_template_parallel.h:246:             for (uint i = 1; i <= xsize; i++) {
	test	r9d, r9d	# _209
	je	.L456	#,
# include/stencil_template_parallel.h:247:                 new[ IDX(i, 0) ] = new[ IDX(i, ysize) ]; // propagate south edge on the north halo
	imul	r8d, ebx	# _187, sizex
# include/stencil_template_parallel.h:248:                 new[ IDX(i, ysize+1) ] = new[ IDX(i,1) ]; // propagate north edge on the south halo
	lea	esi, [r8+rbx]	# _202,
	.p2align 4,,10
	.p2align 3
.L459:
# include/stencil_template_parallel.h:247:                 new[ IDX(i, 0) ] = new[ IDX(i, ysize) ]; // propagate south edge on the north halo
	lea	ecx, [r8+rdx]	# tmp428,
# include/stencil_template_parallel.h:247:                 new[ IDX(i, 0) ] = new[ IDX(i, ysize) ]; // propagate south edge on the north halo
	vmovsd	xmm0, QWORD PTR [rax+rcx*8]	# _196, *_192
# include/stencil_template_parallel.h:247:                 new[ IDX(i, 0) ] = new[ IDX(i, ysize) ]; // propagate south edge on the north halo
	mov	ecx, edx	# i, i
# include/stencil_template_parallel.h:247:                 new[ IDX(i, 0) ] = new[ IDX(i, ysize) ]; // propagate south edge on the north halo
	vmovsd	QWORD PTR [rax+rcx*8], xmm0	# *_195, _196
# include/stencil_template_parallel.h:248:                 new[ IDX(i, ysize+1) ] = new[ IDX(i,1) ]; // propagate north edge on the south halo
	lea	ecx, [rbx+rdx]	# tmp431,
# include/stencil_template_parallel.h:248:                 new[ IDX(i, ysize+1) ] = new[ IDX(i,1) ]; // propagate north edge on the south halo
	vmovsd	xmm0, QWORD PTR [rax+rcx*8]	# _207, *_200
# include/stencil_template_parallel.h:248:                 new[ IDX(i, ysize+1) ] = new[ IDX(i,1) ]; // propagate north edge on the south halo
	lea	ecx, [rsi+rdx]	# tmp433,
# include/stencil_template_parallel.h:246:             for (uint i = 1; i <= xsize; i++) {
	add	edx, 1	# i,
# include/stencil_template_parallel.h:246:             for (uint i = 1; i <= xsize; i++) {
	cmp	r9d, edx	# _209, i
# include/stencil_template_parallel.h:248:                 new[ IDX(i, ysize+1) ] = new[ IDX(i,1) ]; // propagate north edge on the south halo
	vmovsd	QWORD PTR [rax+rcx*8], xmm0	# *_206, _207
# include/stencil_template_parallel.h:246:             for (uint i = 1; i <= xsize; i++) {
	jnb	.L459	#,
	jmp	.L456	#
	.p2align 4,,10
	.p2align 3
.L454:
# include/stencil_template_parallel.h:237:             for (uint j = 1; j <= ysize; j++) {
	test	r8d, r8d	# _162
	je	.L457	#,
	mov	edx, ebx	# ivtmp.502, sizex
	lea	edi, 1[r9]	# tmp499,
	.p2align 4,,10
	.p2align 3
.L458:
# include/stencil_template_parallel.h:238:                 new[ IDX(0,j) ] = new[ IDX(xsize,j) ]; // we propagate east edge on west halo
	lea	esi, [r9+rdx]	# tmp420,
# include/stencil_template_parallel.h:237:             for (uint j = 1; j <= ysize; j++) {
	add	ecx, 1	# j,
# include/stencil_template_parallel.h:238:                 new[ IDX(0,j) ] = new[ IDX(xsize,j) ]; // we propagate east edge on west halo
	vmovsd	xmm0, QWORD PTR [rax+rsi*8]	# _175, *_171
# include/stencil_template_parallel.h:238:                 new[ IDX(0,j) ] = new[ IDX(xsize,j) ]; // we propagate east edge on west halo
	mov	esi, edx	# ivtmp.502, ivtmp.502
# include/stencil_template_parallel.h:238:                 new[ IDX(0,j) ] = new[ IDX(xsize,j) ]; // we propagate east edge on west halo
	vmovsd	QWORD PTR [rax+rsi*8], xmm0	# *_174, _175
# include/stencil_template_parallel.h:239:                 new[ IDX(xsize+1,j) ] = new[ IDX(1,j) ]; // we propagate west edge on east halo
	lea	esi, 1[rdx]	# tmp423,
# include/stencil_template_parallel.h:239:                 new[ IDX(xsize+1,j) ] = new[ IDX(1,j) ]; // we propagate west edge on east halo
	vmovsd	xmm0, QWORD PTR [rax+rsi*8]	# _184, *_179
# include/stencil_template_parallel.h:239:                 new[ IDX(xsize+1,j) ] = new[ IDX(1,j) ]; // we propagate west edge on east halo
	lea	esi, [rdi+rdx]	# tmp426,
# include/stencil_template_parallel.h:237:             for (uint j = 1; j <= ysize; j++) {
	add	edx, ebx	# ivtmp.502, sizex
	cmp	r8d, ecx	# _162, j
# include/stencil_template_parallel.h:239:                 new[ IDX(xsize+1,j) ] = new[ IDX(1,j) ]; // we propagate west edge on east halo
	vmovsd	QWORD PTR [rax+rsi*8], xmm0	# *_183, _184
# include/stencil_template_parallel.h:237:             for (uint j = 1; j <= ysize; j++) {
	jnb	.L458	#,
	jmp	.L457	#
	.p2align 4,,10
	.p2align 3
.L469:
# include/stencil_template_parallel.h:157:                 if ( (N[_y_] == 1) )
	mov	QWORD PTR -488[rbp], r11	# %sfp, current
	vmovsd	xmm8, QWORD PTR -576[rbp]	# pretmp_431, %sfp
	.p2align 4,,10
	.p2align 3
.L439:
# include/stencil_template_parallel.h:133:             int y = Sources[s][_y_];
	mov	ecx, DWORD PTR 4[rax]	# _217, MEM[(unsigned int *)_290 + 4B]
# include/stencil_template_parallel.h:136:             data[ IDX(x,y) ] += energy;
	mov	esi, ebx	# _218, sizex
# include/stencil_template_parallel.h:132:             int x = Sources[s][_x_];
	mov	edx, DWORD PTR [rax]	#, MEM[(unsigned int *)_290]
# include/stencil_template_parallel.h:136:             data[ IDX(x,y) ] += energy;
	imul	esi, ecx	# _218, _217
# include/stencil_template_parallel.h:150:                     if ( x == 1 ){ // source on easr edge case
	cmp	edx, 1	# _216,
# include/stencil_template_parallel.h:136:             data[ IDX(x,y) ] += energy;
	lea	r8d, [rdx+rsi]	# tmp364,
	lea	r8, [r12+r8*8]	# _222,
# include/stencil_template_parallel.h:136:             data[ IDX(x,y) ] += energy;
	vaddsd	xmm0, xmm8, QWORD PTR [r8]	# tmp366, pretmp_431, *_222
	vmovsd	QWORD PTR [r8], xmm0	# *_222, tmp366
# include/stencil_template_parallel.h:150:                     if ( x == 1 ){ // source on easr edge case
	je	.L493	#,
.L446:
# include/stencil_template_parallel.h:153:                     if ( x == plane->size[_x_] ) { // source on west edge case
	cmp	r9d, edx	# _209, _216
	jne	.L448	#,
# include/stencil_template_parallel.h:154:                         data[IDX(0, y)] += energy;} // propagate on east halo
	mov	r8d, esi	# _218, _218
	lea	r8, [r12+r8*8]	# _235,
# include/stencil_template_parallel.h:154:                         data[IDX(0, y)] += energy;} // propagate on east halo
	vaddsd	xmm0, xmm8, QWORD PTR [r8]	# tmp376, pretmp_431, *_235
	vmovsd	QWORD PTR [r8], xmm0	# *_235, tmp376
.L448:
# include/stencil_template_parallel.h:157:                 if ( (N[_y_] == 1) )
	cmp	edi, 1	# _238,
	je	.L494	#,
.L447:
# include/stencil_template_parallel.h:129:     for (int s = 0; s < Nsources; s++) {
	add	rax, 8	# ivtmp.523,
	cmp	r13, rax	# _269, ivtmp.523
	jne	.L439	#,
	mov	r11, QWORD PTR -488[rbp]	# current, %sfp
	jmp	.L463	#
	.p2align 4,,10
	.p2align 3
.L494:
# include/stencil_template_parallel.h:167:                     if ( y == 1 ){ // source on easr edge case
	cmp	ecx, 1	# _217,
# include/stencil_template_parallel.h:168:                         data[IDX(x, plane->size[_y_]+1)] += energy;
	mov	r8d, DWORD PTR -180[r10]	# pretmp_430, MEM <struct plane_t[2]> [(struct plane_t *)&planes][current_106].size[1]
# include/stencil_template_parallel.h:167:                     if ( y == 1 ){ // source on easr edge case
	jne	.L450	#,
# include/stencil_template_parallel.h:168:                         data[IDX(x, plane->size[_y_]+1)] += energy;
	lea	r11d, 1[r8]	# tmp382,
	imul	esi, r11d	# tmp383, tmp382
	add	esi, edx	# tmp385, _216
	mov	esi, esi	# tmp385, tmp385
	lea	rsi, [r12+rsi*8]	# _245,
# include/stencil_template_parallel.h:168:                         data[IDX(x, plane->size[_y_]+1)] += energy;
	vaddsd	xmm0, xmm8, QWORD PTR [rsi]	# tmp387, pretmp_431, *_245
	vmovsd	QWORD PTR [rsi], xmm0	# *_245, tmp387
.L450:
# include/stencil_template_parallel.h:171:                     if ( y == plane->size[_y_] ) { // source on west edge case
	cmp	ecx, r8d	# _217, pretmp_430
	jne	.L447	#,
# include/stencil_template_parallel.h:172:                         data[IDX(x, 0)] += energy;
	lea	rdx, [r12+rdx*8]	# _251,
# include/stencil_template_parallel.h:172:                         data[IDX(x, 0)] += energy;
	vaddsd	xmm0, xmm8, QWORD PTR [rdx]	# tmp391, pretmp_431, *_251
	vmovsd	QWORD PTR [rdx], xmm0	# *_251, tmp391
	jmp	.L447	#
.L491:
# src/stencil_template_parallel.c:168:   double ps_inject  = t_inject     / denom;
	vmovsd	xmm4, QWORD PTR -440[rbp]	# t_inject, %sfp
# src/stencil_template_parallel.c:171:   double ps_unpack  = t_unpack     / denom;
	vmovsd	xmm5, QWORD PTR -456[rbp]	# t_unpack, %sfp
# src/stencil_template_parallel.c:167:   double ps_total   = my_total     / denom;       // per-step (rank 0’s own)
	vdivsd	xmm10, xmm9, xmm0	# ps_total, my_total, iftmp.48_70
# src/stencil_template_parallel.c:177:       int P; MPI_Comm_size(myCOMM_WORLD, &P);
	lea	rsi, -304[rbp]	# tmp458,
	mov	rdi, QWORD PTR -328[rbp]	#, myCOMM_WORLD
	vmovsd	QWORD PTR -480[rbp], xmm8	# %sfp, pretmp_431
# src/stencil_template_parallel.c:167:   double ps_total   = my_total     / denom;       // per-step (rank 0’s own)
	vmovsd	QWORD PTR -472[rbp], xmm9	# %sfp, my_total
# src/stencil_template_parallel.c:168:   double ps_inject  = t_inject     / denom;
	vdivsd	xmm2, xmm4, xmm0	# ps_inject, t_inject, iftmp.48_70
# src/stencil_template_parallel.c:169:   double ps_pack    = t_pack       / denom;
	vmovsd	xmm4, QWORD PTR -448[rbp]	# t_pack, %sfp
	vdivsd	xmm3, xmm4, xmm0	# ps_pack, t_pack, iftmp.48_70
# src/stencil_template_parallel.c:170:   double ps_wait    = t_wait       / denom;
	vmovsd	xmm4, QWORD PTR -432[rbp]	# t_wait, %sfp
# src/stencil_template_parallel.c:174:   double ps_other   = ps_total - (ps_inject + ps_comm + ps_compute);  // may be ~0 ± rounding
	vmovsd	QWORD PTR -432[rbp], xmm10	# %sfp, ps_total
# src/stencil_template_parallel.c:174:   double ps_other   = ps_total - (ps_inject + ps_comm + ps_compute);  // may be ~0 ± rounding
	vmovsd	QWORD PTR -392[rbp], xmm2	# %sfp, ps_inject
# src/stencil_template_parallel.c:170:   double ps_wait    = t_wait       / denom;
	vdivsd	xmm4, xmm4, xmm0	# ps_wait, t_wait, iftmp.48_70
# src/stencil_template_parallel.c:172:   double ps_comm    = ps_pack + ps_wait + ps_unpack;
	vmovsd	QWORD PTR -448[rbp], xmm3	# %sfp, ps_pack
# src/stencil_template_parallel.c:171:   double ps_unpack  = t_unpack     / denom;
	vdivsd	xmm5, xmm5, xmm0	# ps_unpack, t_unpack, iftmp.48_70
# src/stencil_template_parallel.c:172:   double ps_comm    = ps_pack + ps_wait + ps_unpack;
	vaddsd	xmm1, xmm3, xmm4	# tmp455, ps_pack, ps_wait
	vmovsd	QWORD PTR -456[rbp], xmm4	# %sfp, ps_wait
# src/stencil_template_parallel.c:173:   double ps_compute = t_compute    / denom;
	vmovsd	xmm4, QWORD PTR -464[rbp]	# t_compute, %sfp
# src/stencil_template_parallel.c:172:   double ps_comm    = ps_pack + ps_wait + ps_unpack;
	vaddsd	xmm3, xmm1, xmm5	# ps_comm, tmp455, ps_unpack
# src/stencil_template_parallel.c:173:   double ps_compute = t_compute    / denom;
	vdivsd	xmm1, xmm4, xmm0	# ps_compute, t_compute, iftmp.48_70
# src/stencil_template_parallel.c:172:   double ps_comm    = ps_pack + ps_wait + ps_unpack;
	vmovsd	QWORD PTR -440[rbp], xmm5	# %sfp, ps_unpack
# src/stencil_template_parallel.c:174:   double ps_other   = ps_total - (ps_inject + ps_comm + ps_compute);  // may be ~0 ± rounding
	vaddsd	xmm0, xmm3, xmm2	# tmp456, ps_comm, ps_inject
	vmovsd	QWORD PTR -408[rbp], xmm3	# %sfp, ps_comm
# src/stencil_template_parallel.c:174:   double ps_other   = ps_total - (ps_inject + ps_comm + ps_compute);  // may be ~0 ± rounding
	vaddsd	xmm0, xmm0, xmm1	# tmp457, tmp456, ps_compute
	vmovsd	QWORD PTR -416[rbp], xmm1	# %sfp, ps_compute
# src/stencil_template_parallel.c:174:   double ps_other   = ps_total - (ps_inject + ps_comm + ps_compute);  // may be ~0 ± rounding
	vsubsd	xmm5, xmm10, xmm0	# ps_other, ps_total, tmp457
	vmovsd	QWORD PTR -424[rbp], xmm5	# %sfp, ps_other
# src/stencil_template_parallel.c:177:       int P; MPI_Comm_size(myCOMM_WORLD, &P);
	call	MPI_Comm_size@PLT	#
# /usr/include/x86_64-linux-gnu/bits/stdio2.h:86:   return __printf_chk (__USE_FORTIFY_LEVEL - 1, __fmt, __va_arg_pack ());
	lea	rsi, .LC30[rip]	# tmp461,
	mov	edi, 2	#,
	push	rax	#
	vmovsd	xmm10, QWORD PTR -432[rbp]	# ps_total, %sfp
	mov	eax, 8	#,
	push	QWORD PTR -424[rbp]	# %sfp
	vmovsd	xmm9, QWORD PTR -472[rbp]	# my_total, %sfp
	mov	edx, DWORD PTR -304[rbp]	#, MEM[(int *)_201]
	vmovsd	xmm1, xmm10, xmm10	#, ps_total
	vmovsd	xmm7, QWORD PTR -392[rbp]	#, %sfp
	vmovsd	xmm6, QWORD PTR -416[rbp]	#, %sfp
	vmovsd	xmm0, xmm9, xmm9	#, my_total
	vmovsd	xmm5, QWORD PTR -440[rbp]	# ps_unpack, %sfp
	vmovsd	xmm4, QWORD PTR -456[rbp]	# ps_wait, %sfp
	vmovsd	xmm3, QWORD PTR -448[rbp]	# ps_pack, %sfp
	vmovsd	xmm2, QWORD PTR -408[rbp]	#, %sfp
	call	__printf_chk@PLT	#
# src/stencil_template_parallel.c:186:       if (ps_total > 0.0) {
	vmovsd	xmm10, QWORD PTR -432[rbp]	# ps_total, %sfp
	pop	rdx	#
	vxorpd	xmm0, xmm0, xmm0	# tmp462
	vmovsd	xmm8, QWORD PTR -480[rbp]	# pretmp_431, %sfp
	pop	rcx	#
	vcomisd	xmm10, xmm0	# ps_total, tmp462
	jbe	.L464	#,
# src/stencil_template_parallel.c:191:                 100.0*ps_other/ps_total);
	vmovsd	xmm0, QWORD PTR .LC31[rip]	# tmp464,
	vmulsd	xmm3, xmm0, QWORD PTR -424[rbp]	# tmp463, tmp464, %sfp
# /usr/include/x86_64-linux-gnu/bits/stdio2.h:86:   return __printf_chk (__USE_FORTIFY_LEVEL - 1, __fmt, __va_arg_pack ());
	lea	rsi, .LC32[rip]	# tmp475,
	mov	edi, 2	#,
# src/stencil_template_parallel.c:190:                 100.0*ps_inject/ps_total,
	vmulsd	xmm2, xmm0, QWORD PTR -392[rbp]	# tmp466, tmp464, %sfp
# /usr/include/x86_64-linux-gnu/bits/stdio2.h:86:   return __printf_chk (__USE_FORTIFY_LEVEL - 1, __fmt, __va_arg_pack ());
	mov	eax, 4	#,
	vmovsd	QWORD PTR -432[rbp], xmm8	# %sfp, pretmp_431
# src/stencil_template_parallel.c:189:                 100.0*ps_compute/ps_total,
	vmulsd	xmm1, xmm0, QWORD PTR -416[rbp]	# tmp469, tmp464, %sfp
# src/stencil_template_parallel.c:188:                 100.0*ps_comm/ps_total,
	vmulsd	xmm0, xmm0, QWORD PTR -408[rbp]	# tmp472, tmp464, %sfp
# /usr/include/x86_64-linux-gnu/bits/stdio2.h:86:   return __printf_chk (__USE_FORTIFY_LEVEL - 1, __fmt, __va_arg_pack ());
	vdivsd	xmm3, xmm3, xmm10	#, tmp463, ps_total
# src/stencil_template_parallel.c:187:           printf("WEAK_PCT comm=%.1f%% compute=%.1f%% inject=%.1f%% other=%.1f%%\n",
	vdivsd	xmm0, xmm0, xmm10	# tmp474, tmp472, ps_total
# /usr/include/x86_64-linux-gnu/bits/stdio2.h:86:   return __printf_chk (__USE_FORTIFY_LEVEL - 1, __fmt, __va_arg_pack ());
	vdivsd	xmm2, xmm2, xmm10	#, tmp466, ps_total
	vdivsd	xmm1, xmm1, xmm10	#, tmp469, ps_total
	call	__printf_chk@PLT	#
	vmovsd	xmm8, QWORD PTR -432[rbp]	# pretmp_431, %sfp
	jmp	.L464	#
.L493:
# include/stencil_template_parallel.h:151:                         data[IDX(plane->size[_x_]+1, y)] += energy;} // propagate on west halo
	lea	r8d, 1[r9+rsi]	# tmp370,
	lea	r8, [r12+r8*8]	# _230,
# include/stencil_template_parallel.h:151:                         data[IDX(plane->size[_x_]+1, y)] += energy;} // propagate on west halo
	vaddsd	xmm0, xmm8, QWORD PTR [r8]	# tmp372, pretmp_431, *_230
	vmovsd	QWORD PTR [r8], xmm0	# *_230, tmp372
	jmp	.L446	#
.L487:
# /usr/include/x86_64-linux-gnu/bits/stdio2.h:86:   return __printf_chk (__USE_FORTIFY_LEVEL - 1, __fmt, __va_arg_pack ());
	mov	edx, DWORD PTR -360[rbp]	#, Rank
	mov	ecx, r15d	#, current
	mov	edi, 2	#,
	xor	eax, eax	#
	lea	rsi, .LC29[rip]	# tmp331,
	call	__printf_chk@PLT	#
# src/stencil_template_parallel.c:81:     MPI_Finalize();
	call	MPI_Finalize@PLT	#
# src/stencil_template_parallel.c:82:     return 0;
	jmp	.L433	#
.L434:
# src/stencil_template_parallel.c:162:   MPI_Barrier(myCOMM_WORLD);                      // end of measured region
	mov	rdi, QWORD PTR -328[rbp]	#, myCOMM_WORLD
	vmovsd	QWORD PTR -392[rbp], xmm8	# %sfp, pretmp_431
	call	MPI_Barrier@PLT	#
# src/stencil_template_parallel.c:164:   double my_total = MPI_Wtime() - t_loop0;
	call	MPI_Wtime@PLT	#
	vmovsd	xmm8, QWORD PTR -392[rbp]	# pretmp_431, %sfp
# src/stencil_template_parallel.c:164:   double my_total = MPI_Wtime() - t_loop0;
	vsubsd	xmm9, xmm0, QWORD PTR -656[rbp]	# my_total, tmp518, %sfp
# src/stencil_template_parallel.c:90:   double t_pack=0.0, t_wait=0.0, t_unpack=0.0, t_inject=0.0, t_compute=0.0;
	mov	QWORD PTR -464[rbp], 0x000000000	# %sfp,
# src/stencil_template_parallel.c:90:   double t_pack=0.0, t_wait=0.0, t_unpack=0.0, t_inject=0.0, t_compute=0.0;
	mov	QWORD PTR -440[rbp], 0x000000000	# %sfp,
# src/stencil_template_parallel.c:165:   double denom    = (M > 0 ? (double)M : 1.0);
	vmovsd	xmm0, QWORD PTR .LC10[rip]	# iftmp.48_70,
# src/stencil_template_parallel.c:90:   double t_pack=0.0, t_wait=0.0, t_unpack=0.0, t_inject=0.0, t_compute=0.0;
	mov	QWORD PTR -456[rbp], 0x000000000	# %sfp,
# src/stencil_template_parallel.c:90:   double t_pack=0.0, t_wait=0.0, t_unpack=0.0, t_inject=0.0, t_compute=0.0;
	mov	QWORD PTR -448[rbp], 0x000000000	# %sfp,
# src/stencil_template_parallel.c:90:   double t_pack=0.0, t_wait=0.0, t_unpack=0.0, t_inject=0.0, t_compute=0.0;
	mov	QWORD PTR -432[rbp], 0x000000000	# %sfp,
	jmp	.L467	#
.L492:
# src/stencil_template_parallel.c:202: }
	call	__stack_chk_fail@PLT	#
.L486:
# /usr/include/x86_64-linux-gnu/bits/stdio2.h:86:   return __printf_chk (__USE_FORTIFY_LEVEL - 1, __fmt, __va_arg_pack ());
	mov	edi, 2	#,
	mov	ecx, 1	#,
	lea	rsi, .LC27[rip]	# tmp304,
	xor	eax, eax	#
	call	__printf_chk@PLT	#
# src/stencil_template_parallel.c:55:     MPI_Finalize();
	call	MPI_Finalize@PLT	#
# src/stencil_template_parallel.c:56:     exit(1); 
	mov	edi, 1	#,
	call	exit@PLT	#
	.cfi_endproc
.LFE51:
	.size	main, .-main
	.section	.rodata.cst8,"aM",@progbits,8
	.align 8
.LC0:
	.long	0
	.long	1070596096
	.align 8
.LC1:
	.long	0
	.long	1071644672
	.align 8
.LC9:
	.long	10000
	.long	10000
	.align 8
.LC10:
	.long	0
	.long	1072693248
	.align 8
.LC31:
	.long	0
	.long	1079574528
	.ident	"GCC: (Ubuntu 13.3.0-6ubuntu2~24.04) 13.3.0"
	.section	.note.GNU-stack,"",@progbits
	.section	.note.gnu.property,"a"
	.align 8
	.long	1f - 0f
	.long	4f - 1f
	.long	5
0:
	.string	"GNU"
1:
	.align 8
	.long	0xc0000002
	.long	3f - 2f
2:
	.long	0x3
3:
	.align 8
4:
