	.file	"SIMD.cpp"
	.text
	.p2align 4
	.globl	_ZN4poly9poly_base10mem_helper4alocEm
	.type	_ZN4poly9poly_base10mem_helper4alocEm, @function
_ZN4poly9poly_base10mem_helper4alocEm:
.LFB8949:
	.cfi_startproc
	endbr64
	movq	_ZN4poly9poly_base10mem_helper3nowE(%rip), %rax
	addq	%rax, %rdi
	movq	%rdi, _ZN4poly9poly_base10mem_helper3nowE(%rip)
	ret
	.cfi_endproc
.LFE8949:
	.size	_ZN4poly9poly_base10mem_helper4alocEm, .-_ZN4poly9poly_base10mem_helper4alocEm
	.p2align 4
	.globl	_ZN4poly9poly_base10mem_helper5alocPEi
	.type	_ZN4poly9poly_base10mem_helper5alocPEi, @function
_ZN4poly9poly_base10mem_helper5alocPEi:
.LFB8950:
	.cfi_startproc
	endbr64
	movq	_ZN4poly9poly_base10mem_helper3nowE(%rip), %rax
	movslq	%edi, %rdi
	addq	$31, %rax
	andq	$-32, %rax
	leaq	(%rax,%rdi,4), %rdx
	movq	%rdx, _ZN4poly9poly_base10mem_helper3nowE(%rip)
	ret
	.cfi_endproc
.LFE8950:
	.size	_ZN4poly9poly_base10mem_helper5alocPEi, .-_ZN4poly9poly_base10mem_helper5alocPEi
	.p2align 4
	.globl	_ZN4poly3dotEPjiPKj
	.type	_ZN4poly3dotEPjiPKj, @function
_ZN4poly3dotEPjiPKj:
.LFB8951:
	.cfi_startproc
	endbr64
	cmpl	$7, %esi
	jle	.L15
	leal	-8(%rsi), %r8d
	vmovdqa	.LC0(%rip), %ymm6
	vmovdqa	.LC1(%rip), %ymm5
	xorl	%eax, %eax
	shrl	$3, %r8d
	movl	%r8d, %ecx
	incq	%rcx
	salq	$5, %rcx
	.p2align 4,,10
	.p2align 3
.L8:
	vmovdqa	(%rdi,%rax), %ymm0
	vmovdqa	(%rdx,%rax), %ymm4
	vpsrlq	$32, %ymm0, %ymm1
	vpsrlq	$32, %ymm4, %ymm2
	vpmuludq	%ymm2, %ymm1, %ymm2
	vpmulld	%ymm0, %ymm4, %ymm1
	vpmuludq	%ymm4, %ymm0, %ymm0
	vpmulld	%ymm6, %ymm1, %ymm1
	vpshufd	$245, %ymm1, %ymm3
	vpmuludq	%ymm5, %ymm1, %ymm1
	vpmuludq	%ymm5, %ymm3, %ymm3
	vpaddq	%ymm1, %ymm0, %ymm0
	vpsrlq	$32, %ymm0, %ymm0
	vpaddq	%ymm3, %ymm2, %ymm1
	vpblendd	$170, %ymm1, %ymm0, %ymm0
	vmovdqa	%ymm0, (%rdi,%rax)
	addq	$32, %rax
	cmpq	%rax, %rcx
	jne	.L8
	leal	8(,%r8,8), %ecx
	vzeroupper
.L9:
	cmpl	%ecx, %esi
	jle	.L16
	movslq	%ecx, %rcx
	.p2align 4,,10
	.p2align 3
.L10:
	movl	(%rdi,%rcx,4), %eax
	movl	(%rdx,%rcx,4), %r8d
	imulq	%rax, %r8
	imull	$104857599, %r8d, %eax
	imulq	$104857601, %rax, %rax
	addq	%r8, %rax
	shrq	$32, %rax
	movl	%eax, (%rdi,%rcx,4)
	incq	%rcx
	cmpl	%ecx, %esi
	jg	.L10
	ret
	.p2align 4,,10
	.p2align 3
.L16:
	ret
	.p2align 4,,10
	.p2align 3
.L15:
	xorl	%ecx, %ecx
	jmp	.L9
	.cfi_endproc
.LFE8951:
	.size	_ZN4poly3dotEPjiPKj, .-_ZN4poly3dotEPjiPKj
	.p2align 4
	.globl	_ZN4poly3dotEPKjS1_iPj
	.type	_ZN4poly3dotEPKjS1_iPj, @function
_ZN4poly3dotEPKjS1_iPj:
.LFB8952:
	.cfi_startproc
	endbr64
	cmpl	$7, %edx
	jle	.L28
	leal	-8(%rdx), %r9d
	vmovdqa	.LC0(%rip), %ymm6
	vmovdqa	.LC1(%rip), %ymm5
	xorl	%eax, %eax
	shrl	$3, %r9d
	movl	%r9d, %r8d
	incq	%r8
	salq	$5, %r8
	.p2align 4,,10
	.p2align 3
.L21:
	vmovdqa	(%rdi,%rax), %ymm0
	vmovdqa	(%rsi,%rax), %ymm4
	vpsrlq	$32, %ymm0, %ymm1
	vpsrlq	$32, %ymm4, %ymm2
	vpmuludq	%ymm2, %ymm1, %ymm2
	vpmulld	%ymm0, %ymm4, %ymm1
	vpmuludq	%ymm4, %ymm0, %ymm0
	vpmulld	%ymm6, %ymm1, %ymm1
	vpshufd	$245, %ymm1, %ymm3
	vpmuludq	%ymm5, %ymm1, %ymm1
	vpmuludq	%ymm5, %ymm3, %ymm3
	vpaddq	%ymm1, %ymm0, %ymm0
	vpsrlq	$32, %ymm0, %ymm0
	vpaddq	%ymm3, %ymm2, %ymm1
	vpblendd	$170, %ymm1, %ymm0, %ymm0
	vmovdqa	%ymm0, (%rcx,%rax)
	addq	$32, %rax
	cmpq	%rax, %r8
	jne	.L21
	leal	8(,%r9,8), %r8d
	vzeroupper
.L22:
	cmpl	%r8d, %edx
	jle	.L29
	movslq	%r8d, %r8
	.p2align 4,,10
	.p2align 3
.L23:
	movl	(%rdi,%r8,4), %eax
	movl	(%rsi,%r8,4), %r9d
	imulq	%rax, %r9
	imull	$104857599, %r9d, %eax
	imulq	$104857601, %rax, %rax
	addq	%r9, %rax
	shrq	$32, %rax
	movl	%eax, (%rcx,%r8,4)
	incq	%r8
	cmpl	%r8d, %edx
	jg	.L23
	ret
	.p2align 4,,10
	.p2align 3
.L29:
	ret
	.p2align 4,,10
	.p2align 3
.L28:
	xorl	%r8d, %r8d
	jmp	.L22
	.cfi_endproc
.LFE8952:
	.size	_ZN4poly3dotEPKjS1_iPj, .-_ZN4poly3dotEPKjS1_iPj
	.section	.rodata.str1.1,"aMS",@progbits,1
.LC4:
	.string	"%d %d %d\n"
	.text
	.p2align 4
	.globl	_ZN4poly7f_n_t_t8dif_funcEPv
	.type	_ZN4poly7f_n_t_t8dif_funcEPv, @function
_ZN4poly7f_n_t_t8dif_funcEPv:
.LFB8965:
	.cfi_startproc
	endbr64
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	pushq	%r15
	pushq	%r14
	.cfi_offset 15, -24
	.cfi_offset 14, -32
	leaq	_ZN4poly7f_n_t_tL4iab4E(%rip), %r14
	pushq	%r13
	pushq	%r12
	pushq	%rbx
	andq	$-32, %rsp
	subq	$224, %rsp
	.cfi_offset 13, -40
	.cfi_offset 12, -48
	.cfi_offset 3, -56
	movl	(%rdi), %r13d
	movq	8(%rdi), %rdx
	movl	20(%rdi), %r15d
	movl	16(%rdi), %r12d
	movl	%r13d, %eax
	notl	%eax
	sarl	$3, %r12d
	leal	0(,%r15,4), %ebx
	movl	%eax, 28(%rsp)
	movslq	%r13d, %rax
	movq	%rax, 8(%rsp)
	leaq	32(%rdx), %rax
	movq	%rax, 16(%rsp)
	.p2align 4,,10
	.p2align 3
.L37:
	xorl	%eax, %eax
	movl	$23, %r8d
	movl	$1, %ecx
	movl	%r15d, %edx
	leaq	.LC4(%rip), %rsi
	movl	$1, %edi
	call	__printf_chk@PLT
	testl	%r12d, %r12d
	vmovdqa	.LC1(%rip), %ymm12
	vmovdqa	.LC6(%rip), %ymm15
	vmovdqa	.LC7(%rip), %ymm14
	jle	.L31
	movl	28(%rsp), %eax
	movslq	%ebx, %r11
	movslq	%r15d, %rsi
	movl	$-2, %r9d
	vmovdqa	.LC2(%rip), %ymm4
	salq	$5, %r11
	salq	$5, %rsi
	xorl	%r8d, %r8d
	leal	(%rax,%r15), %r10d
	movq	8(%rsp), %rax
	vmovdqa	.LC3(%rip), %ymm11
	vmovdqa	%ymm4, 32(%rsp)
	leaq	(%r10,%rax), %rdi
	leal	(%r15,%r15), %eax
	notq	%r10
	leal	(%rax,%r15), %edx
	salq	$5, %rdi
	movslq	%eax, %rcx
	addq	16(%rsp), %rdi
	salq	$5, %rcx
	movslq	%edx, %rdx
	salq	$5, %r10
	salq	$5, %rdx
	.p2align 4,,10
	.p2align 3
.L36:
	vpmulld	%ymm11, %ymm11, %ymm1
	vpsrlq	$32, %ymm11, %ymm4
	vpmulld	.LC0(%rip), %ymm11, %ymm9
	vpmuludq	%ymm11, %ymm11, %ymm8
	vpmuludq	%ymm4, %ymm4, %ymm0
	vmovdqa	%ymm4, 64(%rsp)
	vpmulld	.LC0(%rip), %ymm1, %ymm1
	vpshufd	$245, %ymm1, %ymm2
	vpmuludq	%ymm12, %ymm1, %ymm1
	vmovdqa	%ymm4, %ymm13
	vpmuludq	%ymm12, %ymm2, %ymm2
	vpaddq	%ymm1, %ymm8, %ymm8
	vpaddq	%ymm2, %ymm0, %ymm0
	vpsrlq	$32, %ymm8, %ymm8
	vpblendd	$170, %ymm0, %ymm8, %ymm8
	vpaddd	.LC5(%rip), %ymm8, %ymm0
	vpsrad	$31, %ymm0, %ymm8
	vpand	%ymm12, %ymm8, %ymm8
	vpaddd	%ymm0, %ymm8, %ymm8
	vpmulld	%ymm9, %ymm8, %ymm1
	vpsrlq	$32, %ymm8, %ymm0
	vmovdqa	%ymm0, 192(%rsp)
	vpmuludq	%ymm11, %ymm8, %ymm10
	vpmuludq	%ymm4, %ymm0, %ymm0
	vpshufd	$245, %ymm1, %ymm2
	vpmuludq	%ymm12, %ymm1, %ymm1
	vpmuludq	%ymm12, %ymm2, %ymm2
	vpaddq	%ymm1, %ymm10, %ymm10
	vpaddq	%ymm2, %ymm0, %ymm0
	vpsrlq	$32, %ymm10, %ymm10
	vpblendd	$170, %ymm0, %ymm10, %ymm10
	vpaddd	.LC5(%rip), %ymm10, %ymm0
	vpsrad	$31, %ymm0, %ymm10
	vpand	%ymm12, %ymm10, %ymm10
	vpaddd	%ymm0, %ymm10, %ymm10
	cmpl	%r13d, %r15d
	jle	.L34
	vpmulld	.LC0(%rip), %ymm8, %ymm4
	vmovdqa	%ymm4, 96(%rsp)
	vpmulld	.LC0(%rip), %ymm10, %ymm4
	leaq	(%r10,%rdi), %rax
	vmovdqa	%ymm4, 128(%rsp)
	vpsrlq	$32, %ymm10, %ymm4
	vmovdqa	%ymm4, 160(%rsp)
	.p2align 4,,10
	.p2align 3
.L35:
	vmovdqa	(%rax,%rsi), %ymm0
	vpaddd	(%rax), %ymm15, %ymm2
	addq	$32, %rax
	vpmulld	%ymm9, %ymm0, %ymm3
	vpsrad	$31, %ymm2, %ymm4
	vpand	%ymm14, %ymm4, %ymm4
	vpaddd	%ymm2, %ymm4, %ymm4
	vpsrlq	$32, %ymm0, %ymm2
	vpmuludq	%ymm11, %ymm0, %ymm0
	vpmuludq	%ymm13, %ymm2, %ymm2
	vpshufd	$245, %ymm3, %ymm5
	vpmuludq	%ymm12, %ymm3, %ymm3
	vpmuludq	%ymm12, %ymm5, %ymm5
	vpaddq	%ymm3, %ymm0, %ymm3
	vmovdqa	-32(%rax,%rcx), %ymm0
	vpaddq	%ymm5, %ymm2, %ymm1
	vpsrlq	$32, %ymm3, %ymm3
	vpmulld	96(%rsp), %ymm0, %ymm5
	vpshufd	$245, %ymm5, %ymm6
	vpmuludq	%ymm12, %ymm5, %ymm5
	vpmuludq	%ymm12, %ymm6, %ymm6
	vpblendd	$170, %ymm1, %ymm3, %ymm1
	vpsrlq	$32, %ymm0, %ymm3
	vpmuludq	%ymm8, %ymm0, %ymm0
	vpmuludq	192(%rsp), %ymm3, %ymm3
	vpaddq	%ymm6, %ymm3, %ymm2
	vpaddq	%ymm5, %ymm0, %ymm5
	vmovdqa	-32(%rax,%rdx), %ymm0
	vpmulld	128(%rsp), %ymm0, %ymm6
	vpshufd	$245, %ymm6, %ymm7
	vpsrlq	$32, %ymm5, %ymm5
	vpmuludq	%ymm12, %ymm6, %ymm6
	vpblendd	$170, %ymm2, %ymm5, %ymm2
	vpsrlq	$32, %ymm0, %ymm5
	vpmuludq	160(%rsp), %ymm5, %ymm5
	vpmuludq	%ymm10, %ymm0, %ymm0
	vpmuludq	%ymm12, %ymm7, %ymm7
	vpaddq	%ymm6, %ymm0, %ymm6
	vpaddq	%ymm7, %ymm5, %ymm3
	vpsrlq	$32, %ymm6, %ymm6
	vpaddd	%ymm14, %ymm1, %ymm0
	vpblendd	$170, %ymm3, %ymm6, %ymm3
	vpsubd	%ymm3, %ymm0, %ymm0
	vpaddd	%ymm3, %ymm1, %ymm1
	vpmulld	.LC9(%rip), %ymm0, %ymm6
	vpshufd	$245, %ymm6, %ymm7
	vpmuludq	%ymm12, %ymm6, %ymm6
	vpmuludq	%ymm12, %ymm7, %ymm7
	vpaddd	%ymm15, %ymm1, %ymm1
	vpsrlq	$32, %ymm0, %ymm5
	vpmuludq	.LC10(%rip), %ymm0, %ymm0
	vpmuludq	.LC8(%rip), %ymm5, %ymm5
	vpaddq	%ymm6, %ymm0, %ymm0
	vpaddd	%ymm15, %ymm2, %ymm6
	vpaddd	%ymm4, %ymm6, %ymm6
	vpsrlq	$32, %ymm0, %ymm0
	vpaddq	%ymm7, %ymm5, %ymm5
	vpblendd	$170, %ymm5, %ymm0, %ymm5
	vpsrad	$31, %ymm6, %ymm0
	vpsubd	%ymm2, %ymm4, %ymm2
	vpand	%ymm14, %ymm0, %ymm0
	vpsrad	$31, %ymm2, %ymm3
	vpaddd	%ymm6, %ymm0, %ymm0
	vpsrad	$31, %ymm1, %ymm6
	vpand	%ymm14, %ymm3, %ymm3
	vpand	%ymm14, %ymm6, %ymm6
	vpaddd	%ymm2, %ymm3, %ymm2
	vpaddd	%ymm1, %ymm6, %ymm1
	vpaddd	%ymm0, %ymm1, %ymm3
	vpsubd	%ymm1, %ymm14, %ymm1
	vpaddd	%ymm0, %ymm1, %ymm0
	vmovdqa	%ymm3, -32(%rax)
	vmovdqa	%ymm0, -32(%rax,%rsi)
	vpaddd	%ymm2, %ymm5, %ymm0
	vpsubd	%ymm5, %ymm14, %ymm5
	vpaddd	%ymm2, %ymm5, %ymm2
	vmovdqa	%ymm0, -32(%rax,%rcx)
	vmovdqa	%ymm2, -32(%rax,%rdx)
	cmpq	%rdi, %rax
	jne	.L35
.L34:
	vmovdqa	32(%rsp), %ymm4
	addl	%ebx, %r8d
	addq	%r11, %rdi
	vpmulld	%ymm4, %ymm9, %ymm9
	vpmuludq	%ymm4, %ymm11, %ymm11
	vpsrlq	$32, %ymm4, %ymm0
	vpmuludq	64(%rsp), %ymm0, %ymm0
	vpshufd	$245, %ymm9, %ymm1
	vpmuludq	%ymm12, %ymm9, %ymm9
	vpmuludq	%ymm12, %ymm1, %ymm1
	vpaddq	%ymm9, %ymm11, %ymm11
	vpaddq	%ymm1, %ymm0, %ymm0
	vpsrlq	$32, %ymm11, %ymm11
	vpblendd	$170, %ymm0, %ymm11, %ymm11
	vpaddd	.LC5(%rip), %ymm11, %ymm11
	vpsrad	$31, %ymm11, %ymm0
	vpand	%ymm12, %ymm0, %ymm0
	vpaddd	%ymm11, %ymm0, %ymm11
	cmpl	%r8d, %r12d
	jle	.L31
	xorl	%eax, %eax
	tzcntl	%r9d, %eax
	decl	%r9d
	cltq
	salq	$5, %rax
	vmovdqa	(%r14,%rax), %ymm4
	vmovdqa	%ymm4, 32(%rsp)
	jmp	.L36
	.p2align 4,,10
	.p2align 3
.L31:
	sarl	$2, %r15d
	sarl	$2, %ebx
	vzeroupper
	jmp	.L37
	.cfi_endproc
.LFE8965:
	.size	_ZN4poly7f_n_t_t8dif_funcEPv, .-_ZN4poly7f_n_t_t8dif_funcEPv
	.p2align 4
	.globl	_ZN4poly7f_n_t_t8dit_funcEPv
	.type	_ZN4poly7f_n_t_t8dit_funcEPv, @function
_ZN4poly7f_n_t_t8dit_funcEPv:
.LFB8967:
	.cfi_startproc
	endbr64
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	pushq	%r15
	.cfi_offset 15, -24
	leaq	_ZN4poly7f_n_t_tL4iab4E(%rip), %r15
	pushq	%r14
	pushq	%r13
	pushq	%r12
	pushq	%rbx
	andq	$-32, %rsp
	subq	$192, %rsp
	.cfi_offset 14, -32
	.cfi_offset 13, -40
	.cfi_offset 12, -48
	.cfi_offset 3, -56
	movl	(%rdi), %r14d
	movq	8(%rdi), %rdx
	movl	20(%rdi), %ebx
	movl	16(%rdi), %r13d
	movl	%r14d, %eax
	vmovdqa	.LC1(%rip), %ymm12
	vmovdqa	.LC7(%rip), %ymm15
	notl	%eax
	vmovdqa	.LC6(%rip), %ymm14
	sarl	$3, %r13d
	leal	0(,%rbx,4), %r12d
	movl	%eax, 12(%rsp)
	movslq	%r14d, %rax
	movq	%rax, 24(%rsp)
	leaq	32(%rdx), %rax
	movq	%rax, 16(%rsp)
	.p2align 4,,10
	.p2align 3
.L42:
	testl	%r13d, %r13d
	jle	.L44
	movl	12(%rsp), %eax
	movslq	%r12d, %r11
	movslq	%ebx, %rsi
	movl	$-1, %r9d
	vmovdqa	.LC3(%rip), %ymm11
	salq	$5, %r11
	salq	$5, %rsi
	xorl	%r8d, %r8d
	leal	(%rax,%rbx), %r10d
	movq	24(%rsp), %rax
	leaq	(%r10,%rax), %rdi
	leal	(%rbx,%rbx), %eax
	notq	%r10
	leal	(%rax,%rbx), %edx
	salq	$5, %rdi
	movslq	%eax, %rcx
	addq	16(%rsp), %rdi
	salq	$5, %rcx
	movslq	%edx, %rdx
	salq	$5, %r10
	salq	$5, %rdx
	.p2align 4,,10
	.p2align 3
.L45:
	vpmulld	%ymm11, %ymm11, %ymm1
	vpsrlq	$32, %ymm11, %ymm6
	vpmulld	.LC0(%rip), %ymm11, %ymm9
	vpmuludq	%ymm11, %ymm11, %ymm8
	vpmuludq	%ymm6, %ymm6, %ymm0
	vmovdqa	%ymm6, 32(%rsp)
	vpmulld	.LC0(%rip), %ymm1, %ymm1
	vpshufd	$245, %ymm1, %ymm2
	vpmuludq	%ymm12, %ymm1, %ymm1
	vmovdqa	%ymm6, %ymm13
	vpmuludq	%ymm12, %ymm2, %ymm2
	vpaddq	%ymm1, %ymm8, %ymm8
	vpaddq	%ymm2, %ymm0, %ymm0
	vpsrlq	$32, %ymm8, %ymm8
	vpblendd	$170, %ymm0, %ymm8, %ymm8
	vpaddd	.LC5(%rip), %ymm8, %ymm0
	vpsrad	$31, %ymm0, %ymm8
	vpand	%ymm12, %ymm8, %ymm8
	vpaddd	%ymm0, %ymm8, %ymm8
	vpmulld	%ymm9, %ymm8, %ymm1
	vpsrlq	$32, %ymm8, %ymm0
	vmovdqa	%ymm0, 160(%rsp)
	vpmuludq	%ymm11, %ymm8, %ymm10
	vpmuludq	%ymm6, %ymm0, %ymm0
	vpshufd	$245, %ymm1, %ymm2
	vpmuludq	%ymm12, %ymm1, %ymm1
	vpmuludq	%ymm12, %ymm2, %ymm2
	vpaddq	%ymm1, %ymm10, %ymm10
	vpaddq	%ymm2, %ymm0, %ymm0
	vpsrlq	$32, %ymm10, %ymm10
	vpblendd	$170, %ymm0, %ymm10, %ymm10
	vpaddd	.LC5(%rip), %ymm10, %ymm0
	vpsrad	$31, %ymm0, %ymm10
	vpand	%ymm12, %ymm10, %ymm10
	vpaddd	%ymm0, %ymm10, %ymm10
	cmpl	%ebx, %r14d
	jge	.L46
	vpmulld	.LC0(%rip), %ymm8, %ymm7
	vmovdqa	%ymm7, 96(%rsp)
	vpmulld	.LC0(%rip), %ymm10, %ymm7
	leaq	(%r10,%rdi), %rax
	vmovdqa	%ymm7, 128(%rsp)
	vpsrlq	$32, %ymm10, %ymm7
	vmovdqa	%ymm7, 64(%rsp)
	.p2align 4,,10
	.p2align 3
.L47:
	vpaddd	(%rax,%rcx), %ymm15, %ymm3
	vpsubd	(%rax,%rdx), %ymm3, %ymm3
	addq	$32, %rax
	vpmulld	.LC12(%rip), %ymm3, %ymm4
	vpshufd	$245, %ymm4, %ymm1
	vpmuludq	%ymm12, %ymm4, %ymm4
	vmovdqa	-32(%rax), %ymm7
	vpmuludq	.LC13(%rip), %ymm3, %ymm2
	vpmuludq	%ymm12, %ymm1, %ymm1
	vpsrlq	$32, %ymm3, %ymm0
	vpmuludq	.LC11(%rip), %ymm0, %ymm0
	vpaddq	%ymm4, %ymm2, %ymm2
	vpaddq	%ymm1, %ymm0, %ymm0
	vpaddd	-32(%rax,%rsi), %ymm7, %ymm1
	vpsrlq	$32, %ymm2, %ymm2
	vmovdqa	-32(%rax,%rcx), %ymm7
	vpblendd	$170, %ymm0, %ymm2, %ymm2
	vpaddd	-32(%rax,%rdx), %ymm7, %ymm4
	vmovdqa	-32(%rax), %ymm7
	vpaddd	%ymm14, %ymm1, %ymm0
	vpsubd	-32(%rax,%rsi), %ymm7, %ymm3
	vpsrad	$31, %ymm0, %ymm1
	vpand	%ymm15, %ymm1, %ymm1
	vpaddd	%ymm0, %ymm1, %ymm1
	vpaddd	%ymm14, %ymm4, %ymm0
	vpsrad	$31, %ymm0, %ymm4
	vpand	%ymm15, %ymm4, %ymm4
	vpaddd	%ymm0, %ymm4, %ymm4
	vpsrad	$31, %ymm3, %ymm0
	vpaddd	%ymm1, %ymm4, %ymm5
	vpand	%ymm15, %ymm0, %ymm0
	vpaddd	%ymm15, %ymm1, %ymm1
	vpaddd	%ymm14, %ymm5, %ymm5
	vpaddd	%ymm3, %ymm0, %ymm0
	vpsrad	$31, %ymm5, %ymm3
	vpsubd	%ymm4, %ymm1, %ymm4
	vpand	%ymm15, %ymm3, %ymm3
	vpaddd	%ymm5, %ymm3, %ymm3
	vmovdqa	%ymm3, -32(%rax)
	vpaddd	%ymm0, %ymm2, %ymm3
	vpaddd	%ymm15, %ymm0, %ymm0
	vpmulld	%ymm9, %ymm3, %ymm6
	vpsrlq	$32, %ymm3, %ymm5
	vpsubd	%ymm2, %ymm0, %ymm2
	vpmuludq	%ymm11, %ymm3, %ymm3
	vpmuludq	%ymm13, %ymm5, %ymm5
	vpsrlq	$32, %ymm2, %ymm1
	vpmuludq	64(%rsp), %ymm1, %ymm1
	vpshufd	$245, %ymm6, %ymm7
	vpmuludq	%ymm12, %ymm6, %ymm6
	vpmuludq	%ymm12, %ymm7, %ymm7
	vpaddq	%ymm6, %ymm3, %ymm6
	vpsrlq	$32, %ymm4, %ymm3
	vpmuludq	160(%rsp), %ymm3, %ymm3
	vpsrlq	$32, %ymm6, %ymm6
	vpaddq	%ymm7, %ymm5, %ymm5
	vpblendd	$170, %ymm5, %ymm6, %ymm5
	vmovdqa	%ymm5, -32(%rax,%rsi)
	vpmulld	96(%rsp), %ymm4, %ymm5
	vpmuludq	%ymm8, %ymm4, %ymm4
	vpshufd	$245, %ymm5, %ymm6
	vpmuludq	%ymm12, %ymm5, %ymm5
	vpmuludq	%ymm12, %ymm6, %ymm6
	vpaddq	%ymm5, %ymm4, %ymm4
	vpaddq	%ymm6, %ymm3, %ymm3
	vpsrlq	$32, %ymm4, %ymm4
	vpblendd	$170, %ymm3, %ymm4, %ymm4
	vmovdqa	%ymm4, -32(%rax,%rcx)
	vpmulld	128(%rsp), %ymm2, %ymm4
	vpmuludq	%ymm10, %ymm2, %ymm2
	vpshufd	$245, %ymm4, %ymm3
	vpmuludq	%ymm12, %ymm4, %ymm4
	vpmuludq	%ymm12, %ymm3, %ymm3
	vpaddq	%ymm4, %ymm2, %ymm2
	vpaddq	%ymm3, %ymm1, %ymm1
	vpsrlq	$32, %ymm2, %ymm2
	vpblendd	$170, %ymm1, %ymm2, %ymm2
	vmovdqa	%ymm2, -32(%rax,%rdx)
	cmpq	%rax, %rdi
	jne	.L47
.L46:
	xorl	%eax, %eax
	addl	%r12d, %r8d
	addq	%r11, %rdi
	tzcntl	%r9d, %eax
	decl	%r9d
	cltq
	salq	$5, %rax
	vmovdqa	640(%rax,%r15), %ymm0
	vpmulld	%ymm9, %ymm0, %ymm9
	vpmuludq	%ymm0, %ymm11, %ymm11
	vpsrlq	$32, %ymm0, %ymm1
	vpmuludq	32(%rsp), %ymm1, %ymm1
	vpshufd	$245, %ymm9, %ymm2
	vpmuludq	%ymm12, %ymm9, %ymm9
	vpmuludq	%ymm12, %ymm2, %ymm2
	vpaddq	%ymm9, %ymm11, %ymm11
	vpaddq	%ymm2, %ymm1, %ymm0
	vpsrlq	$32, %ymm11, %ymm11
	vpblendd	$170, %ymm0, %ymm11, %ymm11
	vpaddd	.LC5(%rip), %ymm11, %ymm0
	vpsrad	$31, %ymm0, %ymm11
	vpand	%ymm12, %ymm11, %ymm11
	vpaddd	%ymm0, %ymm11, %ymm11
	cmpl	%r8d, %r13d
	jg	.L45
.L44:
	leaq	_ZN4poly7f_n_t_t10barr_mergeE(%rip), %rdi
	vzeroupper
	call	pthread_barrier_wait@PLT
	sall	$2, %ebx
	vmovdqa	.LC1(%rip), %ymm12
	vmovdqa	.LC6(%rip), %ymm14
	vmovdqa	.LC7(%rip), %ymm15
	sall	$2, %r12d
	jmp	.L42
	.cfi_endproc
.LFE8967:
	.size	_ZN4poly7f_n_t_t8dit_funcEPv, .-_ZN4poly7f_n_t_t8dit_funcEPv
	.section	.rodata._ZStplIcSt11char_traitsIcESaIcEENSt7__cxx1112basic_stringIT_T0_T1_EEOS8_PKS5_.str1.1,"aMS",@progbits,1
.LC14:
	.string	"basic_string::append"
	.section	.text._ZStplIcSt11char_traitsIcESaIcEENSt7__cxx1112basic_stringIT_T0_T1_EEOS8_PKS5_,"axG",@progbits,_ZStplIcSt11char_traitsIcESaIcEENSt7__cxx1112basic_stringIT_T0_T1_EEOS8_PKS5_,comdat
	.p2align 4
	.weak	_ZStplIcSt11char_traitsIcESaIcEENSt7__cxx1112basic_stringIT_T0_T1_EEOS8_PKS5_
	.type	_ZStplIcSt11char_traitsIcESaIcEENSt7__cxx1112basic_stringIT_T0_T1_EEOS8_PKS5_, @function
_ZStplIcSt11char_traitsIcESaIcEENSt7__cxx1112basic_stringIT_T0_T1_EEOS8_PKS5_:
.LFB9288:
	.cfi_startproc
	endbr64
	pushq	%r13
	.cfi_def_cfa_offset 16
	.cfi_offset 13, -16
	movq	%rdx, %r13
	pushq	%r12
	.cfi_def_cfa_offset 24
	.cfi_offset 12, -24
	movq	%rdi, %r12
	movq	%rdx, %rdi
	pushq	%rbp
	.cfi_def_cfa_offset 32
	.cfi_offset 6, -32
	movq	%rsi, %rbp
	call	strlen@PLT
	movq	%rax, %rdx
	movabsq	$4611686018427387903, %rax
	subq	8(%rbp), %rax
	cmpq	%rax, %rdx
	ja	.L54
	movq	%r13, %rsi
	movq	%rbp, %rdi
	call	_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEE9_M_appendEPKcm@PLT
	leaq	16(%r12), %rdx
	movq	%rdx, (%r12)
	movq	(%rax), %rcx
	leaq	16(%rax), %rdx
	cmpq	%rdx, %rcx
	je	.L55
	movq	%rcx, (%r12)
	movq	16(%rax), %rcx
	movq	%rcx, 16(%r12)
.L52:
	movq	8(%rax), %rcx
	movq	%rdx, (%rax)
	movq	$0, 8(%rax)
	movq	%rcx, 8(%r12)
	movb	$0, 16(%rax)
	movq	%r12, %rax
	popq	%rbp
	.cfi_remember_state
	.cfi_def_cfa_offset 24
	popq	%r12
	.cfi_def_cfa_offset 16
	popq	%r13
	.cfi_def_cfa_offset 8
	ret
	.p2align 4,,10
	.p2align 3
.L55:
	.cfi_restore_state
	vmovdqu	16(%rax), %xmm0
	vmovups	%xmm0, 16(%r12)
	jmp	.L52
.L54:
	leaq	.LC14(%rip), %rdi
	call	_ZSt20__throw_length_errorPKc@PLT
	.cfi_endproc
.LFE9288:
	.size	_ZStplIcSt11char_traitsIcESaIcEENSt7__cxx1112basic_stringIT_T0_T1_EEOS8_PKS5_, .-_ZStplIcSt11char_traitsIcESaIcEENSt7__cxx1112basic_stringIT_T0_T1_EEOS8_PKS5_
	.section	.text._ZN4poly7f_n_t_t11dif_base4x8ILb0ELi1EEEvPji,"axG",@progbits,_ZN4poly7f_n_t_t11dif_base4x8ILb0ELi1EEEvPji,comdat
	.p2align 4
	.weak	_ZN4poly7f_n_t_t11dif_base4x8ILb0ELi1EEEvPji
	.type	_ZN4poly7f_n_t_t11dif_base4x8ILb0ELi1EEEvPji, @function
_ZN4poly7f_n_t_t11dif_base4x8ILb0ELi1EEEvPji:
.LFB9383:
	.cfi_startproc
	endbr64
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	pushq	%r13
	.cfi_offset 13, -24
	movq	%rdi, %r13
	pushq	%r12
	.cfi_offset 12, -32
	movl	%esi, %r12d
	pushq	%rbx
	.cfi_offset 3, -40
	movl	%esi, %ebx
	sarl	$4, %ebx
	andq	$-32, %rsp
	subq	$32, %rsp
	movq	%fs:40, %rax
	movq	%rax, 24(%rsp)
	xorl	%eax, %eax
	movl	%esi, %eax
	sarl	$3, %eax
	tzcntl	%eax, %eax
	testb	$1, %al
	je	.L57
	testl	%ebx, %ebx
	jle	.L58
	leal	-1(%rbx), %edx
	vmovdqa	.LC7(%rip), %ymm3
	movq	%rdi, %rax
	salq	$5, %rdx
	leaq	32(%rdi,%rdx), %rcx
	movslq	%ebx, %rdx
	salq	$5, %rdx
	.p2align 4,,10
	.p2align 3
.L59:
	vmovdqa	(%rax), %ymm0
	vmovdqa	(%rax,%rdx), %ymm1
	addq	$32, %rax
	vpaddd	%ymm1, %ymm0, %ymm2
	vpaddd	%ymm3, %ymm0, %ymm0
	vpsubd	%ymm1, %ymm0, %ymm0
	vmovdqa	%ymm2, -32(%rax)
	vmovdqa	%ymm0, -32(%rax,%rdx)
	cmpq	%rcx, %rax
	jne	.L59
	vzeroupper
.L58:
	movl	%r12d, %ebx
	sarl	$5, %ebx
.L57:
	leaq	_ZN4poly7f_n_t_t10barr_mergeE(%rip), %rdi
	movl	$1, %edx
	xorl	%esi, %esi
	sarl	%ebx
	call	pthread_barrier_init@PLT
	movq	%rsp, %rdi
	movq	%r13, 8(%rsp)
	movl	$0, (%rsp)
	movl	%r12d, 16(%rsp)
	movl	%ebx, 20(%rsp)
	call	_ZN4poly7f_n_t_t8dif_funcEPv
	.cfi_endproc
.LFE9383:
	.size	_ZN4poly7f_n_t_t11dif_base4x8ILb0ELi1EEEvPji, .-_ZN4poly7f_n_t_t11dif_base4x8ILb0ELi1EEEvPji
	.section	.rodata._ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEE12_M_constructIPcEEvT_S7_St20forward_iterator_tag.str1.8,"aMS",@progbits,1
	.align 8
.LC15:
	.string	"basic_string::_M_construct null not valid"
	.section	.text._ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEE12_M_constructIPcEEvT_S7_St20forward_iterator_tag,"axG",@progbits,_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEE12_M_constructIPcEEvT_S7_St20forward_iterator_tag,comdat
	.align 2
	.p2align 4
	.weak	_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEE12_M_constructIPcEEvT_S7_St20forward_iterator_tag
	.type	_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEE12_M_constructIPcEEvT_S7_St20forward_iterator_tag, @function
_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEE12_M_constructIPcEEvT_S7_St20forward_iterator_tag:
.LFB9495:
	.cfi_startproc
	endbr64
	pushq	%r12
	.cfi_def_cfa_offset 16
	.cfi_offset 12, -16
	movq	%rdx, %r12
	pushq	%rbp
	.cfi_def_cfa_offset 24
	.cfi_offset 6, -24
	movq	%rsi, %rbp
	pushq	%rbx
	.cfi_def_cfa_offset 32
	.cfi_offset 3, -32
	movq	%rdi, %rbx
	subq	$16, %rsp
	.cfi_def_cfa_offset 48
	movq	%fs:40, %rax
	movq	%rax, 8(%rsp)
	xorl	%eax, %eax
	testq	%rdx, %rdx
	je	.L67
	testq	%rsi, %rsi
	je	.L84
.L67:
	subq	%rbp, %r12
	movq	%r12, (%rsp)
	cmpq	$15, %r12
	ja	.L85
	movq	(%rbx), %rdi
	cmpq	$1, %r12
	jne	.L70
	movzbl	0(%rbp), %eax
	movb	%al, (%rdi)
	movq	(%rsp), %r12
	movq	(%rbx), %rdi
.L71:
	movq	%r12, 8(%rbx)
	movb	$0, (%rdi,%r12)
	movq	8(%rsp), %rax
	xorq	%fs:40, %rax
	jne	.L86
	addq	$16, %rsp
	.cfi_remember_state
	.cfi_def_cfa_offset 32
	popq	%rbx
	.cfi_def_cfa_offset 24
	popq	%rbp
	.cfi_def_cfa_offset 16
	popq	%r12
	.cfi_def_cfa_offset 8
	ret
	.p2align 4,,10
	.p2align 3
.L70:
	.cfi_restore_state
	testq	%r12, %r12
	je	.L71
	jmp	.L69
	.p2align 4,,10
	.p2align 3
.L85:
	movq	%rbx, %rdi
	movq	%rsp, %rsi
	xorl	%edx, %edx
	call	_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEE9_M_createERmm@PLT
	movq	%rax, (%rbx)
	movq	%rax, %rdi
	movq	(%rsp), %rax
	movq	%rax, 16(%rbx)
.L69:
	movq	%r12, %rdx
	movq	%rbp, %rsi
	call	memcpy@PLT
	movq	(%rsp), %r12
	movq	(%rbx), %rdi
	jmp	.L71
.L84:
	leaq	.LC15(%rip), %rdi
	call	_ZSt19__throw_logic_errorPKc@PLT
.L86:
	call	__stack_chk_fail@PLT
	.cfi_endproc
.LFE9495:
	.size	_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEE12_M_constructIPcEEvT_S7_St20forward_iterator_tag, .-_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEE12_M_constructIPcEEvT_S7_St20forward_iterator_tag
	.section	.text._ZN9__gnu_cxx12__to_xstringINSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEEcEET_PFiPT0_mPKS8_P13__va_list_tagEmSB_z,"axG",@progbits,_ZN9__gnu_cxx12__to_xstringINSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEEcEET_PFiPT0_mPKS8_P13__va_list_tagEmSB_z,comdat
	.p2align 4
	.weak	_ZN9__gnu_cxx12__to_xstringINSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEEcEET_PFiPT0_mPKS8_P13__va_list_tagEmSB_z
	.type	_ZN9__gnu_cxx12__to_xstringINSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEEcEET_PFiPT0_mPKS8_P13__va_list_tagEmSB_z, @function
_ZN9__gnu_cxx12__to_xstringINSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEEcEET_PFiPT0_mPKS8_P13__va_list_tagEmSB_z:
.LFB9097:
	.cfi_startproc
	endbr64
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsi, %r10
	movq	%rdx, %rsi
	movq	%rcx, %rdx
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	pushq	%r13
	pushq	%r12
	.cfi_offset 13, -24
	.cfi_offset 12, -32
	movq	%rdi, %r12
	subq	$208, %rsp
	movq	%r8, -160(%rbp)
	movq	%r9, -152(%rbp)
	testb	%al, %al
	je	.L88
	vmovaps	%xmm0, -144(%rbp)
	vmovaps	%xmm1, -128(%rbp)
	vmovaps	%xmm2, -112(%rbp)
	vmovaps	%xmm3, -96(%rbp)
	vmovaps	%xmm4, -80(%rbp)
	vmovaps	%xmm5, -64(%rbp)
	vmovaps	%xmm6, -48(%rbp)
	vmovaps	%xmm7, -32(%rbp)
.L88:
	movq	%fs:40, %rax
	movq	%rax, -200(%rbp)
	xorl	%eax, %eax
	leaq	39(%rsi), %rax
	movq	%rsp, %rdi
	movq	%rax, %rcx
	andq	$-4096, %rax
	subq	%rax, %rdi
	andq	$-16, %rcx
	movq	%rdi, %rax
	cmpq	%rax, %rsp
	je	.L90
.L97:
	subq	$4096, %rsp
	orq	$0, 4088(%rsp)
	cmpq	%rax, %rsp
	jne	.L97
.L90:
	andl	$4095, %ecx
	subq	%rcx, %rsp
	testq	%rcx, %rcx
	jne	.L98
.L91:
	leaq	16(%rbp), %rax
	leaq	31(%rsp), %r13
	movl	$32, -224(%rbp)
	andq	$-32, %r13
	movq	%rax, -216(%rbp)
	leaq	-192(%rbp), %rax
	leaq	-224(%rbp), %rcx
	movq	%r13, %rdi
	movq	%rax, -208(%rbp)
	movl	$48, -220(%rbp)
	call	*%r10
	leaq	16(%r12), %rdx
	movq	%r13, %rsi
	movq	%r12, %rdi
	movq	%rdx, (%r12)
	movslq	%eax, %rdx
	addq	%r13, %rdx
	call	_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEE12_M_constructIPcEEvT_S7_St20forward_iterator_tag
	movq	-200(%rbp), %rax
	xorq	%fs:40, %rax
	jne	.L99
	leaq	-16(%rbp), %rsp
	movq	%r12, %rax
	popq	%r12
	popq	%r13
	popq	%rbp
	.cfi_remember_state
	.cfi_def_cfa 7, 8
	ret
	.p2align 4,,10
	.p2align 3
.L98:
	.cfi_restore_state
	orq	$0, -8(%rsp,%rcx)
	jmp	.L91
.L99:
	call	__stack_chk_fail@PLT
	.cfi_endproc
.LFE9097:
	.size	_ZN9__gnu_cxx12__to_xstringINSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEEcEET_PFiPT0_mPKS8_P13__va_list_tagEmSB_z, .-_ZN9__gnu_cxx12__to_xstringINSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEEcEET_PFiPT0_mPKS8_P13__va_list_tagEmSB_z
	.section	.rodata.str1.1
.LC16:
	.string	"%d"
.LC17:
	.string	"_"
.LC18:
	.string	".in"
.LC19:
	.string	"r"
	.section	.text.unlikely,"ax",@progbits
.LCOLDB21:
	.section	.text.startup,"ax",@progbits
.LHOTB21:
	.p2align 4
	.globl	main
	.type	main, @function
main:
.LFB8998:
	.cfi_startproc
	.cfi_personality 0x9b,DW.ref.__gxx_personality_v0
	.cfi_lsda 0x1b,.LLSDA8998
	endbr64
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	xorl	%edx, %edx
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	pushq	%r15
	pushq	%r14
	leaq	-288(%rbp), %rdi
	leaq	-296(%rbp), %rsi
	pushq	%r13
	pushq	%r12
	pushq	%rbx
	subq	$312, %rsp
	.cfi_offset 15, -24
	.cfi_offset 14, -32
	.cfi_offset 13, -40
	.cfi_offset 12, -48
	.cfi_offset 3, -56
	movq	%fs:40, %rax
	movq	%rax, -56(%rbp)
	xorl	%eax, %eax
	leaq	-272(%rbp), %rax
	movq	$18, -296(%rbp)
	movq	%rax, -344(%rbp)
	movq	%rax, -288(%rbp)
.LEHB0:
	call	_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEE9_M_createERmm@PLT
.LEHE0:
	movq	-296(%rbp), %rdx
	vmovdqa	.LC20(%rip), %xmm0
	leaq	-256(%rbp), %rdi
	movq	%rax, -288(%rbp)
	movq	vsnprintf@GOTPCREL(%rip), %r12
	leaq	.LC16(%rip), %rcx
	movq	%rdx, -272(%rbp)
	movl	$21588, %edx
	vmovups	%xmm0, (%rax)
	movq	%r12, %rsi
	movw	%dx, 16(%rax)
	movq	-296(%rbp), %rax
	movq	-288(%rbp), %rdx
	movq	%rax, -280(%rbp)
	movb	$0, (%rdx,%rax)
	movl	nn(%rip), %r8d
	movl	$16, %edx
	xorl	%eax, %eax
.LEHB1:
	call	_ZN9__gnu_cxx12__to_xstringINSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEEcEET_PFiPT0_mPKS8_P13__va_list_tagEmSB_z
.LEHE1:
	movl	qq(%rip), %r8d
	movl	$16, %edx
	movq	%r12, %rsi
	xorl	%eax, %eax
	leaq	-224(%rbp), %rdi
	leaq	.LC16(%rip), %rcx
.LEHB2:
	call	_ZN9__gnu_cxx12__to_xstringINSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEEcEET_PFiPT0_mPKS8_P13__va_list_tagEmSB_z
.LEHE2:
	movq	-288(%rbp), %rsi
	movq	-280(%rbp), %rdx
	leaq	-160(%rbp), %r12
	leaq	-144(%rbp), %rbx
	movq	%r12, %rdi
	movq	%rbx, -160(%rbp)
	addq	%rsi, %rdx
.LEHB3:
	call	_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEE12_M_constructIPcEEvT_S7_St20forward_iterator_tag
.LEHE3:
	movq	-248(%rbp), %rdx
	movq	-256(%rbp), %rsi
	movq	%r12, %rdi
.LEHB4:
	call	_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEE9_M_appendEPKcm@PLT
.LEHE4:
	leaq	-128(%rbp), %r13
	leaq	.LC17(%rip), %rdx
	movq	%r12, %rsi
	movq	%r13, %rdi
.LEHB5:
	call	_ZStplIcSt11char_traitsIcESaIcEENSt7__cxx1112basic_stringIT_T0_T1_EEOS8_PKS5_
.LEHE5:
	movq	-216(%rbp), %rdx
	movq	-224(%rbp), %rsi
	movq	%r13, %rdi
.LEHB6:
	call	_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEE9_M_appendEPKcm@PLT
.LEHE6:
	leaq	-80(%rbp), %r12
	leaq	16(%rax), %rdx
	movq	%r12, -96(%rbp)
	movq	(%rax), %rcx
	cmpq	%rdx, %rcx
	je	.L158
	movq	%rcx, -96(%rbp)
	movq	16(%rax), %rcx
	movq	%rcx, -80(%rbp)
.L106:
	movq	8(%rax), %rcx
	movb	$0, 16(%rax)
	leaq	-192(%rbp), %rdi
	leaq	-96(%rbp), %rsi
	movq	%rcx, -88(%rbp)
	movq	%rdx, (%rax)
	leaq	.LC18(%rip), %rdx
	movq	$0, 8(%rax)
.LEHB7:
	call	_ZStplIcSt11char_traitsIcESaIcEENSt7__cxx1112basic_stringIT_T0_T1_EEOS8_PKS5_
.LEHE7:
	movq	-96(%rbp), %rdi
	cmpq	%r12, %rdi
	je	.L107
	call	_ZdlPv@PLT
.L107:
	movq	-128(%rbp), %rdi
	leaq	-112(%rbp), %rax
	cmpq	%rax, %rdi
	je	.L108
	call	_ZdlPv@PLT
.L108:
	movq	-160(%rbp), %rdi
	cmpq	%rbx, %rdi
	je	.L109
	call	_ZdlPv@PLT
.L109:
	movq	-184(%rbp), %rbx
	movq	%rsp, %rcx
	leaq	16(%rbx), %rax
	movq	%rax, %rdx
	andq	$-4096, %rax
	subq	%rax, %rcx
	andq	$-16, %rdx
	movq	%rcx, %rax
	cmpq	%rax, %rsp
	je	.L111
.L159:
	subq	$4096, %rsp
	orq	$0, 4088(%rsp)
	cmpq	%rax, %rsp
	jne	.L159
.L111:
	andl	$4095, %edx
	subq	%rdx, %rsp
	testq	%rdx, %rdx
	jne	.L160
.L112:
	movq	%rsp, %r12
	testq	%rbx, %rbx
	je	.L113
	movq	-192(%rbp), %rsi
	movq	%rbx, %rdx
	movq	%rsp, %rdi
	call	memmove@PLT
.L113:
	movb	$0, (%r12,%rbx)
	movq	stdin(%rip), %rdx
	leaq	.LC19(%rip), %rsi
	movq	%r12, %rdi
.LEHB8:
	call	freopen@PLT
	movl	nn(%rip), %ecx
	movl	$32, %eax
	leal	-2(%rcx,%rcx), %edx
	movl	%ecx, -324(%rbp)
	lzcntl	%edx, %edx
	subl	%edx, %eax
	movl	$1, %edx
	shlx	%eax, %edx, %esi
	movq	_ZN4poly9poly_base10mem_helper3nowE(%rip), %rax
	movslq	%esi, %r12
	movl	%esi, -328(%rbp)
	salq	$2, %r12
	addq	$31, %rax
	andq	$-32, %rax
	leaq	(%rax,%r12), %r13
	movq	%rax, -320(%rbp)
	leaq	31(%r13), %rax
	andq	$-32, %rax
	addq	%rax, %r12
	movq	%rax, -312(%rbp)
	movq	%r12, _ZN4poly9poly_base10mem_helper3nowE(%rip)
	testl	%ecx, %ecx
	jle	.L114
	leal	-1(%rcx), %r14d
	leaq	f(%rip), %rbx
	leaq	0(,%r14,4), %rax
	leaq	4(%rbx), %r14
	movq	%rax, -336(%rbp)
	addq	%rax, %r14
	leaq	_ZSt3cin(%rip), %r15
	.p2align 4,,10
	.p2align 3
.L115:
	movq	%rbx, %rsi
	movq	%r15, %rdi
	call	_ZNSirsERi@PLT
	addq	$4, %rbx
	cmpq	%rbx, %r14
	jne	.L115
	leaq	g(%rip), %rbx
	leaq	_ZSt3cin(%rip), %r15
	leaq	4(%rbx), %r14
	addq	-336(%rbp), %r14
	.p2align 4,,10
	.p2align 3
.L116:
	movq	%rbx, %rsi
	movq	%r15, %rdi
	call	_ZNSirsERi@PLT
	addq	$4, %rbx
	cmpq	%rbx, %r14
	jne	.L116
.L114:
	movq	-320(%rbp), %rbx
	movl	$1600020, %edx
	leaq	f(%rip), %rsi
	movq	%rbx, %rdi
	call	memcpy@PLT
	movq	-312(%rbp), %rdi
	movl	$1600020, %edx
	leaq	g(%rip), %rsi
	call	memcpy@PLT
	call	_ZNSt6chrono3_V212system_clock3nowEv@PLT
	movslq	-324(%rbp), %rax
	salq	$2, %rax
	leaq	(%rbx,%rax), %rdx
	cmpq	%r13, %rdx
	je	.L117
	.p2align 4,,10
	.p2align 3
.L118:
	movl	$0, (%rdx)
	addq	$4, %rdx
	cmpq	%rdx, %r13
	jne	.L118
.L117:
	addq	-312(%rbp), %rax
	cmpq	%r12, %rax
	je	.L119
	.p2align 4,,10
	.p2align 3
.L120:
	movl	$0, (%rax)
	addq	$4, %rax
	cmpq	%rax, %r12
	jne	.L120
.L119:
	movl	-328(%rbp), %esi
	movq	-320(%rbp), %rdi
	call	_ZN4poly7f_n_t_t11dif_base4x8ILb0ELi1EEEvPji
.LEHE8:
.L158:
	vmovdqu	16(%rax), %xmm1
	vmovaps	%xmm1, -80(%rbp)
	jmp	.L106
.L160:
	orq	$0, -8(%rsp,%rdx)
	jmp	.L112
.L135:
	endbr64
	movq	%rax, %r12
	vzeroupper
	jmp	.L132
.L141:
	endbr64
	movq	%rax, %rbx
	jmp	.L127
.L138:
	endbr64
	movq	%rax, %r12
	vzeroupper
	jmp	.L125
.L142:
	endbr64
	movq	%rax, %r12
	jmp	.L102
.L140:
	endbr64
	movq	%rax, %r13
	jmp	.L121
.L139:
	endbr64
	movq	%rax, %r12
	vzeroupper
	jmp	.L123
.L137:
	endbr64
	movq	%rax, %rbx
	vzeroupper
	jmp	.L104
.L136:
	endbr64
	movq	%rax, %rbx
	vzeroupper
	jmp	.L130
	.globl	__gxx_personality_v0
	.section	.gcc_except_table,"a",@progbits
.LLSDA8998:
	.byte	0xff
	.byte	0xff
	.byte	0x1
	.uleb128 .LLSDACSE8998-.LLSDACSB8998
.LLSDACSB8998:
	.uleb128 .LEHB0-.LFB8998
	.uleb128 .LEHE0-.LEHB0
	.uleb128 0
	.uleb128 0
	.uleb128 .LEHB1-.LFB8998
	.uleb128 .LEHE1-.LEHB1
	.uleb128 .L135-.LFB8998
	.uleb128 0
	.uleb128 .LEHB2-.LFB8998
	.uleb128 .LEHE2-.LEHB2
	.uleb128 .L136-.LFB8998
	.uleb128 0
	.uleb128 .LEHB3-.LFB8998
	.uleb128 .LEHE3-.LEHB3
	.uleb128 .L137-.LFB8998
	.uleb128 0
	.uleb128 .LEHB4-.LFB8998
	.uleb128 .LEHE4-.LEHB4
	.uleb128 .L142-.LFB8998
	.uleb128 0
	.uleb128 .LEHB5-.LFB8998
	.uleb128 .LEHE5-.LEHB5
	.uleb128 .L138-.LFB8998
	.uleb128 0
	.uleb128 .LEHB6-.LFB8998
	.uleb128 .LEHE6-.LEHB6
	.uleb128 .L139-.LFB8998
	.uleb128 0
	.uleb128 .LEHB7-.LFB8998
	.uleb128 .LEHE7-.LEHB7
	.uleb128 .L140-.LFB8998
	.uleb128 0
	.uleb128 .LEHB8-.LFB8998
	.uleb128 .LEHE8-.LEHB8
	.uleb128 .L141-.LFB8998
	.uleb128 0
.LLSDACSE8998:
	.section	.text.startup
	.cfi_endproc
	.section	.text.unlikely
	.cfi_startproc
	.cfi_personality 0x9b,DW.ref.__gxx_personality_v0
	.cfi_lsda 0x1b,.LLSDAC8998
	.type	main.cold, @function
main.cold:
.LFSB8998:
.L102:
	.cfi_def_cfa 6, 16
	.cfi_offset 3, -56
	.cfi_offset 6, -16
	.cfi_offset 12, -48
	.cfi_offset 13, -40
	.cfi_offset 14, -32
	.cfi_offset 15, -24
	movq	-160(%rbp), %rdi
	cmpq	%rbx, %rdi
	je	.L154
	vzeroupper
	call	_ZdlPv@PLT
.L103:
	movq	%r12, %rbx
.L104:
	movq	-224(%rbp), %rdi
	leaq	-208(%rbp), %rax
	cmpq	%rax, %rdi
	je	.L130
	call	_ZdlPv@PLT
.L130:
	movq	-256(%rbp), %rdi
	leaq	-240(%rbp), %rax
	cmpq	%rax, %rdi
	je	.L131
	call	_ZdlPv@PLT
.L131:
	movq	%rbx, %r12
.L132:
	movq	-288(%rbp), %rdi
	cmpq	-344(%rbp), %rdi
	je	.L133
	call	_ZdlPv@PLT
.L133:
	movq	%r12, %rdi
.LEHB9:
	call	_Unwind_Resume@PLT
.LEHE9:
.L127:
	movq	-192(%rbp), %rdi
	leaq	-176(%rbp), %rax
	cmpq	%rax, %rdi
	je	.L156
	vzeroupper
	call	_ZdlPv@PLT
	jmp	.L104
.L121:
	movq	-96(%rbp), %rdi
	cmpq	%r12, %rdi
	je	.L155
	vzeroupper
	call	_ZdlPv@PLT
.L122:
	movq	%r13, %r12
.L123:
	movq	-128(%rbp), %rdi
	leaq	-112(%rbp), %rax
	cmpq	%rax, %rdi
	je	.L125
	call	_ZdlPv@PLT
.L125:
	movq	-160(%rbp), %rdi
	cmpq	%rbx, %rdi
	je	.L126
	call	_ZdlPv@PLT
.L126:
	movq	%r12, %rbx
	jmp	.L104
.L155:
	vzeroupper
	jmp	.L122
.L154:
	vzeroupper
	jmp	.L103
.L156:
	vzeroupper
	jmp	.L104
	.cfi_endproc
.LFE8998:
	.section	.gcc_except_table
.LLSDAC8998:
	.byte	0xff
	.byte	0xff
	.byte	0x1
	.uleb128 .LLSDACSEC8998-.LLSDACSBC8998
.LLSDACSBC8998:
	.uleb128 .LEHB9-.LCOLDB21
	.uleb128 .LEHE9-.LEHB9
	.uleb128 0
	.uleb128 0
.LLSDACSEC8998:
	.section	.text.unlikely
	.section	.text.startup
	.size	main, .-main
	.section	.text.unlikely
	.size	main.cold, .-main.cold
.LCOLDE21:
	.section	.text.startup
.LHOTE21:
	.p2align 4
	.type	_GLOBAL__sub_I__ZN4poly9poly_base10mem_helper4_memE, @function
_GLOBAL__sub_I__ZN4poly9poly_base10mem_helper4_memE:
.LFB9565:
	.cfi_startproc
	endbr64
	subq	$8, %rsp
	.cfi_def_cfa_offset 16
	leaq	_ZStL8__ioinit(%rip), %rdi
	call	_ZNSt8ios_base4InitC1Ev@PLT
	movq	_ZNSt8ios_base4InitD1Ev@GOTPCREL(%rip), %rdi
	addq	$8, %rsp
	.cfi_def_cfa_offset 8
	leaq	__dso_handle(%rip), %rdx
	leaq	_ZStL8__ioinit(%rip), %rsi
	jmp	__cxa_atexit@PLT
	.cfi_endproc
.LFE9565:
	.size	_GLOBAL__sub_I__ZN4poly9poly_base10mem_helper4_memE, .-_GLOBAL__sub_I__ZN4poly9poly_base10mem_helper4_memE
	.section	.init_array,"aw"
	.align 8
	.quad	_GLOBAL__sub_I__ZN4poly9poly_base10mem_helper4_memE
	.globl	g
	.bss
	.align 32
	.type	g, @object
	.size	g, 1600020
g:
	.zero	1600020
	.globl	f
	.align 32
	.type	f, @object
	.size	f, 1600020
f:
	.zero	1600020
	.globl	qq
	.data
	.align 16
	.type	qq, @object
	.size	qq, 20
qq:
	.long	1409
	.long	3329
	.long	7681
	.long	12289
	.zero	4
	.globl	nn
	.align 32
	.type	nn, @object
	.size	nn, 44
nn:
	.long	256
	.long	512
	.long	1024
	.long	2048
	.long	4096
	.long	8192
	.long	16384
	.long	32768
	.long	65536
	.long	131072
	.zero	4
	.globl	_ZN4poly7f_n_t_t10barr_mergeE
	.bss
	.align 32
	.type	_ZN4poly7f_n_t_t10barr_mergeE, @object
	.size	_ZN4poly7f_n_t_t10barr_mergeE, 32
_ZN4poly7f_n_t_t10barr_mergeE:
	.zero	32
	.section	.rodata
	.align 32
	.type	_ZN4poly7f_n_t_tL4iab4E, @object
	.size	_ZN4poly7f_n_t_tL4iab4E, 2624
_ZN4poly7f_n_t_tL4iab4E:
	.long	81453865
	.long	81453865
	.long	81453865
	.long	81453865
	.long	81453865
	.long	81453865
	.long	81453865
	.long	81453865
	.long	1975273
	.long	1975273
	.long	1975273
	.long	1975273
	.long	1975273
	.long	1975273
	.long	1975273
	.long	1975273
	.long	29996215
	.long	29996215
	.long	29996215
	.long	29996215
	.long	29996215
	.long	29996215
	.long	29996215
	.long	29996215
	.long	104605285
	.long	104605285
	.long	104605285
	.long	104605285
	.long	104605285
	.long	104605285
	.long	104605285
	.long	104605285
	.long	72830425
	.long	72830425
	.long	72830425
	.long	72830425
	.long	72830425
	.long	72830425
	.long	72830425
	.long	72830425
	.long	99602582
	.long	99602582
	.long	99602582
	.long	99602582
	.long	99602582
	.long	99602582
	.long	99602582
	.long	99602582
	.long	83501717
	.long	83501717
	.long	83501717
	.long	83501717
	.long	83501717
	.long	83501717
	.long	83501717
	.long	83501717
	.long	11851883
	.long	11851883
	.long	11851883
	.long	11851883
	.long	11851883
	.long	11851883
	.long	11851883
	.long	11851883
	.long	35821440
	.long	35821440
	.long	35821440
	.long	35821440
	.long	35821440
	.long	35821440
	.long	35821440
	.long	35821440
	.long	58009840
	.long	58009840
	.long	58009840
	.long	58009840
	.long	58009840
	.long	58009840
	.long	58009840
	.long	58009840
	.long	65601194
	.long	65601194
	.long	65601194
	.long	65601194
	.long	65601194
	.long	65601194
	.long	65601194
	.long	65601194
	.long	67066948
	.long	67066948
	.long	67066948
	.long	67066948
	.long	67066948
	.long	67066948
	.long	67066948
	.long	67066948
	.long	31926659
	.long	31926659
	.long	31926659
	.long	31926659
	.long	31926659
	.long	31926659
	.long	31926659
	.long	31926659
	.long	60997406
	.long	60997406
	.long	60997406
	.long	60997406
	.long	60997406
	.long	60997406
	.long	60997406
	.long	60997406
	.long	36020067
	.long	36020067
	.long	36020067
	.long	36020067
	.long	36020067
	.long	36020067
	.long	36020067
	.long	36020067
	.long	74132622
	.long	74132622
	.long	74132622
	.long	74132622
	.long	74132622
	.long	74132622
	.long	74132622
	.long	74132622
	.long	104432633
	.long	104432633
	.long	104432633
	.long	104432633
	.long	104432633
	.long	104432633
	.long	104432633
	.long	104432633
	.long	53788560
	.long	53788560
	.long	53788560
	.long	53788560
	.long	53788560
	.long	53788560
	.long	53788560
	.long	53788560
	.long	60091236
	.long	60091236
	.long	60091236
	.long	60091236
	.long	60091236
	.long	60091236
	.long	60091236
	.long	60091236
	.long	70190262
	.long	70190262
	.long	70190262
	.long	70190262
	.long	70190262
	.long	70190262
	.long	70190262
	.long	70190262
	.long	50219246
	.long	50219246
	.long	50219246
	.long	50219246
	.long	50219246
	.long	50219246
	.long	50219246
	.long	50219246
	.long	79977691
	.long	79977691
	.long	79977691
	.long	79977691
	.long	79977691
	.long	79977691
	.long	79977691
	.long	79977691
	.long	76300401
	.long	76300401
	.long	76300401
	.long	76300401
	.long	76300401
	.long	76300401
	.long	76300401
	.long	76300401
	.long	90888097
	.long	90888097
	.long	90888097
	.long	90888097
	.long	90888097
	.long	90888097
	.long	90888097
	.long	90888097
	.long	104083916
	.long	104083916
	.long	104083916
	.long	104083916
	.long	104083916
	.long	104083916
	.long	104083916
	.long	104083916
	.long	11909896
	.long	11909896
	.long	11909896
	.long	11909896
	.long	11909896
	.long	11909896
	.long	11909896
	.long	11909896
	.long	42634566
	.long	42634566
	.long	42634566
	.long	42634566
	.long	42634566
	.long	42634566
	.long	42634566
	.long	42634566
	.long	68478099
	.long	68478099
	.long	68478099
	.long	68478099
	.long	68478099
	.long	68478099
	.long	68478099
	.long	68478099
	.long	98926185
	.long	98926185
	.long	98926185
	.long	98926185
	.long	98926185
	.long	98926185
	.long	98926185
	.long	98926185
	.long	13535815
	.long	13535815
	.long	13535815
	.long	13535815
	.long	13535815
	.long	13535815
	.long	13535815
	.long	13535815
	.long	69699293
	.long	69699293
	.long	69699293
	.long	69699293
	.long	69699293
	.long	69699293
	.long	69699293
	.long	69699293
	.long	85408463
	.long	85408463
	.long	85408463
	.long	85408463
	.long	85408463
	.long	85408463
	.long	85408463
	.long	85408463
	.long	66243509
	.long	66243509
	.long	66243509
	.long	66243509
	.long	66243509
	.long	66243509
	.long	66243509
	.long	66243509
	.long	18321273
	.long	18321273
	.long	18321273
	.long	18321273
	.long	18321273
	.long	18321273
	.long	18321273
	.long	18321273
	.long	39369937
	.long	39369937
	.long	39369937
	.long	39369937
	.long	39369937
	.long	39369937
	.long	39369937
	.long	39369937
	.long	87560486
	.long	87560486
	.long	87560486
	.long	87560486
	.long	87560486
	.long	87560486
	.long	87560486
	.long	87560486
	.long	92941318
	.long	92941318
	.long	92941318
	.long	92941318
	.long	92941318
	.long	92941318
	.long	92941318
	.long	92941318
	.long	81319803
	.long	81319803
	.long	81319803
	.long	81319803
	.long	81319803
	.long	81319803
	.long	81319803
	.long	81319803
	.long	62050207
	.long	62050207
	.long	62050207
	.long	62050207
	.long	62050207
	.long	62050207
	.long	62050207
	.long	62050207
	.long	34529057
	.long	34529057
	.long	34529057
	.long	34529057
	.long	34529057
	.long	34529057
	.long	34529057
	.long	34529057
	.long	100663256
	.long	79977691
	.long	81453865
	.long	10721473
	.long	63333991
	.long	71165571
	.long	54638355
	.long	102882328
	.long	100663256
	.long	56834138
	.long	1975273
	.long	22121189
	.long	50219246
	.long	29996215
	.long	33692030
	.long	90131399
	.long	100663256
	.long	79629103
	.long	29996215
	.long	58276067
	.long	94136128
	.long	96403364
	.long	93779359
	.long	21953196
	.long	100663256
	.long	45392609
	.long	104605285
	.long	51515697
	.long	33328271
	.long	54128192
	.long	58045831
	.long	11877514
	.long	100663256
	.long	10245096
	.long	72830425
	.long	18362680
	.long	37724185
	.long	13278476
	.long	90205732
	.long	68603196
	.long	100663256
	.long	3894551
	.long	99602582
	.long	1553345
	.long	36293688
	.long	40741404
	.long	11204615
	.long	69507338
	.long	100663256
	.long	96385508
	.long	83501717
	.long	93479544
	.long	85412354
	.long	33575408
	.long	45144822
	.long	34467861
	.long	100663256
	.long	89612531
	.long	11851883
	.long	69907200
	.long	48703526
	.long	65369202
	.long	82582648
	.long	73166213
	.long	100663256
	.long	20171564
	.long	35821440
	.long	2388175
	.long	43037563
	.long	100678316
	.long	20612882
	.long	24774428
	.long	100663256
	.long	45318479
	.long	58009840
	.long	30678380
	.long	19657302
	.long	25885624
	.long	18955697
	.long	80335646
	.long	100663256
	.long	39111804
	.long	65601194
	.long	21252078
	.long	2451935
	.long	76756754
	.long	79711751
	.long	104509294
	.long	100663256
	.long	100287395
	.long	67066948
	.long	85064342
	.long	38434554
	.long	90521152
	.long	24537050
	.long	72268654
	.long	100663256
	.long	18942077
	.long	31926659
	.long	1653825
	.long	53118571
	.long	44058305
	.long	34150184
	.long	42208922
	.long	100663256
	.long	13007420
	.long	60997406
	.long	75376589
	.long	87845843
	.long	17575122
	.long	94820295
	.long	29430512
	.long	100663256
	.long	78781843
	.long	36020067
	.long	46025930
	.long	81565884
	.long	84361750
	.long	101305130
	.long	489553
	.long	100663256
	.long	51532674
	.long	74132622
	.long	42318809
	.long	61303363
	.long	53705648
	.long	2438185
	.long	62212006
	.long	100663256
	.long	24803510
	.long	104432633
	.long	68221432
	.long	53875641
	.long	80857149
	.long	65358231
	.long	56916487
	.long	100663256
	.long	23367834
	.long	53788560
	.long	96247472
	.long	52346922
	.long	53534049
	.long	98467689
	.long	43531406
	.long	100663256
	.long	49585336
	.long	60091236
	.long	80242354
	.long	82733948
	.long	83264256
	.long	38231435
	.long	5036682
	.long	100663256
	.long	1975273
	.long	50219246
	.long	33692030
	.long	41523610
	.long	94136128
	.long	23403736
	.long	24879910
	.long	100663256
	.long	11078242
	.long	79977691
	.long	71529330
	.long	81453865
	.long	76300401
	.long	10721473
	.long	82970031
	.long	100663256
	.long	75520719
	.long	76300401
	.long	64087498
	.long	71165571
	.long	58045831
	.long	48023463
	.long	37724185
	.long	100663256
	.long	43431937
	.long	90888097
	.long	28572589
	.long	82736412
	.long	72717323
	.long	96403364
	.long	40877461
	.long	100663256
	.long	60967849
	.long	104083916
	.long	70642868
	.long	21953196
	.long	55025433
	.long	56297331
	.long	57610977
	.long	100663256
	.long	20570464
	.long	11909896
	.long	8939973
	.long	58214325
	.long	82022859
	.long	45581499
	.long	45637303
	.long	100663256
	.long	57709732
	.long	42634566
	.long	96804909
	.long	96912524
	.long	33163214
	.long	17315851
	.long	103256248
	.long	100663256
	.long	88979018
	.long	68478099
	.long	45550536
	.long	49094724
	.long	42003225
	.long	62545813
	.long	101179671
	.long	100663256
	.long	100570232
	.long	98926185
	.long	76111707
	.long	71901728
	.long	104182046
	.long	67430922
	.long	18974972
	.long	100663256
	.long	65346906
	.long	13535815
	.long	57841695
	.long	25148861
	.long	10225650
	.long	95088050
	.long	6261576
	.long	100663256
	.long	58744418
	.long	69699293
	.long	15022471
	.long	15002922
	.long	53977665
	.long	49941512
	.long	8818206
	.long	100663256
	.long	76508494
	.long	85408463
	.long	44838138
	.long	44929687
	.long	35965225
	.long	53372271
	.long	54268881
	.long	100663256
	.long	67056758
	.long	66243509
	.long	81281193
	.long	34588821
	.long	57910077
	.long	88672730
	.long	94439121
	.long	100663256
	.long	30987822
	.long	18321273
	.long	49381901
	.long	95146310
	.long	65454282
	.long	84825046
	.long	39844078
	.long	100663256
	.long	39669200
	.long	39369937
	.long	79491187
	.long	85270270
	.long	65604247
	.long	26661899
	.long	87768506
	.long	100663256
	.long	14279339
	.long	87560486
	.long	91918308
	.long	29320965
	.long	16773918
	.long	82990362
	.long	96721039
	.long	100663256
	.long	100181362
	.long	92941318
	.long	97720730
	.long	17969511
	.long	89210527
	.long	83561876
	.long	37615606
	.long	100663256
	.long	43369141
	.long	81319803
	.long	69041041
	.long	73347957
	.long	18290970
	.long	5924775
	.long	86109683
	.long	100663256
	.long	53724180
	.long	62050207
	.long	23875436
	.long	64284422
	.long	56630283
	.long	85425788
	.long	54571250
	.long	100663256
	.long	100663256
	.long	100663256
	.long	63333991
	.long	100663256
	.long	100663256
	.long	100663256
	.long	63333991
	.long	100663256
	.long	100663256
	.long	100663256
	.long	100663256
	.long	100663256
	.long	81453865
	.long	63333991
	.long	54638355
	.long	100663256
	.long	100663256
	.long	100663256
	.long	41523610
	.long	100663256
	.long	100663256
	.long	100663256
	.long	41523610
	.long	100663256
	.long	100663256
	.long	100663256
	.long	100663256
	.long	100663256
	.long	50219246
	.long	41523610
	.long	23403736
	.globl	_ZN4poly9poly_base10mem_helper3nowE
	.section	.data.rel.local,"aw"
	.align 8
	.type	_ZN4poly9poly_base10mem_helper3nowE, @object
	.size	_ZN4poly9poly_base10mem_helper3nowE, 8
_ZN4poly9poly_base10mem_helper3nowE:
	.quad	_ZN4poly9poly_base10mem_helper4_memE
	.globl	_ZN4poly9poly_base10mem_helper4_memE
	.bss
	.align 32
	.type	_ZN4poly9poly_base10mem_helper4_memE, @object
	.size	_ZN4poly9poly_base10mem_helper4_memE, 67108864
_ZN4poly9poly_base10mem_helper4_memE:
	.zero	67108864
	.local	_ZStL8__ioinit
	.comm	_ZStL8__ioinit,1,1
	.section	.rodata.cst32,"aM",@progbits,32
	.align 32
.LC0:
	.long	104857599
	.long	104857599
	.long	104857599
	.long	104857599
	.long	104857599
	.long	104857599
	.long	104857599
	.long	104857599
	.align 32
.LC1:
	.long	104857601
	.long	104857601
	.long	104857601
	.long	104857601
	.long	104857601
	.long	104857601
	.long	104857601
	.long	104857601
	.align 32
.LC2:
	.long	81453865
	.long	81453865
	.long	81453865
	.long	81453865
	.long	81453865
	.long	81453865
	.long	81453865
	.long	81453865
	.align 32
.LC3:
	.long	100663256
	.long	100663256
	.long	100663256
	.long	100663256
	.long	100663256
	.long	100663256
	.long	100663256
	.long	100663256
	.align 32
.LC5:
	.long	-104857601
	.long	-104857601
	.long	-104857601
	.long	-104857601
	.long	-104857601
	.long	-104857601
	.long	-104857601
	.long	-104857601
	.align 32
.LC6:
	.long	-209715202
	.long	-209715202
	.long	-209715202
	.long	-209715202
	.long	-209715202
	.long	-209715202
	.long	-209715202
	.long	-209715202
	.align 32
.LC7:
	.long	209715202
	.long	209715202
	.long	209715202
	.long	209715202
	.long	209715202
	.long	209715202
	.long	209715202
	.long	209715202
	.align 32
.LC8:
	.long	63333991
	.long	0
	.long	63333991
	.long	0
	.long	63333991
	.long	0
	.long	63333991
	.long	0
	.align 32
.LC9:
	.long	-419431
	.long	-419431
	.long	-419431
	.long	-419431
	.long	-419431
	.long	-419431
	.long	-419431
	.long	-419431
	.align 32
.LC10:
	.long	63333991
	.long	63333991
	.long	63333991
	.long	63333991
	.long	63333991
	.long	63333991
	.long	63333991
	.long	63333991
	.align 32
.LC11:
	.long	41523610
	.long	0
	.long	41523610
	.long	0
	.long	41523610
	.long	0
	.long	41523610
	.long	0
	.align 32
.LC12:
	.long	419430
	.long	419430
	.long	419430
	.long	419430
	.long	419430
	.long	419430
	.long	419430
	.long	419430
	.align 32
.LC13:
	.long	41523610
	.long	41523610
	.long	41523610
	.long	41523610
	.long	41523610
	.long	41523610
	.long	41523610
	.long	41523610
	.section	.rodata.cst16,"aM",@progbits,16
	.align 16
.LC20:
	.quad	3414137954165009966
	.quad	5633843057749418340
	.hidden	DW.ref.__gxx_personality_v0
	.weak	DW.ref.__gxx_personality_v0
	.section	.data.rel.local.DW.ref.__gxx_personality_v0,"awG",@progbits,DW.ref.__gxx_personality_v0,comdat
	.align 8
	.type	DW.ref.__gxx_personality_v0, @object
	.size	DW.ref.__gxx_personality_v0, 8
DW.ref.__gxx_personality_v0:
	.quad	__gxx_personality_v0
	.hidden	__dso_handle
	.ident	"GCC: (Ubuntu 9.5.0-1ubuntu1~22.04) 9.5.0"
	.section	.note.GNU-stack,"",@progbits
	.section	.note.gnu.property,"a"
	.align 8
	.long	 1f - 0f
	.long	 4f - 1f
	.long	 5
0:
	.string	 "GNU"
1:
	.align 8
	.long	 0xc0000002
	.long	 3f - 2f
2:
	.long	 0x3
3:
	.align 8
4:
