	.file	"spqlios-fft-avx.s"
	.section	.text.unlikely,"ax",@progbits
.LCOLDB0:
	.text
.LHOTB0:
	.p2align 4,,15
	.globl	fft
	.type	fft, @function
fft:
.LFB0:
	.cfi_startproc

	pushq       %r10
	pushq       %r11
	pushq       %r12
	pushq       %r13
	pushq       %r14
	pushq       %rbx

	movq        %rdi, %rax
	movq        %rsi, %rdi

	movq         0(%rax), %rdx
	movq         8(%rax), %r8

	movq	%rdx, %r9
	shr	$2,%r9
	leaq	(%rdi,%r9,8),%rsi

	vmovapd     size4negation0, %ymm15
	vmovapd     size4negation1, %ymm14
	vmovapd     size4negation2, %ymm13
	vmovapd     size4negation3, %ymm12

	movq	$0,%rax
	movq	%rdi,%r10
	movq	%rsi,%r11
fftsize2loop:
	vmovapd (%r10),%ymm0
	vmovapd (%r11),%ymm1
	vshufpd $0,%ymm0,%ymm0,%ymm2
	vshufpd $15,%ymm0,%ymm0,%ymm3
	vshufpd $0,%ymm1,%ymm1,%ymm4
	vshufpd $15,%ymm1,%ymm1,%ymm5
	vfmadd231pd %ymm3,%ymm12,%ymm2
	vfmadd231pd %ymm5,%ymm12,%ymm4
	vmovapd %ymm2,(%r10)
	vmovapd %ymm4,(%r11)

        leaq 32(%r10),%r10
        leaq 32(%r11),%r11
        addq $4,%rax
	cmpq %r9,%rax
	jb fftsize2loop

	movq	$0, %rax
	movq	%rdi,%r10
	movq	%rsi,%r11
fftsize4loop:
	vmovapd (%r10),%ymm0
	vmovapd (%r11),%ymm1
	vperm2f128 $32,%ymm0,%ymm0,%ymm4
	vperm2f128 $32,%ymm1,%ymm1,%ymm5
	vperm2f128 $49,%ymm0,%ymm0,%ymm6
	vperm2f128 $49,%ymm1,%ymm1,%ymm7
	vshufpd $10,%ymm7,%ymm6,%ymm8
	vshufpd $10,%ymm6,%ymm7,%ymm9
	vfmadd231pd %ymm8,%ymm13,%ymm4
	vfmadd231pd %ymm9,%ymm14,%ymm5
	vmovapd %ymm4,(%r10)
  	vmovapd %ymm5,(%r11)

        leaq 32(%r10),%r10
        leaq 32(%r11),%r11
        addq $4,%rax
	cmpq %r9,%rax
	jb fftsize4loop

	movq %r8,%rdx
	movq $4,%rax
ffthalfnnloop:
	movq $0,%rbx
fftblockloop:
	leaq (%rdi,%rbx,8),%r10
	leaq (%rsi,%rbx,8),%r11
	leaq (%r10,%rax,8),%r12
	leaq (%r11,%rax,8),%r13
	movq %rdx,%r14
	movq $0,%rcx
fftoffloop:
	vmovapd (%r10),%ymm0
	vmovapd (%r11),%ymm1
	vmovapd (%r12),%ymm2
	vmovapd (%r13),%ymm3
	vmovapd (%r14),%ymm4
	vmovapd 32(%r14),%ymm5
	vmulpd	%ymm2,%ymm4,%ymm6
	vmulpd	%ymm2,%ymm5,%ymm7
        vfnmadd231pd %ymm3,%ymm5,%ymm6
        vfmadd231pd %ymm3,%ymm4,%ymm7
	vsubpd	%ymm6,%ymm0,%ymm2
	vsubpd	%ymm7,%ymm1,%ymm3
	vaddpd	%ymm6,%ymm0,%ymm0
	vaddpd	%ymm7,%ymm1,%ymm1
	vmovapd %ymm0,(%r10)
	vmovapd %ymm1,(%r11)
	vmovapd %ymm2,(%r12)
	vmovapd %ymm3,(%r13)

    	leaq 	32(%r10),%r10
    	leaq	32(%r11),%r11
    	leaq 	32(%r12),%r12
    	leaq 	32(%r13),%r13
	leaq 	64(%r14),%r14
	addq 	$4,%rcx
	cmpq	%rax,%rcx
	jb 	fftoffloop

	leaq	(%rbx,%rax,2),%rbx
	cmpq	%r9,%rbx
	jb 	fftblockloop

	shlq	$1,%rax
	leaq	(%rdx,%rax,8),%rdx
	cmpq	%r9,%rax
	jb ffthalfnnloop

	movq $0,%rax
	movq %rdi,%r10
	movq %rsi,%r11
fftfinalloop:
	vmovapd	(%r10),%ymm0
	vmovapd	(%r11),%ymm1
	vmovapd (%rdx),%ymm2
	vmovapd 32(%rdx),%ymm3
	vmulpd %ymm0,%ymm2,%ymm4
	vmulpd %ymm0,%ymm3,%ymm5
	vmulpd %ymm1,%ymm2,%ymm6
	vmulpd %ymm1,%ymm3,%ymm7
	vsubpd %ymm7,%ymm4,%ymm0
	vaddpd %ymm6,%ymm5,%ymm1
	vmovapd %ymm0,(%r10)
	vmovapd %ymm1,(%r11)

    	leaq	32(%r10),%r10
    	leaq	32(%r11),%r11
    	leaq	64(%rdx),%rdx
	addq	$4,%rax
	cmpq	%r9,%rax
	jb fftfinalloop

fftend:
	vzeroall
	popq        %rbx
	popq        %r14
	popq        %r13
	popq        %r12
	popq        %r11
	popq        %r10
	retq

.balign 32
size4negation0: .double +1.0, +1.0, +1.0, -1.0
size4negation1: .double +1.0, -1.0, -1.0, +1.0
size4negation2: .double +1.0, +1.0, -1.0, -1.0
size4negation3: .double +1.0, -1.0, +1.0, -1.0

	.cfi_endproc
.LFE0:
	.size	fft, .-fft
	.section	.text.unlikely
.LCOLDE0:
	.text
.LHOTE0:
	.ident	"GCC: (Ubuntu 5.2.1-22ubuntu2) 5.2.1 20151010"
	.section	.note.GNU-stack,"",@progbits
