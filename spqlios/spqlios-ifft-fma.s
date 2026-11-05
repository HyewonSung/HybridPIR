	.file	"spqlios-ifft-avx.s"
	.section	.text.unlikely,"ax",@progbits
.LCOLDB0:
	.text
.LHOTB0:
	.p2align 4,,15
	.globl	ifft
	.type	ifft, @function
ifft:
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

	movq	%rdx, %r10
	shl	$1, %r10
	add	%r10, %rsi

	shr	$3, %r10
	movq    $0, %rcx
	movq	%r8, %r11
firstloop:
        vmovapd (%rdi,%rcx,8), %ymm0
        vmovapd (%rsi,%rcx,8), %ymm1
	vmovapd 0(%r11), %ymm2
	vmovapd 32(%r11), %ymm3
        vmulpd  %ymm0, %ymm2, %ymm4
        vmulpd  %ymm0, %ymm3, %ymm5
	vfnmadd231pd  %ymm1, %ymm3, %ymm4
	vfmadd231pd  %ymm1, %ymm2, %ymm5
	vmovapd %ymm4, (%rdi,%rcx,8)
	vmovapd %ymm5, (%rsi,%rcx,8)

	leaq    64(%r11), %r11
	addq	$4,%rcx
	cmpq	%r10,%rcx
	jb	firstloop

	movq %r10,%r12
nnloop:
	movq %r12,%r13
	shr  $1,%r13
	leaq (%r8,%r12,8),%r8
	leaq (%r8,%r12,8),%r8
	movq $0,%r11
blockloop:
	leaq (%rdi,%r11,8),%rax
	leaq (%rsi,%r11,8),%rbx
	leaq (%rax,%r13,8),%rcx
	leaq (%rbx,%r13,8),%rdx
	movq $0,%r9
	movq %r8,%r14
offloop:
        vmovapd (%rax,%r9,8), %ymm0
        vmovapd (%rbx,%r9,8), %ymm1
        vmovapd (%rcx,%r9,8), %ymm2
        vmovapd (%rdx,%r9,8), %ymm3
	vaddpd %ymm0,%ymm2,%ymm4
	vaddpd %ymm1,%ymm3,%ymm5
	vsubpd %ymm2,%ymm0,%ymm6
	vsubpd %ymm3,%ymm1,%ymm7
	vmovapd %ymm4,(%rax,%r9,8)
        vmovapd %ymm5,(%rbx,%r9,8)
        vmovapd (%r14),%ymm8
	vmovapd 32(%r14),%ymm9
	vmulpd %ymm6,%ymm8,%ymm4
	vfnmadd231pd %ymm7,%ymm9,%ymm4
	vmulpd %ymm6,%ymm9,%ymm5
	vfmadd231pd %ymm7,%ymm8,%ymm5
	vmovapd %ymm4,(%rcx,%r9,8)
        vmovapd %ymm5,(%rdx,%r9,8)

       	leaq 64(%r14),%r14
        addq $4,%r9
	cmpq %r13,%r9
        jb offloop

	addq %r12,%r11
	cmpq %r10,%r11
	jb blockloop

	movq %r13,%r12
	cmpq $8,%r12
	jae nnloop

	vmovapd     size4negation0, %ymm15
	vmovapd     size4negation1, %ymm14
	vmovapd     size4negation2, %ymm13
	vmovapd     size4negation3, %ymm12
	movq $0,%rax
	movq %rdi,%r11
	movq %rsi,%r12
size4loop:
	vmovapd (%r11),%ymm0
	vmovapd (%r12),%ymm1
	vshufpd $10,%ymm1,%ymm0,%ymm2
	vshufpd $10,%ymm0,%ymm1,%ymm3
	vperm2f128 $32,%ymm2,%ymm0,%ymm4
	vperm2f128 $49,%ymm2,%ymm0,%ymm5
	vperm2f128 $32,%ymm3,%ymm1,%ymm6
	vperm2f128 $49,%ymm3,%ymm1,%ymm7
	vmulpd	%ymm4,%ymm15,%ymm4
	vfmadd231pd	%ymm5,%ymm14,%ymm4
	vfmadd231pd	%ymm7,%ymm13,%ymm6
	vmovapd %ymm4,(%r11)
	vmovapd %ymm6,(%r12)

        leaq 32(%r11),%r11
        leaq 32(%r12),%r12
        addq $4,%rax
	cmpq %r10,%rax
	jb size4loop

	movq $0,%rax
	movq %rdi,%r11
	movq %rsi,%r12
size2loop:
	vmovapd (%r11),%ymm0
	vmovapd (%r12),%ymm1
	vshufpd $0,%ymm0,%ymm0,%ymm2
	vshufpd $15,%ymm0,%ymm0,%ymm3
	vshufpd $0,%ymm1,%ymm1,%ymm4
	vshufpd $15,%ymm1,%ymm1,%ymm5
	vfmadd231pd %ymm3,%ymm12,%ymm2
	vfmadd231pd %ymm5,%ymm12,%ymm4
	vmovapd %ymm2,(%r11)
	vmovapd %ymm4,(%r12)

        leaq 32(%r11),%r11
        leaq 32(%r12),%r12
        addq $4,%rax
	cmpq %r10,%rax
	jb size2loop

end:
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
size4negation1: .double +1.0, +1.0, -1.0, +1.0
size4negation2: .double +1.0, +1.0, -1.0, -1.0
size4negation3: .double +1.0, -1.0, +1.0, -1.0

	.cfi_endproc
.LFE0:
	.size	ifft, .-ifft
	.section	.text.unlikely
.LCOLDE0:
	.text
.LHOTE0:
	.ident	"GCC: (Ubuntu 5.2.1-22ubuntu2) 5.2.1 20151010"
	.section	.note.GNU-stack,"",@progbits
