	.file	"lagrangehalfc_impl_asm.cpp"
	.section	.text.unlikely,"ax",@progbits
.LCOLDB0:
	.text
.LHOTB0:
	.p2align 4,,15
	.globl	LagrangeHalfCPolynomialMul
	.type	LagrangeHalfCPolynomialMul, @function
LagrangeHalfCPolynomialMul:
.LFB7134:
	.cfi_startproc
	pushq	%r12
	pushq	%r13
	pushq	%r14
	pushq	%r15

	movq	8(%rdi), %rax
	movq	$0,%rcx
	movl	8(%rax), %ecx
	movq	(%rdi), %r8
	movq	(%rsi), %r10
	movq	(%rdx), %r12
	leaq	(%r8,%rcx,8),%r9
	leaq	(%r10,%rcx,8),%r11
	leaq	(%r12,%rcx,8),%r13
	movq	%r9,%r14
.loopLHCmul:
	vmovupd (%r10),%ymm0
	vmovupd (%r11),%ymm1
	vmovupd (%r12),%ymm2
	vmovupd (%r13),%ymm3
	vmulpd	%ymm1,%ymm3,%ymm4
	vmulpd	%ymm0,%ymm3,%ymm6
	vfmadd231pd	%ymm2,%ymm1,%ymm6
	vfmsub231pd	%ymm0,%ymm2,%ymm4
	vmovupd	%ymm4,(%r8)
	vmovupd	%ymm6,(%r9)

	addq	$32,%r8
	addq	$32,%r9
	addq	$32,%r10
	addq	$32,%r11
	addq	$32,%r12
	addq	$32,%r13
	cmpq	%r14,%r8
	jb .loopLHCmul

	vzeroall
	popq	%r15
	popq	%r14
	popq	%r13
	popq	%r12
.L1:
	rep ret
	.p2align 4,,10
	.p2align 3
.L22:
	rep ret
	.p2align 4,,10
	.p2align 3
.L21:
	rep ret
	.cfi_endproc
.LFE7134:
	.size	LagrangeHalfCPolynomialMul, .-LagrangeHalfCPolynomialMul
	.section	.text.unlikely
.LCOLDE0:
	.text
.LHOTE0:
	.section	.text.unlikely
.LCOLDB1:
	.text
.LHOTB1:
	.p2align 4,,15
	.globl	LagrangeHalfCPolynomialAddMulASM
	.type	LagrangeHalfCPolynomialAddMulASM, @function
LagrangeHalfCPolynomialAddMulASM:
.LFB7135:
	.cfi_startproc
	pushq	%r12
	pushq	%r13
	pushq	%r14
	pushq	%r15

	movq	%rdi, %r8
	movq    %rsi, %r10
	movq	%rdx, %r12
	leaq	(%r8,%rcx,8),%r9
	leaq	(%r10,%rcx,8),%r11
	leaq	(%r12,%rcx,8),%r13
	movq	%r9,%r14
.loopLHCaddmul:
	vmovupd (%r8),%ymm8
	vmovupd (%r9),%ymm9
	vmovupd (%r10),%ymm0
	vmovupd (%r11),%ymm1
	vmovupd (%r12),%ymm2
	vmovupd (%r13),%ymm3
        vfmsub231pd	%ymm3,%ymm1,%ymm8
        vfmsub231pd	%ymm2,%ymm0,%ymm8
	vfmadd231pd	%ymm3,%ymm0,%ymm9
	vfmadd231pd	%ymm2,%ymm1,%ymm9
	vmovupd	%ymm8,(%r8)
	vmovupd	%ymm9,(%r9)

	addq	$32,%r8
	addq	$32,%r9
	addq	$32,%r10
	addq	$32,%r11
	addq	$32,%r12
	addq	$32,%r13
	cmpq	%r14,%r8
	jb .loopLHCaddmul

	vzeroall
	popq	%r15
	popq	%r14
	popq	%r13
	popq	%r12
.L23:
	rep ret
	.p2align 4,,10
	.p2align 3
.L43:
	rep ret
	.p2align 4,,10
	.p2align 3
.L42:
	rep ret
	.cfi_endproc
.LFE7135:
	.size	LagrangeHalfCPolynomialAddMulASM, .-LagrangeHalfCPolynomialAddMulASM
	.section	.text.unlikely
.LCOLDE1:
	.text
.LHOTE1:
	.section	.text.unlikely
.LCOLDB2:
	.text
.LHOTB2:
	.p2align 4,,15
	.globl	LagrangeHalfCPolynomialSubMul
	.type	LagrangeHalfCPolynomialSubMul, @function
LagrangeHalfCPolynomialSubMul:
.LFB7136:
	.cfi_startproc
	pushq	%r12
	pushq	%r13
	pushq	%r14
	pushq	%r15

	movq	8(%rdi), %rax
	movq	$0,%rcx
	movl	8(%rax), %ecx
	movq	(%rdi), %r8
	movq	(%rsi), %r10
	movq	(%rdx), %r12
	leaq	(%r8,%rcx,8),%r9
	leaq	(%r10,%rcx,8),%r11
	leaq	(%r12,%rcx,8),%r13
	movq	%r9,%r14
.loopLHCsubmul:
	vmovupd (%r8),%ymm8
	vmovupd (%r9),%ymm9
	vmovupd (%r10),%ymm0
	vmovupd (%r11),%ymm1
	vmovupd (%r12),%ymm2
	vmovupd (%r13),%ymm3
        vfnmsub231pd	%ymm3,%ymm1,%ymm8
        vfnmsub231pd	%ymm2,%ymm0,%ymm8
	vfmsub231pd	%ymm3,%ymm0,%ymm9
	vfnmsub231pd	%ymm2,%ymm1,%ymm9
	vmovupd	%ymm8,(%r8)
	vmovupd	%ymm9,(%r9)

	addq	$32,%r8
	addq	$32,%r9
	addq	$32,%r10
	addq	$32,%r11
	addq	$32,%r12
	addq	$32,%r13
	cmpq	%r14,%r8
	jb .loopLHCsubmul

	vzeroall
	popq	%r15
	popq	%r14
	popq	%r13
	popq	%r12
.L44:
	rep ret
	.p2align 4,,10
	.p2align 3
.L64:
	rep ret
	.p2align 4,,10
	.p2align 3
.L63:
	rep ret
	.cfi_endproc
.LFE7136:
	.size	LagrangeHalfCPolynomialSubMul, .-LagrangeHalfCPolynomialSubMul
	.section	.text.unlikely
.LCOLDE2:
	.text
.LHOTE2:
	.ident	"GCC: (Ubuntu 5.2.1-22ubuntu2) 5.2.1 20151010"
	.section	.note.GNU-stack,"",@progbits
