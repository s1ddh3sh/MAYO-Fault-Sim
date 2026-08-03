
simpleserial-mayo-CW308_STM32F4.elf:     file format elf32-littlearm


Disassembly of section .text:

08000188 <deregister_tm_clones>:
 8000188:	4803      	ldr	r0, [pc, #12]	@ (8000198 <deregister_tm_clones+0x10>)
 800018a:	4b04      	ldr	r3, [pc, #16]	@ (800019c <deregister_tm_clones+0x14>)
 800018c:	4283      	cmp	r3, r0
 800018e:	d002      	beq.n	8000196 <deregister_tm_clones+0xe>
 8000190:	4b03      	ldr	r3, [pc, #12]	@ (80001a0 <deregister_tm_clones+0x18>)
 8000192:	b103      	cbz	r3, 8000196 <deregister_tm_clones+0xe>
 8000194:	4718      	bx	r3
 8000196:	4770      	bx	lr
 8000198:	20000538 	.word	0x20000538
 800019c:	20000538 	.word	0x20000538
 80001a0:	00000000 	.word	0x00000000

080001a4 <register_tm_clones>:
 80001a4:	4805      	ldr	r0, [pc, #20]	@ (80001bc <register_tm_clones+0x18>)
 80001a6:	4b06      	ldr	r3, [pc, #24]	@ (80001c0 <register_tm_clones+0x1c>)
 80001a8:	1a1b      	subs	r3, r3, r0
 80001aa:	0fd9      	lsrs	r1, r3, #31
 80001ac:	eb01 01a3 	add.w	r1, r1, r3, asr #2
 80001b0:	1049      	asrs	r1, r1, #1
 80001b2:	d002      	beq.n	80001ba <register_tm_clones+0x16>
 80001b4:	4b03      	ldr	r3, [pc, #12]	@ (80001c4 <register_tm_clones+0x20>)
 80001b6:	b103      	cbz	r3, 80001ba <register_tm_clones+0x16>
 80001b8:	4718      	bx	r3
 80001ba:	4770      	bx	lr
 80001bc:	20000538 	.word	0x20000538
 80001c0:	20000538 	.word	0x20000538
 80001c4:	00000000 	.word	0x00000000

080001c8 <__do_global_dtors_aux>:
 80001c8:	b510      	push	{r4, lr}
 80001ca:	4c06      	ldr	r4, [pc, #24]	@ (80001e4 <__do_global_dtors_aux+0x1c>)
 80001cc:	7823      	ldrb	r3, [r4, #0]
 80001ce:	b943      	cbnz	r3, 80001e2 <__do_global_dtors_aux+0x1a>
 80001d0:	f7ff ffda 	bl	8000188 <deregister_tm_clones>
 80001d4:	4b04      	ldr	r3, [pc, #16]	@ (80001e8 <__do_global_dtors_aux+0x20>)
 80001d6:	b113      	cbz	r3, 80001de <__do_global_dtors_aux+0x16>
 80001d8:	4804      	ldr	r0, [pc, #16]	@ (80001ec <__do_global_dtors_aux+0x24>)
 80001da:	f3af 8000 	nop.w
 80001de:	2301      	movs	r3, #1
 80001e0:	7023      	strb	r3, [r4, #0]
 80001e2:	bd10      	pop	{r4, pc}
 80001e4:	20000538 	.word	0x20000538
 80001e8:	00000000 	.word	0x00000000
 80001ec:	08002f10 	.word	0x08002f10

080001f0 <frame_dummy>:
 80001f0:	b508      	push	{r3, lr}
 80001f2:	4b04      	ldr	r3, [pc, #16]	@ (8000204 <frame_dummy+0x14>)
 80001f4:	b11b      	cbz	r3, 80001fe <frame_dummy+0xe>
 80001f6:	4904      	ldr	r1, [pc, #16]	@ (8000208 <frame_dummy+0x18>)
 80001f8:	4804      	ldr	r0, [pc, #16]	@ (800020c <frame_dummy+0x1c>)
 80001fa:	f3af 8000 	nop.w
 80001fe:	e8bd 4008 	ldmia.w	sp!, {r3, lr}
 8000202:	e7cf      	b.n	80001a4 <register_tm_clones>
 8000204:	00000000 	.word	0x00000000
 8000208:	2000053c 	.word	0x2000053c
 800020c:	08002f10 	.word	0x08002f10

08000210 <memset>:
 8000210:	0783      	lsls	r3, r0, #30
 8000212:	b530      	push	{r4, r5, lr}
 8000214:	d047      	beq.n	80002a6 <memset+0x96>
 8000216:	1e54      	subs	r4, r2, #1
 8000218:	2a00      	cmp	r2, #0
 800021a:	d03e      	beq.n	800029a <memset+0x8a>
 800021c:	b2ca      	uxtb	r2, r1
 800021e:	4603      	mov	r3, r0
 8000220:	e001      	b.n	8000226 <memset+0x16>
 8000222:	3c01      	subs	r4, #1
 8000224:	d339      	bcc.n	800029a <memset+0x8a>
 8000226:	f803 2b01 	strb.w	r2, [r3], #1
 800022a:	079d      	lsls	r5, r3, #30
 800022c:	d1f9      	bne.n	8000222 <memset+0x12>
 800022e:	2c03      	cmp	r4, #3
 8000230:	d92c      	bls.n	800028c <memset+0x7c>
 8000232:	b2cd      	uxtb	r5, r1
 8000234:	eb05 2505 	add.w	r5, r5, r5, lsl #8
 8000238:	2c0f      	cmp	r4, #15
 800023a:	eb05 4505 	add.w	r5, r5, r5, lsl #16
 800023e:	d935      	bls.n	80002ac <memset+0x9c>
 8000240:	f1a4 0210 	sub.w	r2, r4, #16
 8000244:	f022 0c0f 	bic.w	ip, r2, #15
 8000248:	f103 0e10 	add.w	lr, r3, #16
 800024c:	44e6      	add	lr, ip
 800024e:	ea4f 1c12 	mov.w	ip, r2, lsr #4
 8000252:	461a      	mov	r2, r3
 8000254:	e9c2 5500 	strd	r5, r5, [r2]
 8000258:	e9c2 5502 	strd	r5, r5, [r2, #8]
 800025c:	3210      	adds	r2, #16
 800025e:	4572      	cmp	r2, lr
 8000260:	d1f8      	bne.n	8000254 <memset+0x44>
 8000262:	f10c 0201 	add.w	r2, ip, #1
 8000266:	f014 0f0c 	tst.w	r4, #12
 800026a:	eb03 1202 	add.w	r2, r3, r2, lsl #4
 800026e:	f004 0c0f 	and.w	ip, r4, #15
 8000272:	d013      	beq.n	800029c <memset+0x8c>
 8000274:	f1ac 0304 	sub.w	r3, ip, #4
 8000278:	f023 0303 	bic.w	r3, r3, #3
 800027c:	3304      	adds	r3, #4
 800027e:	4413      	add	r3, r2
 8000280:	f842 5b04 	str.w	r5, [r2], #4
 8000284:	4293      	cmp	r3, r2
 8000286:	d1fb      	bne.n	8000280 <memset+0x70>
 8000288:	f00c 0403 	and.w	r4, ip, #3
 800028c:	b12c      	cbz	r4, 800029a <memset+0x8a>
 800028e:	b2c9      	uxtb	r1, r1
 8000290:	441c      	add	r4, r3
 8000292:	f803 1b01 	strb.w	r1, [r3], #1
 8000296:	42a3      	cmp	r3, r4
 8000298:	d1fb      	bne.n	8000292 <memset+0x82>
 800029a:	bd30      	pop	{r4, r5, pc}
 800029c:	4664      	mov	r4, ip
 800029e:	4613      	mov	r3, r2
 80002a0:	2c00      	cmp	r4, #0
 80002a2:	d1f4      	bne.n	800028e <memset+0x7e>
 80002a4:	e7f9      	b.n	800029a <memset+0x8a>
 80002a6:	4603      	mov	r3, r0
 80002a8:	4614      	mov	r4, r2
 80002aa:	e7c0      	b.n	800022e <memset+0x1e>
 80002ac:	461a      	mov	r2, r3
 80002ae:	46a4      	mov	ip, r4
 80002b0:	e7e0      	b.n	8000274 <memset+0x64>
 80002b2:	bf00      	nop

080002b4 <calloc>:
 80002b4:	b410      	push	{r4}
 80002b6:	4c04      	ldr	r4, [pc, #16]	@ (80002c8 <calloc+0x14>)
 80002b8:	4603      	mov	r3, r0
 80002ba:	460a      	mov	r2, r1
 80002bc:	6820      	ldr	r0, [r4, #0]
 80002be:	bc10      	pop	{r4}
 80002c0:	4619      	mov	r1, r3
 80002c2:	f000 b809 	b.w	80002d8 <_calloc_r>
 80002c6:	bf00      	nop
 80002c8:	20000000 	.word	0x20000000

080002cc <__errno>:
 80002cc:	4b01      	ldr	r3, [pc, #4]	@ (80002d4 <__errno+0x8>)
 80002ce:	6818      	ldr	r0, [r3, #0]
 80002d0:	4770      	bx	lr
 80002d2:	bf00      	nop
 80002d4:	20000000 	.word	0x20000000

080002d8 <_calloc_r>:
 80002d8:	b538      	push	{r3, r4, r5, lr}
 80002da:	fba1 1402 	umull	r1, r4, r1, r2
 80002de:	bb54      	cbnz	r4, 8000336 <_calloc_r+0x5e>
 80002e0:	f000 f830 	bl	8000344 <_malloc_r>
 80002e4:	4605      	mov	r5, r0
 80002e6:	b350      	cbz	r0, 800033e <_calloc_r+0x66>
 80002e8:	f850 2c04 	ldr.w	r2, [r0, #-4]
 80002ec:	f022 0203 	bic.w	r2, r2, #3
 80002f0:	3a04      	subs	r2, #4
 80002f2:	2a24      	cmp	r2, #36	@ 0x24
 80002f4:	d810      	bhi.n	8000318 <_calloc_r+0x40>
 80002f6:	2a13      	cmp	r2, #19
 80002f8:	d913      	bls.n	8000322 <_calloc_r+0x4a>
 80002fa:	2a1b      	cmp	r2, #27
 80002fc:	e9c0 4400 	strd	r4, r4, [r0]
 8000300:	d916      	bls.n	8000330 <_calloc_r+0x58>
 8000302:	2a24      	cmp	r2, #36	@ 0x24
 8000304:	e9c0 4402 	strd	r4, r4, [r0, #8]
 8000308:	bf0a      	itet	eq
 800030a:	e9c0 4404 	strdeq	r4, r4, [r0, #16]
 800030e:	f100 0210 	addne.w	r2, r0, #16
 8000312:	f100 0218 	addeq.w	r2, r0, #24
 8000316:	e005      	b.n	8000324 <_calloc_r+0x4c>
 8000318:	4621      	mov	r1, r4
 800031a:	f7ff ff79 	bl	8000210 <memset>
 800031e:	4628      	mov	r0, r5
 8000320:	bd38      	pop	{r3, r4, r5, pc}
 8000322:	4602      	mov	r2, r0
 8000324:	2300      	movs	r3, #0
 8000326:	e9c2 3300 	strd	r3, r3, [r2]
 800032a:	6093      	str	r3, [r2, #8]
 800032c:	4628      	mov	r0, r5
 800032e:	bd38      	pop	{r3, r4, r5, pc}
 8000330:	f100 0208 	add.w	r2, r0, #8
 8000334:	e7f6      	b.n	8000324 <_calloc_r+0x4c>
 8000336:	f7ff ffc9 	bl	80002cc <__errno>
 800033a:	230c      	movs	r3, #12
 800033c:	6003      	str	r3, [r0, #0]
 800033e:	2500      	movs	r5, #0
 8000340:	4628      	mov	r0, r5
 8000342:	bd38      	pop	{r3, r4, r5, pc}

08000344 <_malloc_r>:
 8000344:	e92d 4ff0 	stmdb	sp!, {r4, r5, r6, r7, r8, r9, sl, fp, lr}
 8000348:	f101 050b 	add.w	r5, r1, #11
 800034c:	2d16      	cmp	r5, #22
 800034e:	b083      	sub	sp, #12
 8000350:	4606      	mov	r6, r0
 8000352:	d823      	bhi.n	800039c <_malloc_r+0x58>
 8000354:	2910      	cmp	r1, #16
 8000356:	f200 80af 	bhi.w	80004b8 <_malloc_r+0x174>
 800035a:	f000 fcd9 	bl	8000d10 <__malloc_lock>
 800035e:	2510      	movs	r5, #16
 8000360:	2318      	movs	r3, #24
 8000362:	2002      	movs	r0, #2
 8000364:	4fc1      	ldr	r7, [pc, #772]	@ (800066c <_malloc_r+0x328>)
 8000366:	443b      	add	r3, r7
 8000368:	f1a3 0208 	sub.w	r2, r3, #8
 800036c:	685c      	ldr	r4, [r3, #4]
 800036e:	4294      	cmp	r4, r2
 8000370:	f000 8125 	beq.w	80005be <_malloc_r+0x27a>
 8000374:	6863      	ldr	r3, [r4, #4]
 8000376:	68e2      	ldr	r2, [r4, #12]
 8000378:	68a1      	ldr	r1, [r4, #8]
 800037a:	f023 0303 	bic.w	r3, r3, #3
 800037e:	60ca      	str	r2, [r1, #12]
 8000380:	4423      	add	r3, r4
 8000382:	4630      	mov	r0, r6
 8000384:	6091      	str	r1, [r2, #8]
 8000386:	685a      	ldr	r2, [r3, #4]
 8000388:	f042 0201 	orr.w	r2, r2, #1
 800038c:	605a      	str	r2, [r3, #4]
 800038e:	f000 fcc1 	bl	8000d14 <__malloc_unlock>
 8000392:	3408      	adds	r4, #8
 8000394:	4620      	mov	r0, r4
 8000396:	b003      	add	sp, #12
 8000398:	e8bd 8ff0 	ldmia.w	sp!, {r4, r5, r6, r7, r8, r9, sl, fp, pc}
 800039c:	f035 0507 	bics.w	r5, r5, #7
 80003a0:	f100 808a 	bmi.w	80004b8 <_malloc_r+0x174>
 80003a4:	42a9      	cmp	r1, r5
 80003a6:	f200 8087 	bhi.w	80004b8 <_malloc_r+0x174>
 80003aa:	f000 fcb1 	bl	8000d10 <__malloc_lock>
 80003ae:	f5b5 7ffc 	cmp.w	r5, #504	@ 0x1f8
 80003b2:	f0c0 816d 	bcc.w	8000690 <_malloc_r+0x34c>
 80003b6:	0a6b      	lsrs	r3, r5, #9
 80003b8:	f000 8082 	beq.w	80004c0 <_malloc_r+0x17c>
 80003bc:	2b04      	cmp	r3, #4
 80003be:	f200 8124 	bhi.w	800060a <_malloc_r+0x2c6>
 80003c2:	09ab      	lsrs	r3, r5, #6
 80003c4:	f103 0039 	add.w	r0, r3, #57	@ 0x39
 80003c8:	f103 0e38 	add.w	lr, r3, #56	@ 0x38
 80003cc:	00c3      	lsls	r3, r0, #3
 80003ce:	4fa7      	ldr	r7, [pc, #668]	@ (800066c <_malloc_r+0x328>)
 80003d0:	443b      	add	r3, r7
 80003d2:	f1a3 0c08 	sub.w	ip, r3, #8
 80003d6:	685c      	ldr	r4, [r3, #4]
 80003d8:	45a4      	cmp	ip, r4
 80003da:	d107      	bne.n	80003ec <_malloc_r+0xa8>
 80003dc:	e00d      	b.n	80003fa <_malloc_r+0xb6>
 80003de:	2a00      	cmp	r2, #0
 80003e0:	68e1      	ldr	r1, [r4, #12]
 80003e2:	f280 80e6 	bge.w	80005b2 <_malloc_r+0x26e>
 80003e6:	458c      	cmp	ip, r1
 80003e8:	d007      	beq.n	80003fa <_malloc_r+0xb6>
 80003ea:	460c      	mov	r4, r1
 80003ec:	6863      	ldr	r3, [r4, #4]
 80003ee:	f023 0303 	bic.w	r3, r3, #3
 80003f2:	1b5a      	subs	r2, r3, r5
 80003f4:	2a0f      	cmp	r2, #15
 80003f6:	ddf2      	ble.n	80003de <_malloc_r+0x9a>
 80003f8:	4670      	mov	r0, lr
 80003fa:	f8df 8274 	ldr.w	r8, [pc, #628]	@ 8000670 <_malloc_r+0x32c>
 80003fe:	693c      	ldr	r4, [r7, #16]
 8000400:	4544      	cmp	r4, r8
 8000402:	f000 80c3 	beq.w	800058c <_malloc_r+0x248>
 8000406:	6863      	ldr	r3, [r4, #4]
 8000408:	f023 0c03 	bic.w	ip, r3, #3
 800040c:	ebac 0305 	sub.w	r3, ip, r5
 8000410:	2b0f      	cmp	r3, #15
 8000412:	f300 8141 	bgt.w	8000698 <_malloc_r+0x354>
 8000416:	2b00      	cmp	r3, #0
 8000418:	e9c7 8804 	strd	r8, r8, [r7, #16]
 800041c:	f280 812c 	bge.w	8000678 <_malloc_r+0x334>
 8000420:	f5bc 7f00 	cmp.w	ip, #512	@ 0x200
 8000424:	f8d7 e004 	ldr.w	lr, [r7, #4]
 8000428:	f080 80cf 	bcs.w	80005ca <_malloc_r+0x286>
 800042c:	ea4f 01dc 	mov.w	r1, ip, lsr #3
 8000430:	3101      	adds	r1, #1
 8000432:	ea4f 1c5c 	mov.w	ip, ip, lsr #5
 8000436:	2301      	movs	r3, #1
 8000438:	fa03 f30c 	lsl.w	r3, r3, ip
 800043c:	f857 2031 	ldr.w	r2, [r7, r1, lsl #3]
 8000440:	60a2      	str	r2, [r4, #8]
 8000442:	ea4e 0e03 	orr.w	lr, lr, r3
 8000446:	eb07 03c1 	add.w	r3, r7, r1, lsl #3
 800044a:	3b08      	subs	r3, #8
 800044c:	60e3      	str	r3, [r4, #12]
 800044e:	f8c7 e004 	str.w	lr, [r7, #4]
 8000452:	f847 4031 	str.w	r4, [r7, r1, lsl #3]
 8000456:	60d4      	str	r4, [r2, #12]
 8000458:	1083      	asrs	r3, r0, #2
 800045a:	f04f 0c01 	mov.w	ip, #1
 800045e:	fa0c fc03 	lsl.w	ip, ip, r3
 8000462:	45f4      	cmp	ip, lr
 8000464:	d832      	bhi.n	80004cc <_malloc_r+0x188>
 8000466:	ea1c 0f0e 	tst.w	ip, lr
 800046a:	d108      	bne.n	800047e <_malloc_r+0x13a>
 800046c:	f020 0003 	bic.w	r0, r0, #3
 8000470:	ea4f 0c4c 	mov.w	ip, ip, lsl #1
 8000474:	ea1c 0f0e 	tst.w	ip, lr
 8000478:	f100 0004 	add.w	r0, r0, #4
 800047c:	d0f8      	beq.n	8000470 <_malloc_r+0x12c>
 800047e:	eb07 0ac0 	add.w	sl, r7, r0, lsl #3
 8000482:	46d6      	mov	lr, sl
 8000484:	4681      	mov	r9, r0
 8000486:	f8de 300c 	ldr.w	r3, [lr, #12]
 800048a:	e00b      	b.n	80004a4 <_malloc_r+0x160>
 800048c:	6859      	ldr	r1, [r3, #4]
 800048e:	f021 0103 	bic.w	r1, r1, #3
 8000492:	1b4a      	subs	r2, r1, r5
 8000494:	2a0f      	cmp	r2, #15
 8000496:	461c      	mov	r4, r3
 8000498:	68db      	ldr	r3, [r3, #12]
 800049a:	f300 80c2 	bgt.w	8000622 <_malloc_r+0x2de>
 800049e:	2a00      	cmp	r2, #0
 80004a0:	f280 80d6 	bge.w	8000650 <_malloc_r+0x30c>
 80004a4:	459e      	cmp	lr, r3
 80004a6:	d1f1      	bne.n	800048c <_malloc_r+0x148>
 80004a8:	f109 0901 	add.w	r9, r9, #1
 80004ac:	f019 0f03 	tst.w	r9, #3
 80004b0:	f10e 0e08 	add.w	lr, lr, #8
 80004b4:	d1e7      	bne.n	8000486 <_malloc_r+0x142>
 80004b6:	e118      	b.n	80006ea <_malloc_r+0x3a6>
 80004b8:	230c      	movs	r3, #12
 80004ba:	6033      	str	r3, [r6, #0]
 80004bc:	2400      	movs	r4, #0
 80004be:	e769      	b.n	8000394 <_malloc_r+0x50>
 80004c0:	f44f 7300 	mov.w	r3, #512	@ 0x200
 80004c4:	2040      	movs	r0, #64	@ 0x40
 80004c6:	f04f 0e3f 	mov.w	lr, #63	@ 0x3f
 80004ca:	e780      	b.n	80003ce <_malloc_r+0x8a>
 80004cc:	68bc      	ldr	r4, [r7, #8]
 80004ce:	6863      	ldr	r3, [r4, #4]
 80004d0:	f023 0903 	bic.w	r9, r3, #3
 80004d4:	45a9      	cmp	r9, r5
 80004d6:	d303      	bcc.n	80004e0 <_malloc_r+0x19c>
 80004d8:	eba9 0305 	sub.w	r3, r9, r5
 80004dc:	2b0f      	cmp	r3, #15
 80004de:	dc58      	bgt.n	8000592 <_malloc_r+0x24e>
 80004e0:	f8df b190 	ldr.w	fp, [pc, #400]	@ 8000674 <_malloc_r+0x330>
 80004e4:	f8db 2000 	ldr.w	r2, [fp]
 80004e8:	eb04 0309 	add.w	r3, r4, r9
 80004ec:	3210      	adds	r2, #16
 80004ee:	2008      	movs	r0, #8
 80004f0:	eb02 0805 	add.w	r8, r2, r5
 80004f4:	9300      	str	r3, [sp, #0]
 80004f6:	f000 fbd5 	bl	8000ca4 <sysconf>
 80004fa:	f8d7 1408 	ldr.w	r1, [r7, #1032]	@ 0x408
 80004fe:	3101      	adds	r1, #1
 8000500:	4602      	mov	r2, r0
 8000502:	d005      	beq.n	8000510 <_malloc_r+0x1cc>
 8000504:	f108 38ff 	add.w	r8, r8, #4294967295	@ 0xffffffff
 8000508:	4480      	add	r8, r0
 800050a:	4241      	negs	r1, r0
 800050c:	ea01 0808 	and.w	r8, r1, r8
 8000510:	4641      	mov	r1, r8
 8000512:	4630      	mov	r0, r6
 8000514:	9200      	str	r2, [sp, #0]
 8000516:	f000 fc13 	bl	8000d40 <_sbrk_r>
 800051a:	f1b0 3fff 	cmp.w	r0, #4294967295	@ 0xffffffff
 800051e:	9a00      	ldr	r2, [sp, #0]
 8000520:	4682      	mov	sl, r0
 8000522:	f000 80d4 	beq.w	80006ce <_malloc_r+0x38a>
 8000526:	eb04 0309 	add.w	r3, r4, r9
 800052a:	4283      	cmp	r3, r0
 800052c:	f200 80cd 	bhi.w	80006ca <_malloc_r+0x386>
 8000530:	f8db 1004 	ldr.w	r1, [fp, #4]
 8000534:	4441      	add	r1, r8
 8000536:	f8cb 1004 	str.w	r1, [fp, #4]
 800053a:	4608      	mov	r0, r1
 800053c:	f102 3cff 	add.w	ip, r2, #4294967295	@ 0xffffffff
 8000540:	f040 80f5 	bne.w	800072e <_malloc_r+0x3ea>
 8000544:	ea1a 0f0c 	tst.w	sl, ip
 8000548:	f040 80f1 	bne.w	800072e <_malloc_r+0x3ea>
 800054c:	f8d7 a008 	ldr.w	sl, [r7, #8]
 8000550:	44c8      	add	r8, r9
 8000552:	f048 0001 	orr.w	r0, r8, #1
 8000556:	f8ca 0004 	str.w	r0, [sl, #4]
 800055a:	f8db 202c 	ldr.w	r2, [fp, #44]	@ 0x2c
 800055e:	428a      	cmp	r2, r1
 8000560:	f8db 2030 	ldr.w	r2, [fp, #48]	@ 0x30
 8000564:	bf38      	it	cc
 8000566:	f8cb 102c 	strcc.w	r1, [fp, #44]	@ 0x2c
 800056a:	428a      	cmp	r2, r1
 800056c:	bf38      	it	cc
 800056e:	f8cb 1030 	strcc.w	r1, [fp, #48]	@ 0x30
 8000572:	4654      	mov	r4, sl
 8000574:	f020 0803 	bic.w	r8, r0, #3
 8000578:	45a8      	cmp	r8, r5
 800057a:	eba8 0305 	sub.w	r3, r8, r5
 800057e:	d301      	bcc.n	8000584 <_malloc_r+0x240>
 8000580:	2b0f      	cmp	r3, #15
 8000582:	dc06      	bgt.n	8000592 <_malloc_r+0x24e>
 8000584:	4630      	mov	r0, r6
 8000586:	f000 fbc5 	bl	8000d14 <__malloc_unlock>
 800058a:	e797      	b.n	80004bc <_malloc_r+0x178>
 800058c:	f8d7 e004 	ldr.w	lr, [r7, #4]
 8000590:	e762      	b.n	8000458 <_malloc_r+0x114>
 8000592:	1962      	adds	r2, r4, r5
 8000594:	f043 0301 	orr.w	r3, r3, #1
 8000598:	f045 0501 	orr.w	r5, r5, #1
 800059c:	6065      	str	r5, [r4, #4]
 800059e:	4630      	mov	r0, r6
 80005a0:	60ba      	str	r2, [r7, #8]
 80005a2:	6053      	str	r3, [r2, #4]
 80005a4:	f000 fbb6 	bl	8000d14 <__malloc_unlock>
 80005a8:	3408      	adds	r4, #8
 80005aa:	4620      	mov	r0, r4
 80005ac:	b003      	add	sp, #12
 80005ae:	e8bd 8ff0 	ldmia.w	sp!, {r4, r5, r6, r7, r8, r9, sl, fp, pc}
 80005b2:	68a2      	ldr	r2, [r4, #8]
 80005b4:	4423      	add	r3, r4
 80005b6:	60d1      	str	r1, [r2, #12]
 80005b8:	4630      	mov	r0, r6
 80005ba:	608a      	str	r2, [r1, #8]
 80005bc:	e6e3      	b.n	8000386 <_malloc_r+0x42>
 80005be:	68dc      	ldr	r4, [r3, #12]
 80005c0:	42a3      	cmp	r3, r4
 80005c2:	f47f aed7 	bne.w	8000374 <_malloc_r+0x30>
 80005c6:	3002      	adds	r0, #2
 80005c8:	e717      	b.n	80003fa <_malloc_r+0xb6>
 80005ca:	f5bc 6f20 	cmp.w	ip, #2560	@ 0xa00
 80005ce:	ea4f 225c 	mov.w	r2, ip, lsr #9
 80005d2:	d373      	bcc.n	80006bc <_malloc_r+0x378>
 80005d4:	2a14      	cmp	r2, #20
 80005d6:	f200 810a 	bhi.w	80007ee <_malloc_r+0x4aa>
 80005da:	f102 035c 	add.w	r3, r2, #92	@ 0x5c
 80005de:	00db      	lsls	r3, r3, #3
 80005e0:	325b      	adds	r2, #91	@ 0x5b
 80005e2:	18f9      	adds	r1, r7, r3
 80005e4:	58fb      	ldr	r3, [r7, r3]
 80005e6:	3908      	subs	r1, #8
 80005e8:	4299      	cmp	r1, r3
 80005ea:	d103      	bne.n	80005f4 <_malloc_r+0x2b0>
 80005ec:	e0e6      	b.n	80007bc <_malloc_r+0x478>
 80005ee:	689b      	ldr	r3, [r3, #8]
 80005f0:	4299      	cmp	r1, r3
 80005f2:	d004      	beq.n	80005fe <_malloc_r+0x2ba>
 80005f4:	685a      	ldr	r2, [r3, #4]
 80005f6:	f022 0203 	bic.w	r2, r2, #3
 80005fa:	4562      	cmp	r2, ip
 80005fc:	d8f7      	bhi.n	80005ee <_malloc_r+0x2aa>
 80005fe:	68d9      	ldr	r1, [r3, #12]
 8000600:	e9c4 3102 	strd	r3, r1, [r4, #8]
 8000604:	608c      	str	r4, [r1, #8]
 8000606:	60dc      	str	r4, [r3, #12]
 8000608:	e726      	b.n	8000458 <_malloc_r+0x114>
 800060a:	2b14      	cmp	r3, #20
 800060c:	d962      	bls.n	80006d4 <_malloc_r+0x390>
 800060e:	2b54      	cmp	r3, #84	@ 0x54
 8000610:	f200 80f6 	bhi.w	8000800 <_malloc_r+0x4bc>
 8000614:	0b2b      	lsrs	r3, r5, #12
 8000616:	f103 006f 	add.w	r0, r3, #111	@ 0x6f
 800061a:	f103 0e6e 	add.w	lr, r3, #110	@ 0x6e
 800061e:	00c3      	lsls	r3, r0, #3
 8000620:	e6d5      	b.n	80003ce <_malloc_r+0x8a>
 8000622:	f8d4 c008 	ldr.w	ip, [r4, #8]
 8000626:	4630      	mov	r0, r6
 8000628:	1966      	adds	r6, r4, r5
 800062a:	f045 0501 	orr.w	r5, r5, #1
 800062e:	6065      	str	r5, [r4, #4]
 8000630:	f8cc 300c 	str.w	r3, [ip, #12]
 8000634:	f8c3 c008 	str.w	ip, [r3, #8]
 8000638:	f042 0301 	orr.w	r3, r2, #1
 800063c:	e9c7 6604 	strd	r6, r6, [r7, #16]
 8000640:	e9c6 8802 	strd	r8, r8, [r6, #8]
 8000644:	6073      	str	r3, [r6, #4]
 8000646:	5062      	str	r2, [r4, r1]
 8000648:	f000 fb64 	bl	8000d14 <__malloc_unlock>
 800064c:	3408      	adds	r4, #8
 800064e:	e6a1      	b.n	8000394 <_malloc_r+0x50>
 8000650:	4421      	add	r1, r4
 8000652:	4630      	mov	r0, r6
 8000654:	684a      	ldr	r2, [r1, #4]
 8000656:	f042 0201 	orr.w	r2, r2, #1
 800065a:	604a      	str	r2, [r1, #4]
 800065c:	f854 2f08 	ldr.w	r2, [r4, #8]!
 8000660:	60d3      	str	r3, [r2, #12]
 8000662:	609a      	str	r2, [r3, #8]
 8000664:	f000 fb56 	bl	8000d14 <__malloc_unlock>
 8000668:	e694      	b.n	8000394 <_malloc_r+0x50>
 800066a:	bf00      	nop
 800066c:	20000128 	.word	0x20000128
 8000670:	20000130 	.word	0x20000130
 8000674:	20006980 	.word	0x20006980
 8000678:	44a4      	add	ip, r4
 800067a:	4630      	mov	r0, r6
 800067c:	f8dc 3004 	ldr.w	r3, [ip, #4]
 8000680:	f043 0301 	orr.w	r3, r3, #1
 8000684:	f8cc 3004 	str.w	r3, [ip, #4]
 8000688:	f000 fb44 	bl	8000d14 <__malloc_unlock>
 800068c:	3408      	adds	r4, #8
 800068e:	e681      	b.n	8000394 <_malloc_r+0x50>
 8000690:	08e8      	lsrs	r0, r5, #3
 8000692:	f105 0308 	add.w	r3, r5, #8
 8000696:	e665      	b.n	8000364 <_malloc_r+0x20>
 8000698:	1962      	adds	r2, r4, r5
 800069a:	f043 0101 	orr.w	r1, r3, #1
 800069e:	f045 0501 	orr.w	r5, r5, #1
 80006a2:	6065      	str	r5, [r4, #4]
 80006a4:	4630      	mov	r0, r6
 80006a6:	e9c7 2204 	strd	r2, r2, [r7, #16]
 80006aa:	e9c2 8802 	strd	r8, r8, [r2, #8]
 80006ae:	6051      	str	r1, [r2, #4]
 80006b0:	f844 300c 	str.w	r3, [r4, ip]
 80006b4:	f000 fb2e 	bl	8000d14 <__malloc_unlock>
 80006b8:	3408      	adds	r4, #8
 80006ba:	e66b      	b.n	8000394 <_malloc_r+0x50>
 80006bc:	ea4f 129c 	mov.w	r2, ip, lsr #6
 80006c0:	f102 0339 	add.w	r3, r2, #57	@ 0x39
 80006c4:	00db      	lsls	r3, r3, #3
 80006c6:	3238      	adds	r2, #56	@ 0x38
 80006c8:	e78b      	b.n	80005e2 <_malloc_r+0x29e>
 80006ca:	42bc      	cmp	r4, r7
 80006cc:	d028      	beq.n	8000720 <_malloc_r+0x3dc>
 80006ce:	68bc      	ldr	r4, [r7, #8]
 80006d0:	6860      	ldr	r0, [r4, #4]
 80006d2:	e74f      	b.n	8000574 <_malloc_r+0x230>
 80006d4:	f103 005c 	add.w	r0, r3, #92	@ 0x5c
 80006d8:	f103 0e5b 	add.w	lr, r3, #91	@ 0x5b
 80006dc:	00c3      	lsls	r3, r0, #3
 80006de:	e676      	b.n	80003ce <_malloc_r+0x8a>
 80006e0:	f85a 3908 	ldr.w	r3, [sl], #-8
 80006e4:	4553      	cmp	r3, sl
 80006e6:	f040 80d8 	bne.w	800089a <_malloc_r+0x556>
 80006ea:	f010 0f03 	tst.w	r0, #3
 80006ee:	f100 30ff 	add.w	r0, r0, #4294967295	@ 0xffffffff
 80006f2:	d1f5      	bne.n	80006e0 <_malloc_r+0x39c>
 80006f4:	687b      	ldr	r3, [r7, #4]
 80006f6:	ea23 030c 	bic.w	r3, r3, ip
 80006fa:	607b      	str	r3, [r7, #4]
 80006fc:	ea4f 0c4c 	mov.w	ip, ip, lsl #1
 8000700:	459c      	cmp	ip, r3
 8000702:	f63f aee3 	bhi.w	80004cc <_malloc_r+0x188>
 8000706:	f1bc 0f00 	cmp.w	ip, #0
 800070a:	d104      	bne.n	8000716 <_malloc_r+0x3d2>
 800070c:	e6de      	b.n	80004cc <_malloc_r+0x188>
 800070e:	ea4f 0c4c 	mov.w	ip, ip, lsl #1
 8000712:	f109 0904 	add.w	r9, r9, #4
 8000716:	ea1c 0f03 	tst.w	ip, r3
 800071a:	d0f8      	beq.n	800070e <_malloc_r+0x3ca>
 800071c:	4648      	mov	r0, r9
 800071e:	e6ae      	b.n	800047e <_malloc_r+0x13a>
 8000720:	f8db 0004 	ldr.w	r0, [fp, #4]
 8000724:	4440      	add	r0, r8
 8000726:	f8cb 0004 	str.w	r0, [fp, #4]
 800072a:	f102 3cff 	add.w	ip, r2, #4294967295	@ 0xffffffff
 800072e:	f8d7 1408 	ldr.w	r1, [r7, #1032]	@ 0x408
 8000732:	3101      	adds	r1, #1
 8000734:	d06e      	beq.n	8000814 <_malloc_r+0x4d0>
 8000736:	eb04 0309 	add.w	r3, r4, r9
 800073a:	ebaa 0303 	sub.w	r3, sl, r3
 800073e:	4418      	add	r0, r3
 8000740:	f8cb 0004 	str.w	r0, [fp, #4]
 8000744:	f01a 0307 	ands.w	r3, sl, #7
 8000748:	9300      	str	r3, [sp, #0]
 800074a:	d041      	beq.n	80007d0 <_malloc_r+0x48c>
 800074c:	f1c3 0108 	rsb	r1, r3, #8
 8000750:	448a      	add	sl, r1
 8000752:	44d0      	add	r8, sl
 8000754:	440a      	add	r2, r1
 8000756:	ea08 010c 	and.w	r1, r8, ip
 800075a:	1a52      	subs	r2, r2, r1
 800075c:	ea02 010c 	and.w	r1, r2, ip
 8000760:	4630      	mov	r0, r6
 8000762:	9101      	str	r1, [sp, #4]
 8000764:	f000 faec 	bl	8000d40 <_sbrk_r>
 8000768:	1c42      	adds	r2, r0, #1
 800076a:	d06f      	beq.n	800084c <_malloc_r+0x508>
 800076c:	9901      	ldr	r1, [sp, #4]
 800076e:	eba0 000a 	sub.w	r0, r0, sl
 8000772:	eb00 0801 	add.w	r8, r0, r1
 8000776:	f8db 2004 	ldr.w	r2, [fp, #4]
 800077a:	f8c7 a008 	str.w	sl, [r7, #8]
 800077e:	f048 0001 	orr.w	r0, r8, #1
 8000782:	4411      	add	r1, r2
 8000784:	42bc      	cmp	r4, r7
 8000786:	f8ca 0004 	str.w	r0, [sl, #4]
 800078a:	f8cb 1004 	str.w	r1, [fp, #4]
 800078e:	f43f aee4 	beq.w	800055a <_malloc_r+0x216>
 8000792:	f1b9 0f0f 	cmp.w	r9, #15
 8000796:	d940      	bls.n	800081a <_malloc_r+0x4d6>
 8000798:	6862      	ldr	r2, [r4, #4]
 800079a:	f1a9 000c 	sub.w	r0, r9, #12
 800079e:	f020 0007 	bic.w	r0, r0, #7
 80007a2:	f002 0201 	and.w	r2, r2, #1
 80007a6:	4302      	orrs	r2, r0
 80007a8:	6062      	str	r2, [r4, #4]
 80007aa:	2305      	movs	r3, #5
 80007ac:	1822      	adds	r2, r4, r0
 80007ae:	280f      	cmp	r0, #15
 80007b0:	e9c2 3301 	strd	r3, r3, [r2, #4]
 80007b4:	d852      	bhi.n	800085c <_malloc_r+0x518>
 80007b6:	f8da 0004 	ldr.w	r0, [sl, #4]
 80007ba:	e6ce      	b.n	800055a <_malloc_r+0x216>
 80007bc:	1092      	asrs	r2, r2, #2
 80007be:	f04f 0c01 	mov.w	ip, #1
 80007c2:	fa0c f202 	lsl.w	r2, ip, r2
 80007c6:	ea4e 0e02 	orr.w	lr, lr, r2
 80007ca:	f8c7 e004 	str.w	lr, [r7, #4]
 80007ce:	e717      	b.n	8000600 <_malloc_r+0x2bc>
 80007d0:	eb0a 0108 	add.w	r1, sl, r8
 80007d4:	ea01 010c 	and.w	r1, r1, ip
 80007d8:	1a52      	subs	r2, r2, r1
 80007da:	ea02 010c 	and.w	r1, r2, ip
 80007de:	4630      	mov	r0, r6
 80007e0:	9101      	str	r1, [sp, #4]
 80007e2:	f000 faad 	bl	8000d40 <_sbrk_r>
 80007e6:	1c43      	adds	r3, r0, #1
 80007e8:	d1c0      	bne.n	800076c <_malloc_r+0x428>
 80007ea:	9900      	ldr	r1, [sp, #0]
 80007ec:	e7c3      	b.n	8000776 <_malloc_r+0x432>
 80007ee:	2a54      	cmp	r2, #84	@ 0x54
 80007f0:	d817      	bhi.n	8000822 <_malloc_r+0x4de>
 80007f2:	ea4f 321c 	mov.w	r2, ip, lsr #12
 80007f6:	f102 036f 	add.w	r3, r2, #111	@ 0x6f
 80007fa:	00db      	lsls	r3, r3, #3
 80007fc:	326e      	adds	r2, #110	@ 0x6e
 80007fe:	e6f0      	b.n	80005e2 <_malloc_r+0x29e>
 8000800:	f5b3 7faa 	cmp.w	r3, #340	@ 0x154
 8000804:	d817      	bhi.n	8000836 <_malloc_r+0x4f2>
 8000806:	0beb      	lsrs	r3, r5, #15
 8000808:	f103 0078 	add.w	r0, r3, #120	@ 0x78
 800080c:	f103 0e77 	add.w	lr, r3, #119	@ 0x77
 8000810:	00c3      	lsls	r3, r0, #3
 8000812:	e5dc      	b.n	80003ce <_malloc_r+0x8a>
 8000814:	f8c7 a408 	str.w	sl, [r7, #1032]	@ 0x408
 8000818:	e794      	b.n	8000744 <_malloc_r+0x400>
 800081a:	2301      	movs	r3, #1
 800081c:	f8ca 3004 	str.w	r3, [sl, #4]
 8000820:	e6b0      	b.n	8000584 <_malloc_r+0x240>
 8000822:	f5b2 7faa 	cmp.w	r2, #340	@ 0x154
 8000826:	d823      	bhi.n	8000870 <_malloc_r+0x52c>
 8000828:	ea4f 32dc 	mov.w	r2, ip, lsr #15
 800082c:	f102 0378 	add.w	r3, r2, #120	@ 0x78
 8000830:	00db      	lsls	r3, r3, #3
 8000832:	3277      	adds	r2, #119	@ 0x77
 8000834:	e6d5      	b.n	80005e2 <_malloc_r+0x29e>
 8000836:	f240 5254 	movw	r2, #1364	@ 0x554
 800083a:	4293      	cmp	r3, r2
 800083c:	d823      	bhi.n	8000886 <_malloc_r+0x542>
 800083e:	0cab      	lsrs	r3, r5, #18
 8000840:	f103 007d 	add.w	r0, r3, #125	@ 0x7d
 8000844:	f103 0e7c 	add.w	lr, r3, #124	@ 0x7c
 8000848:	00c3      	lsls	r3, r0, #3
 800084a:	e5c0      	b.n	80003ce <_malloc_r+0x8a>
 800084c:	9b00      	ldr	r3, [sp, #0]
 800084e:	f1a3 0208 	sub.w	r2, r3, #8
 8000852:	4490      	add	r8, r2
 8000854:	eba8 080a 	sub.w	r8, r8, sl
 8000858:	2100      	movs	r1, #0
 800085a:	e78c      	b.n	8000776 <_malloc_r+0x432>
 800085c:	f104 0108 	add.w	r1, r4, #8
 8000860:	4630      	mov	r0, r6
 8000862:	f000 f879 	bl	8000958 <_free_r>
 8000866:	f8db 1004 	ldr.w	r1, [fp, #4]
 800086a:	f8d7 a008 	ldr.w	sl, [r7, #8]
 800086e:	e7a2      	b.n	80007b6 <_malloc_r+0x472>
 8000870:	f240 5354 	movw	r3, #1364	@ 0x554
 8000874:	429a      	cmp	r2, r3
 8000876:	d80c      	bhi.n	8000892 <_malloc_r+0x54e>
 8000878:	ea4f 429c 	mov.w	r2, ip, lsr #18
 800087c:	f102 037d 	add.w	r3, r2, #125	@ 0x7d
 8000880:	00db      	lsls	r3, r3, #3
 8000882:	327c      	adds	r2, #124	@ 0x7c
 8000884:	e6ad      	b.n	80005e2 <_malloc_r+0x29e>
 8000886:	f44f 737e 	mov.w	r3, #1016	@ 0x3f8
 800088a:	207f      	movs	r0, #127	@ 0x7f
 800088c:	f04f 0e7e 	mov.w	lr, #126	@ 0x7e
 8000890:	e59d      	b.n	80003ce <_malloc_r+0x8a>
 8000892:	f44f 737e 	mov.w	r3, #1016	@ 0x3f8
 8000896:	227e      	movs	r2, #126	@ 0x7e
 8000898:	e6a3      	b.n	80005e2 <_malloc_r+0x29e>
 800089a:	687b      	ldr	r3, [r7, #4]
 800089c:	e72e      	b.n	80006fc <_malloc_r+0x3b8>
 800089e:	bf00      	nop

080008a0 <_malloc_trim_r>:
 80008a0:	e92d 43f8 	stmdb	sp!, {r3, r4, r5, r6, r7, r8, r9, lr}
 80008a4:	4606      	mov	r6, r0
 80008a6:	2008      	movs	r0, #8
 80008a8:	4689      	mov	r9, r1
 80008aa:	f000 f9fb 	bl	8000ca4 <sysconf>
 80008ae:	f8df 809c 	ldr.w	r8, [pc, #156]	@ 800094c <_malloc_trim_r+0xac>
 80008b2:	4605      	mov	r5, r0
 80008b4:	4630      	mov	r0, r6
 80008b6:	f000 fa2b 	bl	8000d10 <__malloc_lock>
 80008ba:	f8d8 3008 	ldr.w	r3, [r8, #8]
 80008be:	685f      	ldr	r7, [r3, #4]
 80008c0:	f027 0703 	bic.w	r7, r7, #3
 80008c4:	f1a7 0411 	sub.w	r4, r7, #17
 80008c8:	eba4 0409 	sub.w	r4, r4, r9
 80008cc:	442c      	add	r4, r5
 80008ce:	fbb4 f4f5 	udiv	r4, r4, r5
 80008d2:	3c01      	subs	r4, #1
 80008d4:	fb05 f404 	mul.w	r4, r5, r4
 80008d8:	42a5      	cmp	r5, r4
 80008da:	dc08      	bgt.n	80008ee <_malloc_trim_r+0x4e>
 80008dc:	2100      	movs	r1, #0
 80008de:	4630      	mov	r0, r6
 80008e0:	f000 fa2e 	bl	8000d40 <_sbrk_r>
 80008e4:	f8d8 3008 	ldr.w	r3, [r8, #8]
 80008e8:	443b      	add	r3, r7
 80008ea:	4298      	cmp	r0, r3
 80008ec:	d005      	beq.n	80008fa <_malloc_trim_r+0x5a>
 80008ee:	4630      	mov	r0, r6
 80008f0:	f000 fa10 	bl	8000d14 <__malloc_unlock>
 80008f4:	2000      	movs	r0, #0
 80008f6:	e8bd 83f8 	ldmia.w	sp!, {r3, r4, r5, r6, r7, r8, r9, pc}
 80008fa:	4261      	negs	r1, r4
 80008fc:	4630      	mov	r0, r6
 80008fe:	f000 fa1f 	bl	8000d40 <_sbrk_r>
 8000902:	3001      	adds	r0, #1
 8000904:	d00f      	beq.n	8000926 <_malloc_trim_r+0x86>
 8000906:	4a12      	ldr	r2, [pc, #72]	@ (8000950 <_malloc_trim_r+0xb0>)
 8000908:	f8d8 3008 	ldr.w	r3, [r8, #8]
 800090c:	1b3f      	subs	r7, r7, r4
 800090e:	f047 0701 	orr.w	r7, r7, #1
 8000912:	605f      	str	r7, [r3, #4]
 8000914:	6813      	ldr	r3, [r2, #0]
 8000916:	4630      	mov	r0, r6
 8000918:	1b1b      	subs	r3, r3, r4
 800091a:	6013      	str	r3, [r2, #0]
 800091c:	f000 f9fa 	bl	8000d14 <__malloc_unlock>
 8000920:	2001      	movs	r0, #1
 8000922:	e8bd 83f8 	ldmia.w	sp!, {r3, r4, r5, r6, r7, r8, r9, pc}
 8000926:	2100      	movs	r1, #0
 8000928:	4630      	mov	r0, r6
 800092a:	f000 fa09 	bl	8000d40 <_sbrk_r>
 800092e:	f8d8 2008 	ldr.w	r2, [r8, #8]
 8000932:	1a83      	subs	r3, r0, r2
 8000934:	2b0f      	cmp	r3, #15
 8000936:	ddda      	ble.n	80008ee <_malloc_trim_r+0x4e>
 8000938:	f043 0301 	orr.w	r3, r3, #1
 800093c:	6053      	str	r3, [r2, #4]
 800093e:	4b05      	ldr	r3, [pc, #20]	@ (8000954 <_malloc_trim_r+0xb4>)
 8000940:	4903      	ldr	r1, [pc, #12]	@ (8000950 <_malloc_trim_r+0xb0>)
 8000942:	681b      	ldr	r3, [r3, #0]
 8000944:	1ac0      	subs	r0, r0, r3
 8000946:	6008      	str	r0, [r1, #0]
 8000948:	e7d1      	b.n	80008ee <_malloc_trim_r+0x4e>
 800094a:	bf00      	nop
 800094c:	20000128 	.word	0x20000128
 8000950:	20006984 	.word	0x20006984
 8000954:	20000530 	.word	0x20000530

08000958 <_free_r>:
 8000958:	2900      	cmp	r1, #0
 800095a:	d07c      	beq.n	8000a56 <_free_r+0xfe>
 800095c:	b5f8      	push	{r3, r4, r5, r6, r7, lr}
 800095e:	460c      	mov	r4, r1
 8000960:	4606      	mov	r6, r0
 8000962:	f000 f9d5 	bl	8000d10 <__malloc_lock>
 8000966:	f854 3c04 	ldr.w	r3, [r4, #-4]
 800096a:	4f76      	ldr	r7, [pc, #472]	@ (8000b44 <_free_r+0x1ec>)
 800096c:	f1a4 0508 	sub.w	r5, r4, #8
 8000970:	f023 0101 	bic.w	r1, r3, #1
 8000974:	1868      	adds	r0, r5, r1
 8000976:	f8d7 c008 	ldr.w	ip, [r7, #8]
 800097a:	6842      	ldr	r2, [r0, #4]
 800097c:	4584      	cmp	ip, r0
 800097e:	f022 0203 	bic.w	r2, r2, #3
 8000982:	f000 8083 	beq.w	8000a8c <_free_r+0x134>
 8000986:	07db      	lsls	r3, r3, #31
 8000988:	6042      	str	r2, [r0, #4]
 800098a:	eb00 0c02 	add.w	ip, r0, r2
 800098e:	d433      	bmi.n	80009f8 <_free_r+0xa0>
 8000990:	f854 4c08 	ldr.w	r4, [r4, #-8]
 8000994:	f8dc 3004 	ldr.w	r3, [ip, #4]
 8000998:	1b2d      	subs	r5, r5, r4
 800099a:	4421      	add	r1, r4
 800099c:	68ac      	ldr	r4, [r5, #8]
 800099e:	f107 0c08 	add.w	ip, r7, #8
 80009a2:	4564      	cmp	r4, ip
 80009a4:	f003 0301 	and.w	r3, r3, #1
 80009a8:	d064      	beq.n	8000a74 <_free_r+0x11c>
 80009aa:	f8d5 e00c 	ldr.w	lr, [r5, #12]
 80009ae:	f8c4 e00c 	str.w	lr, [r4, #12]
 80009b2:	f8ce 4008 	str.w	r4, [lr, #8]
 80009b6:	2b00      	cmp	r3, #0
 80009b8:	f000 8081 	beq.w	8000abe <_free_r+0x166>
 80009bc:	f041 0301 	orr.w	r3, r1, #1
 80009c0:	606b      	str	r3, [r5, #4]
 80009c2:	6001      	str	r1, [r0, #0]
 80009c4:	f5b1 7f00 	cmp.w	r1, #512	@ 0x200
 80009c8:	d222      	bcs.n	8000a10 <_free_r+0xb8>
 80009ca:	6878      	ldr	r0, [r7, #4]
 80009cc:	08cb      	lsrs	r3, r1, #3
 80009ce:	2201      	movs	r2, #1
 80009d0:	0949      	lsrs	r1, r1, #5
 80009d2:	3301      	adds	r3, #1
 80009d4:	408a      	lsls	r2, r1
 80009d6:	4302      	orrs	r2, r0
 80009d8:	f857 1033 	ldr.w	r1, [r7, r3, lsl #3]
 80009dc:	607a      	str	r2, [r7, #4]
 80009de:	eb07 02c3 	add.w	r2, r7, r3, lsl #3
 80009e2:	3a08      	subs	r2, #8
 80009e4:	e9c5 1202 	strd	r1, r2, [r5, #8]
 80009e8:	f847 5033 	str.w	r5, [r7, r3, lsl #3]
 80009ec:	60cd      	str	r5, [r1, #12]
 80009ee:	4630      	mov	r0, r6
 80009f0:	e8bd 40f8 	ldmia.w	sp!, {r3, r4, r5, r6, r7, lr}
 80009f4:	f000 b98e 	b.w	8000d14 <__malloc_unlock>
 80009f8:	f8dc 3004 	ldr.w	r3, [ip, #4]
 80009fc:	07db      	lsls	r3, r3, #31
 80009fe:	d52b      	bpl.n	8000a58 <_free_r+0x100>
 8000a00:	f041 0301 	orr.w	r3, r1, #1
 8000a04:	f5b1 7f00 	cmp.w	r1, #512	@ 0x200
 8000a08:	f844 3c04 	str.w	r3, [r4, #-4]
 8000a0c:	6001      	str	r1, [r0, #0]
 8000a0e:	d3dc      	bcc.n	80009ca <_free_r+0x72>
 8000a10:	f5b1 6f20 	cmp.w	r1, #2560	@ 0xa00
 8000a14:	ea4f 2351 	mov.w	r3, r1, lsr #9
 8000a18:	d253      	bcs.n	8000ac2 <_free_r+0x16a>
 8000a1a:	098b      	lsrs	r3, r1, #6
 8000a1c:	f103 0039 	add.w	r0, r3, #57	@ 0x39
 8000a20:	f103 0238 	add.w	r2, r3, #56	@ 0x38
 8000a24:	00c3      	lsls	r3, r0, #3
 8000a26:	18f8      	adds	r0, r7, r3
 8000a28:	58fb      	ldr	r3, [r7, r3]
 8000a2a:	3808      	subs	r0, #8
 8000a2c:	4298      	cmp	r0, r3
 8000a2e:	d103      	bne.n	8000a38 <_free_r+0xe0>
 8000a30:	e061      	b.n	8000af6 <_free_r+0x19e>
 8000a32:	689b      	ldr	r3, [r3, #8]
 8000a34:	4298      	cmp	r0, r3
 8000a36:	d004      	beq.n	8000a42 <_free_r+0xea>
 8000a38:	685a      	ldr	r2, [r3, #4]
 8000a3a:	f022 0203 	bic.w	r2, r2, #3
 8000a3e:	428a      	cmp	r2, r1
 8000a40:	d8f7      	bhi.n	8000a32 <_free_r+0xda>
 8000a42:	68d8      	ldr	r0, [r3, #12]
 8000a44:	e9c5 3002 	strd	r3, r0, [r5, #8]
 8000a48:	6085      	str	r5, [r0, #8]
 8000a4a:	60dd      	str	r5, [r3, #12]
 8000a4c:	4630      	mov	r0, r6
 8000a4e:	e8bd 40f8 	ldmia.w	sp!, {r3, r4, r5, r6, r7, lr}
 8000a52:	f000 b95f 	b.w	8000d14 <__malloc_unlock>
 8000a56:	4770      	bx	lr
 8000a58:	4411      	add	r1, r2
 8000a5a:	f107 0c08 	add.w	ip, r7, #8
 8000a5e:	6883      	ldr	r3, [r0, #8]
 8000a60:	4563      	cmp	r3, ip
 8000a62:	d03f      	beq.n	8000ae4 <_free_r+0x18c>
 8000a64:	68c2      	ldr	r2, [r0, #12]
 8000a66:	60da      	str	r2, [r3, #12]
 8000a68:	6093      	str	r3, [r2, #8]
 8000a6a:	f041 0301 	orr.w	r3, r1, #1
 8000a6e:	606b      	str	r3, [r5, #4]
 8000a70:	5069      	str	r1, [r5, r1]
 8000a72:	e7a7      	b.n	80009c4 <_free_r+0x6c>
 8000a74:	2b00      	cmp	r3, #0
 8000a76:	d15f      	bne.n	8000b38 <_free_r+0x1e0>
 8000a78:	440a      	add	r2, r1
 8000a7a:	e9d0 1302 	ldrd	r1, r3, [r0, #8]
 8000a7e:	60cb      	str	r3, [r1, #12]
 8000a80:	6099      	str	r1, [r3, #8]
 8000a82:	f042 0301 	orr.w	r3, r2, #1
 8000a86:	606b      	str	r3, [r5, #4]
 8000a88:	50aa      	str	r2, [r5, r2]
 8000a8a:	e7b0      	b.n	80009ee <_free_r+0x96>
 8000a8c:	440a      	add	r2, r1
 8000a8e:	07d9      	lsls	r1, r3, #31
 8000a90:	d407      	bmi.n	8000aa2 <_free_r+0x14a>
 8000a92:	f854 3c08 	ldr.w	r3, [r4, #-8]
 8000a96:	1aed      	subs	r5, r5, r3
 8000a98:	441a      	add	r2, r3
 8000a9a:	e9d5 1302 	ldrd	r1, r3, [r5, #8]
 8000a9e:	60cb      	str	r3, [r1, #12]
 8000aa0:	6099      	str	r1, [r3, #8]
 8000aa2:	f042 0301 	orr.w	r3, r2, #1
 8000aa6:	606b      	str	r3, [r5, #4]
 8000aa8:	4b27      	ldr	r3, [pc, #156]	@ (8000b48 <_free_r+0x1f0>)
 8000aaa:	60bd      	str	r5, [r7, #8]
 8000aac:	681b      	ldr	r3, [r3, #0]
 8000aae:	4293      	cmp	r3, r2
 8000ab0:	d89d      	bhi.n	80009ee <_free_r+0x96>
 8000ab2:	4b26      	ldr	r3, [pc, #152]	@ (8000b4c <_free_r+0x1f4>)
 8000ab4:	4630      	mov	r0, r6
 8000ab6:	6819      	ldr	r1, [r3, #0]
 8000ab8:	f7ff fef2 	bl	80008a0 <_malloc_trim_r>
 8000abc:	e797      	b.n	80009ee <_free_r+0x96>
 8000abe:	4411      	add	r1, r2
 8000ac0:	e7cd      	b.n	8000a5e <_free_r+0x106>
 8000ac2:	2b14      	cmp	r3, #20
 8000ac4:	d908      	bls.n	8000ad8 <_free_r+0x180>
 8000ac6:	2b54      	cmp	r3, #84	@ 0x54
 8000ac8:	d81d      	bhi.n	8000b06 <_free_r+0x1ae>
 8000aca:	0b0b      	lsrs	r3, r1, #12
 8000acc:	f103 006f 	add.w	r0, r3, #111	@ 0x6f
 8000ad0:	f103 026e 	add.w	r2, r3, #110	@ 0x6e
 8000ad4:	00c3      	lsls	r3, r0, #3
 8000ad6:	e7a6      	b.n	8000a26 <_free_r+0xce>
 8000ad8:	f103 005c 	add.w	r0, r3, #92	@ 0x5c
 8000adc:	f103 025b 	add.w	r2, r3, #91	@ 0x5b
 8000ae0:	00c3      	lsls	r3, r0, #3
 8000ae2:	e7a0      	b.n	8000a26 <_free_r+0xce>
 8000ae4:	f041 0301 	orr.w	r3, r1, #1
 8000ae8:	e9c7 5504 	strd	r5, r5, [r7, #16]
 8000aec:	e9c5 cc02 	strd	ip, ip, [r5, #8]
 8000af0:	606b      	str	r3, [r5, #4]
 8000af2:	5069      	str	r1, [r5, r1]
 8000af4:	e77b      	b.n	80009ee <_free_r+0x96>
 8000af6:	6879      	ldr	r1, [r7, #4]
 8000af8:	1092      	asrs	r2, r2, #2
 8000afa:	2401      	movs	r4, #1
 8000afc:	fa04 f202 	lsl.w	r2, r4, r2
 8000b00:	430a      	orrs	r2, r1
 8000b02:	607a      	str	r2, [r7, #4]
 8000b04:	e79e      	b.n	8000a44 <_free_r+0xec>
 8000b06:	f5b3 7faa 	cmp.w	r3, #340	@ 0x154
 8000b0a:	d806      	bhi.n	8000b1a <_free_r+0x1c2>
 8000b0c:	0bcb      	lsrs	r3, r1, #15
 8000b0e:	f103 0078 	add.w	r0, r3, #120	@ 0x78
 8000b12:	f103 0277 	add.w	r2, r3, #119	@ 0x77
 8000b16:	00c3      	lsls	r3, r0, #3
 8000b18:	e785      	b.n	8000a26 <_free_r+0xce>
 8000b1a:	f240 5254 	movw	r2, #1364	@ 0x554
 8000b1e:	4293      	cmp	r3, r2
 8000b20:	d806      	bhi.n	8000b30 <_free_r+0x1d8>
 8000b22:	0c8b      	lsrs	r3, r1, #18
 8000b24:	f103 007d 	add.w	r0, r3, #125	@ 0x7d
 8000b28:	f103 027c 	add.w	r2, r3, #124	@ 0x7c
 8000b2c:	00c3      	lsls	r3, r0, #3
 8000b2e:	e77a      	b.n	8000a26 <_free_r+0xce>
 8000b30:	f44f 737e 	mov.w	r3, #1016	@ 0x3f8
 8000b34:	227e      	movs	r2, #126	@ 0x7e
 8000b36:	e776      	b.n	8000a26 <_free_r+0xce>
 8000b38:	f041 0301 	orr.w	r3, r1, #1
 8000b3c:	606b      	str	r3, [r5, #4]
 8000b3e:	6001      	str	r1, [r0, #0]
 8000b40:	e755      	b.n	80009ee <_free_r+0x96>
 8000b42:	bf00      	nop
 8000b44:	20000128 	.word	0x20000128
 8000b48:	20000534 	.word	0x20000534
 8000b4c:	20006980 	.word	0x20006980

08000b50 <malloc>:
 8000b50:	4b02      	ldr	r3, [pc, #8]	@ (8000b5c <malloc+0xc>)
 8000b52:	4601      	mov	r1, r0
 8000b54:	6818      	ldr	r0, [r3, #0]
 8000b56:	f7ff bbf5 	b.w	8000344 <_malloc_r>
 8000b5a:	bf00      	nop
 8000b5c:	20000000 	.word	0x20000000

08000b60 <free>:
 8000b60:	4b02      	ldr	r3, [pc, #8]	@ (8000b6c <free+0xc>)
 8000b62:	4601      	mov	r1, r0
 8000b64:	6818      	ldr	r0, [r3, #0]
 8000b66:	f7ff bef7 	b.w	8000958 <_free_r>
 8000b6a:	bf00      	nop
 8000b6c:	20000000 	.word	0x20000000

08000b70 <memcpy>:
 8000b70:	4684      	mov	ip, r0
 8000b72:	ea41 0300 	orr.w	r3, r1, r0
 8000b76:	f013 0303 	ands.w	r3, r3, #3
 8000b7a:	d16d      	bne.n	8000c58 <memcpy+0xe8>
 8000b7c:	3a40      	subs	r2, #64	@ 0x40
 8000b7e:	d341      	bcc.n	8000c04 <memcpy+0x94>
 8000b80:	f851 3b04 	ldr.w	r3, [r1], #4
 8000b84:	f840 3b04 	str.w	r3, [r0], #4
 8000b88:	f851 3b04 	ldr.w	r3, [r1], #4
 8000b8c:	f840 3b04 	str.w	r3, [r0], #4
 8000b90:	f851 3b04 	ldr.w	r3, [r1], #4
 8000b94:	f840 3b04 	str.w	r3, [r0], #4
 8000b98:	f851 3b04 	ldr.w	r3, [r1], #4
 8000b9c:	f840 3b04 	str.w	r3, [r0], #4
 8000ba0:	f851 3b04 	ldr.w	r3, [r1], #4
 8000ba4:	f840 3b04 	str.w	r3, [r0], #4
 8000ba8:	f851 3b04 	ldr.w	r3, [r1], #4
 8000bac:	f840 3b04 	str.w	r3, [r0], #4
 8000bb0:	f851 3b04 	ldr.w	r3, [r1], #4
 8000bb4:	f840 3b04 	str.w	r3, [r0], #4
 8000bb8:	f851 3b04 	ldr.w	r3, [r1], #4
 8000bbc:	f840 3b04 	str.w	r3, [r0], #4
 8000bc0:	f851 3b04 	ldr.w	r3, [r1], #4
 8000bc4:	f840 3b04 	str.w	r3, [r0], #4
 8000bc8:	f851 3b04 	ldr.w	r3, [r1], #4
 8000bcc:	f840 3b04 	str.w	r3, [r0], #4
 8000bd0:	f851 3b04 	ldr.w	r3, [r1], #4
 8000bd4:	f840 3b04 	str.w	r3, [r0], #4
 8000bd8:	f851 3b04 	ldr.w	r3, [r1], #4
 8000bdc:	f840 3b04 	str.w	r3, [r0], #4
 8000be0:	f851 3b04 	ldr.w	r3, [r1], #4
 8000be4:	f840 3b04 	str.w	r3, [r0], #4
 8000be8:	f851 3b04 	ldr.w	r3, [r1], #4
 8000bec:	f840 3b04 	str.w	r3, [r0], #4
 8000bf0:	f851 3b04 	ldr.w	r3, [r1], #4
 8000bf4:	f840 3b04 	str.w	r3, [r0], #4
 8000bf8:	f851 3b04 	ldr.w	r3, [r1], #4
 8000bfc:	f840 3b04 	str.w	r3, [r0], #4
 8000c00:	3a40      	subs	r2, #64	@ 0x40
 8000c02:	d2bd      	bcs.n	8000b80 <memcpy+0x10>
 8000c04:	3230      	adds	r2, #48	@ 0x30
 8000c06:	d311      	bcc.n	8000c2c <memcpy+0xbc>
 8000c08:	f851 3b04 	ldr.w	r3, [r1], #4
 8000c0c:	f840 3b04 	str.w	r3, [r0], #4
 8000c10:	f851 3b04 	ldr.w	r3, [r1], #4
 8000c14:	f840 3b04 	str.w	r3, [r0], #4
 8000c18:	f851 3b04 	ldr.w	r3, [r1], #4
 8000c1c:	f840 3b04 	str.w	r3, [r0], #4
 8000c20:	f851 3b04 	ldr.w	r3, [r1], #4
 8000c24:	f840 3b04 	str.w	r3, [r0], #4
 8000c28:	3a10      	subs	r2, #16
 8000c2a:	d2ed      	bcs.n	8000c08 <memcpy+0x98>
 8000c2c:	320c      	adds	r2, #12
 8000c2e:	d305      	bcc.n	8000c3c <memcpy+0xcc>
 8000c30:	f851 3b04 	ldr.w	r3, [r1], #4
 8000c34:	f840 3b04 	str.w	r3, [r0], #4
 8000c38:	3a04      	subs	r2, #4
 8000c3a:	d2f9      	bcs.n	8000c30 <memcpy+0xc0>
 8000c3c:	3204      	adds	r2, #4
 8000c3e:	d008      	beq.n	8000c52 <memcpy+0xe2>
 8000c40:	07d2      	lsls	r2, r2, #31
 8000c42:	bf1c      	itt	ne
 8000c44:	f811 3b01 	ldrbne.w	r3, [r1], #1
 8000c48:	f800 3b01 	strbne.w	r3, [r0], #1
 8000c4c:	d301      	bcc.n	8000c52 <memcpy+0xe2>
 8000c4e:	880b      	ldrh	r3, [r1, #0]
 8000c50:	8003      	strh	r3, [r0, #0]
 8000c52:	4660      	mov	r0, ip
 8000c54:	4770      	bx	lr
 8000c56:	bf00      	nop
 8000c58:	2a08      	cmp	r2, #8
 8000c5a:	d313      	bcc.n	8000c84 <memcpy+0x114>
 8000c5c:	078b      	lsls	r3, r1, #30
 8000c5e:	d08d      	beq.n	8000b7c <memcpy+0xc>
 8000c60:	f010 0303 	ands.w	r3, r0, #3
 8000c64:	d08a      	beq.n	8000b7c <memcpy+0xc>
 8000c66:	f1c3 0304 	rsb	r3, r3, #4
 8000c6a:	1ad2      	subs	r2, r2, r3
 8000c6c:	07db      	lsls	r3, r3, #31
 8000c6e:	bf1c      	itt	ne
 8000c70:	f811 3b01 	ldrbne.w	r3, [r1], #1
 8000c74:	f800 3b01 	strbne.w	r3, [r0], #1
 8000c78:	d380      	bcc.n	8000b7c <memcpy+0xc>
 8000c7a:	f831 3b02 	ldrh.w	r3, [r1], #2
 8000c7e:	f820 3b02 	strh.w	r3, [r0], #2
 8000c82:	e77b      	b.n	8000b7c <memcpy+0xc>
 8000c84:	3a04      	subs	r2, #4
 8000c86:	d3d9      	bcc.n	8000c3c <memcpy+0xcc>
 8000c88:	3a01      	subs	r2, #1
 8000c8a:	f811 3b01 	ldrb.w	r3, [r1], #1
 8000c8e:	f800 3b01 	strb.w	r3, [r0], #1
 8000c92:	d2f9      	bcs.n	8000c88 <memcpy+0x118>
 8000c94:	780b      	ldrb	r3, [r1, #0]
 8000c96:	7003      	strb	r3, [r0, #0]
 8000c98:	784b      	ldrb	r3, [r1, #1]
 8000c9a:	7043      	strb	r3, [r0, #1]
 8000c9c:	788b      	ldrb	r3, [r1, #2]
 8000c9e:	7083      	strb	r3, [r0, #2]
 8000ca0:	4660      	mov	r0, ip
 8000ca2:	4770      	bx	lr

08000ca4 <sysconf>:
 8000ca4:	2808      	cmp	r0, #8
 8000ca6:	d102      	bne.n	8000cae <sysconf+0xa>
 8000ca8:	f44f 5080 	mov.w	r0, #4096	@ 0x1000
 8000cac:	4770      	bx	lr
 8000cae:	b508      	push	{r3, lr}
 8000cb0:	f7ff fb0c 	bl	80002cc <__errno>
 8000cb4:	2316      	movs	r3, #22
 8000cb6:	6003      	str	r3, [r0, #0]
 8000cb8:	f04f 30ff 	mov.w	r0, #4294967295	@ 0xffffffff
 8000cbc:	bd08      	pop	{r3, pc}
 8000cbe:	bf00      	nop

08000cc0 <__libc_init_array>:
 8000cc0:	b570      	push	{r4, r5, r6, lr}
 8000cc2:	4b0f      	ldr	r3, [pc, #60]	@ (8000d00 <__libc_init_array+0x40>)
 8000cc4:	4d0f      	ldr	r5, [pc, #60]	@ (8000d04 <__libc_init_array+0x44>)
 8000cc6:	42ab      	cmp	r3, r5
 8000cc8:	eba3 0605 	sub.w	r6, r3, r5
 8000ccc:	d007      	beq.n	8000cde <__libc_init_array+0x1e>
 8000cce:	10b6      	asrs	r6, r6, #2
 8000cd0:	2400      	movs	r4, #0
 8000cd2:	f855 3b04 	ldr.w	r3, [r5], #4
 8000cd6:	3401      	adds	r4, #1
 8000cd8:	4798      	blx	r3
 8000cda:	42a6      	cmp	r6, r4
 8000cdc:	d8f9      	bhi.n	8000cd2 <__libc_init_array+0x12>
 8000cde:	f002 f917 	bl	8002f10 <_init>
 8000ce2:	4d09      	ldr	r5, [pc, #36]	@ (8000d08 <__libc_init_array+0x48>)
 8000ce4:	4b09      	ldr	r3, [pc, #36]	@ (8000d0c <__libc_init_array+0x4c>)
 8000ce6:	1b5e      	subs	r6, r3, r5
 8000ce8:	42ab      	cmp	r3, r5
 8000cea:	ea4f 06a6 	mov.w	r6, r6, asr #2
 8000cee:	d006      	beq.n	8000cfe <__libc_init_array+0x3e>
 8000cf0:	2400      	movs	r4, #0
 8000cf2:	f855 3b04 	ldr.w	r3, [r5], #4
 8000cf6:	3401      	adds	r4, #1
 8000cf8:	4798      	blx	r3
 8000cfa:	42a6      	cmp	r6, r4
 8000cfc:	d8f9      	bhi.n	8000cf2 <__libc_init_array+0x32>
 8000cfe:	bd70      	pop	{r4, r5, r6, pc}
 8000d00:	08027310 	.word	0x08027310
 8000d04:	08027310 	.word	0x08027310
 8000d08:	08027310 	.word	0x08027310
 8000d0c:	08027318 	.word	0x08027318

08000d10 <__malloc_lock>:
 8000d10:	4770      	bx	lr
 8000d12:	bf00      	nop

08000d14 <__malloc_unlock>:
 8000d14:	4770      	bx	lr
 8000d16:	bf00      	nop

08000d18 <exit>:
 8000d18:	b508      	push	{r3, lr}
 8000d1a:	2100      	movs	r1, #0
 8000d1c:	4604      	mov	r4, r0
 8000d1e:	f000 f835 	bl	8000d8c <__call_exitprocs>
 8000d22:	4b03      	ldr	r3, [pc, #12]	@ (8000d30 <exit+0x18>)
 8000d24:	681b      	ldr	r3, [r3, #0]
 8000d26:	b103      	cbz	r3, 8000d2a <exit+0x12>
 8000d28:	4798      	blx	r3
 8000d2a:	4620      	mov	r0, r4
 8000d2c:	f000 f8b8 	bl	8000ea0 <_exit>
 8000d30:	20006aec 	.word	0x20006aec

08000d34 <atexit>:
 8000d34:	2300      	movs	r3, #0
 8000d36:	4601      	mov	r1, r0
 8000d38:	461a      	mov	r2, r3
 8000d3a:	4618      	mov	r0, r3
 8000d3c:	f000 b87e 	b.w	8000e3c <__register_exitproc>

08000d40 <_sbrk_r>:
 8000d40:	b538      	push	{r3, r4, r5, lr}
 8000d42:	4d07      	ldr	r5, [pc, #28]	@ (8000d60 <_sbrk_r+0x20>)
 8000d44:	2200      	movs	r2, #0
 8000d46:	4604      	mov	r4, r0
 8000d48:	4608      	mov	r0, r1
 8000d4a:	602a      	str	r2, [r5, #0]
 8000d4c:	f002 f878 	bl	8002e40 <_sbrk>
 8000d50:	1c43      	adds	r3, r0, #1
 8000d52:	d000      	beq.n	8000d56 <_sbrk_r+0x16>
 8000d54:	bd38      	pop	{r3, r4, r5, pc}
 8000d56:	682b      	ldr	r3, [r5, #0]
 8000d58:	2b00      	cmp	r3, #0
 8000d5a:	d0fb      	beq.n	8000d54 <_sbrk_r+0x14>
 8000d5c:	6023      	str	r3, [r4, #0]
 8000d5e:	bd38      	pop	{r3, r4, r5, pc}
 8000d60:	20006af4 	.word	0x20006af4

08000d64 <__libc_fini_array>:
 8000d64:	b538      	push	{r3, r4, r5, lr}
 8000d66:	4d07      	ldr	r5, [pc, #28]	@ (8000d84 <__libc_fini_array+0x20>)
 8000d68:	4c07      	ldr	r4, [pc, #28]	@ (8000d88 <__libc_fini_array+0x24>)
 8000d6a:	1b2c      	subs	r4, r5, r4
 8000d6c:	10a4      	asrs	r4, r4, #2
 8000d6e:	d005      	beq.n	8000d7c <__libc_fini_array+0x18>
 8000d70:	3c01      	subs	r4, #1
 8000d72:	f855 3d04 	ldr.w	r3, [r5, #-4]!
 8000d76:	4798      	blx	r3
 8000d78:	2c00      	cmp	r4, #0
 8000d7a:	d1f9      	bne.n	8000d70 <__libc_fini_array+0xc>
 8000d7c:	e8bd 4038 	ldmia.w	sp!, {r3, r4, r5, lr}
 8000d80:	f002 b8cc 	b.w	8002f1c <_fini>
 8000d84:	0802731c 	.word	0x0802731c
 8000d88:	08027318 	.word	0x08027318

08000d8c <__call_exitprocs>:
 8000d8c:	e92d 4ff0 	stmdb	sp!, {r4, r5, r6, r7, r8, r9, sl, fp, lr}
 8000d90:	4f29      	ldr	r7, [pc, #164]	@ (8000e38 <__call_exitprocs+0xac>)
 8000d92:	b083      	sub	sp, #12
 8000d94:	683e      	ldr	r6, [r7, #0]
 8000d96:	9001      	str	r0, [sp, #4]
 8000d98:	b35e      	cbz	r6, 8000df2 <__call_exitprocs+0x66>
 8000d9a:	468b      	mov	fp, r1
 8000d9c:	f04f 0900 	mov.w	r9, #0
 8000da0:	f04f 0801 	mov.w	r8, #1
 8000da4:	6874      	ldr	r4, [r6, #4]
 8000da6:	1e65      	subs	r5, r4, #1
 8000da8:	d423      	bmi.n	8000df2 <__call_exitprocs+0x66>
 8000daa:	3401      	adds	r4, #1
 8000dac:	eb06 0484 	add.w	r4, r6, r4, lsl #2
 8000db0:	f1bb 0f00 	cmp.w	fp, #0
 8000db4:	d120      	bne.n	8000df8 <__call_exitprocs+0x6c>
 8000db6:	6873      	ldr	r3, [r6, #4]
 8000db8:	6822      	ldr	r2, [r4, #0]
 8000dba:	3b01      	subs	r3, #1
 8000dbc:	42ab      	cmp	r3, r5
 8000dbe:	bf0c      	ite	eq
 8000dc0:	6075      	streq	r5, [r6, #4]
 8000dc2:	f8c4 9000 	strne.w	r9, [r4]
 8000dc6:	b17a      	cbz	r2, 8000de8 <__call_exitprocs+0x5c>
 8000dc8:	f8d6 1188 	ldr.w	r1, [r6, #392]	@ 0x188
 8000dcc:	f8d6 a004 	ldr.w	sl, [r6, #4]
 8000dd0:	fa08 fc05 	lsl.w	ip, r8, r5
 8000dd4:	ea1c 0f01 	tst.w	ip, r1
 8000dd8:	d11a      	bne.n	8000e10 <__call_exitprocs+0x84>
 8000dda:	4790      	blx	r2
 8000ddc:	6871      	ldr	r1, [r6, #4]
 8000dde:	683a      	ldr	r2, [r7, #0]
 8000de0:	4551      	cmp	r1, sl
 8000de2:	d122      	bne.n	8000e2a <__call_exitprocs+0x9e>
 8000de4:	42b2      	cmp	r2, r6
 8000de6:	d120      	bne.n	8000e2a <__call_exitprocs+0x9e>
 8000de8:	3d01      	subs	r5, #1
 8000dea:	1c6b      	adds	r3, r5, #1
 8000dec:	f1a4 0404 	sub.w	r4, r4, #4
 8000df0:	d1de      	bne.n	8000db0 <__call_exitprocs+0x24>
 8000df2:	b003      	add	sp, #12
 8000df4:	e8bd 8ff0 	ldmia.w	sp!, {r4, r5, r6, r7, r8, r9, sl, fp, pc}
 8000df8:	f8d4 3100 	ldr.w	r3, [r4, #256]	@ 0x100
 8000dfc:	455b      	cmp	r3, fp
 8000dfe:	d0da      	beq.n	8000db6 <__call_exitprocs+0x2a>
 8000e00:	3d01      	subs	r5, #1
 8000e02:	1c6a      	adds	r2, r5, #1
 8000e04:	f1a4 0404 	sub.w	r4, r4, #4
 8000e08:	d1f6      	bne.n	8000df8 <__call_exitprocs+0x6c>
 8000e0a:	b003      	add	sp, #12
 8000e0c:	e8bd 8ff0 	ldmia.w	sp!, {r4, r5, r6, r7, r8, r9, sl, fp, pc}
 8000e10:	f8d6 318c 	ldr.w	r3, [r6, #396]	@ 0x18c
 8000e14:	f8d4 1080 	ldr.w	r1, [r4, #128]	@ 0x80
 8000e18:	ea1c 0f03 	tst.w	ip, r3
 8000e1c:	d109      	bne.n	8000e32 <__call_exitprocs+0xa6>
 8000e1e:	9801      	ldr	r0, [sp, #4]
 8000e20:	4790      	blx	r2
 8000e22:	6871      	ldr	r1, [r6, #4]
 8000e24:	683a      	ldr	r2, [r7, #0]
 8000e26:	4551      	cmp	r1, sl
 8000e28:	d0dc      	beq.n	8000de4 <__call_exitprocs+0x58>
 8000e2a:	2a00      	cmp	r2, #0
 8000e2c:	d0e1      	beq.n	8000df2 <__call_exitprocs+0x66>
 8000e2e:	4616      	mov	r6, r2
 8000e30:	e7b8      	b.n	8000da4 <__call_exitprocs+0x18>
 8000e32:	4608      	mov	r0, r1
 8000e34:	4790      	blx	r2
 8000e36:	e7d1      	b.n	8000ddc <__call_exitprocs+0x50>
 8000e38:	20006af0 	.word	0x20006af0

08000e3c <__register_exitproc>:
 8000e3c:	b470      	push	{r4, r5, r6}
 8000e3e:	4d16      	ldr	r5, [pc, #88]	@ (8000e98 <__register_exitproc+0x5c>)
 8000e40:	682c      	ldr	r4, [r5, #0]
 8000e42:	b31c      	cbz	r4, 8000e8c <__register_exitproc+0x50>
 8000e44:	6865      	ldr	r5, [r4, #4]
 8000e46:	2d1f      	cmp	r5, #31
 8000e48:	dc23      	bgt.n	8000e92 <__register_exitproc+0x56>
 8000e4a:	b938      	cbnz	r0, 8000e5c <__register_exitproc+0x20>
 8000e4c:	1cab      	adds	r3, r5, #2
 8000e4e:	3501      	adds	r5, #1
 8000e50:	6065      	str	r5, [r4, #4]
 8000e52:	f844 1023 	str.w	r1, [r4, r3, lsl #2]
 8000e56:	2000      	movs	r0, #0
 8000e58:	bc70      	pop	{r4, r5, r6}
 8000e5a:	4770      	bx	lr
 8000e5c:	eb04 0c85 	add.w	ip, r4, r5, lsl #2
 8000e60:	2802      	cmp	r0, #2
 8000e62:	f8cc 2088 	str.w	r2, [ip, #136]	@ 0x88
 8000e66:	f8d4 6188 	ldr.w	r6, [r4, #392]	@ 0x188
 8000e6a:	f04f 0201 	mov.w	r2, #1
 8000e6e:	fa02 f205 	lsl.w	r2, r2, r5
 8000e72:	ea46 0602 	orr.w	r6, r6, r2
 8000e76:	f8c4 6188 	str.w	r6, [r4, #392]	@ 0x188
 8000e7a:	f8cc 3108 	str.w	r3, [ip, #264]	@ 0x108
 8000e7e:	d1e5      	bne.n	8000e4c <__register_exitproc+0x10>
 8000e80:	f8d4 318c 	ldr.w	r3, [r4, #396]	@ 0x18c
 8000e84:	4313      	orrs	r3, r2
 8000e86:	f8c4 318c 	str.w	r3, [r4, #396]	@ 0x18c
 8000e8a:	e7df      	b.n	8000e4c <__register_exitproc+0x10>
 8000e8c:	4c03      	ldr	r4, [pc, #12]	@ (8000e9c <__register_exitproc+0x60>)
 8000e8e:	602c      	str	r4, [r5, #0]
 8000e90:	e7d8      	b.n	8000e44 <__register_exitproc+0x8>
 8000e92:	f04f 30ff 	mov.w	r0, #4294967295	@ 0xffffffff
 8000e96:	e7df      	b.n	8000e58 <__register_exitproc+0x1c>
 8000e98:	20006af0 	.word	0x20006af0
 8000e9c:	20006af8 	.word	0x20006af8

08000ea0 <_exit>:
 8000ea0:	e7fe      	b.n	8000ea0 <_exit>
 8000ea2:	bf00      	nop

08000ea4 <m_vec_mul_add.constprop.0>:
 8000ea4:	e92d 4ff7 	stmdb	sp!, {r0, r1, r2, r4, r5, r6, r7, r8, r9, sl, fp, lr}
 8000ea8:	eb01 2141 	add.w	r1, r1, r1, lsl #9
 8000eac:	eb01 4181 	add.w	r1, r1, r1, lsl #18
 8000eb0:	f001 33f0 	and.w	r3, r1, #4042322160	@ 0xf0f0f0f0
 8000eb4:	08dc      	lsrs	r4, r3, #3
 8000eb6:	ea84 1413 	eor.w	r4, r4, r3, lsr #4
 8000eba:	404c      	eors	r4, r1
 8000ebc:	f102 0320 	add.w	r3, r2, #32
 8000ec0:	f3c4 2c03 	ubfx	ip, r4, #8, #4
 8000ec4:	f3c4 4703 	ubfx	r7, r4, #16, #4
 8000ec8:	f3c4 6603 	ubfx	r6, r4, #24, #4
 8000ecc:	f1a2 0508 	sub.w	r5, r2, #8
 8000ed0:	3808      	subs	r0, #8
 8000ed2:	9301      	str	r3, [sp, #4]
 8000ed4:	b2e4      	uxtb	r4, r4
 8000ed6:	f850 1f08 	ldr.w	r1, [r0, #8]!
 8000eda:	f8d0 8004 	ldr.w	r8, [r0, #4]
 8000ede:	f001 3311 	and.w	r3, r1, #286331153	@ 0x11111111
 8000ee2:	f008 3911 	and.w	r9, r8, #286331153	@ 0x11111111
 8000ee6:	fba3 e304 	umull	lr, r3, r3, r4
 8000eea:	fb04 3309 	mla	r3, r4, r9, r3
 8000eee:	f855 9f08 	ldr.w	r9, [r5, #8]!
 8000ef2:	ea8e 0e09 	eor.w	lr, lr, r9
 8000ef6:	ea4f 0951 	mov.w	r9, r1, lsr #1
 8000efa:	686a      	ldr	r2, [r5, #4]
 8000efc:	ea4f 0a58 	mov.w	sl, r8, lsr #1
 8000f00:	f009 3911 	and.w	r9, r9, #286331153	@ 0x11111111
 8000f04:	fba9 b90c 	umull	fp, r9, r9, ip
 8000f08:	f00a 3a11 	and.w	sl, sl, #286331153	@ 0x11111111
 8000f0c:	fb0c 990a 	mla	r9, ip, sl, r9
 8000f10:	4053      	eors	r3, r2
 8000f12:	ea83 0309 	eor.w	r3, r3, r9
 8000f16:	ea4f 0991 	mov.w	r9, r1, lsr #2
 8000f1a:	ea4f 0a98 	mov.w	sl, r8, lsr #2
 8000f1e:	f009 3911 	and.w	r9, r9, #286331153	@ 0x11111111
 8000f22:	ea8e 0e0b 	eor.w	lr, lr, fp
 8000f26:	f00a 3a11 	and.w	sl, sl, #286331153	@ 0x11111111
 8000f2a:	fba9 b907 	umull	fp, r9, r9, r7
 8000f2e:	08c9      	lsrs	r1, r1, #3
 8000f30:	fb07 990a 	mla	r9, r7, sl, r9
 8000f34:	ea4f 08d8 	mov.w	r8, r8, lsr #3
 8000f38:	f001 3111 	and.w	r1, r1, #286331153	@ 0x11111111
 8000f3c:	ea83 0309 	eor.w	r3, r3, r9
 8000f40:	f008 3811 	and.w	r8, r8, #286331153	@ 0x11111111
 8000f44:	fba1 1906 	umull	r1, r9, r1, r6
 8000f48:	ea8e 0e0b 	eor.w	lr, lr, fp
 8000f4c:	fb06 9808 	mla	r8, r6, r8, r9
 8000f50:	ea83 0308 	eor.w	r3, r3, r8
 8000f54:	ea8e 0101 	eor.w	r1, lr, r1
 8000f58:	e9c5 1300 	strd	r1, r3, [r5]
 8000f5c:	9b01      	ldr	r3, [sp, #4]
 8000f5e:	42ab      	cmp	r3, r5
 8000f60:	d1b9      	bne.n	8000ed6 <m_vec_mul_add.constprop.0+0x32>
 8000f62:	b003      	add	sp, #12
 8000f64:	e8bd 8ff0 	ldmia.w	sp!, {r4, r5, r6, r7, r8, r9, sl, fp, pc}

08000f68 <pqmayo_MAYO_1_m4_mayo_keypair_compact>:
 8000f68:	e92d 4ff0 	stmdb	sp!, {r4, r5, r6, r7, r8, r9, sl, fp, lr}
 8000f6c:	f5ad 5d87 	sub.w	sp, sp, #4320	@ 0x10e0
 8000f70:	b087      	sub	sp, #28
 8000f72:	4615      	mov	r5, r2
 8000f74:	9000      	str	r0, [sp, #0]
 8000f76:	f44f 6220 	mov.w	r2, #2560	@ 0xa00
 8000f7a:	4688      	mov	r8, r1
 8000f7c:	f50d 60df 	add.w	r0, sp, #1784	@ 0x6f8
 8000f80:	2100      	movs	r1, #0
 8000f82:	461c      	mov	r4, r3
 8000f84:	f7ff f944 	bl	8000210 <memset>
 8000f88:	4b59      	ldr	r3, [pc, #356]	@ (80010f0 <pqmayo_MAYO_1_m4_mayo_keypair_compact+0x188>)
 8000f8a:	462a      	mov	r2, r5
 8000f8c:	f103 0118 	add.w	r1, r3, #24
 8000f90:	f853 0b04 	ldr.w	r0, [r3], #4
 8000f94:	f842 0b04 	str.w	r0, [r2], #4
 8000f98:	428b      	cmp	r3, r1
 8000f9a:	d1f9      	bne.n	8000f90 <pqmayo_MAYO_1_m4_mayo_keypair_compact+0x28>
 8000f9c:	2318      	movs	r3, #24
 8000f9e:	462a      	mov	r2, r5
 8000fa0:	f44f 71a4 	mov.w	r1, #328	@ 0x148
 8000fa4:	a804      	add	r0, sp, #16
 8000fa6:	f001 f8cd 	bl	8002144 <shake256>
 8000faa:	1ca3      	adds	r3, r4, #2
 8000fac:	f10d 021f 	add.w	r2, sp, #31
 8000fb0:	f204 2072 	addw	r0, r4, #626	@ 0x272
 8000fb4:	f812 1f01 	ldrb.w	r1, [r2, #1]!
 8000fb8:	f001 050f 	and.w	r5, r1, #15
 8000fbc:	0909      	lsrs	r1, r1, #4
 8000fbe:	f803 5c02 	strb.w	r5, [r3, #-2]
 8000fc2:	f803 1c01 	strb.w	r1, [r3, #-1]
 8000fc6:	3302      	adds	r3, #2
 8000fc8:	4283      	cmp	r3, r0
 8000fca:	d1f3      	bne.n	8000fb4 <pqmayo_MAYO_1_m4_mayo_keypair_compact+0x4c>
 8000fcc:	f001 faee 	bl	80025ac <trigger_high>
 8000fd0:	2600      	movs	r6, #0
 8000fd2:	4635      	mov	r5, r6
 8000fd4:	4b47      	ldr	r3, [pc, #284]	@ (80010f4 <pqmayo_MAYO_1_m4_mayo_keypair_compact+0x18c>)
 8000fd6:	2228      	movs	r2, #40	@ 0x28
 8000fd8:	fb02 3006 	mla	r0, r2, r6, r3
 8000fdc:	f50d 5389 	add.w	r3, sp, #4384	@ 0x1120
 8000fe0:	f44f 72a0 	mov.w	r2, #320	@ 0x140
 8000fe4:	681b      	ldr	r3, [r3, #0]
 8000fe6:	00ef      	lsls	r7, r5, #3
 8000fe8:	fb02 3b05 	mla	fp, r2, r5, r3
 8000fec:	f104 0a08 	add.w	sl, r4, #8
 8000ff0:	e013      	b.n	800101a <pqmayo_MAYO_1_m4_mayo_keypair_compact+0xb2>
 8000ff2:	f819 1b01 	ldrb.w	r1, [r9], #1
 8000ff6:	9001      	str	r0, [sp, #4]
 8000ff8:	e9cd 2302 	strd	r2, r3, [sp, #8]
 8000ffc:	f7ff ff52 	bl	8000ea4 <m_vec_mul_add.constprop.0>
 8001000:	9b03      	ldr	r3, [sp, #12]
 8001002:	9a02      	ldr	r2, [sp, #8]
 8001004:	9801      	ldr	r0, [sp, #4]
 8001006:	4599      	cmp	r9, r3
 8001008:	f102 0228 	add.w	r2, r2, #40	@ 0x28
 800100c:	d1f1      	bne.n	8000ff2 <pqmayo_MAYO_1_m4_mayo_keypair_compact+0x8a>
 800100e:	3708      	adds	r7, #8
 8001010:	f5b7 7f1c 	cmp.w	r7, #624	@ 0x270
 8001014:	f100 0028 	add.w	r0, r0, #40	@ 0x28
 8001018:	d005      	beq.n	8001026 <pqmayo_MAYO_1_m4_mayo_keypair_compact+0xbe>
 800101a:	eb04 0907 	add.w	r9, r4, r7
 800101e:	465a      	mov	r2, fp
 8001020:	eb0a 0307 	add.w	r3, sl, r7
 8001024:	e7e5      	b.n	8000ff2 <pqmayo_MAYO_1_m4_mayo_keypair_compact+0x8a>
 8001026:	364e      	adds	r6, #78	@ 0x4e
 8001028:	1b76      	subs	r6, r6, r5
 800102a:	3501      	adds	r5, #1
 800102c:	2d4e      	cmp	r5, #78	@ 0x4e
 800102e:	d1d1      	bne.n	8000fd4 <pqmayo_MAYO_1_m4_mayo_keypair_compact+0x6c>
 8001030:	f001 fac4 	bl	80025bc <trigger_low>
 8001034:	f50d 69df 	add.w	r9, sp, #1784	@ 0x6f8
 8001038:	2700      	movs	r7, #0
 800103a:	f04f 0b28 	mov.w	fp, #40	@ 0x28
 800103e:	2500      	movs	r5, #0
 8001040:	eb04 0a07 	add.w	sl, r4, r7
 8001044:	e010      	b.n	8001068 <pqmayo_MAYO_1_m4_mayo_keypair_compact+0x100>
 8001046:	1998      	adds	r0, r3, r6
 8001048:	eb09 0206 	add.w	r2, r9, r6
 800104c:	f81a 1005 	ldrb.w	r1, [sl, r5]
 8001050:	9301      	str	r3, [sp, #4]
 8001052:	3628      	adds	r6, #40	@ 0x28
 8001054:	f7ff ff26 	bl	8000ea4 <m_vec_mul_add.constprop.0>
 8001058:	f5b6 7fa0 	cmp.w	r6, #320	@ 0x140
 800105c:	9b01      	ldr	r3, [sp, #4]
 800105e:	d1f2      	bne.n	8001046 <pqmayo_MAYO_1_m4_mayo_keypair_compact+0xde>
 8001060:	3508      	adds	r5, #8
 8001062:	f5b5 7f1c 	cmp.w	r5, #624	@ 0x270
 8001066:	d006      	beq.n	8001076 <pqmayo_MAYO_1_m4_mayo_keypair_compact+0x10e>
 8001068:	f50d 5389 	add.w	r3, sp, #4384	@ 0x1120
 800106c:	2600      	movs	r6, #0
 800106e:	681b      	ldr	r3, [r3, #0]
 8001070:	fb0b 3305 	mla	r3, fp, r5, r3
 8001074:	e7e7      	b.n	8001046 <pqmayo_MAYO_1_m4_mayo_keypair_compact+0xde>
 8001076:	3701      	adds	r7, #1
 8001078:	2f08      	cmp	r7, #8
 800107a:	f509 79a0 	add.w	r9, r9, #320	@ 0x140
 800107e:	d1de      	bne.n	800103e <pqmayo_MAYO_1_m4_mayo_keypair_compact+0xd6>
 8001080:	ab04      	add	r3, sp, #16
 8001082:	4644      	mov	r4, r8
 8001084:	ad08      	add	r5, sp, #32
 8001086:	461a      	mov	r2, r3
 8001088:	ca03      	ldmia	r2!, {r0, r1}
 800108a:	42aa      	cmp	r2, r5
 800108c:	6020      	str	r0, [r4, #0]
 800108e:	6061      	str	r1, [r4, #4]
 8001090:	4613      	mov	r3, r2
 8001092:	f104 0408 	add.w	r4, r4, #8
 8001096:	d1f6      	bne.n	8001086 <pqmayo_MAYO_1_m4_mayo_keypair_compact+0x11e>
 8001098:	aa56      	add	r2, sp, #344	@ 0x158
 800109a:	9800      	ldr	r0, [sp, #0]
 800109c:	2308      	movs	r3, #8
 800109e:	f50d 61df 	add.w	r1, sp, #1784	@ 0x6f8
 80010a2:	f001 f872 	bl	800218a <pqmayo_MAYO_1_m4_m_upper>
 80010a6:	2600      	movs	r6, #0
 80010a8:	aa56      	add	r2, sp, #344	@ 0x158
 80010aa:	f640 27f8 	movw	r7, #2808	@ 0xaf8
 80010ae:	eb08 0366 	add.w	r3, r8, r6, asr #1
 80010b2:	4614      	mov	r4, r2
 80010b4:	3310      	adds	r3, #16
 80010b6:	f102 0c20 	add.w	ip, r2, #32
 80010ba:	4625      	mov	r5, r4
 80010bc:	cd03      	ldmia	r5!, {r0, r1}
 80010be:	4565      	cmp	r5, ip
 80010c0:	6018      	str	r0, [r3, #0]
 80010c2:	6059      	str	r1, [r3, #4]
 80010c4:	462c      	mov	r4, r5
 80010c6:	f103 0308 	add.w	r3, r3, #8
 80010ca:	d1f6      	bne.n	80010ba <pqmayo_MAYO_1_m4_mayo_keypair_compact+0x152>
 80010cc:	6828      	ldr	r0, [r5, #0]
 80010ce:	6018      	str	r0, [r3, #0]
 80010d0:	364e      	adds	r6, #78	@ 0x4e
 80010d2:	88a8      	ldrh	r0, [r5, #4]
 80010d4:	79a9      	ldrb	r1, [r5, #6]
 80010d6:	7199      	strb	r1, [r3, #6]
 80010d8:	42be      	cmp	r6, r7
 80010da:	8098      	strh	r0, [r3, #4]
 80010dc:	f102 0228 	add.w	r2, r2, #40	@ 0x28
 80010e0:	d1e5      	bne.n	80010ae <pqmayo_MAYO_1_m4_mayo_keypair_compact+0x146>
 80010e2:	2000      	movs	r0, #0
 80010e4:	f50d 5d87 	add.w	sp, sp, #4320	@ 0x10e0
 80010e8:	b007      	add	sp, #28
 80010ea:	e8bd 8ff0 	ldmia.w	sp!, {r4, r5, r6, r7, r8, r9, sl, fp, pc}
 80010ee:	bf00      	nop
 80010f0:	08027220 	.word	0x08027220
 80010f4:	08002f28 	.word	0x08002f28

080010f8 <pqmayo_MAYO_1_m4_mayo_keypair>:
 80010f8:	f7ff bf36 	b.w	8000f68 <pqmayo_MAYO_1_m4_mayo_keypair_compact>

080010fc <example_fault_P3_OtP2.constprop.0>:
 80010fc:	b570      	push	{r4, r5, r6, lr}
 80010fe:	f240 518c 	movw	r1, #1420	@ 0x58c
 8001102:	f5ad 7d1e 	sub.w	sp, sp, #632	@ 0x278
 8001106:	2001      	movs	r0, #1
 8001108:	f7ff f8d4 	bl	80002b4 <calloc>
 800110c:	2118      	movs	r1, #24
 800110e:	4605      	mov	r5, r0
 8001110:	2001      	movs	r0, #1
 8001112:	f7ff f8cf 	bl	80002b4 <calloc>
 8001116:	4b11      	ldr	r3, [pc, #68]	@ (800115c <example_fault_P3_OtP2.constprop.0+0x60>)
 8001118:	4911      	ldr	r1, [pc, #68]	@ (8001160 <example_fault_P3_OtP2.constprop.0+0x64>)
 800111a:	4604      	mov	r4, r0
 800111c:	f44f 42c3 	mov.w	r2, #24960	@ 0x6180
 8001120:	4618      	mov	r0, r3
 8001122:	f7ff fd25 	bl	8000b70 <memcpy>
 8001126:	4622      	mov	r2, r4
 8001128:	9000      	str	r0, [sp, #0]
 800112a:	ab02      	add	r3, sp, #8
 800112c:	4629      	mov	r1, r5
 800112e:	2000      	movs	r0, #0
 8001130:	f7ff ffe2 	bl	80010f8 <pqmayo_MAYO_1_m4_mayo_keypair>
 8001134:	2400      	movs	r4, #0
 8001136:	f240 568c 	movw	r6, #1420	@ 0x58c
 800113a:	1b31      	subs	r1, r6, r4
 800113c:	29c8      	cmp	r1, #200	@ 0xc8
 800113e:	bf28      	it	cs
 8001140:	21c8      	movcs	r1, #200	@ 0xc8
 8001142:	192a      	adds	r2, r5, r4
 8001144:	b2c9      	uxtb	r1, r1
 8001146:	2070      	movs	r0, #112	@ 0x70
 8001148:	34c8      	adds	r4, #200	@ 0xc8
 800114a:	f001 f98d 	bl	8002468 <simpleserial_put>
 800114e:	f5b4 6fc8 	cmp.w	r4, #1600	@ 0x640
 8001152:	d1f2      	bne.n	800113a <example_fault_P3_OtP2.constprop.0+0x3e>
 8001154:	f50d 7d1e 	add.w	sp, sp, #632	@ 0x278
 8001158:	bd70      	pop	{r4, r5, r6, pc}
 800115a:	bf00      	nop
 800115c:	20000558 	.word	0x20000558
 8001160:	08021090 	.word	0x08021090

08001164 <cmd_keygen>:
 8001164:	b508      	push	{r3, lr}
 8001166:	f7ff ffc9 	bl	80010fc <example_fault_P3_OtP2.constprop.0>
 800116a:	2000      	movs	r0, #0
 800116c:	bd08      	pop	{r3, pc}
	...

08001170 <main>:
 8001170:	b508      	push	{r3, lr}
 8001172:	f001 f97d 	bl	8002470 <platform_init>
 8001176:	f001 f9bb 	bl	80024f0 <init_uart>
 800117a:	f001 f9f7 	bl	800256c <trigger_setup>
 800117e:	f001 f8ff 	bl	8002380 <simpleserial_init>
 8001182:	4a04      	ldr	r2, [pc, #16]	@ (8001194 <main+0x24>)
 8001184:	2100      	movs	r1, #0
 8001186:	206b      	movs	r0, #107	@ 0x6b
 8001188:	f001 f8f6 	bl	8002378 <simpleserial_addcmd>
 800118c:	f001 f90e 	bl	80023ac <simpleserial_get>
 8001190:	e7fc      	b.n	800118c <main+0x1c>
 8001192:	bf00      	nop
 8001194:	08001165 	.word	0x08001165

08001198 <KeccakF1600_StatePermute>:
 8001198:	e92d 4ff0 	stmdb	sp!, {r4, r5, r6, r7, r8, r9, sl, fp, lr}
 800119c:	b0bd      	sub	sp, #244	@ 0xf4
 800119e:	6803      	ldr	r3, [r0, #0]
 80011a0:	930a      	str	r3, [sp, #40]	@ 0x28
 80011a2:	6843      	ldr	r3, [r0, #4]
 80011a4:	930b      	str	r3, [sp, #44]	@ 0x2c
 80011a6:	6883      	ldr	r3, [r0, #8]
 80011a8:	930c      	str	r3, [sp, #48]	@ 0x30
 80011aa:	68c3      	ldr	r3, [r0, #12]
 80011ac:	930d      	str	r3, [sp, #52]	@ 0x34
 80011ae:	6903      	ldr	r3, [r0, #16]
 80011b0:	930e      	str	r3, [sp, #56]	@ 0x38
 80011b2:	6943      	ldr	r3, [r0, #20]
 80011b4:	930f      	str	r3, [sp, #60]	@ 0x3c
 80011b6:	6983      	ldr	r3, [r0, #24]
 80011b8:	9310      	str	r3, [sp, #64]	@ 0x40
 80011ba:	69c3      	ldr	r3, [r0, #28]
 80011bc:	9311      	str	r3, [sp, #68]	@ 0x44
 80011be:	6a03      	ldr	r3, [r0, #32]
 80011c0:	9312      	str	r3, [sp, #72]	@ 0x48
 80011c2:	6a43      	ldr	r3, [r0, #36]	@ 0x24
 80011c4:	9313      	str	r3, [sp, #76]	@ 0x4c
 80011c6:	6a83      	ldr	r3, [r0, #40]	@ 0x28
 80011c8:	9314      	str	r3, [sp, #80]	@ 0x50
 80011ca:	6ac3      	ldr	r3, [r0, #44]	@ 0x2c
 80011cc:	9315      	str	r3, [sp, #84]	@ 0x54
 80011ce:	6b03      	ldr	r3, [r0, #48]	@ 0x30
 80011d0:	9316      	str	r3, [sp, #88]	@ 0x58
 80011d2:	6b43      	ldr	r3, [r0, #52]	@ 0x34
 80011d4:	9317      	str	r3, [sp, #92]	@ 0x5c
 80011d6:	6b83      	ldr	r3, [r0, #56]	@ 0x38
 80011d8:	9318      	str	r3, [sp, #96]	@ 0x60
 80011da:	6bc3      	ldr	r3, [r0, #60]	@ 0x3c
 80011dc:	9319      	str	r3, [sp, #100]	@ 0x64
 80011de:	6c03      	ldr	r3, [r0, #64]	@ 0x40
 80011e0:	931a      	str	r3, [sp, #104]	@ 0x68
 80011e2:	6c43      	ldr	r3, [r0, #68]	@ 0x44
 80011e4:	931b      	str	r3, [sp, #108]	@ 0x6c
 80011e6:	6c83      	ldr	r3, [r0, #72]	@ 0x48
 80011e8:	931c      	str	r3, [sp, #112]	@ 0x70
 80011ea:	6cc3      	ldr	r3, [r0, #76]	@ 0x4c
 80011ec:	931d      	str	r3, [sp, #116]	@ 0x74
 80011ee:	6d03      	ldr	r3, [r0, #80]	@ 0x50
 80011f0:	931e      	str	r3, [sp, #120]	@ 0x78
 80011f2:	6d43      	ldr	r3, [r0, #84]	@ 0x54
 80011f4:	931f      	str	r3, [sp, #124]	@ 0x7c
 80011f6:	6d83      	ldr	r3, [r0, #88]	@ 0x58
 80011f8:	9320      	str	r3, [sp, #128]	@ 0x80
 80011fa:	6dc3      	ldr	r3, [r0, #92]	@ 0x5c
 80011fc:	9321      	str	r3, [sp, #132]	@ 0x84
 80011fe:	6e03      	ldr	r3, [r0, #96]	@ 0x60
 8001200:	9322      	str	r3, [sp, #136]	@ 0x88
 8001202:	6e43      	ldr	r3, [r0, #100]	@ 0x64
 8001204:	9323      	str	r3, [sp, #140]	@ 0x8c
 8001206:	6e83      	ldr	r3, [r0, #104]	@ 0x68
 8001208:	9324      	str	r3, [sp, #144]	@ 0x90
 800120a:	6ec3      	ldr	r3, [r0, #108]	@ 0x6c
 800120c:	9325      	str	r3, [sp, #148]	@ 0x94
 800120e:	6f03      	ldr	r3, [r0, #112]	@ 0x70
 8001210:	9326      	str	r3, [sp, #152]	@ 0x98
 8001212:	6f43      	ldr	r3, [r0, #116]	@ 0x74
 8001214:	9327      	str	r3, [sp, #156]	@ 0x9c
 8001216:	6f83      	ldr	r3, [r0, #120]	@ 0x78
 8001218:	9328      	str	r3, [sp, #160]	@ 0xa0
 800121a:	6fc3      	ldr	r3, [r0, #124]	@ 0x7c
 800121c:	9329      	str	r3, [sp, #164]	@ 0xa4
 800121e:	f8d0 3080 	ldr.w	r3, [r0, #128]	@ 0x80
 8001222:	932a      	str	r3, [sp, #168]	@ 0xa8
 8001224:	f8d0 3084 	ldr.w	r3, [r0, #132]	@ 0x84
 8001228:	932b      	str	r3, [sp, #172]	@ 0xac
 800122a:	f8d0 3088 	ldr.w	r3, [r0, #136]	@ 0x88
 800122e:	932c      	str	r3, [sp, #176]	@ 0xb0
 8001230:	f8d0 308c 	ldr.w	r3, [r0, #140]	@ 0x8c
 8001234:	932d      	str	r3, [sp, #180]	@ 0xb4
 8001236:	f8d0 3090 	ldr.w	r3, [r0, #144]	@ 0x90
 800123a:	932e      	str	r3, [sp, #184]	@ 0xb8
 800123c:	f8d0 3094 	ldr.w	r3, [r0, #148]	@ 0x94
 8001240:	932f      	str	r3, [sp, #188]	@ 0xbc
 8001242:	f8d0 3098 	ldr.w	r3, [r0, #152]	@ 0x98
 8001246:	9330      	str	r3, [sp, #192]	@ 0xc0
 8001248:	f8d0 309c 	ldr.w	r3, [r0, #156]	@ 0x9c
 800124c:	9331      	str	r3, [sp, #196]	@ 0xc4
 800124e:	f8d0 30a0 	ldr.w	r3, [r0, #160]	@ 0xa0
 8001252:	9332      	str	r3, [sp, #200]	@ 0xc8
 8001254:	f8d0 30a4 	ldr.w	r3, [r0, #164]	@ 0xa4
 8001258:	9333      	str	r3, [sp, #204]	@ 0xcc
 800125a:	f8d0 30a8 	ldr.w	r3, [r0, #168]	@ 0xa8
 800125e:	9334      	str	r3, [sp, #208]	@ 0xd0
 8001260:	f8d0 30ac 	ldr.w	r3, [r0, #172]	@ 0xac
 8001264:	9335      	str	r3, [sp, #212]	@ 0xd4
 8001266:	f8d0 30b0 	ldr.w	r3, [r0, #176]	@ 0xb0
 800126a:	9336      	str	r3, [sp, #216]	@ 0xd8
 800126c:	f8d0 30b4 	ldr.w	r3, [r0, #180]	@ 0xb4
 8001270:	9337      	str	r3, [sp, #220]	@ 0xdc
 8001272:	f8d0 30b8 	ldr.w	r3, [r0, #184]	@ 0xb8
 8001276:	9338      	str	r3, [sp, #224]	@ 0xe0
 8001278:	f8d0 30bc 	ldr.w	r3, [r0, #188]	@ 0xbc
 800127c:	9339      	str	r3, [sp, #228]	@ 0xe4
 800127e:	f8d0 30c0 	ldr.w	r3, [r0, #192]	@ 0xc0
 8001282:	933a      	str	r3, [sp, #232]	@ 0xe8
 8001284:	f8d0 60c4 	ldr.w	r6, [r0, #196]	@ 0xc4
 8001288:	4b01      	ldr	r3, [pc, #4]	@ (8001290 <KeccakF1600_StatePermute+0xf8>)
 800128a:	9300      	str	r3, [sp, #0]
 800128c:	e002      	b.n	8001294 <KeccakF1600_StatePermute+0xfc>
 800128e:	bf00      	nop
 8001290:	08027238 	.word	0x08027238
 8001294:	9a14      	ldr	r2, [sp, #80]	@ 0x50
 8001296:	9b0a      	ldr	r3, [sp, #40]	@ 0x28
 8001298:	9915      	ldr	r1, [sp, #84]	@ 0x54
 800129a:	4053      	eors	r3, r2
 800129c:	9a0b      	ldr	r2, [sp, #44]	@ 0x2c
 800129e:	ea82 0501 	eor.w	r5, r2, r1
 80012a2:	9a1e      	ldr	r2, [sp, #120]	@ 0x78
 80012a4:	4053      	eors	r3, r2
 80012a6:	9a1f      	ldr	r2, [sp, #124]	@ 0x7c
 80012a8:	4055      	eors	r5, r2
 80012aa:	9a28      	ldr	r2, [sp, #160]	@ 0xa0
 80012ac:	4053      	eors	r3, r2
 80012ae:	9a29      	ldr	r2, [sp, #164]	@ 0xa4
 80012b0:	4055      	eors	r5, r2
 80012b2:	9a32      	ldr	r2, [sp, #200]	@ 0xc8
 80012b4:	4053      	eors	r3, r2
 80012b6:	9a33      	ldr	r2, [sp, #204]	@ 0xcc
 80012b8:	9308      	str	r3, [sp, #32]
 80012ba:	4055      	eors	r5, r2
 80012bc:	9b0c      	ldr	r3, [sp, #48]	@ 0x30
 80012be:	9a16      	ldr	r2, [sp, #88]	@ 0x58
 80012c0:	ea83 0902 	eor.w	r9, r3, r2
 80012c4:	9b0d      	ldr	r3, [sp, #52]	@ 0x34
 80012c6:	9a17      	ldr	r2, [sp, #92]	@ 0x5c
 80012c8:	ea83 0802 	eor.w	r8, r3, r2
 80012cc:	9b20      	ldr	r3, [sp, #128]	@ 0x80
 80012ce:	9a18      	ldr	r2, [sp, #96]	@ 0x60
 80012d0:	ea89 0903 	eor.w	r9, r9, r3
 80012d4:	9b21      	ldr	r3, [sp, #132]	@ 0x84
 80012d6:	ea88 0803 	eor.w	r8, r8, r3
 80012da:	9b2a      	ldr	r3, [sp, #168]	@ 0xa8
 80012dc:	ea89 0903 	eor.w	r9, r9, r3
 80012e0:	9b2b      	ldr	r3, [sp, #172]	@ 0xac
 80012e2:	ea88 0803 	eor.w	r8, r8, r3
 80012e6:	9b34      	ldr	r3, [sp, #208]	@ 0xd0
 80012e8:	ea89 0903 	eor.w	r9, r9, r3
 80012ec:	9b35      	ldr	r3, [sp, #212]	@ 0xd4
 80012ee:	ea88 0803 	eor.w	r8, r8, r3
 80012f2:	9b0e      	ldr	r3, [sp, #56]	@ 0x38
 80012f4:	ea83 0c02 	eor.w	ip, r3, r2
 80012f8:	9b0f      	ldr	r3, [sp, #60]	@ 0x3c
 80012fa:	9a19      	ldr	r2, [sp, #100]	@ 0x64
 80012fc:	ea83 0702 	eor.w	r7, r3, r2
 8001300:	9b22      	ldr	r3, [sp, #136]	@ 0x88
 8001302:	9a23      	ldr	r2, [sp, #140]	@ 0x8c
 8001304:	ea8c 0c03 	eor.w	ip, ip, r3
 8001308:	4057      	eors	r7, r2
 800130a:	9b2c      	ldr	r3, [sp, #176]	@ 0xb0
 800130c:	9a2d      	ldr	r2, [sp, #180]	@ 0xb4
 800130e:	ea8c 0c03 	eor.w	ip, ip, r3
 8001312:	4057      	eors	r7, r2
 8001314:	9b36      	ldr	r3, [sp, #216]	@ 0xd8
 8001316:	9a37      	ldr	r2, [sp, #220]	@ 0xdc
 8001318:	ea8c 0c03 	eor.w	ip, ip, r3
 800131c:	4057      	eors	r7, r2
 800131e:	9b10      	ldr	r3, [sp, #64]	@ 0x40
 8001320:	9a1a      	ldr	r2, [sp, #104]	@ 0x68
 8001322:	f8dd e074 	ldr.w	lr, [sp, #116]	@ 0x74
 8001326:	ea83 0402 	eor.w	r4, r3, r2
 800132a:	9b11      	ldr	r3, [sp, #68]	@ 0x44
 800132c:	9a1b      	ldr	r2, [sp, #108]	@ 0x6c
 800132e:	ea83 0102 	eor.w	r1, r3, r2
 8001332:	9a24      	ldr	r2, [sp, #144]	@ 0x90
 8001334:	9b12      	ldr	r3, [sp, #72]	@ 0x48
 8001336:	4054      	eors	r4, r2
 8001338:	9a25      	ldr	r2, [sp, #148]	@ 0x94
 800133a:	4051      	eors	r1, r2
 800133c:	9a2e      	ldr	r2, [sp, #184]	@ 0xb8
 800133e:	4054      	eors	r4, r2
 8001340:	9a2f      	ldr	r2, [sp, #188]	@ 0xbc
 8001342:	4051      	eors	r1, r2
 8001344:	9a38      	ldr	r2, [sp, #224]	@ 0xe0
 8001346:	4054      	eors	r4, r2
 8001348:	9a39      	ldr	r2, [sp, #228]	@ 0xe4
 800134a:	4051      	eors	r1, r2
 800134c:	9a1c      	ldr	r2, [sp, #112]	@ 0x70
 800134e:	405a      	eors	r2, r3
 8001350:	9b13      	ldr	r3, [sp, #76]	@ 0x4c
 8001352:	ea83 0e0e 	eor.w	lr, r3, lr
 8001356:	9b26      	ldr	r3, [sp, #152]	@ 0x98
 8001358:	405a      	eors	r2, r3
 800135a:	9b27      	ldr	r3, [sp, #156]	@ 0x9c
 800135c:	ea8e 0e03 	eor.w	lr, lr, r3
 8001360:	9b30      	ldr	r3, [sp, #192]	@ 0xc0
 8001362:	405a      	eors	r2, r3
 8001364:	9b31      	ldr	r3, [sp, #196]	@ 0xc4
 8001366:	ea8e 0e03 	eor.w	lr, lr, r3
 800136a:	9b3a      	ldr	r3, [sp, #232]	@ 0xe8
 800136c:	ea4f 0b49 	mov.w	fp, r9, lsl #1
 8001370:	405a      	eors	r2, r3
 8001372:	ea4b 7bd8 	orr.w	fp, fp, r8, lsr #31
 8001376:	ea4f 0a48 	mov.w	sl, r8, lsl #1
 800137a:	ea8e 0e06 	eor.w	lr, lr, r6
 800137e:	ea8b 0302 	eor.w	r3, fp, r2
 8001382:	ea4a 7ad9 	orr.w	sl, sl, r9, lsr #31
 8001386:	9301      	str	r3, [sp, #4]
 8001388:	ea8a 030e 	eor.w	r3, sl, lr
 800138c:	9302      	str	r3, [sp, #8]
 800138e:	ea4f 0b4c 	mov.w	fp, ip, lsl #1
 8001392:	9b08      	ldr	r3, [sp, #32]
 8001394:	ea4b 7bd7 	orr.w	fp, fp, r7, lsr #31
 8001398:	ea4f 0a47 	mov.w	sl, r7, lsl #1
 800139c:	ea4a 7adc 	orr.w	sl, sl, ip, lsr #31
 80013a0:	ea8b 0303 	eor.w	r3, fp, r3
 80013a4:	9303      	str	r3, [sp, #12]
 80013a6:	ea8a 0305 	eor.w	r3, sl, r5
 80013aa:	ea4f 0a41 	mov.w	sl, r1, lsl #1
 80013ae:	ea4a 7ad4 	orr.w	sl, sl, r4, lsr #31
 80013b2:	9304      	str	r3, [sp, #16]
 80013b4:	ea8a 0308 	eor.w	r3, sl, r8
 80013b8:	ea4f 084e 	mov.w	r8, lr, lsl #1
 80013bc:	ea48 78d2 	orr.w	r8, r8, r2, lsr #31
 80013c0:	0052      	lsls	r2, r2, #1
 80013c2:	ea42 72de 	orr.w	r2, r2, lr, lsr #31
 80013c6:	9305      	str	r3, [sp, #20]
 80013c8:	ea82 030c 	eor.w	r3, r2, ip
 80013cc:	9306      	str	r3, [sp, #24]
 80013ce:	ea88 0307 	eor.w	r3, r8, r7
 80013d2:	9307      	str	r3, [sp, #28]
 80013d4:	9b08      	ldr	r3, [sp, #32]
 80013d6:	006a      	lsls	r2, r5, #1
 80013d8:	ea42 72d3 	orr.w	r2, r2, r3, lsr #31
 80013dc:	005b      	lsls	r3, r3, #1
 80013de:	ea43 73d5 	orr.w	r3, r3, r5, lsr #31
 80013e2:	4063      	eors	r3, r4
 80013e4:	9308      	str	r3, [sp, #32]
 80013e6:	ea82 0301 	eor.w	r3, r2, r1
 80013ea:	ea4f 0b44 	mov.w	fp, r4, lsl #1
 80013ee:	9a0a      	ldr	r2, [sp, #40]	@ 0x28
 80013f0:	9309      	str	r3, [sp, #36]	@ 0x24
 80013f2:	9b01      	ldr	r3, [sp, #4]
 80013f4:	ea4b 7bd1 	orr.w	fp, fp, r1, lsr #31
 80013f8:	ea8b 0b09 	eor.w	fp, fp, r9
 80013fc:	ea83 0902 	eor.w	r9, r3, r2
 8001400:	9b02      	ldr	r3, [sp, #8]
 8001402:	9a0b      	ldr	r2, [sp, #44]	@ 0x2c
 8001404:	9917      	ldr	r1, [sp, #92]	@ 0x5c
 8001406:	9c23      	ldr	r4, [sp, #140]	@ 0x8c
 8001408:	9f2f      	ldr	r7, [sp, #188]	@ 0xbc
 800140a:	ea83 0a02 	eor.w	sl, r3, r2
 800140e:	9b03      	ldr	r3, [sp, #12]
 8001410:	9a16      	ldr	r2, [sp, #88]	@ 0x58
 8001412:	405a      	eors	r2, r3
 8001414:	9b04      	ldr	r3, [sp, #16]
 8001416:	0d15      	lsrs	r5, r2, #20
 8001418:	404b      	eors	r3, r1
 800141a:	ea45 3503 	orr.w	r5, r5, r3, lsl #12
 800141e:	0d1b      	lsrs	r3, r3, #20
 8001420:	ea43 3302 	orr.w	r3, r3, r2, lsl #12
 8001424:	9a22      	ldr	r2, [sp, #136]	@ 0x88
 8001426:	9905      	ldr	r1, [sp, #20]
 8001428:	ea8b 0202 	eor.w	r2, fp, r2
 800142c:	404c      	eors	r4, r1
 800142e:	ea4f 5852 	mov.w	r8, r2, lsr #21
 8001432:	ea48 28c4 	orr.w	r8, r8, r4, lsl #11
 8001436:	0d64      	lsrs	r4, r4, #21
 8001438:	ea44 24c2 	orr.w	r4, r4, r2, lsl #11
 800143c:	992e      	ldr	r1, [sp, #184]	@ 0xb8
 800143e:	9a06      	ldr	r2, [sp, #24]
 8001440:	4051      	eors	r1, r2
 8001442:	9a07      	ldr	r2, [sp, #28]
 8001444:	407a      	eors	r2, r7
 8001446:	ea4f 5e42 	mov.w	lr, r2, lsl #21
 800144a:	ea4e 2ed1 	orr.w	lr, lr, r1, lsr #11
 800144e:	0549      	lsls	r1, r1, #21
 8001450:	9f3a      	ldr	r7, [sp, #232]	@ 0xe8
 8001452:	ea41 21d2 	orr.w	r1, r1, r2, lsr #11
 8001456:	9a08      	ldr	r2, [sp, #32]
 8001458:	407a      	eors	r2, r7
 800145a:	9f09      	ldr	r7, [sp, #36]	@ 0x24
 800145c:	407e      	eors	r6, r7
 800145e:	ea4f 3c86 	mov.w	ip, r6, lsl #14
 8001462:	ea4c 4c92 	orr.w	ip, ip, r2, lsr #18
 8001466:	0392      	lsls	r2, r2, #14
 8001468:	ea42 4296 	orr.w	r2, r2, r6, lsr #18
 800146c:	9e00      	ldr	r6, [sp, #0]
 800146e:	e9d6 7600 	ldrd	r7, r6, [r6]
 8001472:	ea8a 0606 	eor.w	r6, sl, r6
 8001476:	960b      	str	r6, [sp, #44]	@ 0x2c
 8001478:	ea89 0707 	eor.w	r7, r9, r7
 800147c:	ea28 0605 	bic.w	r6, r8, r5
 8001480:	407e      	eors	r6, r7
 8001482:	9f0b      	ldr	r7, [sp, #44]	@ 0x2c
 8001484:	960a      	str	r6, [sp, #40]	@ 0x28
 8001486:	ea24 0603 	bic.w	r6, r4, r3
 800148a:	4077      	eors	r7, r6
 800148c:	ea2e 0604 	bic.w	r6, lr, r4
 8001490:	405e      	eors	r6, r3
 8001492:	970b      	str	r7, [sp, #44]	@ 0x2c
 8001494:	ea21 0708 	bic.w	r7, r1, r8
 8001498:	406f      	eors	r7, r5
 800149a:	9623      	str	r6, [sp, #140]	@ 0x8c
 800149c:	ea25 0509 	bic.w	r5, r5, r9
 80014a0:	ea2c 060e 	bic.w	r6, ip, lr
 80014a4:	ea23 030a 	bic.w	r3, r3, sl
 80014a8:	4074      	eors	r4, r6
 80014aa:	ea83 030c 	eor.w	r3, r3, ip
 80014ae:	ea29 0602 	bic.w	r6, r9, r2
 80014b2:	9722      	str	r7, [sp, #136]	@ 0x88
 80014b4:	ea22 0701 	bic.w	r7, r2, r1
 80014b8:	406a      	eors	r2, r5
 80014ba:	4071      	eors	r1, r6
 80014bc:	943b      	str	r4, [sp, #236]	@ 0xec
 80014be:	922e      	str	r2, [sp, #184]	@ 0xb8
 80014c0:	ea2a 040c 	bic.w	r4, sl, ip
 80014c4:	9a10      	ldr	r2, [sp, #64]	@ 0x40
 80014c6:	932f      	str	r3, [sp, #188]	@ 0xbc
 80014c8:	9b06      	ldr	r3, [sp, #24]
 80014ca:	9116      	str	r1, [sp, #88]	@ 0x58
 80014cc:	ea87 0708 	eor.w	r7, r7, r8
 80014d0:	ea84 010e 	eor.w	r1, r4, lr
 80014d4:	973a      	str	r7, [sp, #232]	@ 0xe8
 80014d6:	9117      	str	r1, [sp, #92]	@ 0x5c
 80014d8:	ea83 0102 	eor.w	r1, r3, r2
 80014dc:	9b07      	ldr	r3, [sp, #28]
 80014de:	9a11      	ldr	r2, [sp, #68]	@ 0x44
 80014e0:	9c1d      	ldr	r4, [sp, #116]	@ 0x74
 80014e2:	9f2b      	ldr	r7, [sp, #172]	@ 0xac
 80014e4:	4053      	eors	r3, r2
 80014e6:	ea4f 7c03 	mov.w	ip, r3, lsl #28
 80014ea:	ea4c 1c11 	orr.w	ip, ip, r1, lsr #4
 80014ee:	0709      	lsls	r1, r1, #28
 80014f0:	ea41 1113 	orr.w	r1, r1, r3, lsr #4
 80014f4:	9a1c      	ldr	r2, [sp, #112]	@ 0x70
 80014f6:	9b08      	ldr	r3, [sp, #32]
 80014f8:	4053      	eors	r3, r2
 80014fa:	9a09      	ldr	r2, [sp, #36]	@ 0x24
 80014fc:	4062      	eors	r2, r4
 80014fe:	0516      	lsls	r6, r2, #20
 8001500:	ea46 3613 	orr.w	r6, r6, r3, lsr #12
 8001504:	051b      	lsls	r3, r3, #20
 8001506:	9c1e      	ldr	r4, [sp, #120]	@ 0x78
 8001508:	ea43 3312 	orr.w	r3, r3, r2, lsr #12
 800150c:	9a01      	ldr	r2, [sp, #4]
 800150e:	ea82 0504 	eor.w	r5, r2, r4
 8001512:	9c1f      	ldr	r4, [sp, #124]	@ 0x7c
 8001514:	9a02      	ldr	r2, [sp, #8]
 8001516:	4062      	eors	r2, r4
 8001518:	ea4f 08c2 	mov.w	r8, r2, lsl #3
 800151c:	ea48 7855 	orr.w	r8, r8, r5, lsr #29
 8001520:	00ed      	lsls	r5, r5, #3
 8001522:	9c2a      	ldr	r4, [sp, #168]	@ 0xa8
 8001524:	ea45 7552 	orr.w	r5, r5, r2, lsr #29
 8001528:	9a03      	ldr	r2, [sp, #12]
 800152a:	4062      	eors	r2, r4
 800152c:	9c04      	ldr	r4, [sp, #16]
 800152e:	ea4f 4ed2 	mov.w	lr, r2, lsr #19
 8001532:	407c      	eors	r4, r7
 8001534:	ea4e 3e44 	orr.w	lr, lr, r4, lsl #13
 8001538:	0ce4      	lsrs	r4, r4, #19
 800153a:	ea44 3442 	orr.w	r4, r4, r2, lsl #13
 800153e:	9a36      	ldr	r2, [sp, #216]	@ 0xd8
 8001540:	9f37      	ldr	r7, [sp, #220]	@ 0xdc
 8001542:	ea8b 0902 	eor.w	r9, fp, r2
 8001546:	9a05      	ldr	r2, [sp, #20]
 8001548:	407a      	eors	r2, r7
 800154a:	ea4f 07d9 	mov.w	r7, r9, lsr #3
 800154e:	ea47 7742 	orr.w	r7, r7, r2, lsl #29
 8001552:	ea25 0a03 	bic.w	sl, r5, r3
 8001556:	08d2      	lsrs	r2, r2, #3
 8001558:	ea42 7249 	orr.w	r2, r2, r9, lsl #29
 800155c:	ea8a 0a01 	eor.w	sl, sl, r1
 8001560:	ea28 0906 	bic.w	r9, r8, r6
 8001564:	ea89 090c 	eor.w	r9, r9, ip
 8001568:	f8cd a0a8 	str.w	sl, [sp, #168]	@ 0xa8
 800156c:	ea2e 0a05 	bic.w	sl, lr, r5
 8001570:	ea8a 0a03 	eor.w	sl, sl, r3
 8001574:	f8cd 90ac 	str.w	r9, [sp, #172]	@ 0xac
 8001578:	ea24 0908 	bic.w	r9, r4, r8
 800157c:	ea89 0906 	eor.w	r9, r9, r6
 8001580:	f8cd a040 	str.w	sl, [sp, #64]	@ 0x40
 8001584:	ea23 0301 	bic.w	r3, r3, r1
 8001588:	ea27 0a0e 	bic.w	sl, r7, lr
 800158c:	ea8a 0505 	eor.w	r5, sl, r5
 8001590:	ea26 060c 	bic.w	r6, r6, ip
 8001594:	407b      	eors	r3, r7
 8001596:	f8cd 9044 	str.w	r9, [sp, #68]	@ 0x44
 800159a:	ea22 0904 	bic.w	r9, r2, r4
 800159e:	951e      	str	r5, [sp, #120]	@ 0x78
 80015a0:	931c      	str	r3, [sp, #112]	@ 0x70
 80015a2:	ea89 0508 	eor.w	r5, r9, r8
 80015a6:	ea86 0302 	eor.w	r3, r6, r2
 80015aa:	951f      	str	r5, [sp, #124]	@ 0x7c
 80015ac:	931d      	str	r3, [sp, #116]	@ 0x74
 80015ae:	ea2c 0502 	bic.w	r5, ip, r2
 80015b2:	9b03      	ldr	r3, [sp, #12]
 80015b4:	9a0c      	ldr	r2, [sp, #48]	@ 0x30
 80015b6:	ea21 0807 	bic.w	r8, r1, r7
 80015ba:	ea83 0102 	eor.w	r1, r3, r2
 80015be:	9b04      	ldr	r3, [sp, #16]
 80015c0:	9a0d      	ldr	r2, [sp, #52]	@ 0x34
 80015c2:	4053      	eors	r3, r2
 80015c4:	ea4f 0c43 	mov.w	ip, r3, lsl #1
 80015c8:	ea4c 7cd1 	orr.w	ip, ip, r1, lsr #31
 80015cc:	ea88 0e0e 	eor.w	lr, r8, lr
 80015d0:	406c      	eors	r4, r5
 80015d2:	0049      	lsls	r1, r1, #1
 80015d4:	ea41 71d3 	orr.w	r1, r1, r3, lsr #31
 80015d8:	f8cd e0d8 	str.w	lr, [sp, #216]	@ 0xd8
 80015dc:	9b18      	ldr	r3, [sp, #96]	@ 0x60
 80015de:	9437      	str	r4, [sp, #220]	@ 0xdc
 80015e0:	9a05      	ldr	r2, [sp, #20]
 80015e2:	9c19      	ldr	r4, [sp, #100]	@ 0x64
 80015e4:	9f31      	ldr	r7, [sp, #196]	@ 0xc4
 80015e6:	f8dd 90cc 	ldr.w	r9, [sp, #204]	@ 0xcc
 80015ea:	4062      	eors	r2, r4
 80015ec:	ea8b 0303 	eor.w	r3, fp, r3
 80015f0:	0196      	lsls	r6, r2, #6
 80015f2:	ea46 6693 	orr.w	r6, r6, r3, lsr #26
 80015f6:	019b      	lsls	r3, r3, #6
 80015f8:	9c24      	ldr	r4, [sp, #144]	@ 0x90
 80015fa:	ea43 6392 	orr.w	r3, r3, r2, lsr #26
 80015fe:	9a06      	ldr	r2, [sp, #24]
 8001600:	ea82 0504 	eor.w	r5, r2, r4
 8001604:	9c25      	ldr	r4, [sp, #148]	@ 0x94
 8001606:	9a07      	ldr	r2, [sp, #28]
 8001608:	4062      	eors	r2, r4
 800160a:	ea4f 6842 	mov.w	r8, r2, lsl #25
 800160e:	ea48 18d5 	orr.w	r8, r8, r5, lsr #7
 8001612:	066d      	lsls	r5, r5, #25
 8001614:	ea45 15d2 	orr.w	r5, r5, r2, lsr #7
 8001618:	9c30      	ldr	r4, [sp, #192]	@ 0xc0
 800161a:	9a08      	ldr	r2, [sp, #32]
 800161c:	4054      	eors	r4, r2
 800161e:	9a09      	ldr	r2, [sp, #36]	@ 0x24
 8001620:	407a      	eors	r2, r7
 8001622:	ea4f 2e02 	mov.w	lr, r2, lsl #8
 8001626:	ea4e 6e14 	orr.w	lr, lr, r4, lsr #24
 800162a:	0224      	lsls	r4, r4, #8
 800162c:	9f32      	ldr	r7, [sp, #200]	@ 0xc8
 800162e:	ea44 6412 	orr.w	r4, r4, r2, lsr #24
 8001632:	9a01      	ldr	r2, [sp, #4]
 8001634:	407a      	eors	r2, r7
 8001636:	9f02      	ldr	r7, [sp, #8]
 8001638:	ea87 0909 	eor.w	r9, r7, r9
 800163c:	ea4f 4789 	mov.w	r7, r9, lsl #18
 8001640:	ea47 3792 	orr.w	r7, r7, r2, lsr #14
 8001644:	ea25 0a03 	bic.w	sl, r5, r3
 8001648:	0492      	lsls	r2, r2, #18
 800164a:	ea42 3299 	orr.w	r2, r2, r9, lsr #14
 800164e:	ea8a 0a01 	eor.w	sl, sl, r1
 8001652:	ea28 0906 	bic.w	r9, r8, r6
 8001656:	ea89 090c 	eor.w	r9, r9, ip
 800165a:	f8cd a060 	str.w	sl, [sp, #96]	@ 0x60
 800165e:	ea24 0a05 	bic.w	sl, r4, r5
 8001662:	ea8a 0a03 	eor.w	sl, sl, r3
 8001666:	f8cd 9064 	str.w	r9, [sp, #100]	@ 0x64
 800166a:	ea23 0301 	bic.w	r3, r3, r1
 800166e:	ea2e 0908 	bic.w	r9, lr, r8
 8001672:	ea89 0906 	eor.w	r9, r9, r6
 8001676:	4053      	eors	r3, r2
 8001678:	ea26 060c 	bic.w	r6, r6, ip
 800167c:	f8cd a0c0 	str.w	sl, [sp, #192]	@ 0xc0
 8001680:	ea22 0a04 	bic.w	sl, r2, r4
 8001684:	ea8a 0505 	eor.w	r5, sl, r5
 8001688:	f8cd 90c4 	str.w	r9, [sp, #196]	@ 0xc4
 800168c:	9332      	str	r3, [sp, #200]	@ 0xc8
 800168e:	ea27 090e 	bic.w	r9, r7, lr
 8001692:	ea86 0307 	eor.w	r3, r6, r7
 8001696:	950c      	str	r5, [sp, #48]	@ 0x30
 8001698:	9333      	str	r3, [sp, #204]	@ 0xcc
 800169a:	ea89 0508 	eor.w	r5, r9, r8
 800169e:	9b08      	ldr	r3, [sp, #32]
 80016a0:	950d      	str	r5, [sp, #52]	@ 0x34
 80016a2:	ea21 0802 	bic.w	r8, r1, r2
 80016a6:	9a12      	ldr	r2, [sp, #72]	@ 0x48
 80016a8:	ea83 0102 	eor.w	r1, r3, r2
 80016ac:	9a13      	ldr	r2, [sp, #76]	@ 0x4c
 80016ae:	9b09      	ldr	r3, [sp, #36]	@ 0x24
 80016b0:	4053      	eors	r3, r2
 80016b2:	ea2c 0507 	bic.w	r5, ip, r7
 80016b6:	ea4f 6cc3 	mov.w	ip, r3, lsl #27
 80016ba:	ea4c 1c51 	orr.w	ip, ip, r1, lsr #5
 80016be:	06c9      	lsls	r1, r1, #27
 80016c0:	ea41 1153 	orr.w	r1, r1, r3, lsr #5
 80016c4:	ea88 0404 	eor.w	r4, r8, r4
 80016c8:	9b01      	ldr	r3, [sp, #4]
 80016ca:	9a14      	ldr	r2, [sp, #80]	@ 0x50
 80016cc:	9424      	str	r4, [sp, #144]	@ 0x90
 80016ce:	ea85 040e 	eor.w	r4, r5, lr
 80016d2:	405a      	eors	r2, r3
 80016d4:	9425      	str	r4, [sp, #148]	@ 0x94
 80016d6:	9b02      	ldr	r3, [sp, #8]
 80016d8:	9c15      	ldr	r4, [sp, #84]	@ 0x54
 80016da:	0f16      	lsrs	r6, r2, #28
 80016dc:	4063      	eors	r3, r4
 80016de:	ea46 1603 	orr.w	r6, r6, r3, lsl #4
 80016e2:	0f1b      	lsrs	r3, r3, #28
 80016e4:	ea43 1302 	orr.w	r3, r3, r2, lsl #4
 80016e8:	9a03      	ldr	r2, [sp, #12]
 80016ea:	9c20      	ldr	r4, [sp, #128]	@ 0x80
 80016ec:	9f2d      	ldr	r7, [sp, #180]	@ 0xb4
 80016ee:	ea82 0504 	eor.w	r5, r2, r4
 80016f2:	9c21      	ldr	r4, [sp, #132]	@ 0x84
 80016f4:	9a04      	ldr	r2, [sp, #16]
 80016f6:	4062      	eors	r2, r4
 80016f8:	ea4f 2882 	mov.w	r8, r2, lsl #10
 80016fc:	ea48 5895 	orr.w	r8, r8, r5, lsr #22
 8001700:	02ad      	lsls	r5, r5, #10
 8001702:	ea45 5592 	orr.w	r5, r5, r2, lsr #22
 8001706:	9a2c      	ldr	r2, [sp, #176]	@ 0xb0
 8001708:	ea8b 0402 	eor.w	r4, fp, r2
 800170c:	9a05      	ldr	r2, [sp, #20]
 800170e:	407a      	eors	r2, r7
 8001710:	ea4f 3ec2 	mov.w	lr, r2, lsl #15
 8001714:	ea4e 4e54 	orr.w	lr, lr, r4, lsr #17
 8001718:	03e4      	lsls	r4, r4, #15
 800171a:	9f38      	ldr	r7, [sp, #224]	@ 0xe0
 800171c:	ea44 4452 	orr.w	r4, r4, r2, lsr #17
 8001720:	9a06      	ldr	r2, [sp, #24]
 8001722:	ea82 0907 	eor.w	r9, r2, r7
 8001726:	9f39      	ldr	r7, [sp, #228]	@ 0xe4
 8001728:	9a07      	ldr	r2, [sp, #28]
 800172a:	407a      	eors	r2, r7
 800172c:	ea4f 2719 	mov.w	r7, r9, lsr #8
 8001730:	ea47 6702 	orr.w	r7, r7, r2, lsl #24
 8001734:	ea25 0a06 	bic.w	sl, r5, r6
 8001738:	0a12      	lsrs	r2, r2, #8
 800173a:	ea42 6209 	orr.w	r2, r2, r9, lsl #24
 800173e:	ea8a 0a01 	eor.w	sl, sl, r1
 8001742:	ea28 0903 	bic.w	r9, r8, r3
 8001746:	ea89 090c 	eor.w	r9, r9, ip
 800174a:	f8cd a0e0 	str.w	sl, [sp, #224]	@ 0xe0
 800174e:	ea24 0a05 	bic.w	sl, r4, r5
 8001752:	ea8a 0a06 	eor.w	sl, sl, r6
 8001756:	f8cd 90e4 	str.w	r9, [sp, #228]	@ 0xe4
 800175a:	ea2e 0908 	bic.w	r9, lr, r8
 800175e:	ea89 0903 	eor.w	r9, r9, r3
 8001762:	f8cd a050 	str.w	sl, [sp, #80]	@ 0x50
 8001766:	ea23 030c 	bic.w	r3, r3, ip
 800176a:	ea27 0a04 	bic.w	sl, r7, r4
 800176e:	ea8a 0505 	eor.w	r5, sl, r5
 8001772:	4053      	eors	r3, r2
 8001774:	f8cd 9054 	str.w	r9, [sp, #84]	@ 0x54
 8001778:	ea26 0601 	bic.w	r6, r6, r1
 800177c:	ea22 090e 	bic.w	r9, r2, lr
 8001780:	952c      	str	r5, [sp, #176]	@ 0xb0
 8001782:	9321      	str	r3, [sp, #132]	@ 0x84
 8001784:	ea89 0508 	eor.w	r5, r9, r8
 8001788:	9b0e      	ldr	r3, [sp, #56]	@ 0x38
 800178a:	952d      	str	r5, [sp, #180]	@ 0xb4
 800178c:	ea21 0807 	bic.w	r8, r1, r7
 8001790:	ea86 0107 	eor.w	r1, r6, r7
 8001794:	ea2c 0502 	bic.w	r5, ip, r2
 8001798:	ea88 0404 	eor.w	r4, r8, r4
 800179c:	9a05      	ldr	r2, [sp, #20]
 800179e:	9120      	str	r1, [sp, #128]	@ 0x80
 80017a0:	990f      	ldr	r1, [sp, #60]	@ 0x3c
 80017a2:	9412      	str	r4, [sp, #72]	@ 0x48
 80017a4:	ea8b 0303 	eor.w	r3, fp, r3
 80017a8:	ea85 040e 	eor.w	r4, r5, lr
 80017ac:	9413      	str	r4, [sp, #76]	@ 0x4c
 80017ae:	089d      	lsrs	r5, r3, #2
 80017b0:	ea82 0401 	eor.w	r4, r2, r1
 80017b4:	ea45 7584 	orr.w	r5, r5, r4, lsl #30
 80017b8:	08a4      	lsrs	r4, r4, #2
 80017ba:	ea44 7483 	orr.w	r4, r4, r3, lsl #30
 80017be:	9a1a      	ldr	r2, [sp, #104]	@ 0x68
 80017c0:	9b06      	ldr	r3, [sp, #24]
 80017c2:	9e26      	ldr	r6, [sp, #152]	@ 0x98
 80017c4:	9f27      	ldr	r7, [sp, #156]	@ 0x9c
 80017c6:	ea83 0102 	eor.w	r1, r3, r2
 80017ca:	9b07      	ldr	r3, [sp, #28]
 80017cc:	9a1b      	ldr	r2, [sp, #108]	@ 0x6c
 80017ce:	405a      	eors	r2, r3
 80017d0:	0a4b      	lsrs	r3, r1, #9
 80017d2:	ea43 53c2 	orr.w	r3, r3, r2, lsl #23
 80017d6:	0a52      	lsrs	r2, r2, #9
 80017d8:	ea42 52c1 	orr.w	r2, r2, r1, lsl #23
 80017dc:	9908      	ldr	r1, [sp, #32]
 80017de:	4071      	eors	r1, r6
 80017e0:	9e09      	ldr	r6, [sp, #36]	@ 0x24
 80017e2:	ea4f 6c51 	mov.w	ip, r1, lsr #25
 80017e6:	4077      	eors	r7, r6
 80017e8:	ea4c 1cc7 	orr.w	ip, ip, r7, lsl #7
 80017ec:	0e7f      	lsrs	r7, r7, #25
 80017ee:	ea47 17c1 	orr.w	r7, r7, r1, lsl #7
 80017f2:	9901      	ldr	r1, [sp, #4]
 80017f4:	9e28      	ldr	r6, [sp, #160]	@ 0xa0
 80017f6:	f8dd e0a4 	ldr.w	lr, [sp, #164]	@ 0xa4
 80017fa:	f8dd 90d4 	ldr.w	r9, [sp, #212]	@ 0xd4
 80017fe:	4071      	eors	r1, r6
 8001800:	9e02      	ldr	r6, [sp, #8]
 8001802:	ea4f 58d1 	mov.w	r8, r1, lsr #23
 8001806:	ea86 060e 	eor.w	r6, r6, lr
 800180a:	ea48 2846 	orr.w	r8, r8, r6, lsl #9
 800180e:	0df6      	lsrs	r6, r6, #23
 8001810:	f8dd e0d0 	ldr.w	lr, [sp, #208]	@ 0xd0
 8001814:	ea46 2641 	orr.w	r6, r6, r1, lsl #9
 8001818:	9903      	ldr	r1, [sp, #12]
 800181a:	ea81 010e 	eor.w	r1, r1, lr
 800181e:	f8dd e010 	ldr.w	lr, [sp, #16]
 8001822:	ea8e 0909 	eor.w	r9, lr, r9
 8001826:	ea4f 0e89 	mov.w	lr, r9, lsl #2
 800182a:	ea4e 7e91 	orr.w	lr, lr, r1, lsr #30
 800182e:	ea2c 0a03 	bic.w	sl, ip, r3
 8001832:	0089      	lsls	r1, r1, #2
 8001834:	ea41 7199 	orr.w	r1, r1, r9, lsr #30
 8001838:	ea8a 0a05 	eor.w	sl, sl, r5
 800183c:	ea27 0902 	bic.w	r9, r7, r2
 8001840:	f8cd a098 	str.w	sl, [sp, #152]	@ 0x98
 8001844:	ea89 0904 	eor.w	r9, r9, r4
 8001848:	ea28 0a0c 	bic.w	sl, r8, ip
 800184c:	f8cd 909c 	str.w	r9, [sp, #156]	@ 0x9c
 8001850:	ea8a 0a03 	eor.w	sl, sl, r3
 8001854:	ea26 0907 	bic.w	r9, r6, r7
 8001858:	ea23 0305 	bic.w	r3, r3, r5
 800185c:	ea89 0902 	eor.w	r9, r9, r2
 8001860:	404b      	eors	r3, r1
 8001862:	ea22 0204 	bic.w	r2, r2, r4
 8001866:	930e      	str	r3, [sp, #56]	@ 0x38
 8001868:	ea82 030e 	eor.w	r3, r2, lr
 800186c:	f8cd a0d0 	str.w	sl, [sp, #208]	@ 0xd0
 8001870:	9a2a      	ldr	r2, [sp, #168]	@ 0xa8
 8001872:	f8cd 90d4 	str.w	r9, [sp, #212]	@ 0xd4
 8001876:	ea21 0a08 	bic.w	sl, r1, r8
 800187a:	ea2e 0906 	bic.w	r9, lr, r6
 800187e:	930f      	str	r3, [sp, #60]	@ 0x3c
 8001880:	9b0a      	ldr	r3, [sp, #40]	@ 0x28
 8001882:	ea8a 0c0c 	eor.w	ip, sl, ip
 8001886:	ea89 0707 	eor.w	r7, r9, r7
 800188a:	f8cd c068 	str.w	ip, [sp, #104]	@ 0x68
 800188e:	971b      	str	r7, [sp, #108]	@ 0x6c
 8001890:	ea25 0c01 	bic.w	ip, r5, r1
 8001894:	ea24 070e 	bic.w	r7, r4, lr
 8001898:	992b      	ldr	r1, [sp, #172]	@ 0xac
 800189a:	405a      	eors	r2, r3
 800189c:	9b0b      	ldr	r3, [sp, #44]	@ 0x2c
 800189e:	407e      	eors	r6, r7
 80018a0:	9629      	str	r6, [sp, #164]	@ 0xa4
 80018a2:	ea83 0601 	eor.w	r6, r3, r1
 80018a6:	9918      	ldr	r1, [sp, #96]	@ 0x60
 80018a8:	9b22      	ldr	r3, [sp, #136]	@ 0x88
 80018aa:	404a      	eors	r2, r1
 80018ac:	9919      	ldr	r1, [sp, #100]	@ 0x64
 80018ae:	404e      	eors	r6, r1
 80018b0:	9938      	ldr	r1, [sp, #224]	@ 0xe0
 80018b2:	404a      	eors	r2, r1
 80018b4:	9939      	ldr	r1, [sp, #228]	@ 0xe4
 80018b6:	404e      	eors	r6, r1
 80018b8:	9926      	ldr	r1, [sp, #152]	@ 0x98
 80018ba:	404a      	eors	r2, r1
 80018bc:	9927      	ldr	r1, [sp, #156]	@ 0x9c
 80018be:	404e      	eors	r6, r1
 80018c0:	9910      	ldr	r1, [sp, #64]	@ 0x40
 80018c2:	ea83 0901 	eor.w	r9, r3, r1
 80018c6:	9911      	ldr	r1, [sp, #68]	@ 0x44
 80018c8:	9b23      	ldr	r3, [sp, #140]	@ 0x8c
 80018ca:	ea8c 0c08 	eor.w	ip, ip, r8
 80018ce:	f8cd c0a0 	str.w	ip, [sp, #160]	@ 0xa0
 80018d2:	ea83 0801 	eor.w	r8, r3, r1
 80018d6:	9b30      	ldr	r3, [sp, #192]	@ 0xc0
 80018d8:	991e      	ldr	r1, [sp, #120]	@ 0x78
 80018da:	f8dd e0cc 	ldr.w	lr, [sp, #204]	@ 0xcc
 80018de:	ea89 0903 	eor.w	r9, r9, r3
 80018e2:	9b31      	ldr	r3, [sp, #196]	@ 0xc4
 80018e4:	ea88 0803 	eor.w	r8, r8, r3
 80018e8:	9b14      	ldr	r3, [sp, #80]	@ 0x50
 80018ea:	ea89 0903 	eor.w	r9, r9, r3
 80018ee:	9b15      	ldr	r3, [sp, #84]	@ 0x54
 80018f0:	ea88 0803 	eor.w	r8, r8, r3
 80018f4:	9b34      	ldr	r3, [sp, #208]	@ 0xd0
 80018f6:	ea89 0903 	eor.w	r9, r9, r3
 80018fa:	9b35      	ldr	r3, [sp, #212]	@ 0xd4
 80018fc:	ea88 0803 	eor.w	r8, r8, r3
 8001900:	9b3a      	ldr	r3, [sp, #232]	@ 0xe8
 8001902:	ea83 0c01 	eor.w	ip, r3, r1
 8001906:	9b3b      	ldr	r3, [sp, #236]	@ 0xec
 8001908:	991f      	ldr	r1, [sp, #124]	@ 0x7c
 800190a:	ea83 0701 	eor.w	r7, r3, r1
 800190e:	9b0c      	ldr	r3, [sp, #48]	@ 0x30
 8001910:	990d      	ldr	r1, [sp, #52]	@ 0x34
 8001912:	ea8c 0c03 	eor.w	ip, ip, r3
 8001916:	404f      	eors	r7, r1
 8001918:	9b2c      	ldr	r3, [sp, #176]	@ 0xb0
 800191a:	992d      	ldr	r1, [sp, #180]	@ 0xb4
 800191c:	ea8c 0c03 	eor.w	ip, ip, r3
 8001920:	404f      	eors	r7, r1
 8001922:	9b1a      	ldr	r3, [sp, #104]	@ 0x68
 8001924:	991b      	ldr	r1, [sp, #108]	@ 0x6c
 8001926:	ea8c 0c03 	eor.w	ip, ip, r3
 800192a:	404f      	eors	r7, r1
 800192c:	9b16      	ldr	r3, [sp, #88]	@ 0x58
 800192e:	9936      	ldr	r1, [sp, #216]	@ 0xd8
 8001930:	ea83 0501 	eor.w	r5, r3, r1
 8001934:	9b17      	ldr	r3, [sp, #92]	@ 0x5c
 8001936:	9937      	ldr	r1, [sp, #220]	@ 0xdc
 8001938:	ea83 0401 	eor.w	r4, r3, r1
 800193c:	9924      	ldr	r1, [sp, #144]	@ 0x90
 800193e:	9b1c      	ldr	r3, [sp, #112]	@ 0x70
 8001940:	404d      	eors	r5, r1
 8001942:	9925      	ldr	r1, [sp, #148]	@ 0x94
 8001944:	404c      	eors	r4, r1
 8001946:	9912      	ldr	r1, [sp, #72]	@ 0x48
 8001948:	404d      	eors	r5, r1
 800194a:	9913      	ldr	r1, [sp, #76]	@ 0x4c
 800194c:	404c      	eors	r4, r1
 800194e:	9928      	ldr	r1, [sp, #160]	@ 0xa0
 8001950:	404d      	eors	r5, r1
 8001952:	9929      	ldr	r1, [sp, #164]	@ 0xa4
 8001954:	404c      	eors	r4, r1
 8001956:	9932      	ldr	r1, [sp, #200]	@ 0xc8
 8001958:	4059      	eors	r1, r3
 800195a:	9b1d      	ldr	r3, [sp, #116]	@ 0x74
 800195c:	ea83 0e0e 	eor.w	lr, r3, lr
 8001960:	9b2e      	ldr	r3, [sp, #184]	@ 0xb8
 8001962:	4059      	eors	r1, r3
 8001964:	9b2f      	ldr	r3, [sp, #188]	@ 0xbc
 8001966:	ea8e 0e03 	eor.w	lr, lr, r3
 800196a:	9b20      	ldr	r3, [sp, #128]	@ 0x80
 800196c:	4059      	eors	r1, r3
 800196e:	9b21      	ldr	r3, [sp, #132]	@ 0x84
 8001970:	ea8e 0e03 	eor.w	lr, lr, r3
 8001974:	9b0e      	ldr	r3, [sp, #56]	@ 0x38
 8001976:	ea4f 0b49 	mov.w	fp, r9, lsl #1
 800197a:	4059      	eors	r1, r3
 800197c:	9b0f      	ldr	r3, [sp, #60]	@ 0x3c
 800197e:	ea4b 7bd8 	orr.w	fp, fp, r8, lsr #31
 8001982:	ea4f 0a48 	mov.w	sl, r8, lsl #1
 8001986:	ea8e 0e03 	eor.w	lr, lr, r3
 800198a:	ea4a 7ad9 	orr.w	sl, sl, r9, lsr #31
 800198e:	ea8b 0301 	eor.w	r3, fp, r1
 8001992:	ea4f 0b4c 	mov.w	fp, ip, lsl #1
 8001996:	9301      	str	r3, [sp, #4]
 8001998:	ea4b 7bd7 	orr.w	fp, fp, r7, lsr #31
 800199c:	ea8a 030e 	eor.w	r3, sl, lr
 80019a0:	ea4f 0a47 	mov.w	sl, r7, lsl #1
 80019a4:	9302      	str	r3, [sp, #8]
 80019a6:	ea4a 7adc 	orr.w	sl, sl, ip, lsr #31
 80019aa:	ea8b 0302 	eor.w	r3, fp, r2
 80019ae:	9303      	str	r3, [sp, #12]
 80019b0:	ea8a 0306 	eor.w	r3, sl, r6
 80019b4:	ea4f 0a44 	mov.w	sl, r4, lsl #1
 80019b8:	ea4a 7ad5 	orr.w	sl, sl, r5, lsr #31
 80019bc:	9304      	str	r3, [sp, #16]
 80019be:	ea8a 0308 	eor.w	r3, sl, r8
 80019c2:	ea4f 084e 	mov.w	r8, lr, lsl #1
 80019c6:	ea48 78d1 	orr.w	r8, r8, r1, lsr #31
 80019ca:	0049      	lsls	r1, r1, #1
 80019cc:	ea41 71de 	orr.w	r1, r1, lr, lsr #31
 80019d0:	9305      	str	r3, [sp, #20]
 80019d2:	ea81 030c 	eor.w	r3, r1, ip
 80019d6:	0071      	lsls	r1, r6, #1
 80019d8:	ea41 71d2 	orr.w	r1, r1, r2, lsr #31
 80019dc:	0052      	lsls	r2, r2, #1
 80019de:	9306      	str	r3, [sp, #24]
 80019e0:	ea42 72d6 	orr.w	r2, r2, r6, lsr #31
 80019e4:	ea88 0307 	eor.w	r3, r8, r7
 80019e8:	9307      	str	r3, [sp, #28]
 80019ea:	ea82 0305 	eor.w	r3, r2, r5
 80019ee:	9308      	str	r3, [sp, #32]
 80019f0:	ea81 0304 	eor.w	r3, r1, r4
 80019f4:	9a01      	ldr	r2, [sp, #4]
 80019f6:	9309      	str	r3, [sp, #36]	@ 0x24
 80019f8:	ea4f 0b45 	mov.w	fp, r5, lsl #1
 80019fc:	9b0a      	ldr	r3, [sp, #40]	@ 0x28
 80019fe:	9e07      	ldr	r6, [sp, #28]
 8001a00:	9f09      	ldr	r7, [sp, #36]	@ 0x24
 8001a02:	ea4b 7bd4 	orr.w	fp, fp, r4, lsr #31
 8001a06:	ea8b 0b09 	eor.w	fp, fp, r9
 8001a0a:	ea83 0902 	eor.w	r9, r3, r2
 8001a0e:	9b0b      	ldr	r3, [sp, #44]	@ 0x2c
 8001a10:	9a02      	ldr	r2, [sp, #8]
 8001a12:	9c05      	ldr	r4, [sp, #20]
 8001a14:	ea83 0a02 	eor.w	sl, r3, r2
 8001a18:	9b10      	ldr	r3, [sp, #64]	@ 0x40
 8001a1a:	9a03      	ldr	r2, [sp, #12]
 8001a1c:	ea83 0102 	eor.w	r1, r3, r2
 8001a20:	9b11      	ldr	r3, [sp, #68]	@ 0x44
 8001a22:	9a04      	ldr	r2, [sp, #16]
 8001a24:	0d0d      	lsrs	r5, r1, #20
 8001a26:	405a      	eors	r2, r3
 8001a28:	9b0c      	ldr	r3, [sp, #48]	@ 0x30
 8001a2a:	ea45 3502 	orr.w	r5, r5, r2, lsl #12
 8001a2e:	0d12      	lsrs	r2, r2, #20
 8001a30:	ea42 3201 	orr.w	r2, r2, r1, lsl #12
 8001a34:	ea83 010b 	eor.w	r1, r3, fp
 8001a38:	9b0d      	ldr	r3, [sp, #52]	@ 0x34
 8001a3a:	ea4f 5851 	mov.w	r8, r1, lsr #21
 8001a3e:	405c      	eors	r4, r3
 8001a40:	ea48 28c4 	orr.w	r8, r8, r4, lsl #11
 8001a44:	0d64      	lsrs	r4, r4, #21
 8001a46:	ea44 24c1 	orr.w	r4, r4, r1, lsl #11
 8001a4a:	9b12      	ldr	r3, [sp, #72]	@ 0x48
 8001a4c:	9906      	ldr	r1, [sp, #24]
 8001a4e:	4059      	eors	r1, r3
 8001a50:	9b13      	ldr	r3, [sp, #76]	@ 0x4c
 8001a52:	405e      	eors	r6, r3
 8001a54:	ea4f 5e46 	mov.w	lr, r6, lsl #21
 8001a58:	ea4e 2ed1 	orr.w	lr, lr, r1, lsr #11
 8001a5c:	0549      	lsls	r1, r1, #21
 8001a5e:	ea41 21d6 	orr.w	r1, r1, r6, lsr #11
 8001a62:	9b08      	ldr	r3, [sp, #32]
 8001a64:	9e0e      	ldr	r6, [sp, #56]	@ 0x38
 8001a66:	405e      	eors	r6, r3
 8001a68:	4633      	mov	r3, r6
 8001a6a:	9e0f      	ldr	r6, [sp, #60]	@ 0x3c
 8001a6c:	407e      	eors	r6, r7
 8001a6e:	ea4f 3c86 	mov.w	ip, r6, lsl #14
 8001a72:	ea4c 4c93 	orr.w	ip, ip, r3, lsr #18
 8001a76:	039b      	lsls	r3, r3, #14
 8001a78:	ea43 4396 	orr.w	r3, r3, r6, lsr #18
 8001a7c:	9e00      	ldr	r6, [sp, #0]
 8001a7e:	e9d6 7602 	ldrd	r7, r6, [r6, #8]
 8001a82:	ea8a 0606 	eor.w	r6, sl, r6
 8001a86:	960b      	str	r6, [sp, #44]	@ 0x2c
 8001a88:	ea89 0707 	eor.w	r7, r9, r7
 8001a8c:	ea28 0605 	bic.w	r6, r8, r5
 8001a90:	407e      	eors	r6, r7
 8001a92:	9f0b      	ldr	r7, [sp, #44]	@ 0x2c
 8001a94:	960a      	str	r6, [sp, #40]	@ 0x28
 8001a96:	ea24 0602 	bic.w	r6, r4, r2
 8001a9a:	4077      	eors	r7, r6
 8001a9c:	ea2e 0604 	bic.w	r6, lr, r4
 8001aa0:	4056      	eors	r6, r2
 8001aa2:	970b      	str	r7, [sp, #44]	@ 0x2c
 8001aa4:	ea21 0708 	bic.w	r7, r1, r8
 8001aa8:	406f      	eors	r7, r5
 8001aaa:	960d      	str	r6, [sp, #52]	@ 0x34
 8001aac:	ea25 0509 	bic.w	r5, r5, r9
 8001ab0:	ea2c 060e 	bic.w	r6, ip, lr
 8001ab4:	4074      	eors	r4, r6
 8001ab6:	ea22 020a 	bic.w	r2, r2, sl
 8001aba:	ea29 0603 	bic.w	r6, r9, r3
 8001abe:	970c      	str	r7, [sp, #48]	@ 0x30
 8001ac0:	ea23 0701 	bic.w	r7, r3, r1
 8001ac4:	406b      	eors	r3, r5
 8001ac6:	9312      	str	r3, [sp, #72]	@ 0x48
 8001ac8:	ea82 030c 	eor.w	r3, r2, ip
 8001acc:	4071      	eors	r1, r6
 8001ace:	9a06      	ldr	r2, [sp, #24]
 8001ad0:	940f      	str	r4, [sp, #60]	@ 0x3c
 8001ad2:	9313      	str	r3, [sp, #76]	@ 0x4c
 8001ad4:	ea2a 040c 	bic.w	r4, sl, ip
 8001ad8:	9b16      	ldr	r3, [sp, #88]	@ 0x58
 8001ada:	9110      	str	r1, [sp, #64]	@ 0x40
 8001adc:	ea84 010e 	eor.w	r1, r4, lr
 8001ae0:	9111      	str	r1, [sp, #68]	@ 0x44
 8001ae2:	ea83 0102 	eor.w	r1, r3, r2
 8001ae6:	9a07      	ldr	r2, [sp, #28]
 8001ae8:	9b17      	ldr	r3, [sp, #92]	@ 0x5c
 8001aea:	9c09      	ldr	r4, [sp, #36]	@ 0x24
 8001aec:	4053      	eors	r3, r2
 8001aee:	ea4f 7c03 	mov.w	ip, r3, lsl #28
 8001af2:	ea4c 1c11 	orr.w	ip, ip, r1, lsr #4
 8001af6:	0709      	lsls	r1, r1, #28
 8001af8:	ea41 1113 	orr.w	r1, r1, r3, lsr #4
 8001afc:	9a08      	ldr	r2, [sp, #32]
 8001afe:	9b1c      	ldr	r3, [sp, #112]	@ 0x70
 8001b00:	4053      	eors	r3, r2
 8001b02:	9a1d      	ldr	r2, [sp, #116]	@ 0x74
 8001b04:	4062      	eors	r2, r4
 8001b06:	0516      	lsls	r6, r2, #20
 8001b08:	ea46 3613 	orr.w	r6, r6, r3, lsr #12
 8001b0c:	051b      	lsls	r3, r3, #20
 8001b0e:	ea43 3312 	orr.w	r3, r3, r2, lsr #12
 8001b12:	9c01      	ldr	r4, [sp, #4]
 8001b14:	9a18      	ldr	r2, [sp, #96]	@ 0x60
 8001b16:	ea82 0504 	eor.w	r5, r2, r4
 8001b1a:	9c02      	ldr	r4, [sp, #8]
 8001b1c:	9a19      	ldr	r2, [sp, #100]	@ 0x64
 8001b1e:	4062      	eors	r2, r4
 8001b20:	ea87 0708 	eor.w	r7, r7, r8
 8001b24:	ea4f 08c2 	mov.w	r8, r2, lsl #3
 8001b28:	ea48 7855 	orr.w	r8, r8, r5, lsr #29
 8001b2c:	00ed      	lsls	r5, r5, #3
 8001b2e:	ea45 7552 	orr.w	r5, r5, r2, lsr #29
 8001b32:	9c03      	ldr	r4, [sp, #12]
 8001b34:	9a14      	ldr	r2, [sp, #80]	@ 0x50
 8001b36:	970e      	str	r7, [sp, #56]	@ 0x38
 8001b38:	4062      	eors	r2, r4
 8001b3a:	9f04      	ldr	r7, [sp, #16]
 8001b3c:	9c15      	ldr	r4, [sp, #84]	@ 0x54
 8001b3e:	ea4f 4ed2 	mov.w	lr, r2, lsr #19
 8001b42:	407c      	eors	r4, r7
 8001b44:	ea4e 3e44 	orr.w	lr, lr, r4, lsl #13
 8001b48:	0ce4      	lsrs	r4, r4, #19
 8001b4a:	ea44 3442 	orr.w	r4, r4, r2, lsl #13
 8001b4e:	9a1a      	ldr	r2, [sp, #104]	@ 0x68
 8001b50:	ea82 090b 	eor.w	r9, r2, fp
 8001b54:	9a1b      	ldr	r2, [sp, #108]	@ 0x6c
 8001b56:	9f05      	ldr	r7, [sp, #20]
 8001b58:	407a      	eors	r2, r7
 8001b5a:	ea4f 07d9 	mov.w	r7, r9, lsr #3
 8001b5e:	ea47 7742 	orr.w	r7, r7, r2, lsl #29
 8001b62:	ea25 0a03 	bic.w	sl, r5, r3
 8001b66:	08d2      	lsrs	r2, r2, #3
 8001b68:	ea42 7249 	orr.w	r2, r2, r9, lsl #29
 8001b6c:	ea8a 0a01 	eor.w	sl, sl, r1
 8001b70:	ea28 0906 	bic.w	r9, r8, r6
 8001b74:	ea89 090c 	eor.w	r9, r9, ip
 8001b78:	f8cd a050 	str.w	sl, [sp, #80]	@ 0x50
 8001b7c:	ea2e 0a05 	bic.w	sl, lr, r5
 8001b80:	ea8a 0a03 	eor.w	sl, sl, r3
 8001b84:	f8cd 9054 	str.w	r9, [sp, #84]	@ 0x54
 8001b88:	ea24 0908 	bic.w	r9, r4, r8
 8001b8c:	ea89 0906 	eor.w	r9, r9, r6
 8001b90:	f8cd a058 	str.w	sl, [sp, #88]	@ 0x58
 8001b94:	ea23 0301 	bic.w	r3, r3, r1
 8001b98:	ea27 0a0e 	bic.w	sl, r7, lr
 8001b9c:	ea8a 0505 	eor.w	r5, sl, r5
 8001ba0:	ea26 060c 	bic.w	r6, r6, ip
 8001ba4:	407b      	eors	r3, r7
 8001ba6:	f8cd 905c 	str.w	r9, [sp, #92]	@ 0x5c
 8001baa:	ea22 0904 	bic.w	r9, r2, r4
 8001bae:	9518      	str	r5, [sp, #96]	@ 0x60
 8001bb0:	931c      	str	r3, [sp, #112]	@ 0x70
 8001bb2:	ea89 0508 	eor.w	r5, r9, r8
 8001bb6:	ea86 0302 	eor.w	r3, r6, r2
 8001bba:	9519      	str	r5, [sp, #100]	@ 0x64
 8001bbc:	931d      	str	r3, [sp, #116]	@ 0x74
 8001bbe:	ea2c 0502 	bic.w	r5, ip, r2
 8001bc2:	9b22      	ldr	r3, [sp, #136]	@ 0x88
 8001bc4:	9a03      	ldr	r2, [sp, #12]
 8001bc6:	f8dd 9008 	ldr.w	r9, [sp, #8]
 8001bca:	ea21 0807 	bic.w	r8, r1, r7
 8001bce:	ea83 0102 	eor.w	r1, r3, r2
 8001bd2:	9a04      	ldr	r2, [sp, #16]
 8001bd4:	9b23      	ldr	r3, [sp, #140]	@ 0x8c
 8001bd6:	9f09      	ldr	r7, [sp, #36]	@ 0x24
 8001bd8:	4053      	eors	r3, r2
 8001bda:	406c      	eors	r4, r5
 8001bdc:	ea4f 0c43 	mov.w	ip, r3, lsl #1
 8001be0:	ea4c 7cd1 	orr.w	ip, ip, r1, lsr #31
 8001be4:	9a1f      	ldr	r2, [sp, #124]	@ 0x7c
 8001be6:	941b      	str	r4, [sp, #108]	@ 0x6c
 8001be8:	0049      	lsls	r1, r1, #1
 8001bea:	9c05      	ldr	r4, [sp, #20]
 8001bec:	ea41 71d3 	orr.w	r1, r1, r3, lsr #31
 8001bf0:	9b1e      	ldr	r3, [sp, #120]	@ 0x78
 8001bf2:	4062      	eors	r2, r4
 8001bf4:	ea83 030b 	eor.w	r3, r3, fp
 8001bf8:	0196      	lsls	r6, r2, #6
 8001bfa:	ea46 6693 	orr.w	r6, r6, r3, lsr #26
 8001bfe:	019b      	lsls	r3, r3, #6
 8001c00:	ea43 6392 	orr.w	r3, r3, r2, lsr #26
 8001c04:	9c06      	ldr	r4, [sp, #24]
 8001c06:	9a24      	ldr	r2, [sp, #144]	@ 0x90
 8001c08:	ea82 0504 	eor.w	r5, r2, r4
 8001c0c:	9c07      	ldr	r4, [sp, #28]
 8001c0e:	9a25      	ldr	r2, [sp, #148]	@ 0x94
 8001c10:	4062      	eors	r2, r4
 8001c12:	ea88 0e0e 	eor.w	lr, r8, lr
 8001c16:	ea4f 6842 	mov.w	r8, r2, lsl #25
 8001c1a:	ea48 18d5 	orr.w	r8, r8, r5, lsr #7
 8001c1e:	066d      	lsls	r5, r5, #25
 8001c20:	ea45 15d2 	orr.w	r5, r5, r2, lsr #7
 8001c24:	9c08      	ldr	r4, [sp, #32]
 8001c26:	9a20      	ldr	r2, [sp, #128]	@ 0x80
 8001c28:	f8cd e068 	str.w	lr, [sp, #104]	@ 0x68
 8001c2c:	4054      	eors	r4, r2
 8001c2e:	9a21      	ldr	r2, [sp, #132]	@ 0x84
 8001c30:	407a      	eors	r2, r7
 8001c32:	ea4f 2e02 	mov.w	lr, r2, lsl #8
 8001c36:	ea4e 6e14 	orr.w	lr, lr, r4, lsr #24
 8001c3a:	0224      	lsls	r4, r4, #8
 8001c3c:	ea44 6412 	orr.w	r4, r4, r2, lsr #24
 8001c40:	9f01      	ldr	r7, [sp, #4]
 8001c42:	9a26      	ldr	r2, [sp, #152]	@ 0x98
 8001c44:	407a      	eors	r2, r7
 8001c46:	9f27      	ldr	r7, [sp, #156]	@ 0x9c
 8001c48:	ea87 0909 	eor.w	r9, r7, r9
 8001c4c:	ea4f 4789 	mov.w	r7, r9, lsl #18
 8001c50:	ea47 3792 	orr.w	r7, r7, r2, lsr #14
 8001c54:	ea25 0a03 	bic.w	sl, r5, r3
 8001c58:	0492      	lsls	r2, r2, #18
 8001c5a:	ea42 3299 	orr.w	r2, r2, r9, lsr #14
 8001c5e:	ea8a 0a01 	eor.w	sl, sl, r1
 8001c62:	ea28 0906 	bic.w	r9, r8, r6
 8001c66:	ea89 090c 	eor.w	r9, r9, ip
 8001c6a:	f8cd a078 	str.w	sl, [sp, #120]	@ 0x78
 8001c6e:	ea24 0a05 	bic.w	sl, r4, r5
 8001c72:	f8cd 907c 	str.w	r9, [sp, #124]	@ 0x7c
 8001c76:	ea8a 0a03 	eor.w	sl, sl, r3
 8001c7a:	ea2e 0908 	bic.w	r9, lr, r8
 8001c7e:	ea23 0301 	bic.w	r3, r3, r1
 8001c82:	ea89 0906 	eor.w	r9, r9, r6
 8001c86:	4053      	eors	r3, r2
 8001c88:	ea26 060c 	bic.w	r6, r6, ip
 8001c8c:	f8cd a080 	str.w	sl, [sp, #128]	@ 0x80
 8001c90:	ea22 0a04 	bic.w	sl, r2, r4
 8001c94:	ea8a 0505 	eor.w	r5, sl, r5
 8001c98:	f8cd 9084 	str.w	r9, [sp, #132]	@ 0x84
 8001c9c:	9326      	str	r3, [sp, #152]	@ 0x98
 8001c9e:	ea27 090e 	bic.w	r9, r7, lr
 8001ca2:	ea86 0307 	eor.w	r3, r6, r7
 8001ca6:	9522      	str	r5, [sp, #136]	@ 0x88
 8001ca8:	9327      	str	r3, [sp, #156]	@ 0x9c
 8001caa:	ea89 0508 	eor.w	r5, r9, r8
 8001cae:	9b2e      	ldr	r3, [sp, #184]	@ 0xb8
 8001cb0:	9523      	str	r5, [sp, #140]	@ 0x8c
 8001cb2:	ea21 0802 	bic.w	r8, r1, r2
 8001cb6:	9a08      	ldr	r2, [sp, #32]
 8001cb8:	ea83 0102 	eor.w	r1, r3, r2
 8001cbc:	9a09      	ldr	r2, [sp, #36]	@ 0x24
 8001cbe:	9b2f      	ldr	r3, [sp, #188]	@ 0xbc
 8001cc0:	4053      	eors	r3, r2
 8001cc2:	ea2c 0507 	bic.w	r5, ip, r7
 8001cc6:	ea4f 6cc3 	mov.w	ip, r3, lsl #27
 8001cca:	ea4c 1c51 	orr.w	ip, ip, r1, lsr #5
 8001cce:	06c9      	lsls	r1, r1, #27
 8001cd0:	ea41 1153 	orr.w	r1, r1, r3, lsr #5
 8001cd4:	ea88 0404 	eor.w	r4, r8, r4
 8001cd8:	9b2a      	ldr	r3, [sp, #168]	@ 0xa8
 8001cda:	9a01      	ldr	r2, [sp, #4]
 8001cdc:	9424      	str	r4, [sp, #144]	@ 0x90
 8001cde:	ea85 040e 	eor.w	r4, r5, lr
 8001ce2:	405a      	eors	r2, r3
 8001ce4:	9425      	str	r4, [sp, #148]	@ 0x94
 8001ce6:	9b2b      	ldr	r3, [sp, #172]	@ 0xac
 8001ce8:	9c02      	ldr	r4, [sp, #8]
 8001cea:	9f05      	ldr	r7, [sp, #20]
 8001cec:	4063      	eors	r3, r4
 8001cee:	0f16      	lsrs	r6, r2, #28
 8001cf0:	ea46 1603 	orr.w	r6, r6, r3, lsl #4
 8001cf4:	0f1b      	lsrs	r3, r3, #28
 8001cf6:	ea43 1302 	orr.w	r3, r3, r2, lsl #4
 8001cfa:	9c03      	ldr	r4, [sp, #12]
 8001cfc:	9a30      	ldr	r2, [sp, #192]	@ 0xc0
 8001cfe:	ea82 0504 	eor.w	r5, r2, r4
 8001d02:	9c04      	ldr	r4, [sp, #16]
 8001d04:	9a31      	ldr	r2, [sp, #196]	@ 0xc4
 8001d06:	4062      	eors	r2, r4
 8001d08:	ea4f 2882 	mov.w	r8, r2, lsl #10
 8001d0c:	ea48 5895 	orr.w	r8, r8, r5, lsr #22
 8001d10:	02ad      	lsls	r5, r5, #10
 8001d12:	ea45 5592 	orr.w	r5, r5, r2, lsr #22
 8001d16:	9a2c      	ldr	r2, [sp, #176]	@ 0xb0
 8001d18:	ea82 040b 	eor.w	r4, r2, fp
 8001d1c:	9a2d      	ldr	r2, [sp, #180]	@ 0xb4
 8001d1e:	407a      	eors	r2, r7
 8001d20:	ea4f 3ec2 	mov.w	lr, r2, lsl #15
 8001d24:	ea4e 4e54 	orr.w	lr, lr, r4, lsr #17
 8001d28:	03e4      	lsls	r4, r4, #15
 8001d2a:	ea44 4452 	orr.w	r4, r4, r2, lsr #17
 8001d2e:	9f06      	ldr	r7, [sp, #24]
 8001d30:	9a28      	ldr	r2, [sp, #160]	@ 0xa0
 8001d32:	ea82 0907 	eor.w	r9, r2, r7
 8001d36:	9a29      	ldr	r2, [sp, #164]	@ 0xa4
 8001d38:	9f07      	ldr	r7, [sp, #28]
 8001d3a:	407a      	eors	r2, r7
 8001d3c:	ea4f 2719 	mov.w	r7, r9, lsr #8
 8001d40:	ea47 6702 	orr.w	r7, r7, r2, lsl #24
 8001d44:	ea25 0a06 	bic.w	sl, r5, r6
 8001d48:	0a12      	lsrs	r2, r2, #8
 8001d4a:	ea42 6209 	orr.w	r2, r2, r9, lsl #24
 8001d4e:	ea8a 0a01 	eor.w	sl, sl, r1
 8001d52:	ea28 0903 	bic.w	r9, r8, r3
 8001d56:	ea89 090c 	eor.w	r9, r9, ip
 8001d5a:	f8cd a0a0 	str.w	sl, [sp, #160]	@ 0xa0
 8001d5e:	ea24 0a05 	bic.w	sl, r4, r5
 8001d62:	ea8a 0a06 	eor.w	sl, sl, r6
 8001d66:	f8cd 90a4 	str.w	r9, [sp, #164]	@ 0xa4
 8001d6a:	ea2e 0908 	bic.w	r9, lr, r8
 8001d6e:	ea89 0903 	eor.w	r9, r9, r3
 8001d72:	f8cd a0a8 	str.w	sl, [sp, #168]	@ 0xa8
 8001d76:	ea27 0a04 	bic.w	sl, r7, r4
 8001d7a:	ea8a 0505 	eor.w	r5, sl, r5
 8001d7e:	f8cd 90ac 	str.w	r9, [sp, #172]	@ 0xac
 8001d82:	ea23 030c 	bic.w	r3, r3, ip
 8001d86:	ea22 090e 	bic.w	r9, r2, lr
 8001d8a:	4053      	eors	r3, r2
 8001d8c:	952c      	str	r5, [sp, #176]	@ 0xb0
 8001d8e:	ea26 0601 	bic.w	r6, r6, r1
 8001d92:	ea89 0508 	eor.w	r5, r9, r8
 8001d96:	952d      	str	r5, [sp, #180]	@ 0xb4
 8001d98:	ea21 0807 	bic.w	r8, r1, r7
 8001d9c:	9331      	str	r3, [sp, #196]	@ 0xc4
 8001d9e:	ea86 0107 	eor.w	r1, r6, r7
 8001da2:	9b3a      	ldr	r3, [sp, #232]	@ 0xe8
 8001da4:	9130      	str	r1, [sp, #192]	@ 0xc0
 8001da6:	ea2c 0502 	bic.w	r5, ip, r2
 8001daa:	9905      	ldr	r1, [sp, #20]
 8001dac:	9a3b      	ldr	r2, [sp, #236]	@ 0xec
 8001dae:	9f02      	ldr	r7, [sp, #8]
 8001db0:	f8dd 9010 	ldr.w	r9, [sp, #16]
 8001db4:	ea83 030b 	eor.w	r3, r3, fp
 8001db8:	404a      	eors	r2, r1
 8001dba:	ea4f 0c93 	mov.w	ip, r3, lsr #2
 8001dbe:	ea4c 7c82 	orr.w	ip, ip, r2, lsl #30
 8001dc2:	0892      	lsrs	r2, r2, #2
 8001dc4:	ea42 7283 	orr.w	r2, r2, r3, lsl #30
 8001dc8:	9906      	ldr	r1, [sp, #24]
 8001dca:	9b36      	ldr	r3, [sp, #216]	@ 0xd8
 8001dcc:	ea88 0404 	eor.w	r4, r8, r4
 8001dd0:	942e      	str	r4, [sp, #184]	@ 0xb8
 8001dd2:	ea85 040e 	eor.w	r4, r5, lr
 8001dd6:	404b      	eors	r3, r1
 8001dd8:	942f      	str	r4, [sp, #188]	@ 0xbc
 8001dda:	9937      	ldr	r1, [sp, #220]	@ 0xdc
 8001ddc:	9c07      	ldr	r4, [sp, #28]
 8001dde:	0a5d      	lsrs	r5, r3, #9
 8001de0:	ea81 0604 	eor.w	r6, r1, r4
 8001de4:	ea45 55c6 	orr.w	r5, r5, r6, lsl #23
 8001de8:	0a76      	lsrs	r6, r6, #9
 8001dea:	ea46 56c3 	orr.w	r6, r6, r3, lsl #23
 8001dee:	9908      	ldr	r1, [sp, #32]
 8001df0:	9b32      	ldr	r3, [sp, #200]	@ 0xc8
 8001df2:	9c09      	ldr	r4, [sp, #36]	@ 0x24
 8001df4:	404b      	eors	r3, r1
 8001df6:	9933      	ldr	r1, [sp, #204]	@ 0xcc
 8001df8:	ea4f 6853 	mov.w	r8, r3, lsr #25
 8001dfc:	404c      	eors	r4, r1
 8001dfe:	ea48 18c4 	orr.w	r8, r8, r4, lsl #7
 8001e02:	0e64      	lsrs	r4, r4, #25
 8001e04:	ea44 14c3 	orr.w	r4, r4, r3, lsl #7
 8001e08:	9901      	ldr	r1, [sp, #4]
 8001e0a:	9b38      	ldr	r3, [sp, #224]	@ 0xe0
 8001e0c:	404b      	eors	r3, r1
 8001e0e:	9939      	ldr	r1, [sp, #228]	@ 0xe4
 8001e10:	ea4f 5ed3 	mov.w	lr, r3, lsr #23
 8001e14:	4079      	eors	r1, r7
 8001e16:	ea4e 2e41 	orr.w	lr, lr, r1, lsl #9
 8001e1a:	0dc9      	lsrs	r1, r1, #23
 8001e1c:	ea41 2143 	orr.w	r1, r1, r3, lsl #9
 8001e20:	9f03      	ldr	r7, [sp, #12]
 8001e22:	9b34      	ldr	r3, [sp, #208]	@ 0xd0
 8001e24:	407b      	eors	r3, r7
 8001e26:	9f35      	ldr	r7, [sp, #212]	@ 0xd4
 8001e28:	ea87 0909 	eor.w	r9, r7, r9
 8001e2c:	ea4f 0789 	mov.w	r7, r9, lsl #2
 8001e30:	ea47 7793 	orr.w	r7, r7, r3, lsr #30
 8001e34:	009b      	lsls	r3, r3, #2
 8001e36:	ea43 7399 	orr.w	r3, r3, r9, lsr #30
 8001e3a:	ea28 0a05 	bic.w	sl, r8, r5
 8001e3e:	ea24 0906 	bic.w	r9, r4, r6
 8001e42:	ea89 0902 	eor.w	r9, r9, r2
 8001e46:	ea8a 0a0c 	eor.w	sl, sl, ip
 8001e4a:	f8cd a0c8 	str.w	sl, [sp, #200]	@ 0xc8
 8001e4e:	f8cd 90cc 	str.w	r9, [sp, #204]	@ 0xcc
 8001e52:	ea2e 0a08 	bic.w	sl, lr, r8
 8001e56:	ea21 0904 	bic.w	r9, r1, r4
 8001e5a:	ea8a 0a05 	eor.w	sl, sl, r5
 8001e5e:	ea89 0906 	eor.w	r9, r9, r6
 8001e62:	f8cd a0d0 	str.w	sl, [sp, #208]	@ 0xd0
 8001e66:	f8cd 90d4 	str.w	r9, [sp, #212]	@ 0xd4
 8001e6a:	ea23 0a0e 	bic.w	sl, r3, lr
 8001e6e:	ea27 0901 	bic.w	r9, r7, r1
 8001e72:	ea8a 0808 	eor.w	r8, sl, r8
 8001e76:	ea89 0404 	eor.w	r4, r9, r4
 8001e7a:	f8cd 80d8 	str.w	r8, [sp, #216]	@ 0xd8
 8001e7e:	9437      	str	r4, [sp, #220]	@ 0xdc
 8001e80:	ea2c 0803 	bic.w	r8, ip, r3
 8001e84:	ea22 0407 	bic.w	r4, r2, r7
 8001e88:	ea25 050c 	bic.w	r5, r5, ip
 8001e8c:	406b      	eors	r3, r5
 8001e8e:	ea88 0e0e 	eor.w	lr, r8, lr
 8001e92:	4061      	eors	r1, r4
 8001e94:	f8cd e0e0 	str.w	lr, [sp, #224]	@ 0xe0
 8001e98:	9139      	str	r1, [sp, #228]	@ 0xe4
 8001e9a:	933a      	str	r3, [sp, #232]	@ 0xe8
 8001e9c:	9b00      	ldr	r3, [sp, #0]
 8001e9e:	3310      	adds	r3, #16
 8001ea0:	9300      	str	r3, [sp, #0]
 8001ea2:	ea26 0602 	bic.w	r6, r6, r2
 8001ea6:	4b3d      	ldr	r3, [pc, #244]	@ (8001f9c <KeccakF1600_StatePermute+0xe04>)
 8001ea8:	9a00      	ldr	r2, [sp, #0]
 8001eaa:	4293      	cmp	r3, r2
 8001eac:	ea86 0607 	eor.w	r6, r6, r7
 8001eb0:	f47f a9f0 	bne.w	8001294 <KeccakF1600_StatePermute+0xfc>
 8001eb4:	9b0a      	ldr	r3, [sp, #40]	@ 0x28
 8001eb6:	6003      	str	r3, [r0, #0]
 8001eb8:	9b0b      	ldr	r3, [sp, #44]	@ 0x2c
 8001eba:	6043      	str	r3, [r0, #4]
 8001ebc:	9b0c      	ldr	r3, [sp, #48]	@ 0x30
 8001ebe:	6083      	str	r3, [r0, #8]
 8001ec0:	9b0d      	ldr	r3, [sp, #52]	@ 0x34
 8001ec2:	60c3      	str	r3, [r0, #12]
 8001ec4:	9b0e      	ldr	r3, [sp, #56]	@ 0x38
 8001ec6:	6103      	str	r3, [r0, #16]
 8001ec8:	9b0f      	ldr	r3, [sp, #60]	@ 0x3c
 8001eca:	6143      	str	r3, [r0, #20]
 8001ecc:	9b10      	ldr	r3, [sp, #64]	@ 0x40
 8001ece:	6183      	str	r3, [r0, #24]
 8001ed0:	9b11      	ldr	r3, [sp, #68]	@ 0x44
 8001ed2:	61c3      	str	r3, [r0, #28]
 8001ed4:	9b12      	ldr	r3, [sp, #72]	@ 0x48
 8001ed6:	6203      	str	r3, [r0, #32]
 8001ed8:	9b13      	ldr	r3, [sp, #76]	@ 0x4c
 8001eda:	6243      	str	r3, [r0, #36]	@ 0x24
 8001edc:	9b14      	ldr	r3, [sp, #80]	@ 0x50
 8001ede:	6283      	str	r3, [r0, #40]	@ 0x28
 8001ee0:	9b15      	ldr	r3, [sp, #84]	@ 0x54
 8001ee2:	62c3      	str	r3, [r0, #44]	@ 0x2c
 8001ee4:	9b16      	ldr	r3, [sp, #88]	@ 0x58
 8001ee6:	6303      	str	r3, [r0, #48]	@ 0x30
 8001ee8:	9b17      	ldr	r3, [sp, #92]	@ 0x5c
 8001eea:	6343      	str	r3, [r0, #52]	@ 0x34
 8001eec:	9b18      	ldr	r3, [sp, #96]	@ 0x60
 8001eee:	6383      	str	r3, [r0, #56]	@ 0x38
 8001ef0:	9b19      	ldr	r3, [sp, #100]	@ 0x64
 8001ef2:	63c3      	str	r3, [r0, #60]	@ 0x3c
 8001ef4:	9b1a      	ldr	r3, [sp, #104]	@ 0x68
 8001ef6:	6403      	str	r3, [r0, #64]	@ 0x40
 8001ef8:	9b1b      	ldr	r3, [sp, #108]	@ 0x6c
 8001efa:	6443      	str	r3, [r0, #68]	@ 0x44
 8001efc:	9b1c      	ldr	r3, [sp, #112]	@ 0x70
 8001efe:	6483      	str	r3, [r0, #72]	@ 0x48
 8001f00:	9b1d      	ldr	r3, [sp, #116]	@ 0x74
 8001f02:	64c3      	str	r3, [r0, #76]	@ 0x4c
 8001f04:	9b1e      	ldr	r3, [sp, #120]	@ 0x78
 8001f06:	6503      	str	r3, [r0, #80]	@ 0x50
 8001f08:	9b1f      	ldr	r3, [sp, #124]	@ 0x7c
 8001f0a:	6543      	str	r3, [r0, #84]	@ 0x54
 8001f0c:	9b20      	ldr	r3, [sp, #128]	@ 0x80
 8001f0e:	6583      	str	r3, [r0, #88]	@ 0x58
 8001f10:	9b21      	ldr	r3, [sp, #132]	@ 0x84
 8001f12:	65c3      	str	r3, [r0, #92]	@ 0x5c
 8001f14:	9b22      	ldr	r3, [sp, #136]	@ 0x88
 8001f16:	6603      	str	r3, [r0, #96]	@ 0x60
 8001f18:	9b23      	ldr	r3, [sp, #140]	@ 0x8c
 8001f1a:	6643      	str	r3, [r0, #100]	@ 0x64
 8001f1c:	9b24      	ldr	r3, [sp, #144]	@ 0x90
 8001f1e:	6683      	str	r3, [r0, #104]	@ 0x68
 8001f20:	9b25      	ldr	r3, [sp, #148]	@ 0x94
 8001f22:	66c3      	str	r3, [r0, #108]	@ 0x6c
 8001f24:	9b26      	ldr	r3, [sp, #152]	@ 0x98
 8001f26:	6703      	str	r3, [r0, #112]	@ 0x70
 8001f28:	9b27      	ldr	r3, [sp, #156]	@ 0x9c
 8001f2a:	6743      	str	r3, [r0, #116]	@ 0x74
 8001f2c:	9b28      	ldr	r3, [sp, #160]	@ 0xa0
 8001f2e:	6783      	str	r3, [r0, #120]	@ 0x78
 8001f30:	9b29      	ldr	r3, [sp, #164]	@ 0xa4
 8001f32:	67c3      	str	r3, [r0, #124]	@ 0x7c
 8001f34:	9b2a      	ldr	r3, [sp, #168]	@ 0xa8
 8001f36:	f8c0 3080 	str.w	r3, [r0, #128]	@ 0x80
 8001f3a:	9b2b      	ldr	r3, [sp, #172]	@ 0xac
 8001f3c:	f8c0 3084 	str.w	r3, [r0, #132]	@ 0x84
 8001f40:	9b2c      	ldr	r3, [sp, #176]	@ 0xb0
 8001f42:	f8c0 3088 	str.w	r3, [r0, #136]	@ 0x88
 8001f46:	9b2d      	ldr	r3, [sp, #180]	@ 0xb4
 8001f48:	f8c0 308c 	str.w	r3, [r0, #140]	@ 0x8c
 8001f4c:	9b2e      	ldr	r3, [sp, #184]	@ 0xb8
 8001f4e:	f8c0 3090 	str.w	r3, [r0, #144]	@ 0x90
 8001f52:	9b2f      	ldr	r3, [sp, #188]	@ 0xbc
 8001f54:	f8c0 3094 	str.w	r3, [r0, #148]	@ 0x94
 8001f58:	9b30      	ldr	r3, [sp, #192]	@ 0xc0
 8001f5a:	f8c0 3098 	str.w	r3, [r0, #152]	@ 0x98
 8001f5e:	9b31      	ldr	r3, [sp, #196]	@ 0xc4
 8001f60:	f8c0 309c 	str.w	r3, [r0, #156]	@ 0x9c
 8001f64:	9b32      	ldr	r3, [sp, #200]	@ 0xc8
 8001f66:	f8c0 30a0 	str.w	r3, [r0, #160]	@ 0xa0
 8001f6a:	9b33      	ldr	r3, [sp, #204]	@ 0xcc
 8001f6c:	f8c0 30a4 	str.w	r3, [r0, #164]	@ 0xa4
 8001f70:	9b34      	ldr	r3, [sp, #208]	@ 0xd0
 8001f72:	f8c0 30a8 	str.w	r3, [r0, #168]	@ 0xa8
 8001f76:	9b35      	ldr	r3, [sp, #212]	@ 0xd4
 8001f78:	f8c0 30ac 	str.w	r3, [r0, #172]	@ 0xac
 8001f7c:	9b36      	ldr	r3, [sp, #216]	@ 0xd8
 8001f7e:	f8c0 30b0 	str.w	r3, [r0, #176]	@ 0xb0
 8001f82:	9b37      	ldr	r3, [sp, #220]	@ 0xdc
 8001f84:	f8c0 30b4 	str.w	r3, [r0, #180]	@ 0xb4
 8001f88:	9b3a      	ldr	r3, [sp, #232]	@ 0xe8
 8001f8a:	f8c0 e0b8 	str.w	lr, [r0, #184]	@ 0xb8
 8001f8e:	f8c0 10bc 	str.w	r1, [r0, #188]	@ 0xbc
 8001f92:	e9c0 3630 	strd	r3, r6, [r0, #192]	@ 0xc0
 8001f96:	b03d      	add	sp, #244	@ 0xf4
 8001f98:	e8bd 8ff0 	ldmia.w	sp!, {r4, r5, r6, r7, r8, r9, sl, fp, pc}
 8001f9c:	080272f8 	.word	0x080272f8

08001fa0 <keccak_squeezeblocks>:
 8001fa0:	e92d 4ff8 	stmdb	sp!, {r3, r4, r5, r6, r7, r8, r9, sl, fp, lr}
 8001fa4:	4606      	mov	r6, r0
 8001fa6:	460d      	mov	r5, r1
 8001fa8:	4610      	mov	r0, r2
 8001faa:	461f      	mov	r7, r3
 8001fac:	ea4f 08d3 	mov.w	r8, r3, lsr #3
 8001fb0:	b90d      	cbnz	r5, 8001fb6 <keccak_squeezeblocks+0x16>
 8001fb2:	e8bd 8ff8 	ldmia.w	sp!, {r3, r4, r5, r6, r7, r8, r9, sl, fp, pc}
 8001fb6:	f7ff f8ef 	bl	8001198 <KeccakF1600_StatePermute>
 8001fba:	f1a0 0208 	sub.w	r2, r0, #8
 8001fbe:	2300      	movs	r3, #0
 8001fc0:	f852 bf08 	ldr.w	fp, [r2, #8]!
 8001fc4:	f8d2 a004 	ldr.w	sl, [r2, #4]
 8001fc8:	eb06 09c3 	add.w	r9, r6, r3, lsl #3
 8001fcc:	2400      	movs	r4, #0
 8001fce:	f1c4 0e20 	rsb	lr, r4, #32
 8001fd2:	f1a4 0c20 	sub.w	ip, r4, #32
 8001fd6:	fa2b f104 	lsr.w	r1, fp, r4
 8001fda:	fa0a fe0e 	lsl.w	lr, sl, lr
 8001fde:	ea41 010e 	orr.w	r1, r1, lr
 8001fe2:	fa2a fc0c 	lsr.w	ip, sl, ip
 8001fe6:	3408      	adds	r4, #8
 8001fe8:	ea41 010c 	orr.w	r1, r1, ip
 8001fec:	2c40      	cmp	r4, #64	@ 0x40
 8001fee:	f809 1b01 	strb.w	r1, [r9], #1
 8001ff2:	d1ec      	bne.n	8001fce <keccak_squeezeblocks+0x2e>
 8001ff4:	3301      	adds	r3, #1
 8001ff6:	4543      	cmp	r3, r8
 8001ff8:	d3e2      	bcc.n	8001fc0 <keccak_squeezeblocks+0x20>
 8001ffa:	443e      	add	r6, r7
 8001ffc:	3d01      	subs	r5, #1
 8001ffe:	e7d7      	b.n	8001fb0 <keccak_squeezeblocks+0x10>

08002000 <keccak_absorb>:
 8002000:	e92d 4ff0 	stmdb	sp!, {r4, r5, r6, r7, r8, r9, sl, fp, lr}
 8002004:	b0b5      	sub	sp, #212	@ 0xd4
 8002006:	460c      	mov	r4, r1
 8002008:	4617      	mov	r7, r2
 800200a:	2100      	movs	r1, #0
 800200c:	22c8      	movs	r2, #200	@ 0xc8
 800200e:	461d      	mov	r5, r3
 8002010:	9001      	str	r0, [sp, #4]
 8002012:	f7fe f8fd 	bl	8000210 <memset>
 8002016:	9b01      	ldr	r3, [sp, #4]
 8002018:	ea4f 08d4 	mov.w	r8, r4, lsr #3
 800201c:	f1a3 0608 	sub.w	r6, r3, #8
 8002020:	42a5      	cmp	r5, r4
 8002022:	d244      	bcs.n	80020ae <keccak_absorb+0xae>
 8002024:	4622      	mov	r2, r4
 8002026:	2100      	movs	r1, #0
 8002028:	a802      	add	r0, sp, #8
 800202a:	f7fe f8f1 	bl	8000210 <memset>
 800202e:	462a      	mov	r2, r5
 8002030:	4639      	mov	r1, r7
 8002032:	a802      	add	r0, sp, #8
 8002034:	f7fe fd9c 	bl	8000b70 <memcpy>
 8002038:	f105 03d0 	add.w	r3, r5, #208	@ 0xd0
 800203c:	eb0d 0503 	add.w	r5, sp, r3
 8002040:	f89d 30f8 	ldrb.w	r3, [sp, #248]	@ 0xf8
 8002044:	f805 3cc8 	strb.w	r3, [r5, #-200]
 8002048:	f104 03cf 	add.w	r3, r4, #207	@ 0xcf
 800204c:	eb0d 0403 	add.w	r4, sp, r3
 8002050:	aa02      	add	r2, sp, #8
 8002052:	f814 3cc8 	ldrb.w	r3, [r4, #-200]
 8002056:	f063 037f 	orn	r3, r3, #127	@ 0x7f
 800205a:	f804 3cc8 	strb.w	r3, [r4, #-200]
 800205e:	2400      	movs	r4, #0
 8002060:	2300      	movs	r3, #0
 8002062:	4696      	mov	lr, r2
 8002064:	461f      	mov	r7, r3
 8002066:	461d      	mov	r5, r3
 8002068:	f81e 0b01 	ldrb.w	r0, [lr], #1
 800206c:	f1a3 0120 	sub.w	r1, r3, #32
 8002070:	f1c3 0c20 	rsb	ip, r3, #32
 8002074:	fa00 f101 	lsl.w	r1, r0, r1
 8002078:	fa20 fc0c 	lsr.w	ip, r0, ip
 800207c:	4098      	lsls	r0, r3
 800207e:	3308      	adds	r3, #8
 8002080:	ea41 010c 	orr.w	r1, r1, ip
 8002084:	2b40      	cmp	r3, #64	@ 0x40
 8002086:	ea47 0700 	orr.w	r7, r7, r0
 800208a:	ea45 0501 	orr.w	r5, r5, r1
 800208e:	d1eb      	bne.n	8002068 <keccak_absorb+0x68>
 8002090:	f856 1f08 	ldr.w	r1, [r6, #8]!
 8002094:	6873      	ldr	r3, [r6, #4]
 8002096:	3401      	adds	r4, #1
 8002098:	4079      	eors	r1, r7
 800209a:	406b      	eors	r3, r5
 800209c:	45a0      	cmp	r8, r4
 800209e:	e9c6 1300 	strd	r1, r3, [r6]
 80020a2:	f102 0208 	add.w	r2, r2, #8
 80020a6:	d8db      	bhi.n	8002060 <keccak_absorb+0x60>
 80020a8:	b035      	add	sp, #212	@ 0xd4
 80020aa:	e8bd 8ff0 	ldmia.w	sp!, {r4, r5, r6, r7, r8, r9, sl, fp, pc}
 80020ae:	4632      	mov	r2, r6
 80020b0:	f04f 0c00 	mov.w	ip, #0
 80020b4:	2300      	movs	r3, #0
 80020b6:	eb07 0bcc 	add.w	fp, r7, ip, lsl #3
 80020ba:	4699      	mov	r9, r3
 80020bc:	469e      	mov	lr, r3
 80020be:	f81b 0b01 	ldrb.w	r0, [fp], #1
 80020c2:	f1a3 0a20 	sub.w	sl, r3, #32
 80020c6:	f1c3 0120 	rsb	r1, r3, #32
 80020ca:	fa00 fa0a 	lsl.w	sl, r0, sl
 80020ce:	fa20 f101 	lsr.w	r1, r0, r1
 80020d2:	4098      	lsls	r0, r3
 80020d4:	3308      	adds	r3, #8
 80020d6:	ea4a 0101 	orr.w	r1, sl, r1
 80020da:	2b40      	cmp	r3, #64	@ 0x40
 80020dc:	ea40 0909 	orr.w	r9, r0, r9
 80020e0:	ea41 0e0e 	orr.w	lr, r1, lr
 80020e4:	d1eb      	bne.n	80020be <keccak_absorb+0xbe>
 80020e6:	f852 1f08 	ldr.w	r1, [r2, #8]!
 80020ea:	6853      	ldr	r3, [r2, #4]
 80020ec:	f10c 0c01 	add.w	ip, ip, #1
 80020f0:	ea81 0109 	eor.w	r1, r1, r9
 80020f4:	ea83 030e 	eor.w	r3, r3, lr
 80020f8:	45e0      	cmp	r8, ip
 80020fa:	e9c2 1300 	strd	r1, r3, [r2]
 80020fe:	d8d9      	bhi.n	80020b4 <keccak_absorb+0xb4>
 8002100:	9801      	ldr	r0, [sp, #4]
 8002102:	f7ff f849 	bl	8001198 <KeccakF1600_StatePermute>
 8002106:	1b2d      	subs	r5, r5, r4
 8002108:	4427      	add	r7, r4
 800210a:	e789      	b.n	8002020 <keccak_absorb+0x20>

0800210c <shake256_absorb>:
 800210c:	b5f7      	push	{r0, r1, r2, r4, r5, r6, r7, lr}
 800210e:	4607      	mov	r7, r0
 8002110:	20c8      	movs	r0, #200	@ 0xc8
 8002112:	460d      	mov	r5, r1
 8002114:	4616      	mov	r6, r2
 8002116:	f7fe fd1b 	bl	8000b50 <malloc>
 800211a:	6038      	str	r0, [r7, #0]
 800211c:	b910      	cbnz	r0, 8002124 <shake256_absorb+0x18>
 800211e:	206f      	movs	r0, #111	@ 0x6f
 8002120:	f7fe fdfa 	bl	8000d18 <exit>
 8002124:	231f      	movs	r3, #31
 8002126:	9300      	str	r3, [sp, #0]
 8002128:	462a      	mov	r2, r5
 800212a:	4633      	mov	r3, r6
 800212c:	2188      	movs	r1, #136	@ 0x88
 800212e:	f7ff ff67 	bl	8002000 <keccak_absorb>
 8002132:	b003      	add	sp, #12
 8002134:	bdf0      	pop	{r4, r5, r6, r7, pc}

08002136 <shake256_squeezeblocks>:
 8002136:	6812      	ldr	r2, [r2, #0]
 8002138:	2388      	movs	r3, #136	@ 0x88
 800213a:	f7ff bf31 	b.w	8001fa0 <keccak_squeezeblocks>

0800213e <shake256_ctx_release>:
 800213e:	6800      	ldr	r0, [r0, #0]
 8002140:	f7fe bd0e 	b.w	8000b60 <free>

08002144 <shake256>:
 8002144:	b5f0      	push	{r4, r5, r6, r7, lr}
 8002146:	b0a5      	sub	sp, #148	@ 0x94
 8002148:	460c      	mov	r4, r1
 800214a:	4606      	mov	r6, r0
 800214c:	4611      	mov	r1, r2
 800214e:	2788      	movs	r7, #136	@ 0x88
 8002150:	461a      	mov	r2, r3
 8002152:	a801      	add	r0, sp, #4
 8002154:	f7ff ffda 	bl	800210c <shake256_absorb>
 8002158:	fbb4 f5f7 	udiv	r5, r4, r7
 800215c:	aa01      	add	r2, sp, #4
 800215e:	4629      	mov	r1, r5
 8002160:	4630      	mov	r0, r6
 8002162:	437d      	muls	r5, r7
 8002164:	f7ff ffe7 	bl	8002136 <shake256_squeezeblocks>
 8002168:	1b64      	subs	r4, r4, r5
 800216a:	d009      	beq.n	8002180 <shake256+0x3c>
 800216c:	aa01      	add	r2, sp, #4
 800216e:	2101      	movs	r1, #1
 8002170:	a802      	add	r0, sp, #8
 8002172:	f7ff ffe0 	bl	8002136 <shake256_squeezeblocks>
 8002176:	4622      	mov	r2, r4
 8002178:	a902      	add	r1, sp, #8
 800217a:	1970      	adds	r0, r6, r5
 800217c:	f7fe fcf8 	bl	8000b70 <memcpy>
 8002180:	a801      	add	r0, sp, #4
 8002182:	f7ff ffdc 	bl	800213e <shake256_ctx_release>
 8002186:	b025      	add	sp, #148	@ 0x94
 8002188:	bdf0      	pop	{r4, r5, r6, r7, pc}

0800218a <pqmayo_MAYO_1_m4_m_upper>:
 800218a:	e92d 4ff0 	stmdb	sp!, {r4, r5, r6, r7, r8, r9, sl, fp, lr}
 800218e:	b085      	sub	sp, #20
 8002190:	f04f 0828 	mov.w	r8, #40	@ 0x28
 8002194:	9203      	str	r2, [sp, #12]
 8002196:	fb03 8208 	mla	r2, r3, r8, r8
 800219a:	2000      	movs	r0, #0
 800219c:	9202      	str	r2, [sp, #8]
 800219e:	3a28      	subs	r2, #40	@ 0x28
 80021a0:	9100      	str	r1, [sp, #0]
 80021a2:	9201      	str	r2, [sp, #4]
 80021a4:	4601      	mov	r1, r0
 80021a6:	4686      	mov	lr, r0
 80021a8:	4299      	cmp	r1, r3
 80021aa:	da34      	bge.n	8002216 <pqmayo_MAYO_1_m4_m_upper+0x8c>
 80021ac:	9a03      	ldr	r2, [sp, #12]
 80021ae:	fb08 220e 	mla	r2, r8, lr, r2
 80021b2:	3a08      	subs	r2, #8
 80021b4:	4684      	mov	ip, r0
 80021b6:	4683      	mov	fp, r0
 80021b8:	460c      	mov	r4, r1
 80021ba:	e000      	b.n	80021be <pqmayo_MAYO_1_m4_m_upper+0x34>
 80021bc:	464a      	mov	r2, r9
 80021be:	9d00      	ldr	r5, [sp, #0]
 80021c0:	f102 0928 	add.w	r9, r2, #40	@ 0x28
 80021c4:	eb0c 0a05 	add.w	sl, ip, r5
 80021c8:	4615      	mov	r5, r2
 80021ca:	e8fa 6702 	ldrd	r6, r7, [sl], #8
 80021ce:	e9e5 6702 	strd	r6, r7, [r5, #8]!
 80021d2:	454d      	cmp	r5, r9
 80021d4:	d1f9      	bne.n	80021ca <pqmayo_MAYO_1_m4_m_upper+0x40>
 80021d6:	42a1      	cmp	r1, r4
 80021d8:	d010      	beq.n	80021fc <pqmayo_MAYO_1_m4_m_upper+0x72>
 80021da:	9d00      	ldr	r5, [sp, #0]
 80021dc:	445d      	add	r5, fp
 80021de:	e9d5 a600 	ldrd	sl, r6, [r5]
 80021e2:	f852 7f08 	ldr.w	r7, [r2, #8]!
 80021e6:	ea87 0a0a 	eor.w	sl, r7, sl
 80021ea:	6857      	ldr	r7, [r2, #4]
 80021ec:	454a      	cmp	r2, r9
 80021ee:	ea86 0607 	eor.w	r6, r6, r7
 80021f2:	e9c2 a600 	strd	sl, r6, [r2]
 80021f6:	f105 0508 	add.w	r5, r5, #8
 80021fa:	d1f0      	bne.n	80021de <pqmayo_MAYO_1_m4_m_upper+0x54>
 80021fc:	9a01      	ldr	r2, [sp, #4]
 80021fe:	3401      	adds	r4, #1
 8002200:	42a3      	cmp	r3, r4
 8002202:	4493      	add	fp, r2
 8002204:	f10c 0c28 	add.w	ip, ip, #40	@ 0x28
 8002208:	d1d8      	bne.n	80021bc <pqmayo_MAYO_1_m4_m_upper+0x32>
 800220a:	1a5a      	subs	r2, r3, r1
 800220c:	4496      	add	lr, r2
 800220e:	9a02      	ldr	r2, [sp, #8]
 8002210:	3101      	adds	r1, #1
 8002212:	4410      	add	r0, r2
 8002214:	e7c8      	b.n	80021a8 <pqmayo_MAYO_1_m4_m_upper+0x1e>
 8002216:	b005      	add	sp, #20
 8002218:	e8bd 8ff0 	ldmia.w	sp!, {r4, r5, r6, r7, r8, r9, sl, fp, pc}

0800221c <check_version>:
 800221c:	2001      	movs	r0, #1
 800221e:	4770      	bx	lr

08002220 <simpleserial_put.part.0>:
 8002220:	b5f8      	push	{r3, r4, r5, r6, r7, lr}
 8002222:	4614      	mov	r4, r2
 8002224:	460e      	mov	r6, r1
 8002226:	f000 f9e3 	bl	80025f0 <putch>
 800222a:	4f0c      	ldr	r7, [pc, #48]	@ (800225c <simpleserial_put.part.0+0x3c>)
 800222c:	1e65      	subs	r5, r4, #1
 800222e:	f1c4 0401 	rsb	r4, r4, #1
 8002232:	1963      	adds	r3, r4, r5
 8002234:	429e      	cmp	r6, r3
 8002236:	dc04      	bgt.n	8002242 <simpleserial_put.part.0+0x22>
 8002238:	e8bd 40f8 	ldmia.w	sp!, {r3, r4, r5, r6, r7, lr}
 800223c:	200a      	movs	r0, #10
 800223e:	f000 b9d7 	b.w	80025f0 <putch>
 8002242:	f815 3f01 	ldrb.w	r3, [r5, #1]!
 8002246:	091b      	lsrs	r3, r3, #4
 8002248:	5cf8      	ldrb	r0, [r7, r3]
 800224a:	f000 f9d1 	bl	80025f0 <putch>
 800224e:	782b      	ldrb	r3, [r5, #0]
 8002250:	f003 030f 	and.w	r3, r3, #15
 8002254:	5cf8      	ldrb	r0, [r7, r3]
 8002256:	f000 f9cb 	bl	80025f0 <putch>
 800225a:	e7ea      	b.n	8002232 <simpleserial_put.part.0+0x12>
 800225c:	080272f8 	.word	0x080272f8

08002260 <ss_num_commands>:
 8002260:	b507      	push	{r0, r1, r2, lr}
 8002262:	4b07      	ldr	r3, [pc, #28]	@ (8002280 <ss_num_commands+0x20>)
 8002264:	f10d 0207 	add.w	r2, sp, #7
 8002268:	681b      	ldr	r3, [r3, #0]
 800226a:	f88d 3007 	strb.w	r3, [sp, #7]
 800226e:	2101      	movs	r1, #1
 8002270:	2072      	movs	r0, #114	@ 0x72
 8002272:	f7ff ffd5 	bl	8002220 <simpleserial_put.part.0>
 8002276:	2000      	movs	r0, #0
 8002278:	b003      	add	sp, #12
 800227a:	f85d fb04 	ldr.w	pc, [sp], #4
 800227e:	bf00      	nop
 8002280:	200066d8 	.word	0x200066d8

08002284 <ss_get_commands>:
 8002284:	b570      	push	{r4, r5, r6, lr}
 8002286:	4c12      	ldr	r4, [pc, #72]	@ (80022d0 <ss_get_commands+0x4c>)
 8002288:	6821      	ldr	r1, [r4, #0]
 800228a:	b098      	sub	sp, #96	@ 0x60
 800228c:	b2cd      	uxtb	r5, r1
 800228e:	2000      	movs	r0, #0
 8002290:	b2c3      	uxtb	r3, r0
 8002292:	42ab      	cmp	r3, r5
 8002294:	f100 0001 	add.w	r0, r0, #1
 8002298:	db09      	blt.n	80022ae <ss_get_commands+0x2a>
 800229a:	eb01 0141 	add.w	r1, r1, r1, lsl #1
 800229e:	466a      	mov	r2, sp
 80022a0:	b2c9      	uxtb	r1, r1
 80022a2:	2072      	movs	r0, #114	@ 0x72
 80022a4:	f7ff ffbc 	bl	8002220 <simpleserial_put.part.0>
 80022a8:	2000      	movs	r0, #0
 80022aa:	b018      	add	sp, #96	@ 0x60
 80022ac:	bd70      	pop	{r4, r5, r6, pc}
 80022ae:	eb03 0243 	add.w	r2, r3, r3, lsl #1
 80022b2:	3260      	adds	r2, #96	@ 0x60
 80022b4:	eb04 1303 	add.w	r3, r4, r3, lsl #4
 80022b8:	446a      	add	r2, sp
 80022ba:	791e      	ldrb	r6, [r3, #4]
 80022bc:	f802 6c60 	strb.w	r6, [r2, #-96]
 80022c0:	689e      	ldr	r6, [r3, #8]
 80022c2:	7c1b      	ldrb	r3, [r3, #16]
 80022c4:	f802 6c5f 	strb.w	r6, [r2, #-95]
 80022c8:	f802 3c5e 	strb.w	r3, [r2, #-94]
 80022cc:	e7e0      	b.n	8002290 <ss_get_commands+0xc>
 80022ce:	bf00      	nop
 80022d0:	200066d8 	.word	0x200066d8

080022d4 <hex_decode>:
 80022d4:	b5f0      	push	{r4, r5, r6, r7, lr}
 80022d6:	2500      	movs	r5, #0
 80022d8:	1c4f      	adds	r7, r1, #1
 80022da:	4285      	cmp	r5, r0
 80022dc:	db01      	blt.n	80022e2 <hex_decode+0xe>
 80022de:	2000      	movs	r0, #0
 80022e0:	e021      	b.n	8002326 <hex_decode+0x52>
 80022e2:	f817 4015 	ldrb.w	r4, [r7, r5, lsl #1]
 80022e6:	f811 3015 	ldrb.w	r3, [r1, r5, lsl #1]
 80022ea:	f1a4 0630 	sub.w	r6, r4, #48	@ 0x30
 80022ee:	b2f6      	uxtb	r6, r6
 80022f0:	2e09      	cmp	r6, #9
 80022f2:	d80c      	bhi.n	800230e <hex_decode+0x3a>
 80022f4:	7016      	strb	r6, [r2, #0]
 80022f6:	f1a3 0430 	sub.w	r4, r3, #48	@ 0x30
 80022fa:	b2e4      	uxtb	r4, r4
 80022fc:	2c09      	cmp	r4, #9
 80022fe:	d815      	bhi.n	800232c <hex_decode+0x58>
 8002300:	7813      	ldrb	r3, [r2, #0]
 8002302:	ea43 1304 	orr.w	r3, r3, r4, lsl #4
 8002306:	7013      	strb	r3, [r2, #0]
 8002308:	3501      	adds	r5, #1
 800230a:	3201      	adds	r2, #1
 800230c:	e7e5      	b.n	80022da <hex_decode+0x6>
 800230e:	f1a4 0641 	sub.w	r6, r4, #65	@ 0x41
 8002312:	2e05      	cmp	r6, #5
 8002314:	d802      	bhi.n	800231c <hex_decode+0x48>
 8002316:	3c37      	subs	r4, #55	@ 0x37
 8002318:	7014      	strb	r4, [r2, #0]
 800231a:	e7ec      	b.n	80022f6 <hex_decode+0x22>
 800231c:	f1a4 0661 	sub.w	r6, r4, #97	@ 0x61
 8002320:	2e05      	cmp	r6, #5
 8002322:	d901      	bls.n	8002328 <hex_decode+0x54>
 8002324:	2001      	movs	r0, #1
 8002326:	bdf0      	pop	{r4, r5, r6, r7, pc}
 8002328:	3c57      	subs	r4, #87	@ 0x57
 800232a:	e7f5      	b.n	8002318 <hex_decode+0x44>
 800232c:	f1a3 0441 	sub.w	r4, r3, #65	@ 0x41
 8002330:	2c05      	cmp	r4, #5
 8002332:	d802      	bhi.n	800233a <hex_decode+0x66>
 8002334:	f1a3 0437 	sub.w	r4, r3, #55	@ 0x37
 8002338:	e7e2      	b.n	8002300 <hex_decode+0x2c>
 800233a:	f1a3 0461 	sub.w	r4, r3, #97	@ 0x61
 800233e:	2c05      	cmp	r4, #5
 8002340:	d8f0      	bhi.n	8002324 <hex_decode+0x50>
 8002342:	7814      	ldrb	r4, [r2, #0]
 8002344:	3b57      	subs	r3, #87	@ 0x57
 8002346:	ea44 1303 	orr.w	r3, r4, r3, lsl #4
 800234a:	e7dc      	b.n	8002306 <hex_decode+0x32>

0800234c <simpleserial_addcmd_flags>:
 800234c:	b570      	push	{r4, r5, r6, lr}
 800234e:	4e09      	ldr	r6, [pc, #36]	@ (8002374 <simpleserial_addcmd_flags+0x28>)
 8002350:	6834      	ldr	r4, [r6, #0]
 8002352:	2c1f      	cmp	r4, #31
 8002354:	dc0b      	bgt.n	800236e <simpleserial_addcmd_flags+0x22>
 8002356:	293f      	cmp	r1, #63	@ 0x3f
 8002358:	d809      	bhi.n	800236e <simpleserial_addcmd_flags+0x22>
 800235a:	eb06 1504 	add.w	r5, r6, r4, lsl #4
 800235e:	e9c5 1202 	strd	r1, r2, [r5, #8]
 8002362:	3401      	adds	r4, #1
 8002364:	7128      	strb	r0, [r5, #4]
 8002366:	742b      	strb	r3, [r5, #16]
 8002368:	6034      	str	r4, [r6, #0]
 800236a:	2000      	movs	r0, #0
 800236c:	bd70      	pop	{r4, r5, r6, pc}
 800236e:	2001      	movs	r0, #1
 8002370:	e7fc      	b.n	800236c <simpleserial_addcmd_flags+0x20>
 8002372:	bf00      	nop
 8002374:	200066d8 	.word	0x200066d8

08002378 <simpleserial_addcmd>:
 8002378:	2300      	movs	r3, #0
 800237a:	f7ff bfe7 	b.w	800234c <simpleserial_addcmd_flags>
	...

08002380 <simpleserial_init>:
 8002380:	b508      	push	{r3, lr}
 8002382:	4a07      	ldr	r2, [pc, #28]	@ (80023a0 <simpleserial_init+0x20>)
 8002384:	2100      	movs	r1, #0
 8002386:	2076      	movs	r0, #118	@ 0x76
 8002388:	f7ff fff6 	bl	8002378 <simpleserial_addcmd>
 800238c:	4a05      	ldr	r2, [pc, #20]	@ (80023a4 <simpleserial_init+0x24>)
 800238e:	2077      	movs	r0, #119	@ 0x77
 8002390:	f7ff fff2 	bl	8002378 <simpleserial_addcmd>
 8002394:	e8bd 4008 	ldmia.w	sp!, {r3, lr}
 8002398:	4a03      	ldr	r2, [pc, #12]	@ (80023a8 <simpleserial_init+0x28>)
 800239a:	2079      	movs	r0, #121	@ 0x79
 800239c:	f7ff bfec 	b.w	8002378 <simpleserial_addcmd>
 80023a0:	0800221d 	.word	0x0800221d
 80023a4:	08002285 	.word	0x08002285
 80023a8:	08002261 	.word	0x08002261

080023ac <simpleserial_get>:
 80023ac:	e92d 41f0 	stmdb	sp!, {r4, r5, r6, r7, r8, lr}
 80023b0:	4d2c      	ldr	r5, [pc, #176]	@ (8002464 <simpleserial_get+0xb8>)
 80023b2:	b0b2      	sub	sp, #200	@ 0xc8
 80023b4:	f000 f90a 	bl	80025cc <getch>
 80023b8:	462a      	mov	r2, r5
 80023ba:	2300      	movs	r3, #0
 80023bc:	f852 1b04 	ldr.w	r1, [r2], #4
 80023c0:	4299      	cmp	r1, r3
 80023c2:	dc3f      	bgt.n	8002444 <simpleserial_get+0x98>
 80023c4:	d03b      	beq.n	800243e <simpleserial_get+0x92>
 80023c6:	eb05 1403 	add.w	r4, r5, r3, lsl #4
 80023ca:	011e      	lsls	r6, r3, #4
 80023cc:	7c23      	ldrb	r3, [r4, #16]
 80023ce:	07db      	lsls	r3, r3, #31
 80023d0:	d513      	bpl.n	80023fa <simpleserial_get+0x4e>
 80023d2:	2300      	movs	r3, #0
 80023d4:	f88d 3008 	strb.w	r3, [sp, #8]
 80023d8:	f000 f8f8 	bl	80025cc <getch>
 80023dc:	f88d 0048 	strb.w	r0, [sp, #72]	@ 0x48
 80023e0:	f000 f8f4 	bl	80025cc <getch>
 80023e4:	aa02      	add	r2, sp, #8
 80023e6:	f88d 0049 	strb.w	r0, [sp, #73]	@ 0x49
 80023ea:	a912      	add	r1, sp, #72	@ 0x48
 80023ec:	2001      	movs	r0, #1
 80023ee:	f7ff ff71 	bl	80022d4 <hex_decode>
 80023f2:	bb20      	cbnz	r0, 800243e <simpleserial_get+0x92>
 80023f4:	f89d 3008 	ldrb.w	r3, [sp, #8]
 80023f8:	60a3      	str	r3, [r4, #8]
 80023fa:	af12      	add	r7, sp, #72	@ 0x48
 80023fc:	2400      	movs	r4, #0
 80023fe:	eb05 0806 	add.w	r8, r5, r6
 8002402:	f8d8 3008 	ldr.w	r3, [r8, #8]
 8002406:	ebb4 0f43 	cmp.w	r4, r3, lsl #1
 800240a:	d321      	bcc.n	8002450 <simpleserial_get+0xa4>
 800240c:	f000 f8de 	bl	80025cc <getch>
 8002410:	280a      	cmp	r0, #10
 8002412:	d001      	beq.n	8002418 <simpleserial_get+0x6c>
 8002414:	280d      	cmp	r0, #13
 8002416:	d112      	bne.n	800243e <simpleserial_get+0x92>
 8002418:	4435      	add	r5, r6
 800241a:	aa02      	add	r2, sp, #8
 800241c:	68ac      	ldr	r4, [r5, #8]
 800241e:	a912      	add	r1, sp, #72	@ 0x48
 8002420:	4620      	mov	r0, r4
 8002422:	f7ff ff57 	bl	80022d4 <hex_decode>
 8002426:	b950      	cbnz	r0, 800243e <simpleserial_get+0x92>
 8002428:	b2e1      	uxtb	r1, r4
 800242a:	68eb      	ldr	r3, [r5, #12]
 800242c:	a802      	add	r0, sp, #8
 800242e:	4798      	blx	r3
 8002430:	aa01      	add	r2, sp, #4
 8002432:	f88d 0004 	strb.w	r0, [sp, #4]
 8002436:	2101      	movs	r1, #1
 8002438:	207a      	movs	r0, #122	@ 0x7a
 800243a:	f7ff fef1 	bl	8002220 <simpleserial_put.part.0>
 800243e:	b032      	add	sp, #200	@ 0xc8
 8002440:	e8bd 81f0 	ldmia.w	sp!, {r4, r5, r6, r7, r8, pc}
 8002444:	011c      	lsls	r4, r3, #4
 8002446:	5ca4      	ldrb	r4, [r4, r2]
 8002448:	4284      	cmp	r4, r0
 800244a:	d0bc      	beq.n	80023c6 <simpleserial_get+0x1a>
 800244c:	3301      	adds	r3, #1
 800244e:	e7b7      	b.n	80023c0 <simpleserial_get+0x14>
 8002450:	f000 f8bc 	bl	80025cc <getch>
 8002454:	280a      	cmp	r0, #10
 8002456:	d0f2      	beq.n	800243e <simpleserial_get+0x92>
 8002458:	280d      	cmp	r0, #13
 800245a:	d0f0      	beq.n	800243e <simpleserial_get+0x92>
 800245c:	f807 0b01 	strb.w	r0, [r7], #1
 8002460:	3401      	adds	r4, #1
 8002462:	e7ce      	b.n	8002402 <simpleserial_get+0x56>
 8002464:	200066d8 	.word	0x200066d8

08002468 <simpleserial_put>:
 8002468:	b10a      	cbz	r2, 800246e <simpleserial_put+0x6>
 800246a:	f7ff bed9 	b.w	8002220 <simpleserial_put.part.0>
 800246e:	4770      	bx	lr

08002470 <platform_init>:
 8002470:	b530      	push	{r4, r5, lr}
 8002472:	2203      	movs	r2, #3
 8002474:	b093      	sub	sp, #76	@ 0x4c
 8002476:	f44f 23a0 	mov.w	r3, #327680	@ 0x50000
 800247a:	e9cd 2306 	strd	r2, r3, [sp, #24]
 800247e:	2102      	movs	r1, #2
 8002480:	f44f 0380 	mov.w	r3, #4194304	@ 0x400000
 8002484:	e9cd 130c 	strd	r1, r3, [sp, #48]	@ 0x30
 8002488:	200c      	movs	r0, #12
 800248a:	23c4      	movs	r3, #196	@ 0xc4
 800248c:	e9cd 030e 	strd	r0, r3, [sp, #56]	@ 0x38
 8002490:	2404      	movs	r4, #4
 8002492:	2307      	movs	r3, #7
 8002494:	2501      	movs	r5, #1
 8002496:	a806      	add	r0, sp, #24
 8002498:	e9cd 4310 	strd	r4, r3, [sp, #64]	@ 0x40
 800249c:	9509      	str	r5, [sp, #36]	@ 0x24
 800249e:	f000 f8bf 	bl	8002620 <HAL_RCC_OscConfig>
 80024a2:	4604      	mov	r4, r0
 80024a4:	b100      	cbz	r0, 80024a8 <platform_init+0x38>
 80024a6:	e7fe      	b.n	80024a6 <platform_init+0x36>
 80024a8:	4601      	mov	r1, r0
 80024aa:	e9cd 0003 	strd	r0, r0, [sp, #12]
 80024ae:	230f      	movs	r3, #15
 80024b0:	9005      	str	r0, [sp, #20]
 80024b2:	a801      	add	r0, sp, #4
 80024b4:	e9cd 3501 	strd	r3, r5, [sp, #4]
 80024b8:	f000 fa50 	bl	800295c <HAL_RCC_ClockConfig>
 80024bc:	4b09      	ldr	r3, [pc, #36]	@ (80024e4 <platform_init+0x74>)
 80024be:	9400      	str	r4, [sp, #0]
 80024c0:	6b5a      	ldr	r2, [r3, #52]	@ 0x34
 80024c2:	4809      	ldr	r0, [pc, #36]	@ (80024e8 <platform_init+0x78>)
 80024c4:	f042 0240 	orr.w	r2, r2, #64	@ 0x40
 80024c8:	635a      	str	r2, [r3, #52]	@ 0x34
 80024ca:	6b5b      	ldr	r3, [r3, #52]	@ 0x34
 80024cc:	f003 0340 	and.w	r3, r3, #64	@ 0x40
 80024d0:	9300      	str	r3, [sp, #0]
 80024d2:	9b00      	ldr	r3, [sp, #0]
 80024d4:	4b05      	ldr	r3, [pc, #20]	@ (80024ec <platform_init+0x7c>)
 80024d6:	6003      	str	r3, [r0, #0]
 80024d8:	7144      	strb	r4, [r0, #5]
 80024da:	f000 fccc 	bl	8002e76 <HAL_RNG_Init>
 80024de:	b013      	add	sp, #76	@ 0x4c
 80024e0:	bd30      	pop	{r4, r5, pc}
 80024e2:	bf00      	nop
 80024e4:	40023800 	.word	0x40023800
 80024e8:	200068dc 	.word	0x200068dc
 80024ec:	50060800 	.word	0x50060800

080024f0 <init_uart>:
 80024f0:	b530      	push	{r4, r5, lr}
 80024f2:	2302      	movs	r3, #2
 80024f4:	b089      	sub	sp, #36	@ 0x24
 80024f6:	f44f 62c0 	mov.w	r2, #1536	@ 0x600
 80024fa:	e9cd 2303 	strd	r2, r3, [sp, #12]
 80024fe:	2201      	movs	r2, #1
 8002500:	e9cd 2305 	strd	r2, r3, [sp, #20]
 8002504:	4c15      	ldr	r4, [pc, #84]	@ (800255c <init_uart+0x6c>)
 8002506:	4816      	ldr	r0, [pc, #88]	@ (8002560 <init_uart+0x70>)
 8002508:	2500      	movs	r5, #0
 800250a:	2307      	movs	r3, #7
 800250c:	9501      	str	r5, [sp, #4]
 800250e:	9307      	str	r3, [sp, #28]
 8002510:	6b23      	ldr	r3, [r4, #48]	@ 0x30
 8002512:	4313      	orrs	r3, r2
 8002514:	6323      	str	r3, [r4, #48]	@ 0x30
 8002516:	6b23      	ldr	r3, [r4, #48]	@ 0x30
 8002518:	4013      	ands	r3, r2
 800251a:	9301      	str	r3, [sp, #4]
 800251c:	a903      	add	r1, sp, #12
 800251e:	9b01      	ldr	r3, [sp, #4]
 8002520:	f000 fab4 	bl	8002a8c <HAL_GPIO_Init>
 8002524:	480f      	ldr	r0, [pc, #60]	@ (8002564 <init_uart+0x74>)
 8002526:	4910      	ldr	r1, [pc, #64]	@ (8002568 <init_uart+0x78>)
 8002528:	9502      	str	r5, [sp, #8]
 800252a:	f44f 4316 	mov.w	r3, #38400	@ 0x9600
 800252e:	e9c0 1303 	strd	r1, r3, [r0, #12]
 8002532:	e9c0 5505 	strd	r5, r5, [r0, #20]
 8002536:	230c      	movs	r3, #12
 8002538:	61c5      	str	r5, [r0, #28]
 800253a:	6245      	str	r5, [r0, #36]	@ 0x24
 800253c:	6203      	str	r3, [r0, #32]
 800253e:	6c63      	ldr	r3, [r4, #68]	@ 0x44
 8002540:	f043 0310 	orr.w	r3, r3, #16
 8002544:	6463      	str	r3, [r4, #68]	@ 0x44
 8002546:	6c63      	ldr	r3, [r4, #68]	@ 0x44
 8002548:	f003 0310 	and.w	r3, r3, #16
 800254c:	9302      	str	r3, [sp, #8]
 800254e:	300c      	adds	r0, #12
 8002550:	9b02      	ldr	r3, [sp, #8]
 8002552:	f000 fb79 	bl	8002c48 <HAL_UART_Init>
 8002556:	b009      	add	sp, #36	@ 0x24
 8002558:	bd30      	pop	{r4, r5, pc}
 800255a:	bf00      	nop
 800255c:	40023800 	.word	0x40023800
 8002560:	40020000 	.word	0x40020000
 8002564:	200068dc 	.word	0x200068dc
 8002568:	40011000 	.word	0x40011000

0800256c <trigger_setup>:
 800256c:	b57f      	push	{r0, r1, r2, r3, r4, r5, r6, lr}
 800256e:	4b0d      	ldr	r3, [pc, #52]	@ (80025a4 <trigger_setup+0x38>)
 8002570:	480d      	ldr	r0, [pc, #52]	@ (80025a8 <trigger_setup+0x3c>)
 8002572:	2100      	movs	r1, #0
 8002574:	9100      	str	r1, [sp, #0]
 8002576:	6b1a      	ldr	r2, [r3, #48]	@ 0x30
 8002578:	f042 0201 	orr.w	r2, r2, #1
 800257c:	631a      	str	r2, [r3, #48]	@ 0x30
 800257e:	6b1b      	ldr	r3, [r3, #48]	@ 0x30
 8002580:	9103      	str	r1, [sp, #12]
 8002582:	f003 0301 	and.w	r3, r3, #1
 8002586:	9300      	str	r3, [sp, #0]
 8002588:	9b00      	ldr	r3, [sp, #0]
 800258a:	f44f 5280 	mov.w	r2, #4096	@ 0x1000
 800258e:	2301      	movs	r3, #1
 8002590:	e9cd 2301 	strd	r2, r3, [sp, #4]
 8002594:	a901      	add	r1, sp, #4
 8002596:	2302      	movs	r3, #2
 8002598:	9304      	str	r3, [sp, #16]
 800259a:	f000 fa77 	bl	8002a8c <HAL_GPIO_Init>
 800259e:	b007      	add	sp, #28
 80025a0:	f85d fb04 	ldr.w	pc, [sp], #4
 80025a4:	40023800 	.word	0x40023800
 80025a8:	40020000 	.word	0x40020000

080025ac <trigger_high>:
 80025ac:	4802      	ldr	r0, [pc, #8]	@ (80025b8 <trigger_high+0xc>)
 80025ae:	2201      	movs	r2, #1
 80025b0:	f44f 5180 	mov.w	r1, #4096	@ 0x1000
 80025b4:	f000 bb42 	b.w	8002c3c <HAL_GPIO_WritePin>
 80025b8:	40020000 	.word	0x40020000

080025bc <trigger_low>:
 80025bc:	4802      	ldr	r0, [pc, #8]	@ (80025c8 <trigger_low+0xc>)
 80025be:	2200      	movs	r2, #0
 80025c0:	f44f 5180 	mov.w	r1, #4096	@ 0x1000
 80025c4:	f000 bb3a 	b.w	8002c3c <HAL_GPIO_WritePin>
 80025c8:	40020000 	.word	0x40020000

080025cc <getch>:
 80025cc:	b513      	push	{r0, r1, r4, lr}
 80025ce:	4c07      	ldr	r4, [pc, #28]	@ (80025ec <getch+0x20>)
 80025d0:	f241 3388 	movw	r3, #5000	@ 0x1388
 80025d4:	2201      	movs	r2, #1
 80025d6:	f10d 0107 	add.w	r1, sp, #7
 80025da:	4620      	mov	r0, r4
 80025dc:	f000 fbe7 	bl	8002dae <HAL_UART_Receive>
 80025e0:	2800      	cmp	r0, #0
 80025e2:	d1f5      	bne.n	80025d0 <getch+0x4>
 80025e4:	f89d 0007 	ldrb.w	r0, [sp, #7]
 80025e8:	b002      	add	sp, #8
 80025ea:	bd10      	pop	{r4, pc}
 80025ec:	200068e8 	.word	0x200068e8

080025f0 <putch>:
 80025f0:	b507      	push	{r0, r1, r2, lr}
 80025f2:	f241 3388 	movw	r3, #5000	@ 0x1388
 80025f6:	f88d 0007 	strb.w	r0, [sp, #7]
 80025fa:	2201      	movs	r2, #1
 80025fc:	f10d 0107 	add.w	r1, sp, #7
 8002600:	4802      	ldr	r0, [pc, #8]	@ (800260c <putch+0x1c>)
 8002602:	f000 fb8f 	bl	8002d24 <HAL_UART_Transmit>
 8002606:	b003      	add	sp, #12
 8002608:	f85d fb04 	ldr.w	pc, [sp], #4
 800260c:	200068e8 	.word	0x200068e8

08002610 <HAL_GetTick>:
 8002610:	4b02      	ldr	r3, [pc, #8]	@ (800261c <HAL_GetTick+0xc>)
 8002612:	6818      	ldr	r0, [r3, #0]
 8002614:	1c42      	adds	r2, r0, #1
 8002616:	601a      	str	r2, [r3, #0]
 8002618:	4770      	bx	lr
 800261a:	bf00      	nop
 800261c:	20006978 	.word	0x20006978

08002620 <HAL_RCC_OscConfig>:
 8002620:	6803      	ldr	r3, [r0, #0]
 8002622:	b573      	push	{r0, r1, r4, r5, r6, lr}
 8002624:	07de      	lsls	r6, r3, #31
 8002626:	4601      	mov	r1, r0
 8002628:	d43b      	bmi.n	80026a2 <HAL_RCC_OscConfig+0x82>
 800262a:	680b      	ldr	r3, [r1, #0]
 800262c:	079d      	lsls	r5, r3, #30
 800262e:	f100 8088 	bmi.w	8002742 <HAL_RCC_OscConfig+0x122>
 8002632:	680b      	ldr	r3, [r1, #0]
 8002634:	0718      	lsls	r0, r3, #28
 8002636:	f100 80d3 	bmi.w	80027e0 <HAL_RCC_OscConfig+0x1c0>
 800263a:	680b      	ldr	r3, [r1, #0]
 800263c:	075a      	lsls	r2, r3, #29
 800263e:	d52a      	bpl.n	8002696 <HAL_RCC_OscConfig+0x76>
 8002640:	2300      	movs	r3, #0
 8002642:	9301      	str	r3, [sp, #4]
 8002644:	4b90      	ldr	r3, [pc, #576]	@ (8002888 <HAL_RCC_OscConfig+0x268>)
 8002646:	4c91      	ldr	r4, [pc, #580]	@ (800288c <HAL_RCC_OscConfig+0x26c>)
 8002648:	6c1a      	ldr	r2, [r3, #64]	@ 0x40
 800264a:	f042 5280 	orr.w	r2, r2, #268435456	@ 0x10000000
 800264e:	641a      	str	r2, [r3, #64]	@ 0x40
 8002650:	6c1b      	ldr	r3, [r3, #64]	@ 0x40
 8002652:	f003 5380 	and.w	r3, r3, #268435456	@ 0x10000000
 8002656:	9301      	str	r3, [sp, #4]
 8002658:	9b01      	ldr	r3, [sp, #4]
 800265a:	6823      	ldr	r3, [r4, #0]
 800265c:	f443 7380 	orr.w	r3, r3, #256	@ 0x100
 8002660:	6023      	str	r3, [r4, #0]
 8002662:	f7ff ffd5 	bl	8002610 <HAL_GetTick>
 8002666:	4605      	mov	r5, r0
 8002668:	6823      	ldr	r3, [r4, #0]
 800266a:	05d8      	lsls	r0, r3, #23
 800266c:	f140 80dc 	bpl.w	8002828 <HAL_RCC_OscConfig+0x208>
 8002670:	688b      	ldr	r3, [r1, #8]
 8002672:	4c85      	ldr	r4, [pc, #532]	@ (8002888 <HAL_RCC_OscConfig+0x268>)
 8002674:	2b01      	cmp	r3, #1
 8002676:	f040 80de 	bne.w	8002836 <HAL_RCC_OscConfig+0x216>
 800267a:	6f23      	ldr	r3, [r4, #112]	@ 0x70
 800267c:	f043 0301 	orr.w	r3, r3, #1
 8002680:	6723      	str	r3, [r4, #112]	@ 0x70
 8002682:	f7ff ffc5 	bl	8002610 <HAL_GetTick>
 8002686:	4d80      	ldr	r5, [pc, #512]	@ (8002888 <HAL_RCC_OscConfig+0x268>)
 8002688:	4604      	mov	r4, r0
 800268a:	f241 3688 	movw	r6, #5000	@ 0x1388
 800268e:	6f2b      	ldr	r3, [r5, #112]	@ 0x70
 8002690:	079a      	lsls	r2, r3, #30
 8002692:	f140 80f1 	bpl.w	8002878 <HAL_RCC_OscConfig+0x258>
 8002696:	698a      	ldr	r2, [r1, #24]
 8002698:	2a00      	cmp	r2, #0
 800269a:	f040 80fd 	bne.w	8002898 <HAL_RCC_OscConfig+0x278>
 800269e:	2000      	movs	r0, #0
 80026a0:	e015      	b.n	80026ce <HAL_RCC_OscConfig+0xae>
 80026a2:	4b79      	ldr	r3, [pc, #484]	@ (8002888 <HAL_RCC_OscConfig+0x268>)
 80026a4:	689a      	ldr	r2, [r3, #8]
 80026a6:	f002 020c 	and.w	r2, r2, #12
 80026aa:	2a04      	cmp	r2, #4
 80026ac:	d007      	beq.n	80026be <HAL_RCC_OscConfig+0x9e>
 80026ae:	689a      	ldr	r2, [r3, #8]
 80026b0:	f002 020c 	and.w	r2, r2, #12
 80026b4:	2a08      	cmp	r2, #8
 80026b6:	d10c      	bne.n	80026d2 <HAL_RCC_OscConfig+0xb2>
 80026b8:	685b      	ldr	r3, [r3, #4]
 80026ba:	025c      	lsls	r4, r3, #9
 80026bc:	d509      	bpl.n	80026d2 <HAL_RCC_OscConfig+0xb2>
 80026be:	4b72      	ldr	r3, [pc, #456]	@ (8002888 <HAL_RCC_OscConfig+0x268>)
 80026c0:	681b      	ldr	r3, [r3, #0]
 80026c2:	0398      	lsls	r0, r3, #14
 80026c4:	d5b1      	bpl.n	800262a <HAL_RCC_OscConfig+0xa>
 80026c6:	684b      	ldr	r3, [r1, #4]
 80026c8:	2b00      	cmp	r3, #0
 80026ca:	d1ae      	bne.n	800262a <HAL_RCC_OscConfig+0xa>
 80026cc:	2001      	movs	r0, #1
 80026ce:	b002      	add	sp, #8
 80026d0:	bd70      	pop	{r4, r5, r6, pc}
 80026d2:	684b      	ldr	r3, [r1, #4]
 80026d4:	4c6c      	ldr	r4, [pc, #432]	@ (8002888 <HAL_RCC_OscConfig+0x268>)
 80026d6:	f5b3 3f80 	cmp.w	r3, #65536	@ 0x10000
 80026da:	d112      	bne.n	8002702 <HAL_RCC_OscConfig+0xe2>
 80026dc:	6823      	ldr	r3, [r4, #0]
 80026de:	f443 3380 	orr.w	r3, r3, #65536	@ 0x10000
 80026e2:	6023      	str	r3, [r4, #0]
 80026e4:	f7ff ff94 	bl	8002610 <HAL_GetTick>
 80026e8:	4d67      	ldr	r5, [pc, #412]	@ (8002888 <HAL_RCC_OscConfig+0x268>)
 80026ea:	4e69      	ldr	r6, [pc, #420]	@ (8002890 <HAL_RCC_OscConfig+0x270>)
 80026ec:	4604      	mov	r4, r0
 80026ee:	682b      	ldr	r3, [r5, #0]
 80026f0:	039a      	lsls	r2, r3, #14
 80026f2:	d49a      	bmi.n	800262a <HAL_RCC_OscConfig+0xa>
 80026f4:	f7ff ff8c 	bl	8002610 <HAL_GetTick>
 80026f8:	1b00      	subs	r0, r0, r4
 80026fa:	42b0      	cmp	r0, r6
 80026fc:	d9f7      	bls.n	80026ee <HAL_RCC_OscConfig+0xce>
 80026fe:	2003      	movs	r0, #3
 8002700:	e7e5      	b.n	80026ce <HAL_RCC_OscConfig+0xae>
 8002702:	f5b3 2fa0 	cmp.w	r3, #327680	@ 0x50000
 8002706:	d104      	bne.n	8002712 <HAL_RCC_OscConfig+0xf2>
 8002708:	6823      	ldr	r3, [r4, #0]
 800270a:	f443 2380 	orr.w	r3, r3, #262144	@ 0x40000
 800270e:	6023      	str	r3, [r4, #0]
 8002710:	e7e4      	b.n	80026dc <HAL_RCC_OscConfig+0xbc>
 8002712:	6822      	ldr	r2, [r4, #0]
 8002714:	f422 3280 	bic.w	r2, r2, #65536	@ 0x10000
 8002718:	6022      	str	r2, [r4, #0]
 800271a:	6822      	ldr	r2, [r4, #0]
 800271c:	f422 2280 	bic.w	r2, r2, #262144	@ 0x40000
 8002720:	6022      	str	r2, [r4, #0]
 8002722:	2b00      	cmp	r3, #0
 8002724:	d1de      	bne.n	80026e4 <HAL_RCC_OscConfig+0xc4>
 8002726:	f7ff ff73 	bl	8002610 <HAL_GetTick>
 800272a:	4e59      	ldr	r6, [pc, #356]	@ (8002890 <HAL_RCC_OscConfig+0x270>)
 800272c:	4605      	mov	r5, r0
 800272e:	6823      	ldr	r3, [r4, #0]
 8002730:	039b      	lsls	r3, r3, #14
 8002732:	f57f af7a 	bpl.w	800262a <HAL_RCC_OscConfig+0xa>
 8002736:	f7ff ff6b 	bl	8002610 <HAL_GetTick>
 800273a:	1b40      	subs	r0, r0, r5
 800273c:	42b0      	cmp	r0, r6
 800273e:	d9f6      	bls.n	800272e <HAL_RCC_OscConfig+0x10e>
 8002740:	e7dd      	b.n	80026fe <HAL_RCC_OscConfig+0xde>
 8002742:	4b51      	ldr	r3, [pc, #324]	@ (8002888 <HAL_RCC_OscConfig+0x268>)
 8002744:	689a      	ldr	r2, [r3, #8]
 8002746:	f012 0f0c 	tst.w	r2, #12
 800274a:	d007      	beq.n	800275c <HAL_RCC_OscConfig+0x13c>
 800274c:	689a      	ldr	r2, [r3, #8]
 800274e:	f002 020c 	and.w	r2, r2, #12
 8002752:	2a08      	cmp	r2, #8
 8002754:	d116      	bne.n	8002784 <HAL_RCC_OscConfig+0x164>
 8002756:	685b      	ldr	r3, [r3, #4]
 8002758:	0258      	lsls	r0, r3, #9
 800275a:	d413      	bmi.n	8002784 <HAL_RCC_OscConfig+0x164>
 800275c:	484a      	ldr	r0, [pc, #296]	@ (8002888 <HAL_RCC_OscConfig+0x268>)
 800275e:	6803      	ldr	r3, [r0, #0]
 8002760:	079a      	lsls	r2, r3, #30
 8002762:	d502      	bpl.n	800276a <HAL_RCC_OscConfig+0x14a>
 8002764:	68cb      	ldr	r3, [r1, #12]
 8002766:	2b01      	cmp	r3, #1
 8002768:	d1b0      	bne.n	80026cc <HAL_RCC_OscConfig+0xac>
 800276a:	6804      	ldr	r4, [r0, #0]
 800276c:	22f8      	movs	r2, #248	@ 0xf8
 800276e:	fa92 f2a2 	rbit	r2, r2
 8002772:	690b      	ldr	r3, [r1, #16]
 8002774:	fab2 f282 	clz	r2, r2
 8002778:	4093      	lsls	r3, r2
 800277a:	f024 02f8 	bic.w	r2, r4, #248	@ 0xf8
 800277e:	4313      	orrs	r3, r2
 8002780:	6003      	str	r3, [r0, #0]
 8002782:	e756      	b.n	8002632 <HAL_RCC_OscConfig+0x12>
 8002784:	68ca      	ldr	r2, [r1, #12]
 8002786:	4b43      	ldr	r3, [pc, #268]	@ (8002894 <HAL_RCC_OscConfig+0x274>)
 8002788:	b1da      	cbz	r2, 80027c2 <HAL_RCC_OscConfig+0x1a2>
 800278a:	2201      	movs	r2, #1
 800278c:	601a      	str	r2, [r3, #0]
 800278e:	f7ff ff3f 	bl	8002610 <HAL_GetTick>
 8002792:	4c3d      	ldr	r4, [pc, #244]	@ (8002888 <HAL_RCC_OscConfig+0x268>)
 8002794:	4605      	mov	r5, r0
 8002796:	6823      	ldr	r3, [r4, #0]
 8002798:	079b      	lsls	r3, r3, #30
 800279a:	d50c      	bpl.n	80027b6 <HAL_RCC_OscConfig+0x196>
 800279c:	6820      	ldr	r0, [r4, #0]
 800279e:	22f8      	movs	r2, #248	@ 0xf8
 80027a0:	fa92 f2a2 	rbit	r2, r2
 80027a4:	690b      	ldr	r3, [r1, #16]
 80027a6:	fab2 f282 	clz	r2, r2
 80027aa:	4093      	lsls	r3, r2
 80027ac:	f020 02f8 	bic.w	r2, r0, #248	@ 0xf8
 80027b0:	4313      	orrs	r3, r2
 80027b2:	6023      	str	r3, [r4, #0]
 80027b4:	e73d      	b.n	8002632 <HAL_RCC_OscConfig+0x12>
 80027b6:	f7ff ff2b 	bl	8002610 <HAL_GetTick>
 80027ba:	1b40      	subs	r0, r0, r5
 80027bc:	2802      	cmp	r0, #2
 80027be:	d9ea      	bls.n	8002796 <HAL_RCC_OscConfig+0x176>
 80027c0:	e79d      	b.n	80026fe <HAL_RCC_OscConfig+0xde>
 80027c2:	601a      	str	r2, [r3, #0]
 80027c4:	f7ff ff24 	bl	8002610 <HAL_GetTick>
 80027c8:	4d2f      	ldr	r5, [pc, #188]	@ (8002888 <HAL_RCC_OscConfig+0x268>)
 80027ca:	4604      	mov	r4, r0
 80027cc:	682b      	ldr	r3, [r5, #0]
 80027ce:	079e      	lsls	r6, r3, #30
 80027d0:	f57f af2f 	bpl.w	8002632 <HAL_RCC_OscConfig+0x12>
 80027d4:	f7ff ff1c 	bl	8002610 <HAL_GetTick>
 80027d8:	1b00      	subs	r0, r0, r4
 80027da:	2802      	cmp	r0, #2
 80027dc:	d9f6      	bls.n	80027cc <HAL_RCC_OscConfig+0x1ac>
 80027de:	e78e      	b.n	80026fe <HAL_RCC_OscConfig+0xde>
 80027e0:	694a      	ldr	r2, [r1, #20]
 80027e2:	4b2c      	ldr	r3, [pc, #176]	@ (8002894 <HAL_RCC_OscConfig+0x274>)
 80027e4:	b182      	cbz	r2, 8002808 <HAL_RCC_OscConfig+0x1e8>
 80027e6:	2201      	movs	r2, #1
 80027e8:	f8c3 2e80 	str.w	r2, [r3, #3712]	@ 0xe80
 80027ec:	f7ff ff10 	bl	8002610 <HAL_GetTick>
 80027f0:	4d25      	ldr	r5, [pc, #148]	@ (8002888 <HAL_RCC_OscConfig+0x268>)
 80027f2:	4604      	mov	r4, r0
 80027f4:	6f6b      	ldr	r3, [r5, #116]	@ 0x74
 80027f6:	079b      	lsls	r3, r3, #30
 80027f8:	f53f af1f 	bmi.w	800263a <HAL_RCC_OscConfig+0x1a>
 80027fc:	f7ff ff08 	bl	8002610 <HAL_GetTick>
 8002800:	1b00      	subs	r0, r0, r4
 8002802:	2802      	cmp	r0, #2
 8002804:	d9f6      	bls.n	80027f4 <HAL_RCC_OscConfig+0x1d4>
 8002806:	e77a      	b.n	80026fe <HAL_RCC_OscConfig+0xde>
 8002808:	f8c3 2e80 	str.w	r2, [r3, #3712]	@ 0xe80
 800280c:	f7ff ff00 	bl	8002610 <HAL_GetTick>
 8002810:	4d1d      	ldr	r5, [pc, #116]	@ (8002888 <HAL_RCC_OscConfig+0x268>)
 8002812:	4604      	mov	r4, r0
 8002814:	6f6b      	ldr	r3, [r5, #116]	@ 0x74
 8002816:	079e      	lsls	r6, r3, #30
 8002818:	f57f af0f 	bpl.w	800263a <HAL_RCC_OscConfig+0x1a>
 800281c:	f7ff fef8 	bl	8002610 <HAL_GetTick>
 8002820:	1b00      	subs	r0, r0, r4
 8002822:	2802      	cmp	r0, #2
 8002824:	d9f6      	bls.n	8002814 <HAL_RCC_OscConfig+0x1f4>
 8002826:	e76a      	b.n	80026fe <HAL_RCC_OscConfig+0xde>
 8002828:	f7ff fef2 	bl	8002610 <HAL_GetTick>
 800282c:	1b40      	subs	r0, r0, r5
 800282e:	2802      	cmp	r0, #2
 8002830:	f67f af1a 	bls.w	8002668 <HAL_RCC_OscConfig+0x48>
 8002834:	e763      	b.n	80026fe <HAL_RCC_OscConfig+0xde>
 8002836:	2b05      	cmp	r3, #5
 8002838:	d104      	bne.n	8002844 <HAL_RCC_OscConfig+0x224>
 800283a:	6f23      	ldr	r3, [r4, #112]	@ 0x70
 800283c:	f043 0304 	orr.w	r3, r3, #4
 8002840:	6723      	str	r3, [r4, #112]	@ 0x70
 8002842:	e71a      	b.n	800267a <HAL_RCC_OscConfig+0x5a>
 8002844:	6f22      	ldr	r2, [r4, #112]	@ 0x70
 8002846:	f022 0201 	bic.w	r2, r2, #1
 800284a:	6722      	str	r2, [r4, #112]	@ 0x70
 800284c:	6f22      	ldr	r2, [r4, #112]	@ 0x70
 800284e:	f022 0204 	bic.w	r2, r2, #4
 8002852:	6722      	str	r2, [r4, #112]	@ 0x70
 8002854:	2b00      	cmp	r3, #0
 8002856:	f47f af14 	bne.w	8002682 <HAL_RCC_OscConfig+0x62>
 800285a:	f7ff fed9 	bl	8002610 <HAL_GetTick>
 800285e:	f241 3688 	movw	r6, #5000	@ 0x1388
 8002862:	4605      	mov	r5, r0
 8002864:	6f23      	ldr	r3, [r4, #112]	@ 0x70
 8002866:	079b      	lsls	r3, r3, #30
 8002868:	f57f af15 	bpl.w	8002696 <HAL_RCC_OscConfig+0x76>
 800286c:	f7ff fed0 	bl	8002610 <HAL_GetTick>
 8002870:	1b40      	subs	r0, r0, r5
 8002872:	42b0      	cmp	r0, r6
 8002874:	d9f6      	bls.n	8002864 <HAL_RCC_OscConfig+0x244>
 8002876:	e742      	b.n	80026fe <HAL_RCC_OscConfig+0xde>
 8002878:	f7ff feca 	bl	8002610 <HAL_GetTick>
 800287c:	1b00      	subs	r0, r0, r4
 800287e:	42b0      	cmp	r0, r6
 8002880:	f67f af05 	bls.w	800268e <HAL_RCC_OscConfig+0x6e>
 8002884:	e73b      	b.n	80026fe <HAL_RCC_OscConfig+0xde>
 8002886:	bf00      	nop
 8002888:	40023800 	.word	0x40023800
 800288c:	40007000 	.word	0x40007000
 8002890:	05f5e100 	.word	0x05f5e100
 8002894:	42470000 	.word	0x42470000
 8002898:	4c2d      	ldr	r4, [pc, #180]	@ (8002950 <HAL_RCC_OscConfig+0x330>)
 800289a:	68a3      	ldr	r3, [r4, #8]
 800289c:	f003 030c 	and.w	r3, r3, #12
 80028a0:	2b08      	cmp	r3, #8
 80028a2:	f43f af13 	beq.w	80026cc <HAL_RCC_OscConfig+0xac>
 80028a6:	4b2b      	ldr	r3, [pc, #172]	@ (8002954 <HAL_RCC_OscConfig+0x334>)
 80028a8:	2a02      	cmp	r2, #2
 80028aa:	f04f 0200 	mov.w	r2, #0
 80028ae:	661a      	str	r2, [r3, #96]	@ 0x60
 80028b0:	d13f      	bne.n	8002932 <HAL_RCC_OscConfig+0x312>
 80028b2:	f7ff fead 	bl	8002610 <HAL_GetTick>
 80028b6:	4e28      	ldr	r6, [pc, #160]	@ (8002958 <HAL_RCC_OscConfig+0x338>)
 80028b8:	4605      	mov	r5, r0
 80028ba:	6823      	ldr	r3, [r4, #0]
 80028bc:	0198      	lsls	r0, r3, #6
 80028be:	d432      	bmi.n	8002926 <HAL_RCC_OscConfig+0x306>
 80028c0:	f647 72c0 	movw	r2, #32704	@ 0x7fc0
 80028c4:	fa92 f2a2 	rbit	r2, r2
 80028c8:	f44f 3540 	mov.w	r5, #196608	@ 0x30000
 80028cc:	fab2 f282 	clz	r2, r2
 80028d0:	fa95 f5a5 	rbit	r5, r5
 80028d4:	f04f 6070 	mov.w	r0, #251658240	@ 0xf000000
 80028d8:	fab5 f585 	clz	r5, r5
 80028dc:	fa90 f0a0 	rbit	r0, r0
 80028e0:	e9d1 3607 	ldrd	r3, r6, [r1, #28]
 80028e4:	4333      	orrs	r3, r6
 80028e6:	6a4e      	ldr	r6, [r1, #36]	@ 0x24
 80028e8:	4096      	lsls	r6, r2
 80028ea:	6a8a      	ldr	r2, [r1, #40]	@ 0x28
 80028ec:	0852      	lsrs	r2, r2, #1
 80028ee:	3a01      	subs	r2, #1
 80028f0:	40aa      	lsls	r2, r5
 80028f2:	4333      	orrs	r3, r6
 80028f4:	4313      	orrs	r3, r2
 80028f6:	6aca      	ldr	r2, [r1, #44]	@ 0x2c
 80028f8:	4d17      	ldr	r5, [pc, #92]	@ (8002958 <HAL_RCC_OscConfig+0x338>)
 80028fa:	fab0 f080 	clz	r0, r0
 80028fe:	4082      	lsls	r2, r0
 8002900:	4313      	orrs	r3, r2
 8002902:	6063      	str	r3, [r4, #4]
 8002904:	4b13      	ldr	r3, [pc, #76]	@ (8002954 <HAL_RCC_OscConfig+0x334>)
 8002906:	4c12      	ldr	r4, [pc, #72]	@ (8002950 <HAL_RCC_OscConfig+0x330>)
 8002908:	2201      	movs	r2, #1
 800290a:	661a      	str	r2, [r3, #96]	@ 0x60
 800290c:	f7ff fe80 	bl	8002610 <HAL_GetTick>
 8002910:	4601      	mov	r1, r0
 8002912:	6823      	ldr	r3, [r4, #0]
 8002914:	019a      	lsls	r2, r3, #6
 8002916:	f53f aec2 	bmi.w	800269e <HAL_RCC_OscConfig+0x7e>
 800291a:	f7ff fe79 	bl	8002610 <HAL_GetTick>
 800291e:	1a40      	subs	r0, r0, r1
 8002920:	42a8      	cmp	r0, r5
 8002922:	d9f6      	bls.n	8002912 <HAL_RCC_OscConfig+0x2f2>
 8002924:	e6eb      	b.n	80026fe <HAL_RCC_OscConfig+0xde>
 8002926:	f7ff fe73 	bl	8002610 <HAL_GetTick>
 800292a:	1b40      	subs	r0, r0, r5
 800292c:	42b0      	cmp	r0, r6
 800292e:	d9c4      	bls.n	80028ba <HAL_RCC_OscConfig+0x29a>
 8002930:	e6e5      	b.n	80026fe <HAL_RCC_OscConfig+0xde>
 8002932:	f7ff fe6d 	bl	8002610 <HAL_GetTick>
 8002936:	4d08      	ldr	r5, [pc, #32]	@ (8002958 <HAL_RCC_OscConfig+0x338>)
 8002938:	4601      	mov	r1, r0
 800293a:	6823      	ldr	r3, [r4, #0]
 800293c:	019b      	lsls	r3, r3, #6
 800293e:	f57f aeae 	bpl.w	800269e <HAL_RCC_OscConfig+0x7e>
 8002942:	f7ff fe65 	bl	8002610 <HAL_GetTick>
 8002946:	1a40      	subs	r0, r0, r1
 8002948:	42a8      	cmp	r0, r5
 800294a:	d9f6      	bls.n	800293a <HAL_RCC_OscConfig+0x31a>
 800294c:	e6d7      	b.n	80026fe <HAL_RCC_OscConfig+0xde>
 800294e:	bf00      	nop
 8002950:	40023800 	.word	0x40023800
 8002954:	42470000 	.word	0x42470000
 8002958:	001e8480 	.word	0x001e8480

0800295c <HAL_RCC_ClockConfig>:
 800295c:	b5f8      	push	{r3, r4, r5, r6, r7, lr}
 800295e:	4b49      	ldr	r3, [pc, #292]	@ (8002a84 <HAL_RCC_ClockConfig+0x128>)
 8002960:	681a      	ldr	r2, [r3, #0]
 8002962:	f002 020f 	and.w	r2, r2, #15
 8002966:	428a      	cmp	r2, r1
 8002968:	4604      	mov	r4, r0
 800296a:	d311      	bcc.n	8002990 <HAL_RCC_ClockConfig+0x34>
 800296c:	6822      	ldr	r2, [r4, #0]
 800296e:	0795      	lsls	r5, r2, #30
 8002970:	d417      	bmi.n	80029a2 <HAL_RCC_ClockConfig+0x46>
 8002972:	07d0      	lsls	r0, r2, #31
 8002974:	d41d      	bmi.n	80029b2 <HAL_RCC_ClockConfig+0x56>
 8002976:	4b43      	ldr	r3, [pc, #268]	@ (8002a84 <HAL_RCC_ClockConfig+0x128>)
 8002978:	681a      	ldr	r2, [r3, #0]
 800297a:	f002 020f 	and.w	r2, r2, #15
 800297e:	428a      	cmp	r2, r1
 8002980:	d867      	bhi.n	8002a52 <HAL_RCC_ClockConfig+0xf6>
 8002982:	6822      	ldr	r2, [r4, #0]
 8002984:	0751      	lsls	r1, r2, #29
 8002986:	d46c      	bmi.n	8002a62 <HAL_RCC_ClockConfig+0x106>
 8002988:	0713      	lsls	r3, r2, #28
 800298a:	d472      	bmi.n	8002a72 <HAL_RCC_ClockConfig+0x116>
 800298c:	2000      	movs	r0, #0
 800298e:	e007      	b.n	80029a0 <HAL_RCC_ClockConfig+0x44>
 8002990:	b2ca      	uxtb	r2, r1
 8002992:	701a      	strb	r2, [r3, #0]
 8002994:	681b      	ldr	r3, [r3, #0]
 8002996:	f003 030f 	and.w	r3, r3, #15
 800299a:	428b      	cmp	r3, r1
 800299c:	d0e6      	beq.n	800296c <HAL_RCC_ClockConfig+0x10>
 800299e:	2001      	movs	r0, #1
 80029a0:	bdf8      	pop	{r3, r4, r5, r6, r7, pc}
 80029a2:	4839      	ldr	r0, [pc, #228]	@ (8002a88 <HAL_RCC_ClockConfig+0x12c>)
 80029a4:	68a5      	ldr	r5, [r4, #8]
 80029a6:	6883      	ldr	r3, [r0, #8]
 80029a8:	f023 03f0 	bic.w	r3, r3, #240	@ 0xf0
 80029ac:	432b      	orrs	r3, r5
 80029ae:	6083      	str	r3, [r0, #8]
 80029b0:	e7df      	b.n	8002972 <HAL_RCC_ClockConfig+0x16>
 80029b2:	6862      	ldr	r2, [r4, #4]
 80029b4:	4b34      	ldr	r3, [pc, #208]	@ (8002a88 <HAL_RCC_ClockConfig+0x12c>)
 80029b6:	2a01      	cmp	r2, #1
 80029b8:	d11d      	bne.n	80029f6 <HAL_RCC_ClockConfig+0x9a>
 80029ba:	681b      	ldr	r3, [r3, #0]
 80029bc:	f413 3f00 	tst.w	r3, #131072	@ 0x20000
 80029c0:	d0ed      	beq.n	800299e <HAL_RCC_ClockConfig+0x42>
 80029c2:	4d31      	ldr	r5, [pc, #196]	@ (8002a88 <HAL_RCC_ClockConfig+0x12c>)
 80029c4:	68ab      	ldr	r3, [r5, #8]
 80029c6:	f023 0303 	bic.w	r3, r3, #3
 80029ca:	4313      	orrs	r3, r2
 80029cc:	60ab      	str	r3, [r5, #8]
 80029ce:	f7ff fe1f 	bl	8002610 <HAL_GetTick>
 80029d2:	6863      	ldr	r3, [r4, #4]
 80029d4:	2b01      	cmp	r3, #1
 80029d6:	4606      	mov	r6, r0
 80029d8:	f241 3788 	movw	r7, #5000	@ 0x1388
 80029dc:	d115      	bne.n	8002a0a <HAL_RCC_ClockConfig+0xae>
 80029de:	68ab      	ldr	r3, [r5, #8]
 80029e0:	f003 030c 	and.w	r3, r3, #12
 80029e4:	2b04      	cmp	r3, #4
 80029e6:	d0c6      	beq.n	8002976 <HAL_RCC_ClockConfig+0x1a>
 80029e8:	f7ff fe12 	bl	8002610 <HAL_GetTick>
 80029ec:	1b80      	subs	r0, r0, r6
 80029ee:	42b8      	cmp	r0, r7
 80029f0:	d9f5      	bls.n	80029de <HAL_RCC_ClockConfig+0x82>
 80029f2:	2003      	movs	r0, #3
 80029f4:	e7d4      	b.n	80029a0 <HAL_RCC_ClockConfig+0x44>
 80029f6:	1e90      	subs	r0, r2, #2
 80029f8:	2801      	cmp	r0, #1
 80029fa:	681b      	ldr	r3, [r3, #0]
 80029fc:	d802      	bhi.n	8002a04 <HAL_RCC_ClockConfig+0xa8>
 80029fe:	f013 7f00 	tst.w	r3, #33554432	@ 0x2000000
 8002a02:	e7dd      	b.n	80029c0 <HAL_RCC_ClockConfig+0x64>
 8002a04:	f013 0f02 	tst.w	r3, #2
 8002a08:	e7da      	b.n	80029c0 <HAL_RCC_ClockConfig+0x64>
 8002a0a:	2b02      	cmp	r3, #2
 8002a0c:	d10a      	bne.n	8002a24 <HAL_RCC_ClockConfig+0xc8>
 8002a0e:	68ab      	ldr	r3, [r5, #8]
 8002a10:	f003 030c 	and.w	r3, r3, #12
 8002a14:	2b08      	cmp	r3, #8
 8002a16:	d0ae      	beq.n	8002976 <HAL_RCC_ClockConfig+0x1a>
 8002a18:	f7ff fdfa 	bl	8002610 <HAL_GetTick>
 8002a1c:	1b80      	subs	r0, r0, r6
 8002a1e:	42b8      	cmp	r0, r7
 8002a20:	d9f5      	bls.n	8002a0e <HAL_RCC_ClockConfig+0xb2>
 8002a22:	e7e6      	b.n	80029f2 <HAL_RCC_ClockConfig+0x96>
 8002a24:	2b03      	cmp	r3, #3
 8002a26:	d10f      	bne.n	8002a48 <HAL_RCC_ClockConfig+0xec>
 8002a28:	68ab      	ldr	r3, [r5, #8]
 8002a2a:	f003 030c 	and.w	r3, r3, #12
 8002a2e:	2b0c      	cmp	r3, #12
 8002a30:	d0a1      	beq.n	8002976 <HAL_RCC_ClockConfig+0x1a>
 8002a32:	f7ff fded 	bl	8002610 <HAL_GetTick>
 8002a36:	1b80      	subs	r0, r0, r6
 8002a38:	42b8      	cmp	r0, r7
 8002a3a:	d9f5      	bls.n	8002a28 <HAL_RCC_ClockConfig+0xcc>
 8002a3c:	e7d9      	b.n	80029f2 <HAL_RCC_ClockConfig+0x96>
 8002a3e:	f7ff fde7 	bl	8002610 <HAL_GetTick>
 8002a42:	1b80      	subs	r0, r0, r6
 8002a44:	42b8      	cmp	r0, r7
 8002a46:	d8d4      	bhi.n	80029f2 <HAL_RCC_ClockConfig+0x96>
 8002a48:	68ab      	ldr	r3, [r5, #8]
 8002a4a:	f013 0f0c 	tst.w	r3, #12
 8002a4e:	d1f6      	bne.n	8002a3e <HAL_RCC_ClockConfig+0xe2>
 8002a50:	e791      	b.n	8002976 <HAL_RCC_ClockConfig+0x1a>
 8002a52:	b2ca      	uxtb	r2, r1
 8002a54:	701a      	strb	r2, [r3, #0]
 8002a56:	681b      	ldr	r3, [r3, #0]
 8002a58:	f003 030f 	and.w	r3, r3, #15
 8002a5c:	428b      	cmp	r3, r1
 8002a5e:	d19e      	bne.n	800299e <HAL_RCC_ClockConfig+0x42>
 8002a60:	e78f      	b.n	8002982 <HAL_RCC_ClockConfig+0x26>
 8002a62:	4909      	ldr	r1, [pc, #36]	@ (8002a88 <HAL_RCC_ClockConfig+0x12c>)
 8002a64:	68e0      	ldr	r0, [r4, #12]
 8002a66:	688b      	ldr	r3, [r1, #8]
 8002a68:	f423 53e0 	bic.w	r3, r3, #7168	@ 0x1c00
 8002a6c:	4303      	orrs	r3, r0
 8002a6e:	608b      	str	r3, [r1, #8]
 8002a70:	e78a      	b.n	8002988 <HAL_RCC_ClockConfig+0x2c>
 8002a72:	4a05      	ldr	r2, [pc, #20]	@ (8002a88 <HAL_RCC_ClockConfig+0x12c>)
 8002a74:	6921      	ldr	r1, [r4, #16]
 8002a76:	6893      	ldr	r3, [r2, #8]
 8002a78:	f423 4360 	bic.w	r3, r3, #57344	@ 0xe000
 8002a7c:	ea43 03c1 	orr.w	r3, r3, r1, lsl #3
 8002a80:	6093      	str	r3, [r2, #8]
 8002a82:	e783      	b.n	800298c <HAL_RCC_ClockConfig+0x30>
 8002a84:	40023c00 	.word	0x40023c00
 8002a88:	40023800 	.word	0x40023800

08002a8c <HAL_GPIO_Init>:
 8002a8c:	e92d 4ff7 	stmdb	sp!, {r0, r1, r2, r4, r5, r6, r7, r8, r9, sl, fp, lr}
 8002a90:	f8df 819c 	ldr.w	r8, [pc, #412]	@ 8002c30 <HAL_GPIO_Init+0x1a4>
 8002a94:	4a67      	ldr	r2, [pc, #412]	@ (8002c34 <HAL_GPIO_Init+0x1a8>)
 8002a96:	2300      	movs	r3, #0
 8002a98:	f04f 0901 	mov.w	r9, #1
 8002a9c:	680c      	ldr	r4, [r1, #0]
 8002a9e:	fa09 fa03 	lsl.w	sl, r9, r3
 8002aa2:	ea0a 0504 	and.w	r5, sl, r4
 8002aa6:	ea3a 0404 	bics.w	r4, sl, r4
 8002aaa:	f040 80ac 	bne.w	8002c06 <HAL_GPIO_Init+0x17a>
 8002aae:	684c      	ldr	r4, [r1, #4]
 8002ab0:	f024 0e10 	bic.w	lr, r4, #16
 8002ab4:	f1be 0f02 	cmp.w	lr, #2
 8002ab8:	d114      	bne.n	8002ae4 <HAL_GPIO_Init+0x58>
 8002aba:	ea4f 0cd3 	mov.w	ip, r3, lsr #3
 8002abe:	eb00 0c8c 	add.w	ip, r0, ip, lsl #2
 8002ac2:	f003 0b07 	and.w	fp, r3, #7
 8002ac6:	f8dc 6020 	ldr.w	r6, [ip, #32]
 8002aca:	ea4f 0b8b 	mov.w	fp, fp, lsl #2
 8002ace:	270f      	movs	r7, #15
 8002ad0:	fa07 f70b 	lsl.w	r7, r7, fp
 8002ad4:	ea26 0707 	bic.w	r7, r6, r7
 8002ad8:	690e      	ldr	r6, [r1, #16]
 8002ada:	fa06 f60b 	lsl.w	r6, r6, fp
 8002ade:	433e      	orrs	r6, r7
 8002ae0:	f8cc 6020 	str.w	r6, [ip, #32]
 8002ae4:	f8d0 b000 	ldr.w	fp, [r0]
 8002ae8:	ea4f 0c43 	mov.w	ip, r3, lsl #1
 8002aec:	2603      	movs	r6, #3
 8002aee:	fa06 f70c 	lsl.w	r7, r6, ip
 8002af2:	ea2b 0b07 	bic.w	fp, fp, r7
 8002af6:	43fe      	mvns	r6, r7
 8002af8:	f004 0703 	and.w	r7, r4, #3
 8002afc:	fa07 f70c 	lsl.w	r7, r7, ip
 8002b00:	f10e 3eff 	add.w	lr, lr, #4294967295	@ 0xffffffff
 8002b04:	ea47 070b 	orr.w	r7, r7, fp
 8002b08:	f1be 0f01 	cmp.w	lr, #1
 8002b0c:	6007      	str	r7, [r0, #0]
 8002b0e:	d810      	bhi.n	8002b32 <HAL_GPIO_Init+0xa6>
 8002b10:	6887      	ldr	r7, [r0, #8]
 8002b12:	ea06 0e07 	and.w	lr, r6, r7
 8002b16:	68cf      	ldr	r7, [r1, #12]
 8002b18:	fa07 f70c 	lsl.w	r7, r7, ip
 8002b1c:	ea47 070e 	orr.w	r7, r7, lr
 8002b20:	6087      	str	r7, [r0, #8]
 8002b22:	6847      	ldr	r7, [r0, #4]
 8002b24:	ea27 0e0a 	bic.w	lr, r7, sl
 8002b28:	0927      	lsrs	r7, r4, #4
 8002b2a:	409f      	lsls	r7, r3
 8002b2c:	ea47 070e 	orr.w	r7, r7, lr
 8002b30:	6047      	str	r7, [r0, #4]
 8002b32:	68c7      	ldr	r7, [r0, #12]
 8002b34:	4037      	ands	r7, r6
 8002b36:	688e      	ldr	r6, [r1, #8]
 8002b38:	fa06 f60c 	lsl.w	r6, r6, ip
 8002b3c:	433e      	orrs	r6, r7
 8002b3e:	60c6      	str	r6, [r0, #12]
 8002b40:	00e6      	lsls	r6, r4, #3
 8002b42:	d560      	bpl.n	8002c06 <HAL_GPIO_Init+0x17a>
 8002b44:	2600      	movs	r6, #0
 8002b46:	9601      	str	r6, [sp, #4]
 8002b48:	f8d8 6044 	ldr.w	r6, [r8, #68]	@ 0x44
 8002b4c:	f446 4680 	orr.w	r6, r6, #16384	@ 0x4000
 8002b50:	f8c8 6044 	str.w	r6, [r8, #68]	@ 0x44
 8002b54:	f8d8 6044 	ldr.w	r6, [r8, #68]	@ 0x44
 8002b58:	f023 0703 	bic.w	r7, r3, #3
 8002b5c:	f107 4780 	add.w	r7, r7, #1073741824	@ 0x40000000
 8002b60:	f406 4680 	and.w	r6, r6, #16384	@ 0x4000
 8002b64:	f507 379c 	add.w	r7, r7, #79872	@ 0x13800
 8002b68:	9601      	str	r6, [sp, #4]
 8002b6a:	f003 0c03 	and.w	ip, r3, #3
 8002b6e:	9e01      	ldr	r6, [sp, #4]
 8002b70:	f8d7 e008 	ldr.w	lr, [r7, #8]
 8002b74:	ea4f 0c8c 	mov.w	ip, ip, lsl #2
 8002b78:	260f      	movs	r6, #15
 8002b7a:	fa06 f60c 	lsl.w	r6, r6, ip
 8002b7e:	ea2e 0e06 	bic.w	lr, lr, r6
 8002b82:	4e2d      	ldr	r6, [pc, #180]	@ (8002c38 <HAL_GPIO_Init+0x1ac>)
 8002b84:	42b0      	cmp	r0, r6
 8002b86:	d045      	beq.n	8002c14 <HAL_GPIO_Init+0x188>
 8002b88:	f506 6680 	add.w	r6, r6, #1024	@ 0x400
 8002b8c:	42b0      	cmp	r0, r6
 8002b8e:	d043      	beq.n	8002c18 <HAL_GPIO_Init+0x18c>
 8002b90:	f506 6680 	add.w	r6, r6, #1024	@ 0x400
 8002b94:	42b0      	cmp	r0, r6
 8002b96:	d041      	beq.n	8002c1c <HAL_GPIO_Init+0x190>
 8002b98:	f506 6680 	add.w	r6, r6, #1024	@ 0x400
 8002b9c:	42b0      	cmp	r0, r6
 8002b9e:	d03f      	beq.n	8002c20 <HAL_GPIO_Init+0x194>
 8002ba0:	f506 6680 	add.w	r6, r6, #1024	@ 0x400
 8002ba4:	42b0      	cmp	r0, r6
 8002ba6:	d03d      	beq.n	8002c24 <HAL_GPIO_Init+0x198>
 8002ba8:	f506 6680 	add.w	r6, r6, #1024	@ 0x400
 8002bac:	42b0      	cmp	r0, r6
 8002bae:	d03b      	beq.n	8002c28 <HAL_GPIO_Init+0x19c>
 8002bb0:	f506 6680 	add.w	r6, r6, #1024	@ 0x400
 8002bb4:	42b0      	cmp	r0, r6
 8002bb6:	d039      	beq.n	8002c2c <HAL_GPIO_Init+0x1a0>
 8002bb8:	f506 6680 	add.w	r6, r6, #1024	@ 0x400
 8002bbc:	42b0      	cmp	r0, r6
 8002bbe:	bf14      	ite	ne
 8002bc0:	2608      	movne	r6, #8
 8002bc2:	2607      	moveq	r6, #7
 8002bc4:	fa06 f60c 	lsl.w	r6, r6, ip
 8002bc8:	ea46 060e 	orr.w	r6, r6, lr
 8002bcc:	60be      	str	r6, [r7, #8]
 8002bce:	6816      	ldr	r6, [r2, #0]
 8002bd0:	43ef      	mvns	r7, r5
 8002bd2:	f414 3f80 	tst.w	r4, #65536	@ 0x10000
 8002bd6:	bf0c      	ite	eq
 8002bd8:	403e      	andeq	r6, r7
 8002bda:	432e      	orrne	r6, r5
 8002bdc:	6016      	str	r6, [r2, #0]
 8002bde:	6856      	ldr	r6, [r2, #4]
 8002be0:	f414 3f00 	tst.w	r4, #131072	@ 0x20000
 8002be4:	bf0c      	ite	eq
 8002be6:	403e      	andeq	r6, r7
 8002be8:	432e      	orrne	r6, r5
 8002bea:	6056      	str	r6, [r2, #4]
 8002bec:	6896      	ldr	r6, [r2, #8]
 8002bee:	f414 1f80 	tst.w	r4, #1048576	@ 0x100000
 8002bf2:	bf0c      	ite	eq
 8002bf4:	403e      	andeq	r6, r7
 8002bf6:	432e      	orrne	r6, r5
 8002bf8:	6096      	str	r6, [r2, #8]
 8002bfa:	68d6      	ldr	r6, [r2, #12]
 8002bfc:	02a4      	lsls	r4, r4, #10
 8002bfe:	bf54      	ite	pl
 8002c00:	403e      	andpl	r6, r7
 8002c02:	432e      	orrmi	r6, r5
 8002c04:	60d6      	str	r6, [r2, #12]
 8002c06:	3301      	adds	r3, #1
 8002c08:	2b10      	cmp	r3, #16
 8002c0a:	f47f af47 	bne.w	8002a9c <HAL_GPIO_Init+0x10>
 8002c0e:	b003      	add	sp, #12
 8002c10:	e8bd 8ff0 	ldmia.w	sp!, {r4, r5, r6, r7, r8, r9, sl, fp, pc}
 8002c14:	2600      	movs	r6, #0
 8002c16:	e7d5      	b.n	8002bc4 <HAL_GPIO_Init+0x138>
 8002c18:	2601      	movs	r6, #1
 8002c1a:	e7d3      	b.n	8002bc4 <HAL_GPIO_Init+0x138>
 8002c1c:	2602      	movs	r6, #2
 8002c1e:	e7d1      	b.n	8002bc4 <HAL_GPIO_Init+0x138>
 8002c20:	2603      	movs	r6, #3
 8002c22:	e7cf      	b.n	8002bc4 <HAL_GPIO_Init+0x138>
 8002c24:	2604      	movs	r6, #4
 8002c26:	e7cd      	b.n	8002bc4 <HAL_GPIO_Init+0x138>
 8002c28:	2605      	movs	r6, #5
 8002c2a:	e7cb      	b.n	8002bc4 <HAL_GPIO_Init+0x138>
 8002c2c:	2606      	movs	r6, #6
 8002c2e:	e7c9      	b.n	8002bc4 <HAL_GPIO_Init+0x138>
 8002c30:	40023800 	.word	0x40023800
 8002c34:	40013c00 	.word	0x40013c00
 8002c38:	40020000 	.word	0x40020000

08002c3c <HAL_GPIO_WritePin>:
 8002c3c:	b10a      	cbz	r2, 8002c42 <HAL_GPIO_WritePin+0x6>
 8002c3e:	6181      	str	r1, [r0, #24]
 8002c40:	4770      	bx	lr
 8002c42:	0409      	lsls	r1, r1, #16
 8002c44:	e7fb      	b.n	8002c3e <HAL_GPIO_WritePin+0x2>
	...

08002c48 <HAL_UART_Init>:
 8002c48:	b570      	push	{r4, r5, r6, lr}
 8002c4a:	4604      	mov	r4, r0
 8002c4c:	2800      	cmp	r0, #0
 8002c4e:	d065      	beq.n	8002d1c <HAL_UART_Init+0xd4>
 8002c50:	f890 3039 	ldrb.w	r3, [r0, #57]	@ 0x39
 8002c54:	f003 02ff 	and.w	r2, r3, #255	@ 0xff
 8002c58:	b90b      	cbnz	r3, 8002c5e <HAL_UART_Init+0x16>
 8002c5a:	f880 2038 	strb.w	r2, [r0, #56]	@ 0x38
 8002c5e:	6821      	ldr	r1, [r4, #0]
 8002c60:	68e2      	ldr	r2, [r4, #12]
 8002c62:	6925      	ldr	r5, [r4, #16]
 8002c64:	69e0      	ldr	r0, [r4, #28]
 8002c66:	2324      	movs	r3, #36	@ 0x24
 8002c68:	f884 3039 	strb.w	r3, [r4, #57]	@ 0x39
 8002c6c:	68cb      	ldr	r3, [r1, #12]
 8002c6e:	f423 5300 	bic.w	r3, r3, #8192	@ 0x2000
 8002c72:	60cb      	str	r3, [r1, #12]
 8002c74:	690b      	ldr	r3, [r1, #16]
 8002c76:	f423 5340 	bic.w	r3, r3, #12288	@ 0x3000
 8002c7a:	4313      	orrs	r3, r2
 8002c7c:	610b      	str	r3, [r1, #16]
 8002c7e:	68a3      	ldr	r3, [r4, #8]
 8002c80:	68ca      	ldr	r2, [r1, #12]
 8002c82:	432b      	orrs	r3, r5
 8002c84:	6965      	ldr	r5, [r4, #20]
 8002c86:	f422 4216 	bic.w	r2, r2, #38400	@ 0x9600
 8002c8a:	432b      	orrs	r3, r5
 8002c8c:	f022 020c 	bic.w	r2, r2, #12
 8002c90:	4303      	orrs	r3, r0
 8002c92:	4313      	orrs	r3, r2
 8002c94:	60cb      	str	r3, [r1, #12]
 8002c96:	694b      	ldr	r3, [r1, #20]
 8002c98:	69a2      	ldr	r2, [r4, #24]
 8002c9a:	f423 7340 	bic.w	r3, r3, #768	@ 0x300
 8002c9e:	4313      	orrs	r3, r2
 8002ca0:	f5b0 4f00 	cmp.w	r0, #32768	@ 0x8000
 8002ca4:	614b      	str	r3, [r1, #20]
 8002ca6:	4a1e      	ldr	r2, [pc, #120]	@ (8002d20 <HAL_UART_Init+0xd8>)
 8002ca8:	6863      	ldr	r3, [r4, #4]
 8002caa:	f04f 0564 	mov.w	r5, #100	@ 0x64
 8002cae:	d127      	bne.n	8002d00 <HAL_UART_Init+0xb8>
 8002cb0:	005b      	lsls	r3, r3, #1
 8002cb2:	fbb2 f2f3 	udiv	r2, r2, r3
 8002cb6:	fbb2 f6f5 	udiv	r6, r2, r5
 8002cba:	fb05 2216 	mls	r2, r5, r6, r2
 8002cbe:	00d3      	lsls	r3, r2, #3
 8002cc0:	3332      	adds	r3, #50	@ 0x32
 8002cc2:	fbb3 f3f5 	udiv	r3, r3, r5
 8002cc6:	f003 0207 	and.w	r2, r3, #7
 8002cca:	005b      	lsls	r3, r3, #1
 8002ccc:	eb02 1206 	add.w	r2, r2, r6, lsl #4
 8002cd0:	f403 73f8 	and.w	r3, r3, #496	@ 0x1f0
 8002cd4:	4413      	add	r3, r2
 8002cd6:	608b      	str	r3, [r1, #8]
 8002cd8:	690b      	ldr	r3, [r1, #16]
 8002cda:	f423 4390 	bic.w	r3, r3, #18432	@ 0x4800
 8002cde:	610b      	str	r3, [r1, #16]
 8002ce0:	694b      	ldr	r3, [r1, #20]
 8002ce2:	f023 032a 	bic.w	r3, r3, #42	@ 0x2a
 8002ce6:	614b      	str	r3, [r1, #20]
 8002ce8:	68cb      	ldr	r3, [r1, #12]
 8002cea:	f443 5300 	orr.w	r3, r3, #8192	@ 0x2000
 8002cee:	60cb      	str	r3, [r1, #12]
 8002cf0:	2000      	movs	r0, #0
 8002cf2:	2320      	movs	r3, #32
 8002cf4:	63e0      	str	r0, [r4, #60]	@ 0x3c
 8002cf6:	f884 3039 	strb.w	r3, [r4, #57]	@ 0x39
 8002cfa:	f884 303a 	strb.w	r3, [r4, #58]	@ 0x3a
 8002cfe:	bd70      	pop	{r4, r5, r6, pc}
 8002d00:	009b      	lsls	r3, r3, #2
 8002d02:	fbb2 f2f3 	udiv	r2, r2, r3
 8002d06:	fbb2 f6f5 	udiv	r6, r2, r5
 8002d0a:	fb05 2316 	mls	r3, r5, r6, r2
 8002d0e:	011b      	lsls	r3, r3, #4
 8002d10:	3332      	adds	r3, #50	@ 0x32
 8002d12:	fbb3 f3f5 	udiv	r3, r3, r5
 8002d16:	eb03 1306 	add.w	r3, r3, r6, lsl #4
 8002d1a:	e7dc      	b.n	8002cd6 <HAL_UART_Init+0x8e>
 8002d1c:	2001      	movs	r0, #1
 8002d1e:	e7ee      	b.n	8002cfe <HAL_UART_Init+0xb6>
 8002d20:	0afb6e90 	.word	0x0afb6e90

08002d24 <HAL_UART_Transmit>:
 8002d24:	4603      	mov	r3, r0
 8002d26:	f890 0039 	ldrb.w	r0, [r0, #57]	@ 0x39
 8002d2a:	2820      	cmp	r0, #32
 8002d2c:	d13d      	bne.n	8002daa <HAL_UART_Transmit+0x86>
 8002d2e:	2900      	cmp	r1, #0
 8002d30:	d039      	beq.n	8002da6 <HAL_UART_Transmit+0x82>
 8002d32:	2a00      	cmp	r2, #0
 8002d34:	d037      	beq.n	8002da6 <HAL_UART_Transmit+0x82>
 8002d36:	f893 0038 	ldrb.w	r0, [r3, #56]	@ 0x38
 8002d3a:	2801      	cmp	r0, #1
 8002d3c:	d035      	beq.n	8002daa <HAL_UART_Transmit+0x86>
 8002d3e:	2001      	movs	r0, #1
 8002d40:	f883 0038 	strb.w	r0, [r3, #56]	@ 0x38
 8002d44:	2000      	movs	r0, #0
 8002d46:	63d8      	str	r0, [r3, #60]	@ 0x3c
 8002d48:	2021      	movs	r0, #33	@ 0x21
 8002d4a:	f883 0039 	strb.w	r0, [r3, #57]	@ 0x39
 8002d4e:	849a      	strh	r2, [r3, #36]	@ 0x24
 8002d50:	84da      	strh	r2, [r3, #38]	@ 0x26
 8002d52:	681a      	ldr	r2, [r3, #0]
 8002d54:	8cd8      	ldrh	r0, [r3, #38]	@ 0x26
 8002d56:	b280      	uxth	r0, r0
 8002d58:	b948      	cbnz	r0, 8002d6e <HAL_UART_Transmit+0x4a>
 8002d5a:	6811      	ldr	r1, [r2, #0]
 8002d5c:	0649      	lsls	r1, r1, #25
 8002d5e:	d5fc      	bpl.n	8002d5a <HAL_UART_Transmit+0x36>
 8002d60:	2220      	movs	r2, #32
 8002d62:	2000      	movs	r0, #0
 8002d64:	f883 2039 	strb.w	r2, [r3, #57]	@ 0x39
 8002d68:	f883 0038 	strb.w	r0, [r3, #56]	@ 0x38
 8002d6c:	4770      	bx	lr
 8002d6e:	8cd8      	ldrh	r0, [r3, #38]	@ 0x26
 8002d70:	3801      	subs	r0, #1
 8002d72:	b280      	uxth	r0, r0
 8002d74:	84d8      	strh	r0, [r3, #38]	@ 0x26
 8002d76:	6898      	ldr	r0, [r3, #8]
 8002d78:	f5b0 5f80 	cmp.w	r0, #4096	@ 0x1000
 8002d7c:	d10c      	bne.n	8002d98 <HAL_UART_Transmit+0x74>
 8002d7e:	6810      	ldr	r0, [r2, #0]
 8002d80:	0600      	lsls	r0, r0, #24
 8002d82:	d5fc      	bpl.n	8002d7e <HAL_UART_Transmit+0x5a>
 8002d84:	8808      	ldrh	r0, [r1, #0]
 8002d86:	f3c0 0008 	ubfx	r0, r0, #0, #9
 8002d8a:	6050      	str	r0, [r2, #4]
 8002d8c:	6918      	ldr	r0, [r3, #16]
 8002d8e:	b908      	cbnz	r0, 8002d94 <HAL_UART_Transmit+0x70>
 8002d90:	3102      	adds	r1, #2
 8002d92:	e7df      	b.n	8002d54 <HAL_UART_Transmit+0x30>
 8002d94:	3101      	adds	r1, #1
 8002d96:	e7dd      	b.n	8002d54 <HAL_UART_Transmit+0x30>
 8002d98:	6810      	ldr	r0, [r2, #0]
 8002d9a:	0600      	lsls	r0, r0, #24
 8002d9c:	d5fc      	bpl.n	8002d98 <HAL_UART_Transmit+0x74>
 8002d9e:	f811 0b01 	ldrb.w	r0, [r1], #1
 8002da2:	6050      	str	r0, [r2, #4]
 8002da4:	e7d6      	b.n	8002d54 <HAL_UART_Transmit+0x30>
 8002da6:	2001      	movs	r0, #1
 8002da8:	4770      	bx	lr
 8002daa:	2002      	movs	r0, #2
 8002dac:	4770      	bx	lr

08002dae <HAL_UART_Receive>:
 8002dae:	4603      	mov	r3, r0
 8002db0:	f890 003a 	ldrb.w	r0, [r0, #58]	@ 0x3a
 8002db4:	2820      	cmp	r0, #32
 8002db6:	d141      	bne.n	8002e3c <HAL_UART_Receive+0x8e>
 8002db8:	2900      	cmp	r1, #0
 8002dba:	d03d      	beq.n	8002e38 <HAL_UART_Receive+0x8a>
 8002dbc:	2a00      	cmp	r2, #0
 8002dbe:	d03b      	beq.n	8002e38 <HAL_UART_Receive+0x8a>
 8002dc0:	f893 0038 	ldrb.w	r0, [r3, #56]	@ 0x38
 8002dc4:	2801      	cmp	r0, #1
 8002dc6:	d039      	beq.n	8002e3c <HAL_UART_Receive+0x8e>
 8002dc8:	2001      	movs	r0, #1
 8002dca:	f883 0038 	strb.w	r0, [r3, #56]	@ 0x38
 8002dce:	2000      	movs	r0, #0
 8002dd0:	63d8      	str	r0, [r3, #60]	@ 0x3c
 8002dd2:	2022      	movs	r0, #34	@ 0x22
 8002dd4:	f883 003a 	strb.w	r0, [r3, #58]	@ 0x3a
 8002dd8:	859a      	strh	r2, [r3, #44]	@ 0x2c
 8002dda:	85da      	strh	r2, [r3, #46]	@ 0x2e
 8002ddc:	8dd8      	ldrh	r0, [r3, #46]	@ 0x2e
 8002dde:	b280      	uxth	r0, r0
 8002de0:	b928      	cbnz	r0, 8002dee <HAL_UART_Receive+0x40>
 8002de2:	2220      	movs	r2, #32
 8002de4:	f883 203a 	strb.w	r2, [r3, #58]	@ 0x3a
 8002de8:	f883 0038 	strb.w	r0, [r3, #56]	@ 0x38
 8002dec:	4770      	bx	lr
 8002dee:	8dda      	ldrh	r2, [r3, #46]	@ 0x2e
 8002df0:	6898      	ldr	r0, [r3, #8]
 8002df2:	3a01      	subs	r2, #1
 8002df4:	b292      	uxth	r2, r2
 8002df6:	f5b0 5f80 	cmp.w	r0, #4096	@ 0x1000
 8002dfa:	85da      	strh	r2, [r3, #46]	@ 0x2e
 8002dfc:	681a      	ldr	r2, [r3, #0]
 8002dfe:	d10e      	bne.n	8002e1e <HAL_UART_Receive+0x70>
 8002e00:	6810      	ldr	r0, [r2, #0]
 8002e02:	0680      	lsls	r0, r0, #26
 8002e04:	d5fc      	bpl.n	8002e00 <HAL_UART_Receive+0x52>
 8002e06:	6918      	ldr	r0, [r3, #16]
 8002e08:	6852      	ldr	r2, [r2, #4]
 8002e0a:	b920      	cbnz	r0, 8002e16 <HAL_UART_Receive+0x68>
 8002e0c:	f3c2 0208 	ubfx	r2, r2, #0, #9
 8002e10:	f821 2b02 	strh.w	r2, [r1], #2
 8002e14:	e7e2      	b.n	8002ddc <HAL_UART_Receive+0x2e>
 8002e16:	b2d2      	uxtb	r2, r2
 8002e18:	f821 2b01 	strh.w	r2, [r1], #1
 8002e1c:	e7de      	b.n	8002ddc <HAL_UART_Receive+0x2e>
 8002e1e:	6810      	ldr	r0, [r2, #0]
 8002e20:	0680      	lsls	r0, r0, #26
 8002e22:	d5fc      	bpl.n	8002e1e <HAL_UART_Receive+0x70>
 8002e24:	6918      	ldr	r0, [r3, #16]
 8002e26:	6852      	ldr	r2, [r2, #4]
 8002e28:	b918      	cbnz	r0, 8002e32 <HAL_UART_Receive+0x84>
 8002e2a:	b2d2      	uxtb	r2, r2
 8002e2c:	f801 2b01 	strb.w	r2, [r1], #1
 8002e30:	e7d4      	b.n	8002ddc <HAL_UART_Receive+0x2e>
 8002e32:	f002 027f 	and.w	r2, r2, #127	@ 0x7f
 8002e36:	e7f9      	b.n	8002e2c <HAL_UART_Receive+0x7e>
 8002e38:	2001      	movs	r0, #1
 8002e3a:	4770      	bx	lr
 8002e3c:	2002      	movs	r0, #2
 8002e3e:	4770      	bx	lr

08002e40 <_sbrk>:
 8002e40:	4a0a      	ldr	r2, [pc, #40]	@ (8002e6c <_sbrk+0x2c>)
 8002e42:	6811      	ldr	r1, [r2, #0]
 8002e44:	b508      	push	{r3, lr}
 8002e46:	4603      	mov	r3, r0
 8002e48:	b909      	cbnz	r1, 8002e4e <_sbrk+0xe>
 8002e4a:	4909      	ldr	r1, [pc, #36]	@ (8002e70 <_sbrk+0x30>)
 8002e4c:	6011      	str	r1, [r2, #0]
 8002e4e:	6810      	ldr	r0, [r2, #0]
 8002e50:	4669      	mov	r1, sp
 8002e52:	4403      	add	r3, r0
 8002e54:	428b      	cmp	r3, r1
 8002e56:	d906      	bls.n	8002e66 <_sbrk+0x26>
 8002e58:	f7fd fa38 	bl	80002cc <__errno>
 8002e5c:	230c      	movs	r3, #12
 8002e5e:	6003      	str	r3, [r0, #0]
 8002e60:	f04f 30ff 	mov.w	r0, #4294967295	@ 0xffffffff
 8002e64:	bd08      	pop	{r3, pc}
 8002e66:	6013      	str	r3, [r2, #0]
 8002e68:	e7fc      	b.n	8002e64 <_sbrk+0x24>
 8002e6a:	bf00      	nop
 8002e6c:	2000697c 	.word	0x2000697c
 8002e70:	20006c88 	.word	0x20006c88

08002e74 <HAL_RNG_MspInit>:
 8002e74:	4770      	bx	lr

08002e76 <HAL_RNG_Init>:
 8002e76:	b510      	push	{r4, lr}
 8002e78:	4604      	mov	r4, r0
 8002e7a:	b1a8      	cbz	r0, 8002ea8 <HAL_RNG_Init+0x32>
 8002e7c:	7903      	ldrb	r3, [r0, #4]
 8002e7e:	2b01      	cmp	r3, #1
 8002e80:	d014      	beq.n	8002eac <HAL_RNG_Init+0x36>
 8002e82:	7943      	ldrb	r3, [r0, #5]
 8002e84:	f003 02ff 	and.w	r2, r3, #255	@ 0xff
 8002e88:	b913      	cbnz	r3, 8002e90 <HAL_RNG_Init+0x1a>
 8002e8a:	7102      	strb	r2, [r0, #4]
 8002e8c:	f7ff fff2 	bl	8002e74 <HAL_RNG_MspInit>
 8002e90:	6822      	ldr	r2, [r4, #0]
 8002e92:	2302      	movs	r3, #2
 8002e94:	7163      	strb	r3, [r4, #5]
 8002e96:	6813      	ldr	r3, [r2, #0]
 8002e98:	f043 0304 	orr.w	r3, r3, #4
 8002e9c:	6013      	str	r3, [r2, #0]
 8002e9e:	2000      	movs	r0, #0
 8002ea0:	2301      	movs	r3, #1
 8002ea2:	7163      	strb	r3, [r4, #5]
 8002ea4:	7120      	strb	r0, [r4, #4]
 8002ea6:	bd10      	pop	{r4, pc}
 8002ea8:	2001      	movs	r0, #1
 8002eaa:	e7fc      	b.n	8002ea6 <HAL_RNG_Init+0x30>
 8002eac:	2002      	movs	r0, #2
 8002eae:	e7fa      	b.n	8002ea6 <HAL_RNG_Init+0x30>

08002eb0 <Reset_Handler>:
 8002eb0:	2100      	movs	r1, #0
 8002eb2:	e003      	b.n	8002ebc <LoopCopyDataInit>

08002eb4 <CopyDataInit>:
 8002eb4:	4b0b      	ldr	r3, [pc, #44]	@ (8002ee4 <LoopForever+0x2>)
 8002eb6:	585b      	ldr	r3, [r3, r1]
 8002eb8:	5043      	str	r3, [r0, r1]
 8002eba:	3104      	adds	r1, #4

08002ebc <LoopCopyDataInit>:
 8002ebc:	480a      	ldr	r0, [pc, #40]	@ (8002ee8 <LoopForever+0x6>)
 8002ebe:	4b0b      	ldr	r3, [pc, #44]	@ (8002eec <LoopForever+0xa>)
 8002ec0:	1842      	adds	r2, r0, r1
 8002ec2:	429a      	cmp	r2, r3
 8002ec4:	d3f6      	bcc.n	8002eb4 <CopyDataInit>
 8002ec6:	4a0a      	ldr	r2, [pc, #40]	@ (8002ef0 <LoopForever+0xe>)
 8002ec8:	e002      	b.n	8002ed0 <LoopFillZerobss>

08002eca <FillZerobss>:
 8002eca:	2300      	movs	r3, #0
 8002ecc:	6013      	str	r3, [r2, #0]
 8002ece:	3204      	adds	r2, #4

08002ed0 <LoopFillZerobss>:
 8002ed0:	4b08      	ldr	r3, [pc, #32]	@ (8002ef4 <LoopForever+0x12>)
 8002ed2:	429a      	cmp	r2, r3
 8002ed4:	d3f9      	bcc.n	8002eca <FillZerobss>
 8002ed6:	f3af 8000 	nop.w
 8002eda:	f7fd fef1 	bl	8000cc0 <__libc_init_array>
 8002ede:	f7fe f947 	bl	8001170 <main>

08002ee2 <LoopForever>:
 8002ee2:	e7fe      	b.n	8002ee2 <LoopForever>
 8002ee4:	0802731c 	.word	0x0802731c
 8002ee8:	20000000 	.word	0x20000000
 8002eec:	20000538 	.word	0x20000538
 8002ef0:	20000538 	.word	0x20000538
 8002ef4:	20006c88 	.word	0x20006c88

08002ef8 <BusFault_Handler>:
 8002ef8:	e7fe      	b.n	8002ef8 <BusFault_Handler>
	...

08002efc <register_fini>:
 8002efc:	4b02      	ldr	r3, [pc, #8]	@ (8002f08 <register_fini+0xc>)
 8002efe:	b113      	cbz	r3, 8002f06 <register_fini+0xa>
 8002f00:	4802      	ldr	r0, [pc, #8]	@ (8002f0c <register_fini+0x10>)
 8002f02:	f7fd bf17 	b.w	8000d34 <atexit>
 8002f06:	4770      	bx	lr
 8002f08:	00000000 	.word	0x00000000
 8002f0c:	08000d65 	.word	0x08000d65

08002f10 <_init>:
 8002f10:	b5f8      	push	{r3, r4, r5, r6, r7, lr}
 8002f12:	bf00      	nop
 8002f14:	bcf8      	pop	{r3, r4, r5, r6, r7}
 8002f16:	bc08      	pop	{r3}
 8002f18:	469e      	mov	lr, r3
 8002f1a:	4770      	bx	lr

08002f1c <_fini>:
 8002f1c:	b5f8      	push	{r3, r4, r5, r6, r7, lr}
 8002f1e:	bf00      	nop
 8002f20:	bcf8      	pop	{r3, r4, r5, r6, r7}
 8002f22:	bc08      	pop	{r3}
 8002f24:	469e      	mov	lr, r3
 8002f26:	4770      	bx	lr
