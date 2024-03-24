.LBB1_4:
        mov     eax, dword ptr [rip + n]
        mov     dword ptr [rbp - 12], eax
.LBB1_5:                                # =>This Loop Header: Depth=1
        cmp     dword ptr [rbp - 12], 1
        jle     .LBB1_12
        mov     dword ptr [rbp - 16], 0
.LBB1_7:                                #   Parent Loop BB1_5 Depth=1
        mov     eax, dword ptr [rbp - 16]
        mov     dword ptr [rbp - 20], eax       # 4-byte Spill
        mov     eax, dword ptr [rbp - 12]
        mov     ecx, 2
        cdq
        idiv    ecx
        mov     ecx, eax
        mov     eax, dword ptr [rbp - 20]       # 4-byte Reload
        cmp     eax, ecx
        jge     .LBB1_10
        mov     eax, dword ptr [rbp - 16]
        shl     eax, 1
        movsxd  rcx, eax
        lea     rax, [rip + a]
        mov     rdx, qword ptr [rax + 8*rcx]
        mov     eax, dword ptr [rbp - 16]
        shl     eax, 1
        add     eax, 1
        movsxd  rcx, eax
        lea     rax, [rip + a]
        add     rdx, qword ptr [rax + 8*rcx]
        movsxd  rcx, dword ptr [rbp - 16]
        lea     rax, [rip + a]
        mov     qword ptr [rax + 8*rcx], rdx
        mov     eax, dword ptr [rbp - 16]
        add     eax, 1
        mov     dword ptr [rbp - 16], eax
        jmp     .LBB1_7
.LBB1_10:                               #   in Loop: Header=BB1_5 Depth=1
        jmp     .LBB1_11
.LBB1_11:                               #   in Loop: Header=BB1_5 Depth=1
        mov     eax, dword ptr [rbp - 12]
        mov     ecx, 2
        cdq
        idiv    ecx
        mov     dword ptr [rbp - 12], eax
        jmp     .LBB1_5
.LBB1_12:
        mov     rax, qword ptr [rip + a]
        mov     qword ptr [rip + sum], rax