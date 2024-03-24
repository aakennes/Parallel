.LBB1_4:
        mov     dword ptr [rbp - 12], 0
.LBB1_5:                                # =>This Inner Loop Header: Depth=1
        mov     eax, dword ptr [rbp - 12]
        cmp     eax, dword ptr [rip + n]
        jge     .LBB1_8
        movsxd  rcx, dword ptr [rbp - 12]
        lea     rax, [rip + a]
        mov     rax, qword ptr [rax + 8*rcx]
        add     rax, qword ptr [rip + sum1]
        mov     qword ptr [rip + sum1], rax
        mov     eax, dword ptr [rbp - 12]
        add     eax, 1
        movsxd  rcx, eax
        lea     rax, [rip + a]
        mov     rax, qword ptr [rax + 8*rcx]
        add     rax, qword ptr [rip + sum2]
        mov     qword ptr [rip + sum2], rax
        mov     eax, dword ptr [rbp - 12]
        add     eax, 2
        mov     dword ptr [rbp - 12], eax
        jmp     .LBB1_5
.LBB1_8:
        mov     rax, qword ptr [rip + sum1]
        add     rax, qword ptr [rip + sum2]
        mov     qword ptr [rip + sum], rax