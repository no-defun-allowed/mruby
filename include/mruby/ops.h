/* operand types:
   + Z: no operand
   + B: 8bit
   + BB: 8+8bit
   + BBB: 8+8+8bit
   + BS: 8+16bit
   + BSS: 8+16+16bit
   + NN: 4+4bit
   + S: 16bit
   + W: 24bit
*/

/*-----------------------------------------------------------------------
operation code    operands      semantics
------------------------------------------------------------------------*/
OPCODE(NOP,        Z)        /* no operation */
OPCODE(MOVE,       NN)       /* R[a] = R[b] */
OPCODE(LOADL,      NN)       /* R[a] = Pool[b] */
OPCODE(LOADI8,     NN)       /* R[a] = mrb_int(b) */
OPCODE(LOADINEG,   NN)       /* R[a] = mrb_int(-b) */
OPCODE(LOADI16,    BS)       /* R[a] = mrb_int(b) */
OPCODE(LOADI32,    BSS)      /* R[a] = mrb_int((b<<16)+c) */
OPCODE(LOADSYM,    NN)       /* R[a] = Syms[b] */
OPCODE(LOADNIL,    B)        /* R[a] = nil */
OPCODE(LOADTRUE,   B)        /* R[a] = true */
OPCODE(LOADFALSE,  B)        /* R[a] = false */
OPCODE(GETGV,      NN)       /* R[a] = getglobal(Syms[b]) */
OPCODE(SETGV,      NN)       /* setglobal(Syms[b], R[a]) */
OPCODE(GETSV,      NN)       /* R[a] = Special[Syms[b]] */
OPCODE(SETSV,      NN)       /* Special[Syms[b]] = R[a] */
OPCODE(GETIV,      NN)       /* R[a] = ivget(Syms[b]) */
OPCODE(SETIV,      NN)       /* ivset(Syms[b],R[a]) */
OPCODE(GETCV,      NN)       /* R[a] = cvget(Syms[b]) */
OPCODE(SETCV,      NN)       /* cvset(Syms[b],R[a]) */
OPCODE(GETCONST,   NN)       /* R[a] = constget(Syms[b]) */
OPCODE(SETCONST,   NN)       /* constset(Syms[b],R[a]) */
OPCODE(GETMCNST,   NN)       /* R[a] = R[a]::Syms[b] */
OPCODE(SETMCNST,   NN)       /* R[a+1]::Syms[b] = R[a] */
OPCODE(GETUPVAR,   BBB)      /* R[a] = uvget(b,c) */
OPCODE(SETUPVAR,   BBB)      /* uvset(b,c,R[a]) */
OPCODE(GETIDX,     B)        /* R[a] = R[a][R[a+1]] */
OPCODE(SETIDX,     B)        /* R[a][R[a+1]] = R[a+2] */
OPCODE(JMP,        S)        /* pc+=a */
OPCODE(JMPIF,      BS)       /* if R[a] pc+=b */
OPCODE(JMPNOT,     BS)       /* if !R[a] pc+=b */
OPCODE(JMPNIL,     BS)       /* if R[a]==nil pc+=b */
OPCODE(JMPUW,      S)        /* unwind_and_jump_to(a) */
OPCODE(EXCEPT,     B)        /* R[a] = exc */
OPCODE(RESCUE,     BB)       /* R[b] = R[a].isa?(R[b]) */
OPCODE(RAISEIF,    B)        /* raise(R[a]) if R[a] */
OPCODE(MATCHERR,   B)        /* raise NoMatchingPatternError unless R[a] */
OPCODE(SSEND,      BBB)      /* R[a] = self.send(Syms[b],R[a+1]..,R[a+n+1]:R[a+n+2]..) (c=n|k<<4) */
OPCODE(SSEND0,     NN)       /* R[a] = self.send(Syms[b]) (no args) */
OPCODE(SSENDB,     BBB)      /* R[a] = self.send(Syms[b],R[a+1]..,R[a+n+1]:R[a+n+2]..,&R[a+n+2k+1]) */
OPCODE(SEND,       BBB)      /* R[a] = R[a].send(Syms[b],R[a+1]..,R[a+n+1]:R[a+n+2]..) (c=n|k<<4) */
OPCODE(SEND0,      NN)       /* R[a] = R[a].send(Syms[b]) (no args) */
OPCODE(SENDB,      BBB)      /* R[a] = R[a].send(Syms[b],R[a+1]..,R[a+n+1]:R[a+n+2]..,&R[a+n+2k+1]) */
OPCODE(CALL,       Z)        /* self.call(*, **, &) (But overlay the current call frame; tailcall) */
OPCODE(BLKCALL,    BB)       /* R[a] = R[a].call(R[a+1],... ,R[a+b]); direct block call */
OPCODE(SUPER,      BB)       /* R[a] = super(R[a+1],... ,R[a+b+1]) */
OPCODE(ARGARY,     BS)       /* R[a] = argument array (16=m5:r1:m5:d1:lv4) */
OPCODE(ENTER,      W)        /* arg setup according to flags (24=n1:m5:o5:r1:m5:k5:d1:b1) */
OPCODE(KEY_P,      NN)       /* R[a] = kdict.key?(Syms[b]) */
OPCODE(KEYEND,     Z)        /* raise unless kdict.empty? */
OPCODE(KARG,       NN)       /* R[a] = kdict[Syms[b]]; kdict.delete(Syms[b]) */
OPCODE(RETURN,     B)        /* return R[a] (normal) */
OPCODE(RETURN_BLK, B)        /* return R[a] (in-block return) */
OPCODE(BREAK,      B)        /* break R[a] */
OPCODE(BLKPUSH,    BS)       /* R[a] = block (16=m5:r1:m5:d1:lv4) */
OPCODE(ADD,        B)        /* R[a] = R[a]+R[a+1] */
OPCODE(ADDI,       BB)       /* R[a] = R[a]+mrb_int(b) */
OPCODE(SUB,        B)        /* R[a] = R[a]-R[a+1] */
OPCODE(SUBI,       BB)       /* R[a] = R[a]-mrb_int(b) */
OPCODE(ADDILV,     BBB)      /* R[a] = R[a]+mrb_int(c); R[b],R[b+1] for method call */
OPCODE(SUBILV,     BBB)      /* R[a] = R[a]-mrb_int(c); R[b],R[b+1] for method call */
OPCODE(MUL,        B)        /* R[a] = R[a]*R[a+1] */
OPCODE(DIV,        B)        /* R[a] = R[a]/R[a+1] */
OPCODE(EQ,         B)        /* R[a] = R[a]==R[a+1] */
OPCODE(LT,         B)        /* R[a] = R[a]<R[a+1] */
OPCODE(LE,         B)        /* R[a] = R[a]<=R[a+1] */
OPCODE(GT,         B)        /* R[a] = R[a]>R[a+1] */
OPCODE(GE,         B)        /* R[a] = R[a]>=R[a+1] */
OPCODE(ARRAY,      NN)       /* R[a] = ary_new(R[a],R[a+1]..R[a+b]) */
OPCODE(ARRAY2,     BBB)      /* R[a] = ary_new(R[b],R[b+1]..R[b+c]) */
OPCODE(ARYCAT,     B)        /* ary_cat(R[a],R[a+1]) */
OPCODE(ARYPUSH,    BB)       /* ary_push(R[a],R[a+1]..R[a+b]) */
OPCODE(ARYSPLAT,   B)        /* R[a] = ary_splat(R[a]) */
OPCODE(AREF,       BBB)      /* R[a] = R[b][c] */
OPCODE(ASET,       BBB)      /* R[b][c] = R[a] */
OPCODE(APOST,      BBB)      /* *R[a],R[a+1]..R[a+c] = R[a][b..] */
OPCODE(INTERN,     B)        /* R[a] = intern(R[a]) */
OPCODE(SYMBOL,     BB)       /* R[a] = intern(Pool[b]) */
OPCODE(STRING,     NN)       /* R[a] = str_dup(Pool[b]) */
OPCODE(STRCAT,     B)        /* str_cat(R[a],R[a+1]) */
OPCODE(HASH,       BB)       /* R[a] = hash_new(R[a],R[a+1]..R[a+b*2-1]) */
OPCODE(HASHADD,    BB)       /* hash_push(R[a],R[a+1]..R[a+b*2]) */
OPCODE(HASHCAT,    B)        /* R[a] = hash_cat(R[a],R[a+1]) */
OPCODE(LAMBDA,     BB)       /* R[a] = lambda(Irep[b],L_LAMBDA) */
OPCODE(BLOCK,      NN)       /* R[a] = lambda(Irep[b],L_BLOCK) */
OPCODE(METHOD,     BB)       /* R[a] = lambda(Irep[b],L_METHOD) */
OPCODE(RANGE_INC,  B)        /* R[a] = range_new(R[a],R[a+1],FALSE) */
OPCODE(RANGE_EXC,  B)        /* R[a] = range_new(R[a],R[a+1],TRUE) */
OPCODE(OCLASS,     B)        /* R[a] = ::Object */
OPCODE(CLASS,      BB)       /* R[a] = newclass(R[a],Syms[b],R[a+1]) */
OPCODE(MODULE,     BB)       /* R[a] = newmodule(R[a],Syms[b]) */
OPCODE(EXEC,       BB)       /* R[a] = blockexec(R[a],Irep[b]) */
OPCODE(DEF,        BB)       /* R[a].newmethod(Syms[b],R[a+1]); R[a] = Syms[b] */
OPCODE(TDEF,       BBB)      /* target_class.newmethod(Syms[b],Irep[c]); R[a] = Syms[b] */
OPCODE(SDEF,       BBB)      /* R[a].singleton_class.newmethod(Syms[b],Irep[c]); R[a] = Syms[b] */
OPCODE(ALIAS,      BB)       /* alias_method(target_class,Syms[a],Syms[b]) */
OPCODE(UNDEF,      B)        /* undef_method(target_class,Syms[a]) */
OPCODE(SCLASS,     B)        /* R[a] = R[a].singleton_class */
OPCODE(TCLASS,     B)        /* R[a] = target_class */
OPCODE(DEBUG,      BBB)      /* print a,b,c */
OPCODE(ERR,        B)        /* raise(LocalJumpError, Pool[a]) */
OPCODE(EXT1,       Z)        /* make operands at least 8-bit */
OPCODE(EXT2,       Z)        /* make operands 16-bit */
OPCODE(STOP,       Z)        /* stop VM */
OPCODE(ENTER0,     Z)        /* zero-argument entry */
OPCODE(ENTER1,     Z)        /* one-argument entry */
OPCODE(ENTER2,     Z)        /* two-argument entry */
OPCODE(RETURN0,    Z)        /* return R[0] */
OPCODE(RETURN1,    Z)        /* ... */
OPCODE(RETURN2,    Z)
OPCODE(RETURN3,    Z)
OPCODE(RETURN4,    Z)
OPCODE(RETURN5,    Z)
OPCODE(RETURN6,    Z)
OPCODE(RETURN7,    Z)           /* return R[7] */
OPCODE(SSEND1,     NN)          /* R[a] = self.send(Syms[b], R[a+1]) */
OPCODE(SEND1,      NN)          /* R[a] = R[a].send(Syms[b], R[a+1]) */
