// Barely-tested multi-cycle TACKY sample solution
// Written by Henry Dietz, http://aggregate.org/hankd

// basic sizes of things
`define WORD	[15:0]
`define	TAGWORD	[16:0]
`define	OP	[4:0]
`define	STATE	[4:0] // some state numbers are OPs
`define REGSIZE	[7:0]
`define MEMSIZE	[65535:0]

// field placements and values
`define	TAG	[16]    // type tag in registers
`define Op0	[15:11] // opcode field
`define	Reg0	[10:8]  // register number field
`define Op1	[7:3]   // second opcode
`define Reg1	[2:0]   // second register number, also size
`define	Imm8	[7:0]   // also used as size
`define	OpPack	[15:14]	// if not 2'b11, packed inst

// TACKY data type flag values
`define	TFLOAT	1'b1
`define	TINT	1'b0

// opcode values, also state numbers
`define OPa2r	5'b00000
`define OPr2a	5'b00001
`define OPadd	5'b00010
`define OPand	5'b00011
`define OPcvt	5'b00100
`define OPdiv	5'b00101
`define OPmul	5'b00110
`define OPnot	5'b00111
`define OPor	5'b01000
`define OPsh	5'b01001
`define OPslt	5'b01010
`define OPsub	5'b01011
`define OPxor	5'b01100

`define OPjr	5'b10000
`define OPlf	5'b10001
`define OPli	5'b10010
`define OPst	5'b10011

`define OPcf8	5'b11000
`define OPci8	5'b11001
`define OPjnz8	5'b11010
`define OPjz8	5'b11011

`define OPjp8	5'b11100
`define OPpre	5'b11101
`define OPsys	5'b11110

// state numbers
`define Start	`OPa2r
`define Decode	`OPr2a
`define Pack0	`OPadd
`define Pack1	`OPand


// Floating point Verilog modules for CPE480
// Created February 19, 2019 by Henry Dietz, http://aggregate.org/hankd
// Distributed under CC BY 4.0, https://creativecommons.org/licenses/by/4.0/

// Field definitions
`define	WORD	[15:0]	// generic machine word size
`define	INT	signed [15:0]	// integer size
`define FLOAT	[15:0]	// half-precision float size
`define FSIGN	[15]	// sign bit
`define FEXP	[14:7]	// exponent
`define FFRAC	[6:0]	// fractional part (leading 1 implied)

// Constants
`define	FZERO	16'b0	  // float 0
`define F32767  16'h46ff  // closest approx to 32767, actually 32640
`define F32768  16'hc700  // -32768

// Count leading zeros, 16-bit (5-bit result) d=lead0s(s)
module lead0s(d, s);
output wire [4:0] d;
input wire `WORD s;
wire [4:0] t;
wire [7:0] s8;
wire [3:0] s4;
wire [1:0] s2;
assign t[4] = 0;
assign {t[3],s8} = ((|s[15:8]) ? {1'b0,s[15:8]} : {1'b1,s[7:0]});
assign {t[2],s4} = ((|s8[7:4]) ? {1'b0,s8[7:4]} : {1'b1,s8[3:0]});
assign {t[1],s2} = ((|s4[3:2]) ? {1'b0,s4[3:2]} : {1'b1,s4[1:0]});
assign t[0] = !s2[1];
assign d = (s ? t : 16);
endmodule

// Float set-less-than, 16-bit (1-bit result) torf=a<b
module fslt(torf, a, b);
output wire torf;
input wire `FLOAT a, b;
assign torf = (a `FSIGN && !(b `FSIGN)) ||
	      (a `FSIGN && b `FSIGN && (a[14:0] > b[14:0])) ||
	      (!(a `FSIGN) && !(b `FSIGN) && (a[14:0] < b[14:0]));
endmodule

// Floating-point addition, 16-bit r=a+b
module fadd(r, a, b);
output wire `FLOAT r;
input wire `FLOAT a, b;
wire `FLOAT s;
wire [8:0] sexp, sman, sfrac;
wire [7:0] texp, taman, tbman;
wire [4:0] slead;
wire ssign, aegt, amgt, eqsgn;
assign r = ((a == 0) ? b : ((b == 0) ? a : s));
assign aegt = (a `FEXP > b `FEXP);
assign texp = (aegt ? (a `FEXP) : (b `FEXP));
assign taman = (aegt ? {1'b1, (a `FFRAC)} : ({1'b1, (a `FFRAC)} >> (texp - a `FEXP)));
assign tbman = (aegt ? ({1'b1, (b `FFRAC)} >> (texp - b `FEXP)) : {1'b1, (b `FFRAC)});
assign eqsgn = (a `FSIGN == b `FSIGN);
assign amgt = (taman > tbman);
assign sman = (eqsgn ? (taman + tbman) : (amgt ? (taman - tbman) : (tbman - taman)));
lead0s m0(slead, {sman, 7'b0});
assign ssign = (amgt ? (a `FSIGN) : (b `FSIGN));
assign sfrac = sman << slead;
assign sexp = (texp + 1) - slead;
assign s = (sman ? (sexp ? {ssign, sexp[7:0], sfrac[7:1]} : 0) : 0);
endmodule

// Floating-point multiply, 16-bit r=a*b
module fmul(r, a, b);
output wire `FLOAT r;
input wire `FLOAT a, b;
wire [15:0] m; // double the bits in a fraction, we need high bits
wire [7:0] e;
wire s;
assign s = (a `FSIGN ^ b `FSIGN);
assign m = ({1'b1, (a `FFRAC)} * {1'b1, (b `FFRAC)});
assign e = (((a `FEXP) + (b `FEXP)) -127 + m[15]);
assign r = (((a == 0) || (b == 0)) ? 0 : (m[15] ? {s, e, m[14:8]} : {s, e, m[13:7]}));
endmodule

// Floating-point reciprocal, 16-bit r=1.0/a
// Note: requires initialized inverse fraction lookup table
module frecip(r, a);
output wire `FLOAT r;
input wire `FLOAT a;
reg [6:0] look[127:0];
initial $readmemh0(look);
assign r `FSIGN = a `FSIGN;
assign r `FEXP = 253 + (!(a `FFRAC)) - a `FEXP;
assign r `FFRAC = look[a `FFRAC];
endmodule

// Floating-point shift, 16 bit
// Shift +left,-right by integer
module fshift(r, f, i);
output wire `FLOAT r;
input wire `FLOAT f;
input wire `INT i;
assign r `FFRAC = f `FFRAC;
assign r `FSIGN = f `FSIGN;
assign r `FEXP = (f ? (f `FEXP + i) : 0);
endmodule

// Integer to float conversion, 16 bit
module i2f(f, i);
output wire `FLOAT f;
input wire `INT i;
wire [4:0] lead;
wire `WORD pos;
assign pos = (i[15] ? (-i) : i);
lead0s m0(lead, pos);
assign f `FFRAC = (i ? ({pos, 8'b0} >> (16 - lead)) : 0);
assign f `FSIGN = i[15];
assign f `FEXP = (i ? (128 + (14 - lead)) : 0);
endmodule

// Float to integer conversion, 16 bit
// Note: out-of-range values go to -32768 or 32767
module f2i(i, f);
output wire `INT i;
input wire `FLOAT f;
wire `FLOAT ui;
wire tiny, big;
fslt m0(tiny, f, `F32768);
fslt m1(big, `F32767, f);
assign ui = {1'b1, f `FFRAC, 16'b0} >> ((128+22) - f `FEXP);
assign i = (tiny ? 0 : (big ? 32767 : (f `FSIGN ? (-ui) : ui)));
endmodule

// End of float library


module alu(valid, o, op, a, b);
output valid;
output `TAGWORD o;
input `OP op;
input `TAGWORD a, b;
reg `TAGWORD t;
reg ok;
wire sltf;
wire `WORD cvi;
wire `FLOAT addf, cvf, recipf, mulf, shf;
wire signed `WORD sa, sb;

fadd myfadd(addf, a `WORD, ((op == `OPsub) ? (b `WORD ^ 16'h8000) : b `WORD));
i2f myi2f(cvf, b `WORD);
f2i myf2i(cvi, b `WORD);
frecip myfrecip(recipf, b `WORD);
fmul myfmul(mulf, a `WORD, ((op == `OPmul) ? b `WORD : recipf));
fshift myfshift(shf, a `WORD, b `WORD);
fslt myfslt(sltf, a `WORD, b `WORD);
assign sa = a `WORD; // signed version of a
assign sb = b `WORD; // signed version of b
assign o = t;
assign valid = ok;

always @* begin
  t = a;
  ok = 1;
  case ({a `TAG, op})
    {`TINT, `OPr2a},
    {`TFLOAT, `OPr2a}: t = b;
    {`TINT, `OPadd}:   t = {`TINT, (a `WORD + b `WORD)};
    {`TFLOAT, `OPadd}: t = {`TFLOAT,  addf};
    {`TINT, `OPand},
    {`TFLOAT, `OPand}: t `WORD = a `WORD & b `WORD;
    {`TINT, `OPcvt},
    {`TFLOAT, `OPcvt}: t = ((b `TAG == `TFLOAT) ? {`TINT, cvi} : {`TFLOAT, cvf});
    {`TINT, `OPdiv}:   t = {`TINT, (a `WORD / b `WORD)};
    {`TFLOAT, `OPdiv}: t = {`TFLOAT, mulf};
    {`TINT, `OPmul}:   t = {`TINT, (a `WORD * b `WORD)};
    {`TFLOAT, `OPmul}: t = {`TFLOAT, mulf};
    {`TINT, `OPnot},
    {`TFLOAT, `OPnot}: t `WORD = ~(b `WORD);
    {`TINT, `OPor},
    {`TFLOAT, `OPor}:  t `WORD = (a `WORD | b `WORD);
    {`TINT, `OPsh}:    t = {`TINT, ((sb < 0) ? (sa >> -sb) : (sa << sb))};
    {`TFLOAT, `OPsh}:  t = {`TFLOAT, shf};
    {`TINT, `OPslt}:   t = {`TINT, (sa < sb)};
    {`TFLOAT, `OPslt}: t = {`TINT, 15'b0, sltf};
    {`TINT, `OPsub}:   t `WORD = sa - sb;
    {`TFLOAT, `OPsub}: t = {`TFLOAT, addf};
    {`TINT, `OPxor},
    {`TFLOAT, `OPxor}: t `WORD = (a `WORD ^ b `WORD);
    default:           ok = 0;
  endcase
end
endmodule

module processor(halt, reset, clk);
output reg halt;
input reset, clk;
reg `TAGWORD r `REGSIZE;
wire `TAGWORD aluv;
reg `WORD text `MEMSIZE; // instruction memory
reg `WORD data `MEMSIZE; // data memory
reg `WORD pc = 0;
reg `WORD ir;
reg `Imm8 pre, sys;
reg `STATE s, op;
reg `Reg1 rn;
wire valid;

alu myalu(valid, aluv, op, r[s == `Pack1], r[rn]);

always @(posedge reset) begin
  halt <= 0;
  pc <= 0;
  s <= `Start;

  // initialize some stuff...
  // r[0] = 17'h00001; // 1
  // r[1] = 17'h13f80; // 1.0
  // r[2] = 17'h00002; // 2
  // r[3] = 17'h14000; // 2.0
  // text[0] = { `OPadd, 3'd2, `OPadd, 3'd3 };
  // r[0] = 17'h00003, r[1] = 17'h14040
  // text[1] = { `OPsub, 3'd2, `OPsub, 3'd3 };
  // text[2] = { `OPsys, 3'd0, 8'd0 };

  // the better way would be using readmem...
  readmemh1(text);
  readmemh2(data);
end

always @(posedge clk) begin
  case (s)
    `Start:  begin ir <= text[pc]; s <= `Decode; end
    `Decode: begin
               pc <= pc + 1;
	       op <= ir `Op0; rn <= ir `Reg0;
               s <= ((ir `OpPack == 2'b11) ? ir `Op0 : `Pack0);
	     end
    `OPcf8:  begin r[rn] <= {`TFLOAT, pre, ir `Imm8}; s <= `Start; end
    `OPci8:  begin r[rn] <= {`TINT, pre, ir `Imm8}; s <= `Start; end
    `OPjnz8: begin if (r[rn] `WORD) pc <= {pre, ir `Imm8}; s <= `Start; end
    `OPjz8:  begin if (!r[rn] `WORD) pc <= {pre, ir `Imm8}; s <= `Start; end
    `OPjp8:  begin pc <= {pre, ir `Imm8}; s <= `Start; end
    `OPpre:  begin pre <= ir `Imm8; s <= `Start; end
    `Pack0:  begin
               case (op)
                 `OPjr:   pc <= r[rn] `WORD;
		 `OPlf:   r[rn] <= {`TFLOAT, data[r[0] `WORD]};
		 `OPli:   r[rn] <= {`TINT, data[r[0] `WORD]};
		 `OPst:   data[r[rn]] <= r[0] `WORD;
		 `OPa2r:  r[rn] <= r[0];
		 default: if (valid) r[0] <= aluv;
               endcase
	       op <= ir `Op1; rn <= ir `Reg1;
	       s <= `Pack1;
	     end
    `Pack1:  begin
               case (ir `Op1)
                 `OPjr:   pc <= r[rn] `WORD;
		 `OPlf:   r[rn] <= {`TFLOAT, data[r[1] `WORD]};
		 `OPli:   r[rn] <= {`TINT, data[r[1] `WORD]};
		 `OPst:   data[r[rn]] <= r[1] `WORD;
		 `OPa2r:  r[rn] <= r[1];
		 default: if (valid) r[1] <= aluv;
               endcase
	       s <= `Start;
             end
    // default is OPsys and some illegal instructions
    default: begin sys <= ir `Imm8; halt <= 1; end
  endcase
end
endmodule


module testbench;
reg reset = 0;
reg clk = 0;
wire halted;
processor PE(halted, reset, clk);
initial begin
  $dumpfile;
  $dumpvars(1, PE.pc, PE.r[0], PE.r[1]); // would normally trace 0, PE
  #10 reset = 1;
  #10 reset = 0;
  while (!halted) begin
    #10 clk = 1;
    #10 clk = 0;
  end
  $finish;
end
endmodule
