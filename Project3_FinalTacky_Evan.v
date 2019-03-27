// basic sizes of things
`define WORD        [15:0]
`define OPcode1 [15:11]
`define OPcode2 [7:3]
`define REG1        [10:8]
`define REG2        [2:0]
`define IMM8        [7:0]
`define STATE       [5:0]
`define TYPEDREG    [16:0]
`define REGSIZE     [7:0]
`define MEMSIZE     [65535:0]




// opcode values, also state numbers
`define OPno        5'b00000
`define OPpre       5'b00001
`define OPjp8       5'b00011
`define OPsys       5'b00010
`define OPcf8       5'b00110
`define OPci8       5'b00111
`define OPjnz8      5'b00101
`define OPjz8       5'b00100
`define OPa2r       5'b01100
`define OPr2a       5'b01101
`define OPlf        5'b01111
`define OPli        5'b01110
`define OPst        5'b01010
`define OPcvt       5'b01011
`define OPjr        5'b01001
`define OPadd       5'b11000
`define OPsub       5'b11001
`define OPmul       5'b11011
`define OPdiv       5'b11010
`define OPand       5'b11110
`define OPor        5'b11111
`define OPxor       5'b11101
`define OPnot       5'b11100
`define OPsh        5'b10100
`define OPslt       5'b10101


//state numbers only
`define Start       5'b01000
`define Start1      5'b10000
`define Start2      5'b10001


//Module stuff,
`define ALU     		5'b1xxxx

// Floating point Verilog modules for CPE480
// Created February 19, 2019 by Henry Dietz, http://aggregate.org/hankd
// Distributed under CC BY 4.0, https://creativecommons.org/licenses/by/4.0/

// Field definitions
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
initial $readmemh3(look);
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

//ALU MODULE
module ALUmod(out, in1, in2, op, type);
input [4:0] op;
input type; //0->integer arithmetic, 1->floating point arithmetic
//$acc should always be in1 and $r should be in2.
input signed `WORD in1, in2;
output reg `TYPEDREG out;
wire signed `WORD recr, addr, subr, shr, divr, mulr, sltr;
wire signed `WORD outand, outor, outnot, outxor, outslt;
//Assign the bitwise operations.
assign outand = in1 & in2;
assign outor = in1 | in2;
assign outxor = in1 ^ in2;
assign outnot = ~in2;
assign outslt = in1 < in2;
//Instantiate the floating point modules.
fadd fa(addr, in1, in2);
//Use fadd with the negated version of in2.
//Negate by flipping the top bit of the float.
fadd fsu(subr, in1, {~in2[15], in2[14:0]});
fmul fm(mulr, in1, in2);
frecip fr(recr, in2);
fmul fd(divr, in1, recr);
fshift fs(shr, in1, in2);
fslt fsl(sltr, in1, in2);
always @(*) begin
 case(op)
   `OPadd: begin
           case(type)
           0: begin out <= in1 + in2; end
           1: begin out <= {type,addr}; end
           endcase
           end
   `OPsub: begin
           case(type)
           0: begin out <= in1 - in2; end
           1: begin out <= {type,subr}; end
           endcase
           end
   `OPmul: begin
           case(type)
           0: begin out <= in1 * in2; end
           1: begin out <= {type,mulr}; end
           endcase
           end
   `OPdiv: begin
           case(type)
           0: begin out <= in1 / in2; end
           1: begin out <= {type,divr}; end
           endcase
           end
   `OPand: begin out <= {type,outand}; end
   `OPor:  begin out <= {type,outor}; end
   `OPxor: begin out <= {type,outxor}; end
   `OPnot: begin out <= {type,outnot}; end
   //Positive indicates left shift.
   `OPsh:  begin
           case(type)
           0:  begin out <= in1 << in2; end
           1:  begin out <= {type,shr}; end
           endcase
           end
   `OPslt: begin
           case(type)
           0:  begin out = outslt; end
           1:  begin out <= {type,sltr}; end
           endcase
           end
   default: out <= 16'b0;
 endcase
end
endmodule

//PROCSSOR MODULE
module processor(halt, reset, clk);
output reg halt;
input reset, clk;

reg `TYPEDREG regfile `REGSIZE;
reg `IMM8 pre;
reg `WORD mainmem `MEMSIZE;
reg `WORD pc = 0;
reg `WORD ir;
reg `STATE s = `Start;
reg field;
integer a;
integer testcount;
wire `TYPEDREG R0;
wire `TYPEDREG R1;
wire `TYPEDREG R2;
wire `TYPEDREG R3;
wire `OPcode1 op1;
wire `OPcode2 op2;
wire `REG1 reg1;
wire `REG2 reg2;
wire `TYPEDREG memReg;
assign R0 = regfile[0];
assign R1 = regfile[1];
assign R2 = regfile[2];
assign R3 = regfile[3];
assign op1 = ir `OPcode1;
assign op2 = ir `OPcode2;
assign reg1 = ir `REG1;
assign reg2 = ir `REG2;
assign memReg = mainmem[16'h0FAB];
wire `TYPEDREG ALU1out, ALU2out;
wire `WORD cvt2int1, cvt2int2, cvt2float1, cvt2float2;
ALUmod ALU1(ALU1out, regfile[0]`WORD, regfile[ir `REG1]`WORD, ir `OPcode1, regfile[0][16]);
ALUmod ALU2(ALU2out, regfile[1]`WORD, regfile[ir `REG2]`WORD, ir `OPcode2, regfile[1][16]);
f2i cvti1(cvt2int1, regfile[ir `REG1]`WORD);
f2i cvti2(cvt2int2, regfile[ir `REG2]`WORD);
i2f cvtf1(cvt2float1, regfile[ir `REG1]`WORD);
i2f cvtf2(cvt2float2, regfile[ir `REG2]`WORD);

always @(reset) begin
 halt = 0;
 pc = 0;
 s = `Start;
 $readmemh0(regfile);
 $readmemh1(mainmem);
 pre <= 0;
end

always @(posedge clk)
begin
	case (s)
	`Start: begin
	        ir <= mainmem[pc];
	        s <= `Start1;
	         end

	`Start1:   begin
	                pc <= pc + 1;            	   // bump pc
	                s <= ir `OPcode1;   // most instructions, state # is opcode
	                field <= 0;
	                end
	`Start2:    begin
	                s <= ir `OPcode2;
	                field <= 1;
	            end
	`OPpre: 		begin
	                pre <= ir `IMM8;
	                s <= `Start;
	            end
	`OPsys: 		begin
	                halt <= 1;
	                s <= `OPno;
	            end
	`OPcf8: 		begin
	                regfile[ir `REG1] <= {1'b1, pre `IMM8, ir `IMM8};
	                s <= `Start;
	            end
	`OPci8: 		begin
	                regfile[ir `REG1] <= {1'b0, pre, ir `IMM8};
	                s <= `Start;
	            end
	`OPa2r: 		begin
	                regfile[ir `REG1] <= (field ? regfile[ir `REG1] : regfile[field]);
	                regfile[ir `REG2] <= (!field ? regfile[ir `REG2] : regfile[field]);
	                s <= ( field ? `Start :  `Start2);
	            end
	`OPr2a: 		begin
	                regfile[field] <= (field ? regfile[ir `REG2] : regfile[ir `REG1]);
	                s <= ( field ? `Start :  `Start2);
	            end
	`OPlf:      begin
	        				regfile[ir `REG1] <= (field ? regfile[ir `REG1]:{1'b1, mainmem[regfile[field]]`WORD});
	            		regfile[ir `REG2] <= (!field ? regfile[ir `REG2]:{1'b1, mainmem[regfile[field]]`WORD});
	        				s <= ( field ? `Start :  `Start2);
	        		end
	`OPli:      begin
	        				regfile[ir `REG1] <= (field ? regfile[ir `REG1] : {1'b0,mainmem[regfile[field]]`WORD});
	            		regfile[ir `REG2] <= (!field ? regfile[ir `REG2] : {1'b0,mainmem[regfile[field]]`WORD});
	        				s <= ( field ? `Start :  `Start2);
	        		end
	`OPst:      begin
	        		mainmem[regfile[ir `REG1]] <= (field ? regfile[ir `REG1] : regfile[field]);
	        		mainmem[regfile[ir `REG2]] <= (!field ? regfile[ir `REG2] : regfile[field]);
	        		s <= ( field ? `Start :  `Start2);
	        		end
	`OPjr:      begin
	            pc <= regfile[ir `REG1];
	            s <= `Start;
	            end
	`OPjp8: 		begin
	            pc <= {pre,ir `IMM8};
	            s <= `Start;
	            end
	`OPjz8: 		begin
	            pc <= (regfile[ir `REG1]`WORD ? pc : {pre,ir `IMM8});
	            s <= `Start;
	            end
	`OPjnz8:    begin
	            pc <= (regfile[ir `REG1]`WORD ? {pre,ir `IMM8} : pc);
	            s <= `Start;
	            end
	`OPcvt: 		begin
	            regfile[field] <= (field ? (regfile[ir`REG1][16] ? cvt2int1 : cvt2float1) :  (regfile[ir`REG2][16] ? cvt2int2 : cvt2float2));
	            end
	`OPadd:			 begin
	        		regfile[field] <= (field ? ALU2out : ALU1out);
	        		s <= (field ? `Start : `Start2);
	        		end
	`OPsub:			 begin
	        		regfile[field] <= (field ? ALU2out : ALU1out);
	        		s <= (field ? `Start : `Start2);
	        		end
	`OPmul:			 begin
	        		regfile[field] <= (field ? ALU2out : ALU1out);
	        		s <= (field ? `Start : `Start2);
	        		end
	`OPdiv:			 begin
	        		regfile[field] <= (field ? ALU2out : ALU1out);
	        		s <= (field ? `Start : `Start2);
	        		end
	`OPand:			 begin
	        		regfile[field] <= (field ? ALU2out : ALU1out);
	        		s <= (field ? `Start : `Start2);
	        		end
	`OPor :			 begin
	        		regfile[field] <= (field ? ALU2out : ALU1out);
	        		s <= (field ? `Start : `Start2);
	        		end
	`OPxor:			 begin
	        		regfile[field] <= (field ? ALU2out : ALU1out);
	        		s <= (field ? `Start : `Start2);
	        		end
	`OPnot:			 begin
	        		regfile[field] <= (field ? ALU2out : ALU1out);
	        		s <= (field ? `Start : `Start2);
	        		end
	`OPsh :			 begin
	        		regfile[field] <= (field ? ALU2out : ALU1out);
	        		s <= (field ? `Start : `Start2);
	        		end
	`OPslt:     begin
	        		regfile[field] <= (field ? ALU2out : ALU1out);
	        		s <= (field ? `Start : `Start2);
	        		end
	        		default: halt <= 1;
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
  $dumpvars(1, PE);
  #10 reset = 1;
  #10 reset = 0;
  while (!halted) begin
    #10 clk = 1;
    #10 clk = 0;
  end
  $finish;
end
endmodule
