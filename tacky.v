// basic sizes of things
`define WORD        [15:0]
`define HALFWORD    [7:0]
`define OPcode1     [15:11]
`define OPcode2     [7:3]
`define REG1        [10:8]
`define REG2        [2:0]
`define IMM8        [7:0]
`define STATE       [5:0]
`define TYPEDREG    [16:0]
`define REGSIZE     [7:0]
`define MEMSIZE     [65535:0]
`define REGNUM      [2:0]
`define OPcodeID    [4:0]

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
`define OPcvt       5'b10111

//state numbers only
`define Start       5'b01000
`define Start1      5'b10000
`define Start2      5'b10001

//Module stuff
`define ALU     	5'b1xxxx

//Accumulator values
`define Acc0        [0]
`define Acc1        [1]

//Type values
`define Int         1'b0
`define Float       1'b1

//Boolean values
`define true        1'b1
`define false       1'b0

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

//Identifies which registers and accumulators are read from, if any
module RegistersReadFrom(Field1_ACC0, Field1_REG1, Field2_ACC1, Field2_REG2, RR_JumpFlag, RR_SysFlag, RR_inst, clk);
	output reg [3:0] Field1_ACC0;
	output reg [3:0] Field1_REG1;
	output reg [3:0] Field2_ACC1;
	output reg [3:0] Field2_REG2;
	output reg RR_JumpFlag;
	output reg RR_SysFlag;
	
	input clk;
	input `WORD RR_inst;
	
	always @ (posedge clk) begin
		if (RR_inst `OPcode1 == `OPsys) begin
			RR_SysFlag <= `true;
			RR_JumpFlag <= `false;
		end
		else if ((RR_inst `OPcode1 == `OPjz8) || (RR_inst `OPcode1 == `OPjnz8) || (RR_inst `OPcode1 == `OPjr) || (RR_inst `OPcode1 == `OPjp8) || (RR_inst `OPcode2 == `OPjr)) begin  
			RR_JumpFlag <= `true;
			RR_SysFlag <= `false;
		end
		else begin
			RR_JumpFlag <= `false;
			RR_SysFlag <= `false;
			//if the first opcode is a register reader, assign the appropriate register to Field1_REG1
			if ((RR_inst `OPcode1 == `OPst) || (RR_inst `OPcode1 == `OPcvt) || (RR_inst `OPcode1 == `OPr2a) || (RR_inst `OPcode1 == `OPsh) || (RR_inst `OPcode1 == `OPslt) || (RR_inst `OPcode1 == `OPadd) || (RR_inst `OPcode1 == `OPsub) || (RR_inst `OPcode1 == `OPdiv) || (RR_inst `OPcode1 == `OPmul) || (RR_inst `OPcode1 == `OPnot) || (RR_inst `OPcode1 == `OPxor) || (RR_inst `OPcode1 == `OPand) || (RR_inst `OPcode1 == `OPor)) begin
				
				Field1_REG1 <= RR_inst `REG1; 
			end
			else begin
				
				Field1_REG1 <= 4'b1111;
			end
			
			//if the first opcode is an accumulator reader, assign the appropriate accumulator to Field1_ACC0
			if ((RR_inst `OPcode1 == `OPst) || (RR_inst `OPcode1 == `OPa2r) || (RR_inst `OPcode1 == `OPli) || (RR_inst `OPcode1 == `OPlf) || (RR_inst `OPcode1 == `OPsh) || (RR_inst `OPcode1 == `OPslt) || (RR_inst `OPcode1 == `OPxor) || (RR_inst `OPcode1 == `OPand) || (RR_inst `OPcode1 == `OPor)) begin
			
				Field1_ACC0 <= 4'b0000;
			end
			else begin
			
				Field1_ACC0 <= 4'b1111;
			end
			
			//if the instruction can be a two opcode word
			if (RR_inst `OPcode1 > `OPjr) begin	
			
				//if the possible second opcode is a register reader, assign the appropriate register to Field2_REG2
				if ((RR_inst `OPcode2 == `OPst) || (RR_inst `OPcode2 == `OPcvt) || (RR_inst `OPcode2 == `OPr2a) || (RR_inst `OPcode2 == `OPsh) || (RR_inst `OPcode2 == `OPslt) || (RR_inst `OPcode2 == `OPadd) || (RR_inst `OPcode2 == `OPsub) || (RR_inst `OPcode2 == `OPdiv) || (RR_inst `OPcode2 == `OPmul) || (RR_inst `OPcode2 == `OPnot) || (RR_inst `OPcode2 == `OPxor) || (RR_inst `OPcode2 == `OPand) || (RR_inst `OPcode2 == `OPor)) begin
					
					Field2_REG2 <= RR_inst `REG2; 
				end
				else begin 
					
					Field2_REG2 <= 4'b1111;
				end
				
				//if the possible second opcode is an accumulator reader, assign the appropriate accumulator to Field2_ACC1
				if ((RR_inst `OPcode2 == `OPst) || (RR_inst `OPcode2 == `OPa2r) || (RR_inst `OPcode2 == `OPli) || (RR_inst `OPcode2 == `OPlf) || (RR_inst `OPcode2 == `OPsh) || (RR_inst `OPcode2 == `OPslt) || (RR_inst `OPcode2 == `OPxor) || (RR_inst `OPcode2 == `OPand) || (RR_inst `OPcode2 == `OPor)) begin	
					
					Field2_ACC1 <= 4'b0001;
				end
				else begin
					
					Field2_ACC1 <= 4'b1111;
				end
			end
		end
	end
endmodule

//Identifies which registers and accumulators are written to, if any
module RegistersWrittenTo(Write_ACC0, Write_REG1, Write_ACC1, Write_REG2, WR_inst, clk);
	output reg [3:0] Write_ACC0;
	output reg [3:0] Write_REG1;
	output reg [3:0] Write_ACC1;
	output reg [3:0] Write_REG2;
	
	input clk;
	
	input `WORD WR_inst;
	
	always @(posedge clk) begin
		//if the first opcode is a register writer, assign the appropriate register to Write_REG1
		if ((WR_inst `OPcode1 == `OPcf8) || (WR_inst `OPcode1 == `OPci8) || (WR_inst `OPcode1 == `OPa2r) || (WR_inst `OPcode1 == `OPli) || (WR_inst `OPcode1 == `OPlf)) begin
			
			Write_REG1 <= WR_inst `REG1; 
		end
		else begin
			
			Write_REG1 <= 4'b1110;
		end
		
		//if the first opcode is an accumulator writer, assign the appropriate accumulator to Write_ACC0
		if ((WR_inst `OPcode1 == `OPcvt) || (WR_inst `OPcode1 == `OPr2a) || (WR_inst `OPcode1 == `OPsh) || (WR_inst `OPcode1 == `OPslt) || (WR_inst `OPcode1 == `OPadd) || (WR_inst `OPcode1 == `OPsub) || (WR_inst `OPcode1 == `OPdiv) || (WR_inst `OPcode1 == `OPmul) || (WR_inst `OPcode1 == `OPnot) || (WR_inst `OPcode1 == `OPxor) || (WR_inst `OPcode1 == `OPand) || (WR_inst `OPcode1 == `OPor)) begin
		
			Write_ACC0 <= 4'b0000;
		end
		else begin
		
			Write_ACC0 <= 4'b1110;
		end
		
		//if the instruction can be a two opcode word
		if (WR_inst `OPcode1 > `OPjr) begin
			
			//if the possible second opcode is a register writer, assign the appropriate register to Write_REG2
			if ((WR_inst `OPcode2 == `OPcf8) || (WR_inst `OPcode2 == `OPci8) || (WR_inst `OPcode2 == `OPa2r) || (WR_inst `OPcode2 == `OPli) || (WR_inst `OPcode2 == `OPlf)) begin
				
				Write_REG2 <= WR_inst `REG2; 
			end
			else begin 
				
				Write_REG2 <= 4'b1110;
			end
			
			//if the second opcode is an accumulator writer, assign the appropriate accumulator to Write_ACC1
			if ((WR_inst `OPcode1 == `OPcvt) || (WR_inst `OPcode1 == `OPr2a) || (WR_inst `OPcode1 == `OPsh) || (WR_inst `OPcode1 == `OPslt) || (WR_inst `OPcode1 == `OPadd) || (WR_inst `OPcode1 == `OPsub) || (WR_inst `OPcode1 == `OPdiv) || (WR_inst `OPcode1 == `OPmul) || (WR_inst `OPcode1 == `OPnot) || (WR_inst `OPcode1 == `OPxor) || (WR_inst `OPcode1 == `OPand) || (WR_inst `OPcode1 == `OPor)) begin
		
			Write_ACC1 <= 4'b0001;
			end
			else begin
		
			Write_ACC1 <= 4'b1110;
			end
		end
	end
endmodule

// ALU0 - First phase of the ALU
// Outputs
//	 	outVal is the output of the first phase alu
// 		out1 is the value of the first register used in the instruction.
// 		out2 is the value of the second register used in the isntruction.
// 		oreg1 is the first register used in the alu
// 		oreg2 is the second register used in the alu
// 		oop is the opcode of the instruction
//		otyp is the type of the numbers
// Inputs
// 		in1 is the value of the first register used in the instruction.
// 		in2 is the value of the second register used in the isntruction.
// 		ireg1 is the first register used in the alu.
// 		ireg2 is the second register used in the alu.
// 		iop is the opcode of the instruction.
// 		typ is the type of number used in the operation.
module ALU0(outVal, out1, out2, oimm, oinst, rin1, rin2, iimm, iinst, iver, clk);
	input signed `TYPEDREG rin1, rin2;
	input `WORD iinst;
	input iver;
	input `TYPEDREG iimm;
	input clk;
	output signed `TYPEDREG out1, out2;
	output reg `TYPEDREG outVal;
	output `TYPEDREG oimm;
	output `WORD oinst;

	reg `OPcodeID op;
	wire typ;
	wire `REGNUM regNum;
	wire `WORD in1, in2;

	assign typ = rin2[16];
	assign out1 = rin1;
	assign out2 = rin2;
	assign in1 = rin1`WORD;
	assign in2 = rin2`WORD;
	assign oimm = iimm;
	assign oinst = iinst;
	assign over = iver;

	wire signed `WORD recr, addr, subr, shr, mulr;
	wire signed `WORD outand, outor, outnot, outxor, outslt;
	wire signed `WORD cvti, cvtf;
	wire signed sltr;

	//Assign the bitwise operations.
	assign outand = in1 & in2;
	assign outor = in1 | in2;
	assign outxor = in1 ^ in2;
	assign outnot = ~in2;
	assign outslt = in1 < in2;
	fadd fa(addr, in1, in2);
	fadd fsu(subr, in1, {~in2[15], in2[14:0]});
	fmul fm(mulr, in1, in2);
	frecip fr(recr, in2);
	fshift fs(shr, in1, in2);
	fslt fsl(sltr, in1, in2);
	i2f icvt(cvti, in2);
	f2i fcvt(cvtf, in2);
	always @(*) begin
		if (iver == 1'b0)
			op = iinst`OPcode1;
		else
			op = iinst`OPcode2;
		case(op)
			`OPadd: begin
				case(typ)
					0: begin outVal = {typ, in1 + in2}; end
					1: begin outVal = {typ, addr}; end
				endcase
			end
			`OPsub: begin
				case(typ)
					0: begin outVal = {typ, in1 - in2}; end
					1: begin outVal = {typ, subr}; end
				endcase
			end
			`OPmul: begin
				case(typ)
					0: begin outVal = {typ, in1 * in2}; end
					1: begin outVal = {typ, mulr}; end
				endcase
			end
			`OPdiv: begin
				case(typ)
					0: begin outVal = {typ, in1 / in2}; end
					1: begin outVal = {typ, recr}; end // Only perform the first half of the fp division here.
				endcase
			end
			`OPand: begin outVal = {typ, outand}; end
			`OPor:  begin outVal = {typ, outor}; end
			`OPxor: begin outVal = {typ, outxor}; end
			`OPnot: begin outVal = {typ, outnot}; end
		//Positive indicates left shift.
			`OPsh:  begin
				case(typ)
					0:  begin outVal = {typ, in1 << in2}; end
					1:  begin outVal = {typ, shr}; end
				endcase
			end
			`OPslt: begin
				case(typ)
					0:  begin outVal = {typ, outslt}; end
					1:  begin outVal = {typ, 15'b0, sltr}; end
				endcase
			end
			`OPcvt: begin
				case(typ)
					0:  begin outVal = {!typ, cvti}; end
					1:  begin outVal = {!typ, cvtf}; end
				endcase
			end
			`OPjz8: begin
				outVal = {0'b0, in2};
			end
			`OPjnz8: begin
				outVal = {0'b0, in2};
			end
			default: outVal = 17'b0;
		endcase
	end
endmodule

// ALU1 - second phase of the ALU
// Outputs
//	 	outVal is the output of the second phase alu
// 		out1 is the value of the first register used in the instruction.
// 		out2 is the value of the second register used in the isntruction.
// 		oreg1 is the first register used in the alu
// 		oreg2 is the second register used in the alu
// 		oop is the opcode of the instruction
//		otyp is the type of the numbers
// Inputs
// 		inVal is the value received from ALU part 1.
// 		in1 is the value of the first register used in the instruction.
// 		in2 is the value of the second register used in the isntruction.
// 		ireg1 is the first register used in the alu.
// 		ireg2 is the second register used in the alu.
// 		iop is the opcode of the instruction.
// 		typ is the type of number used in the operation.
module ALU1(outVal, out1, out2, oimm, oinst, inVal, rin1, rin2, iimm, iinst, iver, clk);
	input signed `TYPEDREG rin1, rin2;
	input `WORD iinst;
	input `TYPEDREG inVal, iimm;
	input iver;
	input clk;
	output signed `TYPEDREG out1, out2;
	output reg `TYPEDREG outVal;
	output `TYPEDREG oimm;
	output `WORD oinst;

	reg `OPcodeID op;
	wire typ;
	wire `REGNUM regNum;
	wire `WORD in1, in2;

	assign typ = rin2[16];
	assign out1 = rin1;
	assign out2 = rin2;
	assign in1 = rin1`WORD;
	assign in2 = rin2`WORD;
	assign oimm = iimm;
	assign oinst = iinst;

	wire signed `WORD divr;

	fmul fd(divr, in1, inVal`WORD);
	always @(*) begin
		if (iver == 1'b0)
			op = iinst`OPcode1;
		else
			op = iinst`OPcode2;
		case(op)
			`OPdiv: begin
				case(typ)
					0: begin outVal = inVal; end
					1: begin outVal = {typ, divr}; end
				endcase
			end
			`OPlf: begin
				outVal[15] = 1'b1;
			end	
			`OPli: begin
				outVal[15] = 1'b0;
			end	
			default: outVal = inVal;
		endcase
	end
endmodule

module tacky_processor(halt, reset, clk);

output reg halt;
input reset;
input clk;

//stage 0 regs & memory
reg `WORD pc, pc_inc, instruction;
reg `WORD instruction_mem `MEMSIZE;

//stage 1 regs
reg `TYPEDREG regfile `REGSIZE;
reg `HALFWORD pre;
reg `TYPEDREG acc0_val, acc1_val, r1_val, r2_val, imm_to_ALUMEM;
reg `WORD ins_to_ALUMEM;

//stage 2 regs & memory
reg `WORD data_mem `MEMSIZE;
reg `WORD ins_to_ALU2;
reg `TYPEDREG data1_to_ALU2, data2_to_ALU2, imm_to_ALU2;
reg `TYPEDREG r0_to_ALU2, r1_to_ALU2, r2_to_ALU2, r3_to_ALU2;

wire `TYPEDREG alu0_0outVal, alu1_0outVal, alu0_1outVal, alu1_1outVal;
wire `TYPEDREG alu0_0out1, alu1_0out1, alu0_1out1, alu1_1out1;
wire `TYPEDREG alu0_0out2, alu1_0out2, alu0_1out2, alu1_1out2;
wire `TYPEDREG alu0_0oimm, alu1_0oimm, alu0_1oimm, alu1_1oimm;
wire `WORD alu0_0oinst, alu1_0oinst, alu0_1oinst, alu1_1oinst;
wire `TYPEDREG alu1_0inVal, alu1_1inVal;
wire `TYPEDREG alu0_0in1, alu1_0in1, alu0_1in1, alu1_1in1;
wire `TYPEDREG alu0_0in2, alu1_0in2, alu0_1in2, alu1_1in2;
wire `TYPEDREG alu0_0iimm, alu1_0iimm, alu0_1iimm, alu1_1iimm;
wire `WORD alu0_0iinst, alu1_0iinst, alu0_1iinst, alu1_1iinst;

ALU0 alu0_0(alu0_0outVal, alu0_0out1, alu0_0out2, alu0_0oimm, alu0_0oinst, alu0_0in1, alu0_0in2, alu0_0iimm, alu0_0iinst, 1'b0, clk);
ALU0 alu0_1(alu0_1outVal, alu0_1out1, alu0_1out2, alu0_1oimm, alu0_1oinst, alu0_1in1, alu0_1in2, alu0_1iimm, alu0_1iinst, 1'b1, clk);

//stage 3 regs
reg `WORD ins_to_WB;
reg `TYPEDREG data1_to_WB, data2_to_WB, imm_to_WB;
reg `TYPEDREG ALU1_result, ALU2_result;
	
ALU1 alu1_0(alu1_0outVal, alu1_0out1, alu1_0out2, alu1_0oimm, alu1_0oinst, alu1_0inVal, alu1_0in1, alu1_0in2, alu1_0iimm, alu1_0iinst, 1'b0, clk);
ALU1 alu1_1(alu1_1outVal, alu1_1out1, alu1_1out2, alu1_1oimm, alu1_1oinst, alu1_1inVal, alu1_1in1, alu1_1in2, alu1_1iimm, alu1_1iinst, 1'b1, clk);

//stage 4 regs
reg jump_flag;
reg `WORD pc_next;

wire IF_JumpFlag;
wire IF_SysFlag;

//Determines which registers are being read from in stage 0 (1111 if not read from)
RegistersReadFrom RegsRead(R_ACC0, R_REG1, R_ACC1, R_REG2, IF_JumpFlag, IF_SysFlag, instruction, clk);
//Determines which registers are being written to in stage 1 (1110 if not written to)
RegistersWrittenTo RegsWritten1(W1_ACC0, W1_REG1, W1_ACC1, W1_REG2, ins_to_ALUMEM, clk);
//Determines which registers are being written to in stage 2 (1110 if not written to)
RegistersWrittenTo RegsWritten2(W2_ACC0, W2_REG1, W2_ACC1, W2_REG2, ins_to_ALU2, clk);
//Determines which registers are being written to in stage 3 (1110 if not written to)
RegistersWrittenTo RegsWritten3(W3_ACC0, W3_REG1, W3_ACC1, W3_REG2, ins_to_WB, clk);


//Registers read from in stage 0
wire [3:0] R_ACC0;
wire [3:0] R_REG1;
wire [3:0] R_ACC1;
wire [3:0] R_REG2;

//Registers written to in stage 1
wire [3:0] W1_ACC0;
wire [3:0] W1_REG1;
wire [3:0] W1_ACC1;
wire [3:0] W1_REG2;

//Registers written to in stage 2
wire [3:0] W2_ACC0;
wire [3:0] W2_REG1;
wire [3:0] W2_ACC1;
wire [3:0] W2_REG2;

//Registers written to in stage 3
wire [3:0] W3_ACC0;
wire [3:0] W3_REG1;
wire [3:0] W3_ACC1;
wire [3:0] W3_REG2;

//NOPs needed to avoid dependency
//NOPs needed to avoid dependency
reg [2:0] NOPs_win1;
reg [2:0] NOPs_win2;
reg [2:0] NOPs1;
reg [2:0] NOPs2;
reg [2:0] NOPs3;
reg [2:0] NOPs4;
reg [2:0] NOPs, NOP_timer;

always@(posedge reset) begin
    $readmemh0(regfile);
    $readmemh1(instruction_mem);
    $readmemh2(data_mem);
    pc <= 0;
    pc_inc <= 1;
    pre <= 0;
    halt <= 0;
    pre <= 0;
	NOPs <= 0;
	NOP_timer <= 0;
	jump_flag <= 0;
end

//stage 0: instruction fetch
always@(posedge clk) begin
    if(NOP_timer == 0  && !IF_SysFlag) begin
        pc <= (jump_flag) ? pc_next : pc_inc;
    end
    instruction <= instruction_mem[pc];
    pc_inc <= pc + 1; 
	
    //If a sys call, make sure everything finishes by padding with NOPs
    if (IF_SysFlag == `true) begin
		$display("sys-nop");
		NOPs = 6;
    end
	//If a jump, pad NOPs until jump is gone
    else if (IF_JumpFlag == `true) begin
		$display("jump-nop");
		NOPs = 7;
    end
    //Else, check for dependencies
    else begin 
	    //Checks for dependencies on accumulator 0
	    case (R_ACC0)
            W1_ACC0 : NOPs1 = 4;
            W2_ACC0 : NOPs1 = 3;
            W3_ACC0 : NOPs1 = 2;
            default : NOPs1 = 0;
	    endcase

	    //Checks for dependencies on accumulator 1
	    case (R_ACC1)
	    	W1_ACC1 : NOPs2 = 4;
            W2_ACC1 : NOPs2 = 3;
            W3_ACC1 : NOPs2 = 2;
            default : NOPs2 = 0;
        endcase

	    //Checks for dependencies on reg 1 
	    case (R_REG1)
            W1_ACC0 : NOPs3 = 4;
            W1_ACC1 : NOPs3 = 4;
            W1_REG1 : NOPs3 = 4;
            W1_REG2 : NOPs3 = 4;
            W2_ACC0 : NOPs3 = 3;
            W2_ACC1 : NOPs3 = 3;
            W2_REG1 : NOPs3 = 3;
            W2_REG2 : NOPs3 = 3;
            W3_ACC0 : NOPs3 = 2;
            W3_ACC1 : NOPs3 = 2;
            W3_REG1 : NOPs3 = 2;
            W3_REG2 : NOPs3 = 2;
            default : NOPs3 = 0;	
	    endcase

	    //Checks for dependencies on reg 2 
	    case (R_REG2)
            W1_ACC0 : NOPs4 = 4;
            W1_ACC1 : NOPs4 = 4;
            W1_REG1 : NOPs4 = 4;
            W1_REG2 : NOPs4 = 4;
            W2_ACC0 : NOPs4 = 3;
            W2_ACC1 : NOPs4 = 3;
            W2_REG1 : NOPs4 = 3;
            W2_REG2 : NOPs4 = 3;
            W3_ACC0 : NOPs4 = 2;
            W3_ACC1 : NOPs4 = 2;
            W3_REG1 : NOPs4 = 2;
            W3_REG2 : NOPs4 = 2;
            default : NOPs4 = 0;	
	    endcase
		
	    if (NOPs1 > NOPs2) begin
	 	NOPs_win1 = NOPs1;
        end
	    else begin
		NOPs_win1 = NOPs2;
	    end
		
	    if (NOPs3 > NOPs4) begin
		NOPs_win2 = NOPs3;
	    end
	    else begin
		NOPs_win2 = NOPs4;
	    end
		
	    if (NOPs_win1 > NOPs_win2) begin
		NOPs = NOPs_win1;
	    end
	    else begin
		NOPs = NOPs_win2;
	    end
    end
end

//stage 1: register read

always@(posedge clk) begin
    NOP_timer <= (NOPs > 0 && NOP_timer == 0) ? NOPs : NOP_timer;
    if(NOP_timer > 0) begin
        ins_to_ALUMEM <= {`OPno, 3'b000, `OPno, 3'b001 };
        NOP_timer <= NOP_timer - 1; 
    end
    else begin
        if(instruction `OPcode1 == `OPpre) pre <= instruction `IMM8;
        if(instruction `OPcode1 == `OPcf8) imm_to_ALUMEM <= {`Float, pre, instruction `IMM8};
        if(instruction `OPcode1 == `OPci8 || (instruction `OPcode1 >= `OPjp8 && instruction `OPcode1 <= `OPjnz8) ) imm_to_ALUMEM <= {`Int, pre, instruction `IMM8};
        acc0_val <= regfile[0];
        acc1_val <= regfile[1];
        r1_val <= regfile[instruction`REG1];
        r2_val <= regfile[instruction`REG2];
        ins_to_ALUMEM <= instruction;
    end
end

//stage 2: ALU/MEM

assign alu0_0iinst = ins_to_ALUMEM;
assign alu0_1iinst = ins_to_ALUMEM;
assign alu0_0in1 = acc0_val;
assign alu0_0in2 = r1_val;
assign alu0_1in1 = acc1_val;
assign alu0_1in2 = r2_val;
assign alu0_0iimm = imm_to_ALUMEM;
assign alu0_1iimm = imm_to_ALUMEM;

always@(posedge clk) begin
	ins_to_ALU2 <= alu0_0oinst;
	imm_to_ALU2 <= alu0_0oimm;
	r0_to_ALU2 <= alu0_0out1;
	r1_to_ALU2 <= alu0_0out2;
	r2_to_ALU2 <= alu0_1out1;
	r3_to_ALU2 <= alu0_1out2;
	
	case (alu0_0iinst`OPcode1)
		`OPlf: begin
			data1_to_ALU2 <= data_mem[alu0_0out1];
		end	
		`OPli: begin
			data1_to_ALU2 <= data_mem[alu0_0out1];
		end	
		`OPst: begin
			data_mem[alu0_0out2] <= alu0_0out1;
			data1_to_ALU2 <= 17'b0;
		end
		default: data1_to_ALU2 <= alu0_0outVal;
	endcase
	
	case (alu0_0iinst`OPcode2)
		`OPlf: begin
			data2_to_ALU2 <= data_mem[alu0_1out1];
		end	
		`OPli: begin
			data2_to_ALU2 <= data_mem[alu0_1out1];
		end	
		`OPst: begin
			data_mem[alu0_1out2] <= alu0_1out1;
			data1_to_ALU2 <= 17'b0;
		end
		default: data2_to_ALU2 <= alu0_1outVal;
	endcase
end

//stage 3: ALU2

assign alu1_0iinst = ins_to_ALU2;
assign alu1_1iinst = ins_to_ALU2;
assign alu1_0in1 = r0_to_ALU2;
assign alu1_0in2 = r1_to_ALU2;
assign alu1_1in1 = r2_to_ALU2;
assign alu1_1in2 = r3_to_ALU2;
assign alu1_0iimm = imm_to_ALU2;
assign alu1_1iimm = imm_to_ALU2;
assign alu1_0inVal = data1_to_ALU2;
assign alu1_1inVal = data2_to_ALU2;

always@(posedge clk) begin
    ins_to_WB <= alu1_0oinst;
    imm_to_WB <= alu1_0oimm;
    data1_to_WB <= alu1_0outVal;
    data2_to_WB <= alu1_1outVal;
	ALU1_result <= alu1_0outVal;
	ALU2_result <= alu1_1outVal;
end

//stage 4: writeback
always@(posedge clk) begin
    //Opcode1 >= sh -> ALU op
    //Opcode1 >= jr -> 2 OP/Word
    //Opcode1 <= jr -> 1 OP/Word
    //reg1_load <= ins_to_WB `REG1
    //reg2_load <= ins_to_WB `REG2
    
    if(ins_to_WB `OPcode1 == `OPsys) begin
        halt = 1'b1;
    end

    //First instruction logic WB
    if(ins_to_WB `OPcode1 >= `OPsh || ins_to_WB `OPcode1 == `OPr2a) begin
        regfile `Acc0 <= ALU1_result;
    end
    else if (ins_to_WB `OPcode1 == `OPlf || ins_to_WB `OPcode1 == `OPli) begin
        regfile[ins_to_WB `REG1] <= data1_to_WB;
    end
    else if (ins_to_WB `OPcode1 == `OPa2r) begin
        regfile[ins_to_WB `REG1] <= ALU1_result;
    end
    else if(ins_to_WB `OPcode1 == `OPcf8 || ins_to_WB `OPcode1 == `OPci8) begin
        regfile[ins_to_WB `REG1] <= imm_to_WB;
    end
    
    //Second instruction logic WB (if present)
    if(ins_to_WB `OPcode1 >= `OPjr)  begin
        if(ins_to_WB `OPcode2 >= `OPsh || ins_to_WB `OPcode2== `OPr2a) begin
            regfile `Acc1 <= ALU2_result;
        end
        else if (ins_to_WB `OPcode2 == `OPlf || ins_to_WB `OPcode2 == `OPli) begin
            regfile[ins_to_WB `REG2] <= data2_to_WB;
        end
        else if (ins_to_WB `OPcode2 == `OPa2r) begin
            regfile[ins_to_WB `REG2] <= ALU2_result;
        end
    end
    
    //jump instruction logic WB
    if(ins_to_WB `OPcode1 == `OPjp8) begin
        pc_next <= ALU1_result;
        jump_flag <= `true;
    end
    else if (ins_to_WB `OPcode1 == `OPjz8) begin
        pc_next <= imm_to_WB;
        jump_flag <= (ALU1_result == 0) ? `true : `false;
    end
    else if (ins_to_WB `OPcode1 == `OPjnz8) begin
        pc_next <= imm_to_WB;
        jump_flag <= (ALU1_result != 0) ? `true : `false;
    end
    else if(ins_to_WB `OPcode1 == `OPjr) begin
        pc_next <= ALU1_result;
        jump_flag <= `true;
    end
    else if(ins_to_WB `OPcode1 > `OPjr && ins_to_WB `OPcode2 == `OPjr) begin
        pc_next <= ALU1_result;
        jump_flag <= `true;
    end
    else begin
        jump_flag <= `false;
    end 
end
endmodule

module testbench;
reg reset = 0;
reg clk = 0;
wire halted;

tacky_processor PE(halted, reset, clk);

initial begin
  $dumpfile;
  $dumpvars(1, PE.halt, PE.NOPs, PE.NOP_timer, PE.pre, PE.pc, PE.pc_next, PE.pc_inc, PE.jump_flag, PE.instruction, PE.ins_to_ALUMEM, PE.ins_to_ALU2, PE.ins_to_WB, PE.imm_to_ALUMEM, PE.alu0_0outVal, PE.alu1_0outVal, PE.ALU1_result, PE.IF_JumpFlag, PE.IF_SysFlag, clk);
  $dumpvars(2, PE.alu0_0.op, PE.alu0_0.outVal);
  #10 reset = 1;
  #10 reset = 0;
  while (!halted) begin
    #10 clk = 1;
    #10 clk = 0;
  end
  $display("finished");
  $finish;
end

initial begin
    #(10**4) $finish;
end

endmodule
