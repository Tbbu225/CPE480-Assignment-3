

`Defines
// basic sizes of things
`define WORD        [15:0]
`define Opcode1 [15:11]
`define Opcode2 [7:3]
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
`define Start       5’b01000


//Module stuff,
`define ALU     5’b1xxxx


module ALUmod(out, in1, in2, op, type);
input [4:0] op;
input type; //0->integer arithmetic, 1->floating point arithmetic
//$acc should always be in1 and $r should be in2.
input signed `WORD in1, in2;
output reg `WORD out;
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
           1: begin out <= addr; end
           endcase
           end
   `OPsub: begin
           case(type)
           0: begin out <= in1 - in2; end
           1: begin out <= subr; end
           endcase
           end
   `OPmul: begin
           case(type)
           0: begin out <= in1 * in2; end
           1: begin out <= mulr; end
           endcase
           end
   `OPdiv: begin
           case(type)
           0: begin out <= in1 / in2; end
           1: begin out <= divr; end
           endcase
           end
   `OPand: begin out <= outand; end
   `OPor:  begin out <= outor; end
   `OPxor: begin out <= outxor; end
   `OPnot: begin out <= outnot; end
   //Positive indicates left shift.
   `OPsh:  begin
           case(type)
           0:  begin out <= in1 << in2; end
           1:  begin out <= shr; end
           endcase
           end
   `OPslt: begin
           case(type)
           0:  begin out = outslt; end
           1:  begin out <= sltr; end
           endcase
           end
   default: out <= 16'b0;
 endcase
end
endmodule


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

wire `WORD ALU1out, ALU2out;
ALUmod ALU1(ALU1out, regfile[0], regfile[`REG1 ir], `OPcode1 ir, regfile[0][16]);
ALUmod ALU2(ALU2out, regfile[1], regfile[`REG2 ir], `OPcode2 ir, regfile[1][16]);

always @(reset) begin
 halt = 0;
 pc = 0;
 s = `Start;
 $readmemh0(regfile);
 $readmemh1(mainmem);
end





always @(posedge clk) begin
 case (s)
    `Start:        begin
ir <= mainmem[pc];  
s <= `Start1;
end


`Start1:   begin
                   pc <= pc + 1;            	   // bump pc
                   s <= ir `Opcode1;   // most instructions, state # is opcode
field <= 0;
end
    `Start2:    begin
s <= ir `Opcode2;
                field <= 1;
            end

    //Load an 8-bit immediate into the 8-bit “pre” register. Set next state to `Start because the //“pre” instruction will always take up the entire instruction word.
    `OPpre: begin
                pre <= `IMM8 ir;
                s <= `Start;
            end

    //Set the processor halt output to 1, then set next state to no-op(?).
    `OPsys: begin
                halt <= 1;
                s <= `OPno;
            end

    //Loads the bitfield generated by concatenating together a 1 (its a float), the 8-bit “pre” register, and the 8-bit immediate that is part of the instruction into the register specified by the instruction.
    `OPcf8: begin
                regfile[`REG1 ir] <= {1, pre, `IMM8 ir};
                s <= `Start;
            end



    //Loads the bitfield generated by concatenating together a 0 (its an int), the 8-bit “pre” register, and the 8-bit immediate that is part of the instruction into the register specified by the instruction.
    `OPci8: begin
                regfile[`REG1 ir] <= {0, pre, `IMM8 ir};
                s <= `Start;
            end

    //Loads the value of the accumulator associated with the field the instruction was in into the register specified by the instruction. The logic described here only latches the value in the respective accumulator if the instruction was in the correct field.
    `OPa2r: begin
                regfile[`REG1 ir] <= (field ? regfile[`REG1 ir] : regfile[field]);
                regfile[`REG2 ir] <= (!field ? regfile[`REG2 ir] : regfile[field]);
                s <= ( field ? `Start :  `Start2);
            end

    //Loads the value of the specified register into the accumulator that is associated with the field that the instruction is in. Once loaded, the machine will restart at the Start state if the instruction was in the second field. If it was in the first field, it will move to the Start2 state so that any instruction in the second field may be considered.
    `OPr2a: begin
                regfile[field] <= (field ? regfile[field] : regfile[`REG1 ir]);
                regfile[field] <= (!field ? regfile[field] : regfile[`REG2 ir]);
                s <= ( field ? `Start :  `Start2);
            end

    //Loads a floating point value from the accumulator that is associated with the field that the instruction is in into a specified register. To ensure the value is typed as a floating point, a 1 is concatenated to the most significant bit of the accumulator value. Once loaded, the machine will restart at the Start state if the instruction was in the second field. If it was in the first field, it will move to the Start2 state so that any instruction in the second field may be considered.
`OPlf:      begin
        regfile[`REG1 ir] <= (field ? regfile[`REG1 ir]:{1, mainmem[regfile[field]]});
            regfile[`REG2 ir] <= (!field ? regfile[`REG2 ir]:{1, mainmem[regfile[field]]});
        s <= ( field ? `Start :  `Start2);
        end





//Loads an integer value from the accumulator that is associated with the field that the instruction is in into a specified register. To ensure the value typed as an integer, a 1 is concatenated to the most significant bit of the accumulator value. Once loaded, the machine will restart at the Start state if the instruction was in the second field. If it was in the first field, it will move to the Start2 state so that any instruction in the second field may be considered.
    `OPli:      begin
        regfile[`REG1 ir] <= (field ? regfile[`REG1 ir] : {0,mainmem[regfile[field]]});
            regfile[`REG2 ir] <= (!field ? regfile[`REG2 ir] : {0,mainmem[regfile[field]]});
        s <= ( field ? `Start :  `Start2);
        end

//Stores the accumulator that associated with the field that the instruction is in into a memory address that is provided by a specified register. Once loaded, the machine will restart at the Start state if the instruction was in the second field. If it was in the first field, it will move to the Start2 state so that any instruction in the second field may be considered.
`OPst:      begin
        mainmem[`REG1 ir] <= (field ? regfile[`REG1 ir] : regfile[field]);
            mainmem[`REG2 ir] <= (!field ? regfile[`REG2 ir] : regfile[field]);
s <= ( field ? `Start :  `Start2);
end

    //Directs the program to jump to a specified register by loading the register number into the program counter register (pc). The state machine then restarts at the Start state.
    `OPjr:      begin
                pc <= regfile[`REG1 ir];
                s <= `Start;
            end

    //Directs the program to jump to the concatenation of the pre register and an 8-bit immediate value by loading this value into the program counter. The state machine then restarts at the Start state.
`OPjp8: begin
                pc <= {pre,`IMM8 ir};
                s <= `Start;
            end

    //Directs the program to jump to the concatenation of the pre register and an 8-bit immediate value by loading this value into the program counter if the value within a specified register is equal to zero. The state machine then restarts at the Start state. Spans both fields.
    `OPjz8: begin
                pc <= (`REG1 ir ? pc : {pre,`IMM8 ir});
                s <= `Start;
            end
//Directs the program to jump to the concatenation of the pre register and an 8-bit immediate value by loading this value into the program counter if the value within a specified register is NOT equal to zero. The state machine then restarts at the Start state. Spans both fields.
`OPjnz8:    begin
                pc <= (`REG1 ir ? {pre,`IMM8 ir} : pc);
                s <= `Start;
            end

    //Changes the type tag of a specified register by checking the register’s type tag and editing it such that it will be the opposite of what it originally was. The editing is performed by simply reloading the original value with the most significant bit (the type tag) flipped into the specified register.
    `OPcvt: begin
                regfile[field] <= (field ? ([16]regfile`REG1 ? cvt2int : cvt2float)) :  [16]regfile`REG2 ? cvt2int : cvt2float)
            end

//ALUmod takes in two 16-bit operands as in1 and in2, then outputs the result in out. It determines which operation to do by feeding the opcode into the module also. “Type” determines whether the operation should be integer or float, 0 for integer and 1 for float.
`ALU:       begin
regfile[field] <= (field ? ALU2out : ALU1out);
s <= (field ? `Start : `Start2);
end
