#!/usr/bin/env python
# VAMPyRE
# Verilog-A Model Pythonic Rule Enforcer
# version 1.5, 30-Oct-2020
#
# intended for checking for issues like:
# 1) hidden state (variables used before assigned)
# 2) unused parameters or variables
# 3) division by zero for parameters (1/parm where parm's range allows 0)
#    or domain errors for ln() and sqrt()
# 4) integer division (1/2 = 0)
# 5) incorrect ddx() usage
# 6) ports without direction and/or discipline
# 7) incorrect access functions for discipline
# 8) unnamed noise sources
# 9) use of features not appropriate for compact models
#10) various issues of poor coding style
#11) compliance with CMC Verilog-A Code Standards:
#   - use of lower-case identifiers
#   - proper use of tref and dtemp
#   - proper use of the multiplicity attribute
#12) misuse of limexp
#13) various problems with binning equations
#
# Copyright (c) 2020 Analog Devices, Inc.
# 
# Licensed under Educational Community License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License. You may
# obtain a copy of the license at http://opensource.org/licenses/ECL-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.



################################################################################
# setup

import os
import sys
import argparse
from copy import deepcopy

################################################################################
# global variables

gNatures       = {}
gDisciplines   = {}
gModuleName    = ""
gParameters    = {}
gVariables     = {}
gHiddenState   = {}
gUserFunctions = {}
gAccessFuncs   = {}
gPortnames     = {} # terminals
gNodenames     = {} # internal nodes
gBranches      = {}
gBlocknames    = {}
gStatementInCurrentBlock = False
gMissingConstantsFile = ""
gMacros        = {}
gIfDefStatus   = []
gScopeList     = []
gConditions    = [] # if-conditions in play
gCondBiasDep   = [] # whether conditions are bias-dependent:
                    # 0 (no), 1 (from if cond), 2 (yes), 3 (deriv error), 4 (ddx)
gAnalogInitial = []
gCurrentFunc   = 0
gFileName      = []
gLineNo        = []
gIncDir        = []
gDebug         = False
gFixIndent     = False
gMaxNum        = 5
gSpcPerInd     = 4
gStyle         = False
gVerbose       = False
gBinning       = False

# dictionaries to track how many of each type of message
gErrorMsgDict   = {}
gWarningMsgDict = {}
gNoticeMsgDict  = {}
gStyleMsgDict   = {}

# tokens
# (0-255 are ascii and latin-1)
TOKEN_UNUSED     = 256
TOKEN_NUMBER     = 257
TOKEN_IDENTIFIER = 258
TOKEN_STRING     = 259

# keywords
# from Verilog-AMS Language Reference Manual, version 2.4, Annex B
# VAMS LRM available from https://accellera.org/downloads/standards/v-ams
gVAMSkeywords = [
    "above", "abs", "absdelay", "absdelta", "abstol", "access", "acos", "acosh",
    "ac_stim", "aliasparam", "always", "analog", "analysis", "and", "asin",
    "asinh", "assert", "assign", "atan", "atan2", "atanh", "automatic",
    "begin", "branch", "buf", "bufif0", "bufif1",
    "case", "casex", "casez", "ceil", "cell", "cmos", "config", "connect",
    "connectmodule", "connectrules", "continuous", "cos", "cosh", "cross",
    "ddt", "ddt_nature", "ddx", "deassign", "default", "defparam", "design",
    "disable", "discipline", "discrete", "domain", "driver_update",
    "edge", "else", "end", "endcase", "endconfig", "endconnectrules",
    "enddiscipline", "endfunction", "endgenerate", "endmodule", "endnature",
    "endparamset", "endprimitive", "endspecify", "endtable", "endtask",
    "event", "exclude", "exp",
    "final_step", "flicker_noise", "floor", "flow", "for", "force", "forever",
    "fork", "from", "function",
    "generate", "genvar", "ground",
    "highz0", "highz1", "hypot",
    "idt", "idtmod", "idt_nature", "if", "ifnone", "incdir", "include", "inf",
    "initial", "initial_step", "inout", "input", "instance", "integer",
    "join",
    "laplace_nd", "laplace_np", "laplace_zd", "laplace_zp", "large",
    "last_crossing", "liblist", "library", "limexp", "ln", "localparam", "log",
    "macromodule", "max", "medium", "merged", "min", "module",
    "nand", "nature", "negedge", "net_resolution", "nmos", "noise_table",
    "noise_table_log", "nor", "noshowcancelled", "not", "notif0", "notif1",
    "or", "output",
    "parameter", "paramset", "pmos", "posedge", "potential", "pow", "primitive",
    "pull0", "pull1", "pulldown", "pullup", "pulsestyle_onevent",
    "pulsestyle_ondetect",
    "rcmos", "real", "realtime", "reg", "release", "repeat", "resolveto",
    "rnmos", "rpmos", "rtran", "rtranif0", "rtranif1",
    "scalared", "sin", "sinh", "showcancelled", "signed", "slew", "small",
    "specify", "specparam", "split", "sqrt", "string", "strong0", "strong1",
    "supply0", "supply1",
    "table", "tan", "tanh", "task", "time", "timer", "tran", "tranif0",
    "tranif1", "transition", "tri", "tri0", "tri1", "triand", "trior", "trireg",
    "units", "unsigned", "use", "uwire",
    "vectored",
    "wait", "wand", "weak0", "weak1", "while", "white_noise", "wire",
    "wor", "wreal",
    "xnor", "xor",
    "zi_nd", "zi_np", "zi_zd", "zi_zp"]

# compiler directives from VAMS LRM Section 10
gVAMScompdir = [
    "begin_keywords", "celldefine", "default_discipline", "default_nettype",
    "default_transition", "define", "else", "elsif", "end_keywords",
    "endcelldefine", "endif", "ifdef", "ifndef", "include", "line",
    "nounconnected_drive", "pragma", "resetall", "timescale",
    "unconnected_drive", "undef"]

# hierarchical system parameters from VAMS LRM Section 6.3.6
gVAMShiersysparm = ["$xposition", "$yposition", "$angle", "$hflip", "$vflip"]

# math functions and number of arguments
gMathFunctions = { "abs": 1, "limexp": 1, "log": 1, "$log10":1, "min": 2, "max": 2 }
# math functions with both name and $name form
for fn in ["acos", "acosh", "asin", "asinh", "atan", "atanh", "ceil",
           "cos", "cosh", "exp", "floor", "ln", "sin", "sinh", "sqrt",
           "tan", "tanh"]:
    gMathFunctions[fn] = 1
    gMathFunctions["$"+fn] = 1
for fn in ["atan2", "hypot", "pow"]:
    gMathFunctions[fn] = 2
    gMathFunctions["$"+fn] = 2
# not for use in compact models
for fn in ["idtmod", "absdelay", "analysis", "transition", "slew", "last_crossing",
           "laplace_nd", "laplace_np", "laplace_zd", "laplace_zp",
           "zi_nd", "zi_np", "zi_zd", "zi_zp"]:
    gMathFunctions[fn] = 99

# units for multiplicity checking
gUnitsMultiply = ["A", "Amps", "F", "Farads", "C", "coul", "Coulomb", "W", "Watts", "A/V", "S", "mho", "A*V", "A*V/K", "A/K", "C/K", "s*W/K"]
gUnitsDivide   = ["Ohm", "Ohms", "V/A", "K/W"]


################################################################################
# classes

class Nature:
  def __init__(self, name, defn):
    self.name = name
    self.defined = defn # False during pre-processing
    self.units = ""
    self.access = ""
    self.idt_nature = ""

class Discipline:
  def __init__(self, name):
    self.name = name
    self.potential = ""
    self.flow = ""
    self.domain = ""

class Parameter:
  def __init__(self, name, type, file, line):
    self.name = name
    self.type = type         # real, integer, string
    self.declare = [file, line]
    self.defv = 0
    self.units = 0
    self.range = "nb"        # nb, cc, co, oc, oo, cz, oz, sw
    self.max = float("inf")
    self.min = float("-inf")
    self.exclude = []
    self.inst = False        # instance parameter (or model)
    self.used = False        # whether used
    self.is_alias = False    # True for aliasparam (should not be used)

class Variable:
  def __init__(self, name, type, oppt, file, line):
    self.name = name
    self.type = type         # real, integer
    self.range = []
    self.oppt = oppt         # operating-point variable
    self.declare = [file, line]
    self.assign = -1         # line where value assigned
    self.conditions = []     # conditions when assigned
    self.used = False        # whether used
    self.bias_dep = 0        # bias-dependent: 0=no, 1=from if cond, 2=yes, 3=deriv error, 4=ddx
    self.biases = []         # node voltages (for ddx check)
    self.funcNameAndArg = ["", 0] # tracing whether function output arg is used

class Port:
  def __init__(self, name, file, line):
    self.name = name
    self.direction = ""      # inout, input, output
    self.discipline = ""     # electrical, thermal, ...
    self.declare = [file, line]
    self.is_bus = False
    self.msb = 0
    self.lsb = 0

class Branch:
  def __init__(self, name):
    self.name = name
    self.node1 = ""
    self.node2 = ""
    self.discipline = "" # electrical, thermal, ...
    self.lhs_flow = 0    # flow contrib: 0=no, 1=yes, 2=bias-dep condition
    self.lhs_pot  = 0    # potential contrib: 0=no, 1=yes, 2=bias-dep condition
    #self.rhs_flow = 0    # flow probe
    #self.rhs_pot  = 0    # potential probe

class Function:
  def __init__(self, name, type):
    self.name = name
    self.type = type         # real, integer
    self.args = []
    self.inputs = []
    self.outputs = []
    self.used = False        # whether used
    self.outarg_used = []    # whether output args are used in any call to the function

class Expression:
  def __init__(self, type):
    self.type = type         # NAME, STRING, NUMBER, FUNCCALL, or operator
    self.e1 = 0
    self.e2 = 0
    self.e3 = 0
    self.number = 0
    self.is_int = False      # is the expression an integer?
    self.args = []

  def getDependencies(self, assign_context, branch_contrib):
    if self.type == "NAME":
        return [self.e1]
    elif self.type == "NUMBER" or self.type == "STRING":
        return []
    elif self.type == "ARRAY":
        deps = []
        for arg in self.args:
            deps += arg.getDependencies(assign_context, branch_contrib)
        return deps
    elif self.type == "FUNCCALL":
        deps = []
        if self.e1 in gAccessFuncs:
            if len(gAnalogInitial) > 0 and gAnalogInitial[-1]:
                error("Branch access in analog initial block")
            nargs = len(self.args)
            if nargs > 0:
                nname1 = self.args[0].e1
                nname2 = ""
                if nname1 in gPortnames:
                    dname = gPortnames[nname1].discipline
                elif nname1 in gNodenames:
                    dname = gNodenames[nname1].discipline
                elif nname1 in gBranches:
                    dname = gBranches[nname1].discipline
                else:
                    dname = ""
                if dname in gDisciplines:
                    disc = gDisciplines[dname]
                    pname = disc.potential
                    fname = disc.flow
                    # use "potential" and "flow" to indicate bias-dependence
                    # (these keywords cannot collide with variable names)
                    if pname in gNatures and self.e1 == gNatures[pname].access:
                        #deps.append("potential")
                        pass
                    elif fname in gNatures and self.e1 == gNatures[fname].access:
                        #deps.append("flow")
                        pass
                    else:
                        error("Incorrect access function '%s' for %s '%s'" %
                              (self.e1, dname, nname1))
                elif dname == "":
                    error("Cannot determine discipline of '%s'" % nname1)
                if nargs > 2 or (nargs > 1 and (nname1 in gBranches \
                                   or self.args[1].e1 in gBranches)):
                    error("Invalid potential or flow access")
                elif nargs == 2:
                    nname2 = self.args[1].e1
                    dname2 = ""
                    if nname2 in gPortnames:
                        dname2 = gPortnames[nname2].discipline
                    elif nname2 in gNodenames:
                        dname2 = gNodenames[nname2].discipline
                    if dname2 != "" and dname2 != dname:
                        error("Nodes '%s' and '%s' belong to different disciplines" % (nname1, nname2))
                if nname1 != "":
                    if nname1 in gBranches:
                        branch = gBranches[nname1]
                        nname1 = branch.node1
                        nname2 = branch.node2
                    dep = self.e1 + "(" + nname1 + ")"
                    deps.append(dep)
                    if nname2 != "":
                        dep = self.e1 + "(" + nname2 + ")"
                        deps.append(dep)
            elif nargs == 0 and not branch_contrib:
                # if branch_contrib, error was reported already
                error("Missing node/branch name in potential or flow access")
        elif self.e1 == "$port_connected" or self.e1 == "$param_given":
            if len(self.args) == 1:
                nname = self.args[0].e1
                if self.e1 == "$port_connected" and not nname in gPortnames:
                    error("Invalid argument '%s' to %s, must be a port" % (nname, self.e1))
                if self.e1 == "$param_given" and not nname in gParameters:
                    error("Invalid argument '%s' to %s, must be a parameter" % (nname, self.e1))
            else:
                error("Incorrect number of arguments to %s" % self.e1)
        elif self.e1 == "limexp":
            if len(self.args) == 1:
                arg = self.args[0]
                deps = arg.getDependencies(assign_context, branch_contrib)
                [bias_dep, biases] = checkDependencies(deps, "Call to function '%s' depends on" % self.e1, 0, False, False)
                if not bias_dep:
                    warning("Call to %s(%s), but argument is not bias-dependent" % (self.e1, arg.asString()))
            # else: warning issued elsewhere
        elif self.e1 == "$limit":
            if len(self.args) == 0:
                error("Require at least one argument to %s" % self.e1)
            else:
                access = self.args[0]
                deps += access.getDependencies(assign_context, branch_contrib)
                if access.type != "FUNCCALL" or not access.e1 in gAccessFuncs:
                    error("First argument to $limit must be a potential or flow access function")
            if len(self.args) >= 2:
                limfn = self.args[1]
                if limfn.type == "NAME":
                    # identifier: user-defined function
                    fname = limfn.e1
                    if fname in gUserFunctions:
                        gUserFunctions[fname].used = True
                    else:
                        error("Unknown identifier '%s' specified as limiting function" % fname)
                elif limfn.type == "STRING":
                    # string: simulator function
                    fname = limfn.e1
                    if not fname in ["\"pnjlim\"", "\"fetlim\""]:
                        notice("Non-standard limiting function %s" % fname)
                else:
                    warning("Second argument to $limit should be a string or identifier")
        elif self.e1 == "ddx":
            if len(self.args) == 2:
                arg = self.args[0]
                argdeps = arg.getDependencies(assign_context, branch_contrib)
                deps = ["ddx"]
                deps += argdeps
                acc = self.args[1]
                # note: acc does not add a dependency if arg does not already depend on acc
                if acc.type == "FUNCCALL" and acc.e1 in gAccessFuncs:
                    accstr = acc.asString()
                    [bias_dep, biases] = checkDependencies(argdeps, "", 0, False, False)
                    if not accstr in biases:
                        warning("Call to ddx(), but %s does not depend on %s" % (arg.asString(), acc.asString()))
                        if gDebug:
                            print("    Depends on: %s" % biases)
            # else error printed during parsing
        elif branch_contrib and self.e1 in ["white_noise", "flicker_noise", "noise_table", "noise_table_log", \
                                            "laplace_zp", "laplace_zd", "laplace_np", "laplace_nd", \
                                            "zi_zp", "zi_zd", "zi_np", "zi_nd"]:
            # ignore dependencies in noise functions and laplace and z operators
            pass
        else:
            out_arg_pos = []
            if self.e1 in gUserFunctions:
                # find positions of output arguments
                funcdef = gUserFunctions[self.e1]
                funcdef.used = True
                args = funcdef.args
                for i in range(len(args)):
                    if args[i] in funcdef.outputs:
                        out_arg_pos.append(i)
            for i in range(len(self.args)):
                arg = self.args[i]
                if i in out_arg_pos:
                    if not assign_context:
                        error("Function '%s' with output arguments called in unexpected context" % self.e1)
                elif arg.type == "NOTHING":
                    pass # null argument
                else:
                    deps += arg.getDependencies(assign_context, branch_contrib)
            if assign_context and len(out_arg_pos) > 0:
                [bias_dep, biases] = checkDependencies(deps, "Call to function '%s' depends on" % self.e1, 0, False, True)
                for i in range(len(self.args)):
                    if i in out_arg_pos:
                        arg = self.args[i]
                        markVariableAsSet(arg.e1, bias_dep, Expression("NOTHING"), biases, False, True, self.e1, i)
        return deps

    elif self.type in ["!", "~"]:
        return self.e1.getDependencies(assign_context, branch_contrib)
    elif self.type in ["+", "-"]:
        dep1 = self.e1.getDependencies(assign_context, branch_contrib)
        if self.e2:
            dep2 = self.e2.getDependencies(assign_context, branch_contrib)
        else: # unary +/-
            dep2 = []
        return dep1 + dep2
    elif self.type in ["*", "/", "%", "**", \
                       "==", "!=", "<", ">", "<=", ">=", \
                       "&&", "||", "&", "|", "^"]:
        dep1 = self.e1.getDependencies(assign_context, branch_contrib)
        dep2 = self.e2.getDependencies(assign_context, branch_contrib)
        return dep1 + dep2
    elif self.type in ["<<", ">>", "<<<", ">>>", "===", "!==", \
                       "^~", "~^", "~&", "~|"]:
        # not expected in Verilog-A
        dep1 = self.e1.getDependencies(assign_context, branch_contrib)
        dep2 = self.e2.getDependencies(assign_context, branch_contrib)
        return dep1 + dep2
    elif self.type in ["?:"]:
        dep1 = self.e1.getDependencies(assign_context, branch_contrib)
        dep2 = self.e2.getDependencies(assign_context, branch_contrib)
        dep3 = self.e3.getDependencies(assign_context, branch_contrib)
        return dep1 + dep2 + dep3
    else: # pragma: no cover
        fatal("Unhandled expression type '%s' in getDependencies" % self.type)

  def asString(self):
    if self.type == "NAME" or self.type == "STRING":
        return self.e1
    elif self.type == "NUMBER":
        return str(self.number)
    elif self.type == "FUNCCALL":
        ret = self.e1 + "("
        first = True
        for arg in self.args:
            if first:
                first = False
            else:
                ret += ","
            ret += arg.asString()
        ret += ")"
        return ret
    elif self.type in ["!", "~"]:
        return self.type + "(" + self.e1.asString() + ")"
    elif self.type == "NOT":
        return "NOT(" + self.e1.asString() + ")"
    elif self.type in ["+", "-"]:
        if self.e2:
            return "(" + self.e1.asString() + self.type + self.e2.asString() + ")"
        else: # unary +/-
            return self.type + "(" + self.e1.asString() + ")"
    elif self.type == "&&":
        if self.e1.type in ["==", "!=", "<", ">", "<=", ">=", "&&", "||", "!", "NOT"] \
                and self.e2.type in ["==", "!=", "<", ">", "<=", ">=", "&&", "||", "!", "NOT"]:
            # don't need extra ()
            return self.e1.asString() + "&&" + self.e2.asString()
        else:
            return "(" + self.e1.asString() + self.type + self.e2.asString() + ")"
    elif self.type in ["*", "/", "%", "**", \
                       "==", "!=", "<", ">", "<=", ">=", \
                       "||", "&", "|", "^"]:
            return "(" + self.e1.asString() + self.type + self.e2.asString() + ")"
    elif self.type in ["<<", ">>", "<<<", ">>>", "===", "!==", \
                       "^~", "~^", "~&", "~|"]:
        # not expected in Verilog-A
        return "(" + self.e1.asString() + self.type + self.e2.asString() + ")"
    elif self.type in ["?:"]:
        return "(" + self.e1.asString() + ")?(" +  self.e2.asString() \
               + "):(" + self.e3.asString() + ")"
    elif self.type == "EVENT":
        ret = "@("
        first = True
        for arg in self.args:
            if first:
                first = False
            else:
                ret += ","
            ret += arg.asString()
        ret += ")"
        return ret
    #elif self.type == "PORT_FLOW":
    #    return "<" + self.e1 + ">"
    #elif self.type == "NOTHING":
    #    return "TODO"
    else: # pragma: no cover
        fatal("Unhandled expression type '%s' in asString" % self.type)


class Macro:
  def __init__(self, name, text):
    self.name = name
    self.text = text
    self.args = []

class Parser:
  def __init__(self, line):
    self.line = line
    self.cp   = 0
    self.end  = len(line)
    self.token = 0
    self.number = 0
    self.ternary_condition = Expression("NOTHING")
    self.is_int = False
    self.escaped = False
    self.string = ""

  def isNumber(self):
    global TOKEN_NUMBER
    return self.token == TOKEN_NUMBER

  def isIdentifier(self):
    global TOKEN_IDENTIFIER
    return self.token == TOKEN_IDENTIFIER

  def isString(self):
    global TOKEN_STRING
    return self.token == TOKEN_STRING

  def isEscaped(self):
    return self.escaped

  def getNumber(self):
    return self.number

  def getString(self):
    return self.string

  def getRestOfLine(self):
    retstr = self.line[self.cp:]
    self.cp = self.end
    return retstr

  def peekRestOfLine(self):
    return self.line[self.cp:]

  def eatSpace(self):
    while self.peekChar().isspace():
        self.getChar()

  def peekChar(self):
    if self.cp < self.end:
        return self.line[self.cp]
    else:
        return ""

  def getChar(self):
    if self.cp < self.end:
        ch = self.line[self.cp]
        self.cp += 1
        return ch
    else: # pragma: no cover
        return ""

  def ungetChar(self, ch):
    if self.cp > 0:
        self.cp -= 1
        if self.line[self.cp] != ch: # pragma: no cover
            fatal("Parse failure in ungetChar")

  def ungetChars(self, str):
    for i in range(len(str),0,-1):
        ch = str[i-1]
        self.ungetChar(ch)

  def lexNumber(self):
    value = 0
    str = ""
    sci_not = ""
    is_int = True
    while self.peekChar().isdigit():
        str += self.getChar()
    if self.peekChar() == '.':
        str += self.getChar()
        is_int = False
        while self.peekChar().isdigit():
            str += self.getChar()
    if self.peekChar() == '.': # pragma: no cover
        error("Invalid number")
    ch = self.peekChar()
    if ch == 'e' or ch == 'E':
        sci_not = self.getChar()
        ch = self.peekChar()
        last_ch = ""
        if ch == '+' or ch == '-':
            sci_not += self.getChar()
        while self.peekChar().isdigit():
            last_ch = self.getChar()
            sci_not += last_ch
        if last_ch.isdigit():
            str += sci_not
            is_int = False
        else:
            self.ungetChars(sci_not)
            sci_not = ""
    value = float(str)

    # if not sci not, look for SI prefixes
    if sci_not == "":
        ch = self.peekChar()
        if ch == 'T':
            value *= 1e12
            str += self.getChar()
        elif ch == 'G':
            value *= 1e9
            str += self.getChar()
        elif ch == 'M':
            value *= 1e6
            str += self.getChar()
        elif ch == 'K' or ch == 'k':
            value *= 1e3
            str += self.getChar()
        elif ch == 'm':
            value *= 1e-3
            str += self.getChar()
            is_int = False
        elif ch == 'u':
            value *= 1e-6
            str += self.getChar()
            is_int = False
        elif ch == 'n':
            value *= 1e-9
            str += self.getChar()
            is_int = False
        elif ch == 'p':
            value *= 1e-12
            str += self.getChar()
            is_int = False
        elif ch == 'f':
            value *= 1e-15
            str += self.getChar()
            is_int = False
        elif ch == 'a':
            value *= 1e-18
            str += self.getChar()
            is_int = False

    return [value, is_int]
  # end of lexNumber

  def lexName(self):
    str = ""
    while self.peekChar() == '$':
        str += self.getChar()
    while self.peekChar().isalnum() or self.peekChar() == '_':
        str += self.getChar()
    return str

  def lexString(self, qstr):
    str = self.getChar()
    ch = ""
    escaped = False
    while ch != qstr or escaped:
        ch = self.getChar()
        if escaped:
            str += ch
            ch = " "
            escaped = False
        elif ch == "\\":
            escaped = True
        else:
            str += ch
    return str

  def lex(self):
    global TOKEN_NUMBER, TOKEN_IDENTIFIER, TOKEN_STRING
    self.eatSpace()
    ch = self.peekChar()
    #print("ch='%s'" % ch)
    if ch == "":
        self.token = 0
    elif ch.isdigit():
        number = self.lexNumber()
        self.number = number[0]
        self.is_int = number[1]
        self.token = TOKEN_NUMBER
    elif ch.isalpha() or ch == '_' or ch == '$':
        self.string = self.lexName()
        self.token = TOKEN_IDENTIFIER
    elif ch == '\\':
        # escaped identifier
        self.string = ""
        ch = self.getChar()
        ch = self.getChar()
        while ch != "" and not ch.isspace():
            self.string += ch
            ch = self.getChar()
        if self.string != "":
            self.token = TOKEN_IDENTIFIER
            self.escaped = True
        else:
            ch = '\\'
            self.token = ord(ch)
    elif ch == '"':
        self.string = self.lexString(ch)
        self.token = TOKEN_STRING
    else:
        ch = self.getChar()
        self.token = ord(ch)
    return self.token

  # Verilog-A operators: VAMS LRM table 4-1
  #   Arithmetic: + - * / ** %
  #   Relational: > >= < <=
  #   Logical:    ! && || == !=
  #   Bitwise:    ~ & | ^
  #   Ternary:    ? :
  # Not supported in in Verilog-A:
  #   Concatenation, replication: {} {{}}
  #   Case: === !==
  #   Bitwise: ^~ ~^ ~& ~|
  #   Shifts: << >> <<< >>>

  def peekOper(self):
    oper = ""
    ch = self.peekChar()
    if ch in "+-*/%><!&|=~^?:":
        oper = ch
        if self.cp < self.end:
            next = self.line[self.cp+1]
        else:
            next = ""
        if ch == "+" and next == "+":
            oper += next
        elif ch == "-" and next == "-":
            oper += next
        elif ch == "*" and next == "*":
            oper += next
        elif ch in "+-*/" and next == "=": # pragma: no cover
            oper += next
        elif ch == "&" and next == "&":
            oper += next
        elif ch == "|" and next == "|":
            oper += next
        elif ch == "=" and next == "=":
            oper += next
        elif ch == "!" and next == "=":
            oper += next
        elif ch == "<" and next == "<":
            oper += next
        elif ch == "<" and next == "=":
            oper += next
        elif ch == ">" and next == ">":
            oper += next
        elif ch == ">" and next == "=":
            oper += next
        elif ch == "~" and next == "^":
            oper += next
        elif ch == "^" and next == "~":
            oper += next
        elif ch == "~" and next == "&":
            oper += next
        elif ch == "~" and next == "|":
            oper += next
        if len(oper) == 2 and self.cp < self.end-2:
            nextnext = self.line[self.cp+2]
            if oper == "<<" and nextnext == "<":
                oper = "<<<"
            elif oper == ">>" and nextnext == ">":
                oper = "<<<"
            elif oper == "==" and nextnext == "=":
                oper = "==="
            elif oper == "!=" and nextnext == "=":
                oper = "!=="
        # error printed by getOper() if it's used
        #if oper in ["++", "+=", "--", "-=", "*=", "/=", "<<", ">>", "<<<", ">>>", "===", "!==", "^~", "~^", "~&", "~|"]:
        #    error("Operator '%s' not valid in Verilog-A" % oper)
    return oper

  def getOper(self):
    oper = ""
    ch = self.peekChar()
    if ch in "+-*/%><!&|=~^?:":
        oper = self.getChar()
        next = self.peekChar()
        if ch == "+" and next == "+":
            oper += self.getChar()
        elif ch == "-" and next == "-":
            oper += self.getChar()
        elif ch == "*" and next == "*":
            oper += self.getChar()
        elif ch in "+-*/" and next == "=": # pragma: no cover
            oper += self.getChar()
        elif ch == "&" and next == "&":
            oper += self.getChar()
        elif ch == "|" and next == "|":
            oper += self.getChar()
        elif ch == "=" and next == "=":
            oper += self.getChar()
        elif ch == "!" and next == "=":
            oper += self.getChar()
        elif ch == "<" and next == "<":
            oper += self.getChar()
        elif ch == "<" and next == "=":
            oper += self.getChar()
        elif ch == ">" and next == ">":
            oper += self.getChar()
        elif ch == ">" and next == "=":
            oper += self.getChar()
        elif ch == "~" and next == "^":
            oper += self.getChar()
        elif ch == "^" and next == "~":
            oper += self.getChar()
        elif ch == "~" and next == "&":
            oper += self.getChar()
        elif ch == "~" and next == "|":
            oper += self.getChar()
        if len(oper) == 2 and self.cp < self.end:
            nextnext = self.peekChar()
            if oper == "<<" and nextnext == "<":
                oper += self.getChar()
            elif oper == ">>" and nextnext == ">":
                oper += self.getChar()
            elif oper == "==" and nextnext == "=":
                oper += self.getChar()
            elif oper == "!=" and nextnext == "=":
                oper += self.getChar()
        if oper in ["++", "+=", "--", "-=", "*=", "/=", "<<", ">>", "<<<", ">>>", "===", "!==", "^~", "~^", "~&", "~|"]:
            error("Operator '%s' not valid in Verilog-A" % oper)
    self.eatSpace()
    return oper

  # expression parsing - watch for order of operations

  def parseArgList(self, is_cond=False):
    args = []
    expr = self.getExpression(is_cond)
    if expr.type == "NOTHING" and self.peekChar() == ')':
        return []
    args.append(expr)
    self.eatSpace()
    ch = self.peekChar()
    while ch == ',':
        self.getChar()
        expr = self.getExpression()
        args.append(expr)
        self.eatSpace()
        ch = self.peekChar()
    return args

  def parsePrimary(self, is_cond):
    val = self.lex()
    if self.isIdentifier():
        expr = Expression("NAME")
        expr.e1 = self.getString()
    elif self.isNumber():
        expr = Expression("NUMBER")
        expr.number = self.getNumber()
        expr.is_int = self.is_int
    elif self.isString():
        expr = Expression("STRING")
        expr.e1 = self.getString()
    elif val == ord('('):
        expr = self.getExpression(is_cond)
        self.eatSpace()
        if self.peekChar() == ')':
            self.getChar()
            self.eatSpace()
        else:
            error("Missing ')'")
    elif val == ord('\''):
        if self.peekChar() == '{':
            self.getChar()
            expr = Expression("ARRAY")
            expr.args = self.parseArgList(is_cond)
            self.eatSpace()
            if self.peekChar() == '}':
                self.getChar()
                self.eatSpace()
            else:
                error("Missing '}'")
        else: # pragma: no cover
            self.ungetChar(chr(val))
            expr = Expression("NOTHING")
    elif val == ord('{'):
        expr = Expression("ARRAY")
        expr.args = self.parseArgList(is_cond)
        self.eatSpace()
        if self.peekChar() == '}':
            self.getChar()
            self.eatSpace()
        else:
            error("Missing '}'")
    elif val == ord('`'):
        # undefined macro
        self.lex() # drop `
        expr = Expression("NAME")
        expr.e1 = "`" + self.getString()

    else:
        self.ungetChar(chr(val))
        expr = Expression("NOTHING")
    while self.peekChar().isspace():
        self.getChar()
    return expr

  def parsePostfix(self, is_cond):
    expr = self.parsePrimary(is_cond)
    if expr.type == "NAME":
        ch = self.peekChar()
        if ch == '(':
            self.getChar()
            fname = expr.e1
            is_port_flow = False
            ch = self.peekChar()
            if fname in gAccessFuncs and ch == '<':
                self.getChar()
                is_port_flow = True
                pname = ""
                ch = self.peekChar()
                while ch != 0 and ch != '>' and ch != ')':
                    pname += self.getChar()
                    ch = self.peekChar()
                if ch == '>':
                    self.getChar()
                    if not pname in gPortnames:
                        error("%s(<%s>) is not a valid port flow access" % (fname, pname))
                else:
                    error("Missing '>' in port flow access")
                arg = Expression("PORT_FLOW")
                arg.e1 = pname
                args = [arg]
            elif ch == ')':
                args = []
            else:
                args = self.parseArgList(is_cond)
            expr = Expression("FUNCCALL")
            expr.e1 = fname
            expr.args = args
            if fname in gMathFunctions:
                nargs = gMathFunctions[fname]
                if nargs == 99:
                    warning("Function '%s' should not be used in a compact model" % fname)
                elif nargs != len(args):
                    error("Incorrect number of arguments (%d) for function %s (expect %d)"
                          % (len(args), fname, nargs))
                if fname == "log":
                    warning("Found base-10 log(); use ln() for natural log")
                elif (fname == "ln" or fname == "sqrt") and len(args) == 1:
                    arg = args[0]
                    if arg.type == "NAME" and arg.e1 in gParameters:
                        pname = arg.e1
                        check_cond = False
                        allow_zero = False
                        if fname == "ln" and not checkParamRangePos(pname, False):
                            check_cond = True
                        elif fname == "sqrt" and not checkParamRangePos(pname, True):
                            check_cond = True
                            allow_zero = True
                        if check_cond:
                            if not checkConditionsRequirePositive(pname, allow_zero, self.ternary_condition):
                                warning("Possible invalid function call for '%s(%s)'" % (fname, pname))
            elif fname in gUserFunctions:
                if gUserFunctions[fname].type == "integer":
                    expr.is_int = True
            elif fname == "white_noise":
                if len(args) == 1:
                    error("Missing name for white_noise()")
                elif len(args) > 2 or len(args) == 0:
                    error("Incorrect number of arguments (%d) for function %s (expect 1 or 2)"
                          % (len(args), fname))
            elif fname == "flicker_noise":
                if len(args) == 2:
                    error("Missing name for flicker_noise()")
                elif len(args) > 3 or len(args) < 2:
                    error("Incorrect number of arguments (%d) for function %s (expect 2 or 3)"
                           % (len(args), fname))
            elif fname == "ddt":
                if len(args) != 1:
                    error("Incorrect number of arguments (%d) for function %s (expect 1)"
                          % (len(args), fname))
            elif fname == "ddx":
                if len(args) == 2:
                    acc = args[1]
                    if acc.type != "FUNCCALL" or not acc.e1 in gAccessFuncs:
                        error("Second argument to %s must be a branch probe" % fname)
                    elif len(acc.args) > 1:
                        error("Second argument to %s must be a node potential or branch flow" % fname)
                else:
                    error("Incorrect number of arguments (%d) for function %s (expect 2)"
                          % (len(args), fname))
            elif fname == "$simparam":
                if len(args) == 1 or len(args) == 2:
                    name = args[0]
                    if name.type != "STRING":
                        error("First argument to %s must be a string" % fname)
                    elif name.e1 == "\"gmin\"":
                        if len(args) == 1: 
                            warning("$simparam(\"gmin\") call should provide default value 0 for second argument")
                        else:
                            defval = args[1]
                            if defval.type != "NUMBER" or defval.number != 0:
                                warning("$simparam(\"gmin\", %s) should use 0 for default value" % defval.asString())
                else:
                    error("Incorrect number of arguments (%d) for function %s (expect 1 or 2)"
                          % (len(args), fname))

            if self.peekChar() == ')':
                self.getChar()
                self.eatSpace()
            else:
                error("Missing ')' in function call")
        elif ch == '[':
            # bus select
            self.checkBusIndex(expr.e1)
        else:
            if ch == '+' or ch == '-':
                oper = self.peekOper()
                if oper == "++" or oper == "--":
                    # error printed by getOper(); just discard it
                    self.getOper()
            vname = expr.e1
            scope = getCurrentScope()
            found = False
            while not found:
                vn = scope + vname
                if vn in gVariables:
                    found = True
                    if gVariables[vn].type == "integer":
                        expr.is_int = True
                    break
                if scope != "":
                    scopes = scope.split(".")
                    scope = ""
                    if len(scopes) > 2:
                        for sc in scopes[0:-2]:
                            scope += sc + "."
                else:
                    break
            if not found and vname in gParameters:
                if gParameters[vname].type == "integer":
                    expr.is_int = True
    return expr

  def parsePrefix(self, is_cond):
    oper = self.peekOper()
    if oper == "-":               # Unary minus
        oper = self.getOper()
        arg = self.parsePrefix(is_cond)
        if arg.type == "NUMBER":
            expr = arg
            expr.number = -arg.number
        else:
            expr = Expression(oper)
            expr.e1 = arg
            expr.is_int = arg.is_int
    elif oper == "+":             # Unary plus
        oper = self.getOper()
        expr = self.parsePrefix(is_cond)
    elif oper == "!":             # Logical not
        oper = self.getOper()
        arg = self.parsePrefix(is_cond)
        expr = Expression(oper)
        expr.e1 = arg
    elif oper == "~":             # Bitwise not
        oper = self.getOper()
        warning("Unexpected bitwise operator '%s'" % oper)
        arg = self.parsePrefix(is_cond)
        expr = Expression(oper)
        expr.e1 = arg
    else:
        expr = self.parsePostfix(is_cond)
    return expr

  def parsePower(self, is_cond):
    expr = self.parsePrefix(is_cond)
    if( self.peekOper() == "**" ):
        oper = self.getOper()
        lhs = expr
        rhs = self.parsePower(is_cond)
        expr = Expression(oper)
        expr.e1 = lhs
        expr.e2 = rhs
        expr.is_int = lhs.is_int and rhs.is_int
    return expr

  def parseMultiplicative(self, is_cond):
    expr = self.parsePower(is_cond)
    while self.peekOper() in [ "*", "/", "%"]:
        oper = self.getOper()
        lhs = expr
        rhs = self.parsePower(is_cond)
        expr = Expression(oper)
        expr.e1 = lhs
        expr.e2 = rhs
        expr.is_int = lhs.is_int and rhs.is_int
        if oper == "/":
            if expr.is_int:
                warning("Integer divide")
            if expr.e2.type == "NAME" and expr.e2.e1 in gParameters:
                pname = expr.e2.e1
                if not checkParameterRangeExclude(pname, 0):
                    if not checkConditionsExclude(pname, 0, self.ternary_condition):
                        warning("Possible division by zero for parameter %s" % pname)
    return expr

  def parseAdditive(self, is_cond):
    expr = self.parseMultiplicative(is_cond)
    while self.peekOper() in [ "+", "-"]:
        oper = self.getOper()
        lhs = expr
        rhs = self.parseMultiplicative(is_cond)
        expr = Expression(oper)
        expr.e1 = lhs
        expr.e2 = rhs
        expr.is_int = lhs.is_int and rhs.is_int
    return expr

  def parseShift(self, is_cond):
    expr = self.parseAdditive(is_cond)
    while self.peekOper() in ["<<", ">>", "<<<", ">>>"]:
        oper = self.getOper()
        lhs = expr
        rhs = self.parseAdditive(is_cond)
        expr = Expression(oper)
        expr.e1 = lhs
        expr.e2 = rhs
    return expr

  def parseRelation(self, is_cond):
    expr = self.parseShift(is_cond)
    while self.peekOper() in ["<", ">", "<=", ">="]:
        oper = self.getOper()
        lhs = expr
        rhs = self.parseShift(is_cond)
        expr = Expression(oper)
        expr.e1 = lhs
        expr.e2 = rhs
    return expr

  def parseEquality(self, is_cond):
    expr = self.parseRelation(is_cond)
    oper = self.peekOper()
    if is_cond and oper == '=':
         error("Found '=' in condition; should be '=='")
         oper = "=="
    while oper in ["==", "!=", "===", "!=="]:
        oper = self.getOper()
        if is_cond and oper == '=':
            oper = "=="
        lhs = expr
        rhs = self.parseRelation(is_cond)
        expr = Expression(oper)
        expr.e1 = lhs
        expr.e2 = rhs
        oper = self.peekOper()
    return expr

  def parseBitOpers(self, is_cond):
    expr = self.parseEquality(is_cond)
    while self.peekOper() in ["^~", "~^", "~&", "~|"]:
        oper = self.getOper()
        lhs = expr
        rhs = self.parseEquality(is_cond)
        expr = Expression(oper)
        expr.e1 = lhs
        expr.e2 = rhs
    return expr

  def parseBitAnd(self, is_cond):
    expr = self.parseBitOpers(is_cond)
    while self.peekOper() == "&":
        oper = self.getOper()
        lhs = expr
        rhs = self.parseBitOpers(is_cond)
        warning("Unexpected bitwise operator '%s'" % oper)
        expr = Expression(oper)
        expr.e1 = lhs
        expr.e2 = rhs
    return expr

  def parseBitXor(self, is_cond):
    expr = self.parseBitAnd(is_cond)
    while self.peekOper() == "^":
        oper = self.getOper()
        lhs = expr
        rhs = self.parseBitAnd(is_cond)
        warning("Unexpected bitwise operator '%s'" % oper)
        expr = Expression(oper)
        expr.e1 = lhs
        expr.e2 = rhs
    return expr

  def parseBitOr(self, is_cond):
    expr = self.parseBitXor(is_cond)
    while self.peekOper() == "|":
        oper = self.getOper()
        lhs = expr
        rhs = self.parseBitXor(is_cond)
        warning("Unexpected bitwise operator '%s'" % oper)
        expr = Expression(oper)
        expr.e1 = lhs
        expr.e2 = rhs
    return expr

  def parseLogicalAnd(self, is_cond):
    expr = self.parseBitOr(is_cond)
    while self.peekOper() == "&&":
        oper = self.getOper()
        lhs = expr
        rhs = self.parseBitOr(is_cond)
        expr = Expression(oper)
        expr.e1 = lhs
        expr.e2 = rhs
    return expr

  def parseLogicalOr(self, is_cond):
    expr = self.parseLogicalAnd(is_cond)
    while self.peekOper() == "||":
        oper = self.getOper()
        lhs = expr
        rhs = self.parseLogicalAnd(is_cond)
        expr = Expression(oper)
        expr.e1 = lhs
        expr.e2 = rhs
    return expr

  def parseTernary(self, is_cond):
    cond = self.parseLogicalOr(is_cond)
    if( self.peekOper() == '?' ):
        self.ternary_condition = cond
        self.getOper()
        second = self.parseTernary(is_cond)
        has_colon = False
        if( self.peekOper() == ':' ):
            self.getOper()
            has_colon = True
            self.ternary_condition = Expression("!")
            self.ternary_condition.e1 = cond
        else:
            error("Missing ':' in ternary operator")
        third = self.parseTernary(is_cond)
        if third.type == "NOTHING":
            if has_colon:
                error("Missing expression after ':' in ternary operator")
            third = Expression("NUMBER")
            third.number = 0
        expr = Expression("?:")
        expr.e1 = cond
        expr.e2 = second
        expr.e3 = third
        expr.is_int = second.is_int and third.is_int
        self.ternary_condition = Expression("NOTHING")
    else:
        expr = cond
    return expr

  def getExpression(self, is_cond=False):
    self.eatSpace()
    return self.parseTernary(is_cond)

  def checkBusIndex(self, name):
    if self.peekChar() == '[':
        self.token = ord(self.getChar())
    if self.token == ord('['):
        if name in gPortnames:
            port = gPortnames[name]
        elif name in gNodenames:
            port = gNodenames[name]
        else:
            port = False
        if name in gVariables:
            var = gVariables[name]
        else:
            var = False
            sname = getCurrentScope() + name
            if sname in gVariables:
                var = gVariables[sname]
        if not port and not var:
            error("Identifier '%s' is not a vector port or array, cannot index" % name)
        sel = self.getExpression()
        if sel.type == "NUMBER":
            index = sel.number
            if port:
                if port.is_bus:
                    if (index > port.msb and index > port.lsb) \
                          or (index < port.msb and index < port.lsb):
                        error("Invalid index %d, port '%s' range is [%d:%d]" % (index, name, port.msb, port.lsb))
                else:
                    error("Port '%s' is not a bus, cannot index" % name)
            elif var:
                if len(var.range) == 2:
                    msb = var.range[0]
                    lsb = var.range[1]
                    if (index > msb and index > lsb) \
                          or (index < msb and index < lsb):
                        error("Invalid index %d, variable '%s' range is [%d:%d]" % (index, name, msb, lsb))
                elif var.range == []:
                    error("Identifier '%s' is not an array, cannot index" % name)
                else:
                    error("Unexpected range of variable '%s%s'" % (name, var.range))
        elif gVerbose:
            notice("Not validating vector index '%s[%s]'" % (name, sel.asString()))
        if self.peekChar() == ']':
            self.getChar()
            self.eatSpace()
        else:
            error("missing ']' after bus index")
    return True


  def getBusRange(self):
    bus_range = []
    msb = 0
    lsb = 0
    valid = True
    if valid:
        val = self.lex()
        if self.isNumber():
            msb = self.getNumber()
        elif self.isIdentifier():
            name = self.getString()
            if name in gParameters:
                gParameters[name].used = True
                defv = gParameters[name].defv
                if defv.type == "NUMBER":
                    msb = defv.number
                else:
                    warning("Cannot determine range msb from '%s'" % name)
            else:
                error("Expected constant expression for msb after '[', got '%s'" % name)
        else:
            error("Expected constant expression for msb after '['")
            valid = False
    if valid:
        val = self.lex()
        if val != ord(':'):
            error("Expected ':' in bus range")
            valid = False
    if valid:
        val = self.lex()
        if self.isNumber():
            lsb = self.getNumber()
        elif self.isIdentifier():
            name = self.getString()
            if name in gParameters:
                gParameters[name].used = True
                defv = gParameters[name].defv
                if defv.type == "NUMBER":
                    lsb = defv.number
                else:
                    warning("Cannot determine range lsb from '%s'" % name)
            else:
                error("Expected constant expression for lsb after ':', got '%s'" % name)
        else:
            error("Expected constant expression for lsb after ':'")
            if val == ord(']'):
                valid = False
    if valid:
        val = self.lex()
        if val == ord('-') or val == ord('+'):
            # inout [0:width-1] port;
            oper = chr(val)
            val = self.lex()
            if self.isNumber():
                num = self.getNumber()
                if oper == ord('-'):
                    lsb -= num
                elif oper == ord('+'):
                    lsb += num
                val = self.lex()
        if val != ord(']'): # pragma: no cover
            valid = False
    if valid:
        val = self.lex()
        bus_range = [msb,lsb]
    else:
        while val != ord(']'):
            val = self.lex()
            if val == 0: # pragma: no cover
                fatal("Parse error looking for ']'")
    return bus_range

# end of class Parser


################################################################################

# predefined macros
mac1 = Macro("__VAMS_ENABLE__", "1")
gMacros["__VAMS_ENABLE__"]           = mac1
mac2 = Macro("__VAMS_COMPACT_MODELING__", "1")
gMacros["__VAMS_COMPACT_MODELING__"] = mac2

# $temperature, $vt, $mfactor, $abstime
temp = Variable("$temperature", "real", False, "", 0)
temp.assign = 0
gVariables[temp.name] = temp
vt = Variable("$vt", "real", False, "", 0)
vt.assign = 0
gVariables[vt.name] = vt
mfac = Variable("$mfactor", "real", False, "", 0)
mfac.assign = 0
gVariables[mfac.name] = mfac
atime = Variable("$abstime", "real", False, "", 0)
atime.assign = -2
gVariables[atime.name] = atime


################################################################################
# functions

def fatal( message ): # pragma: no cover
    global gFileName, gLineNo
    if len(gLineNo) > 0:
        print("FATAL: File %s, line %d: %s" % (gFileName[-1], gLineNo[-1], message))
    else:
        print("FATAL: %s" % message)
    sys.exit(1)


def error( message ):
    global gErrorMsgDict, gMaxNum, gFileName, gLineNo
    count = 0
    if message in gErrorMsgDict:
        count = gErrorMsgDict[message]
    count += 1
    gErrorMsgDict[message] = count
    if count <= gMaxNum or gMaxNum == 0:
        if len(gLineNo) > 0:
            print("ERROR in file %s, line %d: %s" % (gFileName[-1], gLineNo[-1], message))
        else: # pragma: no cover
            print("ERROR: %s" % message)
        if count == gMaxNum:
            print("    Further errors of this type will be suppressed")


def warning( message, type=None ):
    global gWarningMsgDict, gMaxNum, gFileName, gLineNo
    if type is None:
        type = message
    count = 0
    if type in gWarningMsgDict:
        count = gWarningMsgDict[type]
    count += 1
    gWarningMsgDict[type] = count
    if count <= gMaxNum or gMaxNum == 0:
        if len(gLineNo) > 0:
            print("WARNING in file %s, line %d: %s" % (gFileName[-1], gLineNo[-1], message))
        else: # pragma: no cover
            print("WARNING: %s" % message)
        if count == gMaxNum:
            print("    Further warnings of this type will be suppressed")


def style( message, type=None ):
    global gStyleMsgDict, gMaxNum, gFileName, gLineNo
    if type is None:
        type = message
    count = 0
    if type in gStyleMsgDict:
        count = gStyleMsgDict[type]
    count += 1
    gStyleMsgDict[type] = count
    if count <= gMaxNum or gMaxNum == 0:
        if len(gLineNo) > 0:
            print("STYLE in file %s, line %d: %s" % (gFileName[-1], gLineNo[-1], message))
        else: # pragma: no cover
            print("STYLE: %s" % message)
        if count == gMaxNum:
            print("    Further style comments of this type will be suppressed")


def notice( message ):
    global gNoticeMsgDict, gMaxNum, gFileName, gLineNo
    count = 0
    if message in gNoticeMsgDict:
        count = gNoticeMsgDict[message]
    count += 1
    gNoticeMsgDict[message] = count
    if count <= gMaxNum or gMaxNum == 0:
        if len(gLineNo) > 0:
            print("NOTICE in file %s, line %d: %s" % (gFileName[-1], gLineNo[-1], message))
        else: # pragma: no cover
            print("NOTICE: %s" % message)
        if count == gMaxNum:
            print("    Further notices of this type will be suppressed")


def print_list( keys, start, maxlen ):
    outstr = start
    for item in sorted(keys):
        if outstr == start:
            outstr += item
        elif len(outstr) + len(item) + 2 > maxlen:
            print(outstr)
            outstr = start + item
        else:
            outstr += ", " + item
    print(outstr)


def format_char( chrnum ): # pragma: no cover
    if chrnum > 127:
        retstr = "non-ASCII character"
    elif chrnum < 32 or chrnum == 127:
        retstr = "non-printable character"
    elif chrnum == 34 or chrnum == 39 or chrnum == 96:
        # various quotation marks: " ' `
        retstr = ("character %s" % chr(chrnum))
    else:
        retstr = ("character '%s'" % chr(chrnum))
    return retstr


# check if parameter range excludes the value 'val'
def checkParameterRangeExclude( pname, val ):
    excludes = False
    parm = gParameters[pname]
    if val in parm.exclude:
        excludes = True
    if parm.range == "cc":
        if parm.min.type == "NUMBER" and parm.min.number > val:
            excludes = True
        elif parm.max.type == "NUMBER" and parm.max.number < val:
            excludes = True
    elif parm.range == "co":
        if parm.min.type == "NUMBER" and parm.min.number > val:
            excludes = True
        elif parm.max.type == "NUMBER" and parm.max.number <= val:
            excludes = True
    elif parm.range == "oc":
        if parm.min.type == "NUMBER" and parm.min.number >= val:
            excludes = True
        elif parm.max.type == "NUMBER" and parm.max.number < val:
            excludes = True
    elif parm.range == "oo":
        if parm.min.type == "NUMBER" and parm.min.number >= val:
            excludes = True
        elif parm.max.type == "NUMBER" and parm.max.number <= val:
            excludes = True
    # else "nb"
    return excludes


def checkParamRangePos( pname, allow_zero ):
    parm = gParameters[pname]
    if parm.range != "nb" and parm.min.type == "NUMBER":
        a = parm.min.number
        if parm.range == "cc" or parm.range == "co":
            # from [a:inf), a>0
            if a > 0 or (allow_zero and a == 0):
                return True
        elif parm.range == "oc" or parm.range == "oo":
            # from (a:inf), a>=0
            if a >= 0:
                return True
    return False


def checkParamRangeNeg( pname, allow_zero ):
    parm = gParameters[pname]
    if parm.range != "nb" and parm.max.type == "NUMBER":
        b = parm.max.number
        if parm.range == "cc" or parm.range == "oc":
            # from (-inf:b], b<0
            if b < 0 or (allow_zero and b == 0):
                return True
        elif parm.range == "oc" or parm.range == "oo":
            # from (-inf:b), b<=0
            if b <= 0:
                return True
    return False


# recursive!
def checkConditionsExcludeOper( oper, e1, e2, negated, pname, exval, recurse ):
    excludes = False
    if oper in ["!=", "==", ">=", ">", "<=", "<"]:
        val_eq = False
        val_ne = False
        val_gt = False
        val_ge = False
        val_lt = False
        val_le = False
        val = 0
        if e1.type == "NAME" and pname == e1.e1:
            val = e2
            if   (oper == "==" and not negated) or (oper == "!=" and negated):
                val_ne = True
            elif (oper == "!=" and not negated) or (oper == "==" and negated):
                val_eq = True
            elif (oper == ">"  and not negated) or (oper == "<=" and negated):
                val_ge = True
            elif (oper == "<"  and not negated) or (oper == ">=" and negated):
                val_le = True
            elif (oper == ">=" and not negated) or (oper == "<" and negated):
                val_gt = True
            elif (oper == "<=" and not negated) or (oper == ">" and negated):
                val_lt = True
        elif e2.type == "NAME" and pname == e2.e1:
            val = e1
            if   (oper == "==" and not negated) or (oper == "!=" and negated):
                val_ne = True
            elif (oper == "!=" and not negated) or (oper == "==" and negated):
                val_eq = True
            elif (oper == ">"  and not negated) or (oper == "<=" and negated):
                val_le = True
            elif (oper == "<"  and not negated) or (oper == ">=" and negated):
                val_ge = True
            elif (oper == ">=" and not negated) or (oper == "<" and negated):
                val_lt = True
            elif (oper == "<=" and not negated) or (oper == ">" and negated):
                val_gt = True
        if val:
            if   val_eq:
                if val.type == "NUMBER" and val.number == exval:
                    excludes = True
            elif val_ne:
                if val.type == "NUMBER" and val.number != exval:
                    excludes = True
            elif val_ge:
                if val.type == "NUMBER" and val.number >= exval:
                    excludes = True
                elif val.type == "NAME" and val.e1 in gParameters and exval == 0:
                    if checkParamRangePos(val.e1, True):
                        excludes = True
            elif val_le:
                if val.type == "NUMBER" and val.number <= exval:
                    excludes = True
                elif val.type == "NAME" and val.e1 in gParameters and exval == 0:
                    if checkParamRangeNeg(val.e1, True):
                        excludes = True
            elif val_gt:
                if val.type == "NUMBER" and val.number >  exval:
                    excludes = True
                elif val.type == "NAME" and val.e1 in gParameters and exval == 0:
                    if checkParamRangePos(val.e1, False):
                        excludes = True
            elif val_lt:
                if val.type == "NUMBER" and val.number <  exval:
                    excludes = True
                elif val.type == "NAME" and val.e1 in gParameters and exval == 0:
                    if checkParamRangeNeg(val.e1, False):
                        excludes = True

    elif oper == "&&" and not negated:
        c1 = checkConditionsExcludeOper(e1.type, e1.e1, e1.e2, negated, pname, exval, False)
        c2 = checkConditionsExcludeOper(e2.type, e2.e1, e2.e2, negated, pname, exval, False)
        excludes = c1 or c2
    elif oper == "||" and negated and recurse: # else of (c1 || c2)
        c1 = checkConditionsExcludeOper(e1.type, e1.e1, e1.e2, negated, pname, exval, False)
        c2 = checkConditionsExcludeOper(e2.type, e2.e1, e2.e2, negated, pname, exval, False)
        excludes = c1 or c2
    elif oper == "!" or oper == "NOT":
        excludes = checkConditionsExcludeOper(e1.type, e1.e1, e1.e2, not negated, pname, exval, True)
    return excludes


# check if conditions in effect exclude the value 'exval'
def checkConditionsExclude( pname, exval, extra_cond ):
    global gConditions
    excludes = False
    for cond in gConditions:
        if cond.type != "NOTHING":
            excludes |= checkConditionsExcludeOper(cond.type, cond.e1, cond.e2, False, pname, exval, True)
    if extra_cond.type != "NOTHING":
        # ternary condition
        excludes |= checkConditionsExcludeOper(extra_cond.type, extra_cond.e1, extra_cond.e2, False, pname, exval, True)
    return excludes


# recursive!
def checkConditionsRequirePosOper( oper, e1, e2, negated, pname, allow_zero, recurse ):
    requires = False
    if oper in [">=", ">", "<=", "<"]:
        val_gt = False
        val_ge = False
        val = 0
        if e1.type == "NAME" and pname == e1.e1:
            # PAR > val, PAR >= val
            val = e2
            if   (oper == ">"  and not negated) or (oper == "<=" and negated):
                val_ge = True
            elif (oper == ">=" and not negated) or (oper == "<" and negated):
                val_gt = True
        elif e2.type == "NAME" and pname == e2.e1:
            # val < PAR, val <= PAR
            val = e1
            if   (oper == "<"  and not negated) or (oper == ">=" and negated):
                val_ge = True
            elif (oper == "<=" and not negated) or (oper == ">" and negated):
                val_gt = True
        if val:
            if val_ge:
                if val.type == "NUMBER" and val.number >= 0:
                    requires = True
                elif val.type == "NAME" and val.e1 in gParameters:
                    if checkParamRangePos(val.e1, True):
                        requires = True
            elif val_gt:
                if val.type == "NUMBER" and (val.number > 0 or (allow_zero and val.number == 0)):
                    requires = True
                elif val.type == "NAME" and val.e1 in gParameters:
                    if checkParamRangePos(val.e1, allow_zero):
                        requires = True

    elif oper == "&&" and not negated:
        c1 = checkConditionsRequirePosOper(e1.type, e1.e1, e1.e2, negated, pname, allow_zero, False)
        c2 = checkConditionsRequirePosOper(e2.type, e2.e1, e2.e2, negated, pname, allow_zero, False)
        requires = c1 or c2
    elif oper == "||" and negated and recurse: # else of (c1 || c2)
        c1 = checkConditionsRequirePosOper(e1.type, e1.e1, e1.e2, negated, pname, allow_zero, False)
        c2 = checkConditionsRequirePosOper(e2.type, e2.e1, e2.e2, negated, pname, allow_zero, False)
        requires = c1 or c2
    elif oper == "!" or oper == "NOT":
        requires = checkConditionsRequirePosOper(e1.type, e1.e1, e1.e2, not negated, pname, allow_zero, True)
    return requires


# check if conditions in effect require positive parameter value
def checkConditionsRequirePositive( pname, allow_zero, extra_cond ):
    global gConditions
    requires = False
    for cond in gConditions:
        if cond.type != "NOTHING":
            requires |= checkConditionsRequirePosOper(cond.type, cond.e1, cond.e2, False, pname, allow_zero, True)
    if extra_cond.type != "NOTHING":
        # ternary condition
        requires |= checkConditionsRequirePosOper(extra_cond.type, extra_cond.e1, extra_cond.e2, False, pname, allow_zero, True)
    return requires


def checkIdentifierCollisions( basename, name, escaped, type ):
    valid = True
    if basename in gVAMSkeywords and not escaped:
        error("%s '%s' collides with Verilog-AMS keyword" % (type,name))
        valid = False
    elif name == gModuleName:
        error("%s '%s' collides with module name" % (type,name))
        valid = False
    elif name in gNatures and gNatures[name].defined:
        if type == "Nature":
            error("Duplicate declaration of nature '%s'" % name)
        else:
            error("%s '%s' collides with previously-declared nature" % (type,name))
        valid = False
    elif name in gDisciplines:
        if type == "Discipline":
            error("Duplicate declaration of discipline '%s'" % name)
        else:
            error("%s '%s' collides with previously-declared discipline" % (type,name))
        valid = False
    elif name in gParameters:
        if type == "Parameter":
            error("Duplicate declaration of parameter '%s'" % name)
        else:
            error("%s '%s' collides with previously-declared parameter" % (type,name))
        valid = False
    elif name in gVariables:
        if type == "Variable":
            error("Duplicate declaration of variable '%s'" % name)
        else:
            error("%s '%s' collides with previously-declared variable" % (type,name))
        valid = False
    elif name in gNodenames:
        if type == "Node name":
            error("Duplicate declaration of internal node '%s'" % name)
        else:
            error("%s '%s' collides with previously-declared node" % (type,name))
        valid = False
    elif name in gBranches:
        if type == "Branch":
            error("Duplicate declaration of branch '%s'" % name)
        else:
            error("%s '%s' collides with previously-declared branch" % (type,name))
        valid = False
    elif name in gUserFunctions:
        if type == "Function":
            error("Duplicate user-defined function '%s'" % name)
        else:
            error("%s '%s' collides with previously-defined function" % (type,name))
        valid = False
    elif name in gBlocknames:
        if type == "Block name":
            error("Duplicate named block '%s'" % name)
        else: # pragma: no cover
            error("%s '%s' collides with previously-found named block" % (type,name))
        valid = False
    elif name in gPortnames and type != "Port name":
        error("%s '%s' collides with previously-declared port (terminal)" % (type,name))
        valid = False
    elif type == "Port name" and not name in gPortnames:
        if len(gScopeList) == 1:
            error("Identifier '%s' is not a port (terminal), cannot set direction" % name)
        valid = False
    return valid


def getAttributes( line ):
    retval = []
    do_append = True
    start = line.find("(*")
    while start >= 0:
        stop = line.find("*)")
        semi = line.find(";")
        if semi >= 0 and semi < start:
            do_append = False
            error("Attribute found after ';', please add line break")
        if stop > start:
            attrib = line[start+2:stop]
            attrib = attrib.strip()
            line = line[:start] + line[stop+2:]
            line = line.strip()
            i = 0
            while i < len(attrib):
                while i < len(attrib) and attrib[i].isspace():
                    i += 1
                name = ""
                while i < len(attrib) and (attrib[i].isalpha() or attrib[i] == '_'):
                    name += attrib[i]
                    i += 1
                while i < len(attrib) and (attrib[i].isalnum() or attrib[i] == '_'):
                    name += attrib[i]
                    i += 1
                if do_append:
                    retval.append(name)
                while i < len(attrib) and attrib[i].isspace():
                    i += 1
                value = ""
                if i < len(attrib) and attrib[i] == '=':
                    # (* name = value *)
                    i += 1
                    while i < len(attrib) and attrib[i].isspace():
                        i += 1
                    in_quote = False
                    in_brace = 0
                    while i < len(attrib) and (in_quote or in_brace or attrib[i] != ','):
                        if attrib[i] == '"':
                            in_quote = not in_quote
                        if not in_quote:
                            if attrib[i] == '{':
                                in_brace += 1
                            elif attrib[i] == '}':
                                in_brace -= 1
                        value += attrib[i]
                        i += 1
                    value = value.strip()
                # else:
                #     (* name *)
                if do_append:
                    retval.append(value)
                if i < len(attrib):
                    if attrib[i] == ',':
                        i += 1
                    elif value == "":
                        if attrib[i].isalnum():
                            error("Missing '=' or ',' in attribute")
                        else:
                            error("Unexpected character '%s' in attribute" % attrib[i])
                            i += 1
        else: # pragma: no cover
            fatal("Malformed attribute")
        start = line.find("(*")
    retval.insert(0, line)
    return retval


def parseMacro( line ):
    global gVAMScompdir
    if line.startswith("`define"):
        i = len("`define")
        max = len(line)
        if i < max:
            if not line[i].isspace():
                error("Missing space after `define")
            while i < max and line[i].isspace():
                i += 1
            name = ""
            if i < max and not (line[i].isalpha() or line[i] == '_'):
                error("Macro name must start with a letter or underscore (_)")
            while i < max and (line[i].isalnum() or line[i] == '_'):
                name += line[i]
                i += 1
            if len(name) > 7 and name[0:7] == "__VAMS_":
                error("Macro name may not start with __VAMS_")
            args = []
            if i < max:
                if line[i] == '(':
                    while i < max and line[i] != ')':
                        i += 1
                        while i < max and line[i].isspace():
                            i += 1
                        arg = ""
                        if i < max and not (line[i].isalpha() or line[i] == '_'):
                            error("Macro formal arguments must start with a letter or underscore (_)")
                        while i < max and (line[i].isalnum() or line[i] == '_'):
                            arg += line[i]
                            i += 1
                        args.append(arg)
                        while i < max and line[i].isspace():
                            i += 1
                        if i < max and line[i] != ',' and line[i] != ')':
                            error("Macro formal arguments should be separated by a comma (,)")
                    if i< max and line[i] == ')':
                        i += 1
                    else:
                        error("Missing ')' for macro definition")
                j = i
                while j < max and line[j].isspace():
                    if line[j] == '\r' or line[j] == '\n':
                        i = j+1
                    j += 1
            if i < max:
                value = line[i:]
            else:
                value = ""
            #print("Define %s as %s" % (name, value))
            if name in gVAMScompdir:
                error("Invalid macro name '%s'" % name)
            else:
                if name in gMacros:
                    warning("Redefining macro '%s'" % name)
                mac = Macro(name, value)
                mac.args = args
                gMacros[name] = mac
        else:
            error("Missing macro name after `define")


def parseUndef( line ):
    if line.startswith("`undef"):
        token = "`undef"
        if line.startswith("`undefine"):
            token = "`undefine"
            error("Invalid compiler directive `undefine (use `undef)")
        i = len(token)
        max = len(line)
        if i < max:
            if not line[i].isspace():
                error("Missing space after %s" % token)
            while i < max and line[i].isspace():
                i += 1
            name = ""
            while i < max and (line[i].isalnum() or line[i] == '_'):
                name += line[i]
                i += 1
            if i < max:
                warning("Unexpected characters after macro name")
            if len(name) > 7 and name[0:7] == "__VAMS_":
                warning("`undef has no effect for pre-defined macros")
            elif name in gMacros:
                gMacros.pop(name)
            else:
                error("Macro '%s' is not defined" % name)
        else:
            error("Missing macro name after `undef")


def replaceFormalWithActual( text, args, actuals, check_collision ):
    #print("replaceFormalWithActual")
    #print("    %s" % text)
    #print(args)
    #print(actuals)
    # check if actual argument contains match for later formal argument
    if check_collision:
        collide = False
        for i in range(len(actuals)):
            act = actuals[i]
            for j in range(i+1, len(args)):
                arg = args[j]
                if act.find(arg) >= 0:
                    collide = True
                    #print("actual %d:%s" % (i,act))
                    #print("formal %d:%s" % (j,arg))
                    break
            if collide:
                break
        if collide:
            newargs = []
            for arg in args:
                arg = "ARGUMENT__" + arg
                newargs.append(arg)
            newtext = replaceFormalWithActual(text, args, newargs, False)
            text = newtext
            args = newargs
    for i in range(len(args)):
        arg = args[i]
        #print("Looking for '%s' in text, to replace with %s" % (arg,actuals[i]))
        in_quote = False
        word_bnd = True
        j = 0
        while j < len(text):
            ch = text[j]
            if ch == '"':
                in_quote = not in_quote
            elif not in_quote:
                if ch.isalnum() or ch == '_':
                    if word_bnd:
                        word_bnd = False
                        if len(text) >= j + len(arg) and text[j:j+len(arg)] == arg:
                           if len(text) == j + len(arg) or not (text[j+len(arg)].isalnum() or text[j+len(arg)] == '_'):
                               str = text[0:j]
                               str += actuals[i]
                               str += text[j+len(arg):]
                               text = str
                               j += len(actuals[i]) - 1
                else:
                    word_bnd = True
            j += 1
    #print("    %s" % text)
    return text


def findFirstUnquotedTic( line, start ):
    pos = -1
    i = start
    max = len(line)
    in_quote = False
    escaped = False
    while i < max:
        ch = line[i]
        if escaped:
            escaped = False
        elif ch == '\\':
            escaped = True
        elif ch == '"':
            in_quote = not in_quote
        elif ch == '`' and not in_quote:
            pos = i
            break
        i += 1
    return pos

def expandMacro( line ):
    start = 0
    pos = findFirstUnquotedTic(line, start)
    while pos >= 0:
        name = ""
        i = pos + 1
        while i < len(line) and (line[i].isalnum() or line[i] == '_'):
            name += line[i]
            i += 1
        if name in gMacros:
            text = gMacros[name].text
            args = gMacros[name].args
            if i < len(line) and line[i] == '(' and len(args) > 0:
                nparens = 1
                i += 1
                in_quote = False
                actuals = []
                arg = ""
                while i < len(line) and nparens > 0:
                    ch = line[i]
                    if ch == '"':
                        in_quote = not in_quote
                    if in_quote:
                        if ch == '\n':
                            ch = ' '
                        arg += ch
                    else:
                        if ch == '(':
                            nparens += 1
                            arg += ch
                        elif ch == ')':
                            nparens -= 1
                            if nparens == 0:
                                if arg != "":
                                    actuals.append(arg)
                                    arg = ""
                            else:
                                arg += ch
                        elif ch == ',' and nparens == 1:
                            actuals.append(arg)
                            arg = ""
                        elif ch == "\n":
                            arg += ' '
                        else:
                            arg += ch
                    i += 1
                if arg != "":
                    actuals.append(arg)
                if nparens > 0:
                    error("Missing ')' in macro call")
                if len(args) == len(actuals):
                    text = replaceFormalWithActual(text, args, actuals, True)
                else: # pragma: no cover
                    fatal("Macro `%s called with %d arguments, expected %d" % (name, len(actuals), len(args)))
            elif len(args) > 0:
                error("Macro `%s called without arguments, expected %d" % (name, len(args)))
            str = line[0:pos]
            str += text
            str += line[i:]
            line = str
        elif name in ["ifdef", "ifndef"]:
            ts = i
            while line[ts].isspace():
                ts += 1
            te = ts + 1
            while te < len(line) and not line[te].isspace():
                te += 1
            token = line[ts:te]
            t1 = line.find("`else")
            t2 = line.find("`endif")
            if t2 < 0: # pragma: no cover
                fatal("Found `%s in macro text without corresponding `endif" % name)
            else:
                i = t2 + 6
            rest = line[te:t2]
            if rest.find("`ifdef") > 0 or rest.find("`ifndef") > 0: # pragma: no cover
                fatal("Found nested `ifdef in macro text, not supported")
            if t1 >= 0:
                def_text = line[te:t1]
                not_text = line[t1+5:t2]
            else:
                def_text = line[te:t2]
                not_text = ""
            ifdef_text = ""
            if (token in gMacros and name == "ifdef") or (token not in gMacros and name == "ifndef"):
                ifdef_text = def_text
            else:
                ifdef_text = not_text
            str = line[0:pos]
            str += ifdef_text
            str += line[i:]
            line = str
        elif name in ["else", "endif"]: # pragma: no cover
            fatal("Found `%s in macro text without corresponding `ifdef" % name)
        elif name in gVAMScompdir:
            error("Can't handle compiler directive `%s in macro text" % name)
            start = i
        else:
            if gMissingConstantsFile != "" and name[0:2] in ["M_", "P_"]:
                error("Undefined macro `%s; please check %s" % (name, gMissingConstantsFile))
            else:
                error("Undefined macro `%s" % name)
            start = i
        pos = findFirstUnquotedTic(line,start)
    return line


def verifyEndScope( what ):
    global gScopeList

    endname = "end" + what
    scopename = what.upper()
    if len(gScopeList) == 1:
        if not gScopeList[0].startswith(scopename):
            found = "end" + gScopeList[0].lower()
            colon = found.find(":")
            if colon > 0:
                found = found[0:colon]
            error("Found '%s' but expected '%s'" % (endname, found))
        gScopeList.pop()
    elif len(gScopeList) == 0:
        error("Found '%s' without corresponding '%s'" % (endname, what))
    else: # pragma: no cover
        scope = getCurrentScope()
        error("Unexpected '%s' while still in scope %s" % (endname, scope))


def parseNatureDecl( line, defined ):
    global gScopeList

    parser = Parser(line)
    val = parser.lex()
    if not parser.isIdentifier() or parser.getString() != "nature": # pragma: no cover
        fatal("Parse error looking for 'nature'")

    name = ""
    escaped = False
    val = parser.lex()
    if parser.isIdentifier():
        name = parser.getString()
        escaped = parser.isEscaped()
    else: # pragma: no cover
        fatal("Parse error looking for nature name")

    if name != "":
        if defined:
            valid = checkIdentifierCollisions(name, name, escaped, "Nature")
            if valid:
                nature = Nature(name, defined)
                gNatures[name] = nature
            gScopeList.append("NATURE::" + name)
        elif not name in gNatures:
            nature = Nature(name, defined)
            gNatures[name] = nature


def parseNatureLine( line ):
    global gScopeList

    name = (gScopeList[0])[len("NATURE::"):]

    parser = Parser(line)
    val = parser.lex()
    attrib = ""
    if parser.isIdentifier():
        attrib = parser.getString()
    else: # pragma: no cover
        error("Parse error in nature '%s'" % name)

    if attrib == "endnature":
        if name in gNatures and gNatures[name].defined:
            if gNatures[name].access == "":
                warning("No access function specified for nature '%s'" % name)
            if gNatures[name].units == "":
                warning("No units specified for nature '%s'" % name)
        # else name collided with a different identifier
        verifyEndScope("nature")

    else:
        val = parser.lex()
        if val == ord('='):
            val = parser.lex()
        elif attrib.startswith("end"):
            verifyEndScope(attrib[3:])
            name = ""
        else:
            error("Missing '=' for '%s' in nature '%s'" % (attrib, name))

        if name in gNatures:
            if attrib == "units":
                if parser.isString():
                    if gNatures[name].units != "":
                        error("Duplicate specification for '%s' in nature '%s'" % (attrib, name))
                    gNatures[name].units = parser.getString()
                else:
                    error("Expected string value for units in nature '%s'" % name)

            elif attrib == "access":
                if parser.isIdentifier():
                    if gNatures[name].access != "":
                        error("Duplicate specification for '%s' in nature '%s'" % (attrib, name))
                    access = parser.getString()
                    gNatures[name].access = access
                    gAccessFuncs[access] = 1
                else:
                    error("Expected identifier for access function in nature '%s'" % name)

            elif attrib == "idt_nature":
                if parser.isIdentifier():
                    if gNatures[name].idt_nature != "":
                        error("Duplicate specification for '%s' in nature '%s'" % (attrib, name))
                    idt = parser.getString()
                    gNatures[name].idt_nature = idt
                    if not idt in gNatures:
                        error("In nature '%s', idt_nature references an unknown nature '%s'" % (name, idt))
                else:
                    error("Expected identifier for %s in nature '%s'" % (attrib, name))

        val = parser.lex()
        if name != "" and val != ord(';'):
            error("Missing ';' at end of item in nature '%s'" % name)
# end of parseNatureLine


def parseDisciplineDecl( line ):
    global gScopeList

    parser = Parser(line)
    val = parser.lex()
    if not parser.isIdentifier() or parser.getString() != "discipline": # pragma: no cover
        fatal("Parse error looking for 'discipline'")

    name = ""
    escaped = False
    val = parser.lex()
    if parser.isIdentifier():
        name = parser.getString()
        escaped = parser.isEscaped()
    else: # pragma: no cover
        fatal("Parse error looking for discipline name")

    valid = checkIdentifierCollisions(name, name, escaped, "Discipline")
    if valid:
        disc = Discipline(name)
        gDisciplines[name] = disc
    gScopeList.append("DISCIPLINE::"+name)


def parseDisciplineLine( line ):
    global gScopeList

    disc = (gScopeList[0])[len("DISCIPLINE::"):]

    parser = Parser(line)
    val = parser.lex()
    attrib = ""
    if parser.isIdentifier():
        attrib = parser.getString()
    else:
        error("Parse error in discipline '%s'" % disc)

    if attrib == "enddiscipline":
        if gVerbose:
            if disc in gDisciplines and gDisciplines[disc].domain == "":
                if gDisciplines[disc].potential == "":
                    warning("No potential nature specified for discipline '%s'" % disc)
                if gDisciplines[disc].flow == "":
                    warning("No flow nature specified for discipline '%s'" % disc)
        verifyEndScope("discipline")
    else:
        if attrib == "potential" or attrib == "flow":
            val = parser.lex()
            nature = ""
            if parser.isIdentifier() or parser.isString():
                nature = parser.getString()
                if not nature in gNatures:
                    error("Undefined nature '%s' for discipline '%s'" % (nature, disc))
            elif val == ord(';'):
                error("Missing nature for '%s' in discipline '%s'" % (attrib, disc))
            else: # pragma: no cover
                fatal("Missing nature for '%s' in discipline '%s'" % (attrib, disc))

            if disc in gDisciplines:
                if attrib == "potential":
                    if gDisciplines[disc].potential != "":
                        error("Duplicate specification of '%s' for discipline '%s'" % (attrib, disc))
                    gDisciplines[disc].potential = nature
                elif attrib == "flow":
                    if gDisciplines[disc].flow != "":
                        error("Duplicate specification of '%s' for discipline '%s'" % (attrib, disc))
                    gDisciplines[disc].flow = nature
                    if nature in gNatures:
                        units = gNatures[nature].units.strip('"')
                        if units != "" and not units in gUnitsMultiply:
                            gUnitsMultiply.append(units)
        elif attrib == "domain":
            val = parser.lex()
            if parser.isIdentifier() or parser.isString():
                value = parser.getString()
                if disc in gDisciplines:
                    gDisciplines[disc].domain = value
            else:
                error("Unexpected value for '%s' of discipline '%s'" % (attrib, disc))
        elif attrib.startswith("end"):
            verifyEndScope(attrib[3:])
            disc = ""
        else:
            val = parser.lex()
            if val != ord(';') and val != 0:
                error("Unexpected '%s' in discipline '%s'" % (attrib, disc))

        if val != ord(';'):
            val = parser.lex()
        if disc != "" and val != ord(';'):
            error("Missing ';' at end of item in discipline '%s'" % disc)
# end of parseDisciplineLine


def parseModuleDecl( line ):
    global gModuleName, gFileName, gLineNo
    valid = True

    parser = Parser(line)
    val = parser.lex()
    mtype = ""
    if parser.isIdentifier():
        mtype = parser.getString()
    if mtype != "module" and mtype != "macromodule": # pragma: no cover
        fatal("Parse error looking for 'module'")

    mname = ""
    val = parser.lex()
    if parser.isIdentifier():
        mname = parser.getString()
    else: # pragma: no cover
        fatal("Parse error looking for module name")

    if mname in gVAMSkeywords:
        error("Module name '%s' collides with Verilog-AMS keyword" % mname)
        valid = False
    else:
        if gModuleName != "":
            warning("File contains multiple modules: %s and %s" % (gModuleName, mname))
        gModuleName = mname

    val = parser.lex()
    if val == ord('('):
        val = parser.lex()
        while parser.isIdentifier():
            pname = parser.getString()
            if pname in gPortnames:
                error("Duplicate definition of port (terminal) '%s'" % pname)
            else:
                port = Port(pname, gFileName[-1], gLineNo[-1])
                gPortnames[pname] = port
            val = parser.lex()
            if val == ord(','):
                val = parser.lex()
            elif val != ord(')') and val != ord(';') and not parser.isIdentifier():
                error("Unexpected %s in port list" % format_char(val))
                val = parser.lex()
        if val == ord(')'):
            val = parser.lex()
        else:
            error("Missing ')' after list of ports in module declaration")
    if val != ord(';'):
        error("Missing ';' at end of module declaration")
# end of parseModuleDecl


def parseFunction( line ):
    global gScopeList, gCurrentFunction
    valid = True
    analog = False

    parser = Parser(line)
    val = parser.lex()
    if parser.isIdentifier() and parser.getString() == "analog":
        analog = True
        val = parser.lex()

    if not parser.isIdentifier() or parser.getString() != "function": # pragma: no cover
        fatal("Parse error looking for 'function'")

    if not analog:
        warning("Recommend 'analog' before 'function'")

    fname = ""
    ftype = ""
    escaped = False
    val = parser.lex()
    if parser.isIdentifier():
        ftype = parser.getString()

    val = parser.lex()
    if parser.isIdentifier():
        fname = parser.getString()
        escaped = parser.isEscaped()
    elif val == ord(';'):
        fname = ftype
        ftype = "real"
        warning("Function type not specified for '%s', assuming real" % fname)

    if len(gScopeList) == 0:
        error("Functions may not be defined outside a module")
        # dummy scope for error reporting
        gScopeList.append("FUNCTION::" + fname)
    elif len(gScopeList) > 1:
        error("Functions should be defined at module scope")

    if not ftype in ["real", "integer"]:
        error("Invalid type '%s' for function '%s'" % (ftype, fname))

    gScopeList.append("FUNCTION::" + fname)
    gCurrentFunction = Function(fname, ftype)
    valid = checkIdentifierCollisions(fname, fname, escaped, "Function")
    if valid:
        gUserFunctions[fname] = gCurrentFunction
        # implicit declaration of return variable
        scope = getCurrentScope()
        vname = scope + fname
        retval = Variable(vname, ftype, False, gFileName[-1], gLineNo[-1])
        gVariables[vname] = retval
# end of parseFunction


def parseParamDecl( line, attribs ):
    global gScopeList, gFileName, gLineNo
    in_module = True

    parm_or_alias = ""
    parser = Parser(line)
    val = parser.lex()
    if parser.isIdentifier():
        parm_or_alias = parser.getString()
    if parm_or_alias != "parameter" and parm_or_alias != "aliasparam" and parm_or_alias != "localparam": # pragma: no cover
        fatal("Parse error looking for 'parameter'")

    inst = False
    units = ""
    if len(attribs) > 0:
        for i in range(0, len(attribs)-1, 2):
            if attribs[i] == "type":
                val = attribs[i+1].lower()
                if val == "instance" or val == "\"instance\"":
                    inst = True
                elif val == "both" or val == "\"both\"":
                    inst = True
            elif attribs[i] == "units":
                units = attribs[i+1].strip('"').strip()

    ptype = ""
    if parm_or_alias == "aliasparam":
        ptype = "alias"
    else:
        val = parser.lex()
        if parser.isIdentifier():
            ptype = parser.getString()
        else: # pragma: no cover
            fatal("Parse error looking for parameter type")

    pname = ""
    escaped = False
    val = parser.lex()
    if parser.isIdentifier():
        pname = parser.getString()
        escaped = parser.isEscaped()
    elif val == ord('='):
        if ptype == "alias":
            error("Missing aliasparam name")
        elif ptype in ["real", "integer", "string"]:
            error("Missing parameter name")
        else:
            pname = ptype
            ptype = "real"
            # technically, type is determined by default value
            # but we prefer explicit declaration
            warning("Type not specified for parameter '%s', assuming real" % pname)

    parm = Parameter(pname, ptype, gFileName[-1], gLineNo[-1])
    parm.inst = inst
    parm.units = units

    if val != ord('='):
        val = parser.lex()
    if len(gScopeList) == 0:
        error("Parameter declared outside of module")
        in_module = False
    if parm_or_alias == "aliasparam":
        if parser.isIdentifier():
            if pname in ["real", "integer", "string"]:
                ptype = pname
                pname = parser.getString()
                error("Parameter type should not be specified for aliasparam '%s'" % pname)
            elif ptype != "":
                error("Parse error for aliasparam %s" % ptype)
            val = parser.lex()
    else:
        if not ptype in ["real", "integer", "string"]:
            error("Invalid type '%s' for parameter '%s'" % (ptype, pname))

    if pname != "":
        valid = checkIdentifierCollisions(pname, pname, escaped, "Parameter")
    else:
        valid = False
    # insert into global map
    if valid and in_module:
        gParameters[pname] = parm

    if val != ord('='):
        error("Incomplete %s specification" % parm_or_alias)
        return

    parm.defv = parser.getExpression()
    #print("parameter %s defv type %s" % (pname, parm.defv.type))

    if parm_or_alias == "aliasparam":
        if parm.defv.type == "NAME" and parm.defv.e1 in gParameters:
            parm.type = gParameters[parm.defv.e1].type
            gParameters[parm.defv.e1].used = True
        else:
            error("Default value for aliasparam '%s' must be a parameter" % pname)
        val = parser.lex()
        parm.used = True
        parm.is_alias = True

    else:
        if parm.defv.type == "NAME" and parm.defv.e1 in gParameters:
            dparm = gParameters[parm.defv.e1]
            if dparm.units != units:
                c1 = units.find("^")
                p1 = units.find(pname)
                c2 = dparm.units.find("^")
                p2 = dparm.units.find(dparm.name)
                if c1 > 0 and c2 > 0 and p1 > c1 and p2 > c2:
                    # look for RDDR with units V^(-PRDDR), default RSDR with units V^(-PRSDR)
                    tmp = units[0:p1] + dparm.name + units[p1+len(dparm.name):]
                    if tmp == dparm.units:
                        units = dparm.units
                elif c2 > 0 and p2 > c2 and pname.find(dparm.name) >= 0 and units.find(dparm.name) >= 0:
                    # PTWGLR with units m^PTWGLEXPR, PTWGL with units m^PTWGLEXP
                    units = dparm.units
            if dparm.units != units:
                if units == "":
                    u1 = "no units"
                else:
                    u1 = "units '" + units + "'"
                if dparm.units == "":
                    u2 = "no units"
                else:
                    u2 = "units '" + dparm.units + "'"
                notice("Parameter '%s' with %s has default '%s' with %s" % (pname, u1, dparm.name, u2))
        deps = parm.defv.getDependencies(False, False)
        for dep in deps:
            if dep == pname:
                error("Parameter '%s' default value may not depend on itself" % pname)
            elif dep in gParameters:
                gParameters[dep].used = True
            else:
                if dep in gVariables:
                    error("Parameter '%s' default value may not depend on variable '%s'" % (pname, dep))
                elif dep in gPortnames:
                    error("Parameter '%s' default value may not depend on port '%s'" % (pname, dep))
                elif dep in gNodenames:
                    error("Parameter '%s' default value may not depend on node '%s'" % (pname, dep))
                elif dep in gBranches:
                    error("Parameter '%s' default value may not depend on branch '%s'" % (pname, dep))
                else:
                    error("Parameter '%s' default value depends on unknown identifier '%s'" % (pname, dep))

        val = parser.lex()
        while parser.isIdentifier():
            str = parser.getString()
            if ptype == "string":
                val = parser.lex()
                if val == ord('\''):
                    val = parser.lex()
                if val == ord('{'):
                    val = parser.lex()
                    while parser.isString():
                        val = parser.lex()
                        if val == ord(','):
                            val = parser.lex()
                else:
                    error("Invalid %s in range" % format_char(val))
                    break
                if val == ord('}'):
                    val = parser.lex()
                else:
                    error("Invalid %s in range" % format_char(val))
                    break
            elif str == "from":
                val = parser.lex()
                lend = ""
                if val == ord('['):
                    lend = "c"
                elif val == ord('('):
                    lend = "o"
                elif val == ord(';'):
                    error("Missing range after 'from'")
                    break
                else:
                    error("Invalid %s in range" % format_char(val))
                parm.min = parser.getExpression()
                if parm.min.type == "NOTHING":
                    error("Missing lower bound in range")
                val = parser.lex()
                if val == ord(';'):
                    error("Incomplete range")
                    break
                elif val != ord(':'):
                    error("Invalid %s for range separator" % format_char(val))
                parm.max = parser.getExpression()
                if parm.max.type == "NOTHING":
                    error("Missing upper bound in range")
                val = parser.lex()
                rend = ""
                if val == ord(']'):
                    rend = "c"
                elif val == ord(')'):
                    rend = "o"
                elif val == ord(';'):
                    error("Missing ')' or ']' in range")
                    break
                else:
                    error("Invalid %s in range" % format_char(val))
                parm.range = lend + rend
                val = parser.lex()
            elif str == "exclude":
                exclude = parser.getExpression()
                if exclude.type == "NUMBER":
                    parm.exclude.append(exclude.number)
                elif exclude.type == "NAME":
                    if not exclude.e1 in gParameters:
                        warning("Parameter range 'exclude %s' is not valid" % exclude.e1)
                elif exclude.type == "NOTHING":
                    error("Missing value for 'exclude'")
                else:
                    warning("Unsupported exclusion of type %s" % exclude.type)
                val = parser.lex()
            else: # pragma: no cover
                fatal("Unexpected '%s' in parameter declaration" % str)

    if val != ord(';'):
        error("Missing ';' at end of %s declaration" % parm_or_alias)
# end of parseParamDecl


def parseVariableDecl( line, attribs ):
    global gScopeList, gCurrentFunction, gFileName, gLineNo
    in_module = True

    if len(gScopeList) == 0:
        error("Variable declared outside of module")
        in_module = False
    checkDeclarationContext("Variable")

    oppt = False
    units = ""
    multi = ""
    if len(attribs) > 0:
        for i in range(0, len(attribs), 2):
            if attribs[i] == "units":
                oppt = True
                units = attribs[i+1].strip('"')
            elif attribs[i] == "desc":
                oppt = True
            elif attribs[i] == "multiplicity":
                multi = attribs[i+1]

    vtype = ""
    vname = ""

    parser = Parser(line)
    val = parser.lex()
    while val != 0:
        err_msg = ""
        if vtype == "":
            if parser.isIdentifier():
                vtype = parser.getString()
                if not vtype in ["integer", "real", "genvar", "string"]:
                    err_msg = ("Unexpected variable type '%s'" % vtype)
                val = parser.lex()

        if parser.isIdentifier():
            str = parser.getString()
            escaped = parser.isEscaped()
            scope = getCurrentScope()
            vname = scope + str
            valid = checkIdentifierCollisions(str, vname, escaped, "Variable")
            # insert into global map
            if valid and in_module:
                var = Variable(vname, vtype, oppt, gFileName[-1], gLineNo[-1])
                if oppt:
                    var.used = True
                    if units in gUnitsMultiply:
                        if multi != "\"multiply\"":
                            warning("Operating-point variable %s should have multiplicity=\"multiply\"" % var.name, "multiplicity")
                    elif units in gUnitsDivide:
                        if multi != "\"divide\"":
                            warning("Operating-point variable %s should have multiplicity=\"divide\"" % var.name, "multiplicity")
                    elif multi != "":
                        warning("Unexpected multiplicity=%s for operating-point variable %s with units %s" % (multi, var.name, units), "multiplicity")
                gVariables[vname] = var
                if len(gScopeList) > 0 and gScopeList[-1].startswith("FUNCTION::"):
                    if gCurrentFunction:
                        if str in gCurrentFunction.inputs:
                            # hack, not really where it was assigned a value
                            var.assign = gLineNo[-1]
                        elif str in gCurrentFunction.outputs:
                            # marks the variable within the function as used
                            # (separate question whether the output is used outside the function)
                            var.used = True
            val = parser.lex()
            if val == ord('['):
                bus_range = parser.getBusRange()
                gVariables[vname].range = bus_range
                val = parser.token
            if val == ord('='):
                # real var = 1;
                error("Variable initializers not supported")
                defval = parser.getExpression()
                val = parser.lex()
        else:
            if vtype == "": # pragma: no cover
                fatal("Parse error looking for variable type")
            elif vtype in ["integer", "real", "genvar", "string"]:
                err_msg = ("Missing identifier after '%s' " % vtype)
            elif parser.isNumber():
                num = parser.getNumber()
                err_msg = ("Unexpected number %g in variable declaration" % num)
                val = parser.lex()
            else:
                err_msg = ("Missing variable type before '%s' " % vtype)
            if err_msg != "":
                error(err_msg)

        if val == ord(','):
            val = parser.lex()
        elif val == ord(';'):
            vtype = ""
            # clear info from attributes
            oppt = False
            units = ""
            multi = ""
            val = parser.lex()
        elif parser.isIdentifier():
            str = parser.getString()
            if str in ["inout", "input", "output", "real", "integer", "genvar"]:
                error("Missing ';' after list of variables")
                vtype = str;
                val = parser.lex()
            else:
                error("Missing ',' in list of variables")
        elif val != 0:
            error("Invalid %s in variable declaration" % format_char(val))
            val = parser.lex()
    if line[-1:] != ';':
        error("Missing ';' at end of variable declaration")
# end of parseVariableDecl


def parsePortDirection( line ):
    global gScopeList, gCurrentFunction
    in_module = True

    if len(gScopeList) == 0:
        error("Port declared outside of module")
        in_module = False

    ptype = ""
    pname = ""
    bus_range = []

    parser = Parser(line)
    val = parser.lex()
    while val != 0:
        if ptype == "":
            if parser.isIdentifier():
                ptype = parser.getString()
                if ptype == "input" or ptype == "output":
                    if len(gScopeList) == 1 and gScopeList[-1] == "MODULE":
                        warning("Port direction '%s' should be 'inout'" % ptype)
                elif ptype == "inout":
                    if len(gScopeList) > 0 and gScopeList[-1].startswith("FUNCTION::"):
                        warning("'inout' function arguments are not recommended")
                elif ptype == "real" or ptype == "integer":
                    # input x; real x; OK for functions
                    if gStyle and (len(gScopeList) == 0 or not gScopeList[-1].startswith("FUNCTION::")):
                        style("Variable declaration on same line as port declaration")
                else:
                    error("Invalid port direction '%s', expect 'inout'" % ptype)
                val = parser.lex()

        if val == ord('['):
            bus_range = parser.getBusRange()
            val = parser.token

        if parser.isIdentifier():
            pname = parser.getString()
            escaped = parser.isEscaped()
            if len(gScopeList) > 0 and gScopeList[-1].startswith("FUNCTION::"):
                if gCurrentFunction:
                    if pname in gCurrentFunction.args:
                        if ptype == "real" or ptype == "integer":
                            # input x; real x;
                            scope = getCurrentScope()
                            vname = scope + pname
                            valid = checkIdentifierCollisions(pname, vname, escaped, "Variable")
                            # insert into global map
                            if valid and in_module:
                                var = Variable(vname, ptype, False, gFileName[-1], gLineNo[-1])
                                gVariables[vname] = var
                                if pname in gCurrentFunction.inputs:
                                    # hack, not really where it was assigned a value
                                    var.assign = gLineNo[-1]
                                elif pname in gCurrentFunction.outputs:
                                    # marks the variable within the function as used
                                    # (separate question whether the output is used outside the function)
                                    var.used = True
                        else:
                            error("Function '%s' already had an argument '%s'"
                                  % (gCurrentFunction.name, pname))
                    else:
                        gCurrentFunction.args.append(pname)
                        if ptype == "input" or ptype == "inout":
                            gCurrentFunction.inputs.append(pname)
                            gCurrentFunction.outarg_used.append(True)
                        else:
                            gCurrentFunction.outputs.append(pname)
                            gCurrentFunction.outarg_used.append(False)
                else: # pragma: no cover
                    fatal("Programming error %s" % gScopeList[-1])
            elif ptype == "real" or ptype == "integer":
                scope = getCurrentScope()
                vname = scope + pname
                valid = checkIdentifierCollisions(pname, vname, escaped, "Variable")
                # insert into global map
                if valid and in_module:
                    var = Variable(vname, ptype, False, gFileName[-1], gLineNo[-1])
                    gVariables[vname] = var
            else:
                valid = checkIdentifierCollisions(pname, pname, escaped, "Port name")
                # set direction in global map
                if valid and in_module:
                    port = gPortnames[pname]
                    if port.direction != "":
                        error("Port (terminal) '%s' already had direction '%s'" % (pname, port.direction))
                    else:
                        port.direction = ptype
                    if len(bus_range) == 2:
                        port.is_bus = True
                        port.msb = bus_range[0]
                        port.lsb = bus_range[1]
            val = parser.lex()
        elif parser.isNumber():
            num = parser.getNumber()
            if ptype in ["real", "integer"]:
                error("Unexpected number %g in variable declaration" % num)
            else:
                error("Unexpected number %g in port declaration" % num)
            val = parser.lex()
        else:
            if ptype == "": # pragma: no cover
                fatal("Parse error looking for port direction")
            elif ptype in ["inout", "input", "output", "real", "integer"]:
                error("Missing identifier after '%s' " % ptype)
            else:
                error("Missing port direction before '%s' " % ptype)

        if val == ord(','):
            val = parser.lex()
        elif val == ord(';'):
            ptype = ""
            bus_range = []
            val = parser.lex()
        elif parser.isIdentifier():
            str = parser.getString()
            if str in ["inout", "input", "output", "real", "integer"]:
                error("Missing ';' after list of ports")
                ptype = str
                val = parser.lex()
            else:
                error("Missing ',' in list of ports")
        elif val != 0:
            error("Invalid %s in port declaration" % format_char(val))
            val = parser.lex()
    if line[-1:] != ';':
        error("Missing ';' at end of port direction declaration")
# end of parsePortDirection


def checkDeclarationContext( decl_type ):
    global gStatementInCurrentBlock
    if gStatementInCurrentBlock:
        if len(gScopeList) > 1 and gScopeList[-1] == "":
            error("%s declarations only allowed in named blocks" % decl_type)
        else:
            error("%s declarations must preceed all statements" % decl_type)


def registerContrib( afunc, type, node1, node2, cond_in, bias_dep_in ):
    [conds, bias_dep] = getCurrentConditions(cond_in, bias_dep_in)
    bname = ""
    bprint = node1
    if node2 == "":
        if node1 in gBranches:
            bname = node1
        else:
            bname = "UNNAMED_" + node1 + "_GND"
    else:
        bname = "UNNAMED_" + node1 + "_" + node2
        bprint = node1 + "," + node2
    if bname in gBranches:
        branch = gBranches[bname]
    else:
        branch = Branch(bname)
        branch.node1 = node1
        branch.node2 = node2
        gBranches[bname] = branch
    if type == "flow":
        if bias_dep:
            branch.lhs_flow = 2
        elif branch.lhs_flow == 0:
            branch.lhs_flow = 1
    else:
        if bias_dep:
            branch.lhs_pot = 2
        elif branch.lhs_pot == 0:
            branch.lhs_pot = 1
    if branch.lhs_flow + branch.lhs_pot > 2:
        warning("Switch branch (%s) with bias-dependent condition" % bprint)


def checkPorts():
    global gFileName, gLineNo
    for port in gPortnames.values():
        gFileName.append(port.declare[0])
        gLineNo.append(port.declare[1])
        if port.direction == "":
            error("No direction specified for port (terminal) '%s'" % port.name)
        if port.discipline == "":
            error("No discipline specified for port (terminal) '%s'" % port.name)
        if port.name.lower() != port.name:
            if port.name.upper() != port.name:
                warning("Mixed-case port name '%s'" % port.name)
            elif gStyle:
                style("Port names should be lower-case ('%s' not '%s')" % (port.name.lower(), port.name))
        gLineNo.pop()
        gFileName.pop()


def checkParmsAndVars():
    global gFileName, gLineNo
    temp = gVariables["$temperature"]
    has_tref = False
    has_dtemp = False
    all_lower = []
    all_upper = []
    not_lower = []
    not_upper = []
    param_names_lower = []
    for par in gParameters.values():
        gFileName.append(par.declare[0])
        gLineNo.append(par.declare[1])
        if not par.used:
            warning("Parameter '%s' was never used" % par.name)
        if par.name.lower() in ["tref", "tnom"]:
            has_tref = True
            if par.is_alias:
                defv = par.defv
                if defv.type == "NAME":
                    defv = defv.e1
                    if defv in gParameters:
                        par2 = gParameters[defv]
                        if par2.inst:
                            warning("Reference temperature parameter '%s' (aliased as %s) should not be an instance parameter" % (par2.name, par.name))
            elif par.inst:
                warning("Reference temperature parameter '%s' should not be an instance parameter" % par.name)
        if par.name.lower() in ["dtemp", "trise"]:
            has_dtemp = True
            if par.is_alias:
                defv = par.defv
                if defv.type == "NAME":
                    defv = defv.e1
                    if defv in gParameters:
                        par2 = gParameters[defv]
                        if not par2.inst:
                            warning("Temperature offset parameter '%s' (aliased as %s) should be an instance parameter" % (par2.name, par.name))
            elif not par.inst:
                warning("Temperature offset parameter '%s' should be an instance parameter" % par.name)
        param_names_lower.append(par.name.lower())
        if par.name.lower() == par.name:
            all_lower.append(par.name)
        else:
            not_lower.append(par.name)
        if par.name.upper() == par.name:
            all_upper.append(par.name)
        else:
            not_upper.append(par.name)
        gLineNo.pop()
        gFileName.pop()
    if temp.used > 0:
        if not has_tref:
            warning("Module uses $temperature but does not have 'tref'")
        if not has_dtemp:
            warning("Module uses $temperature but does not have 'dtemp'")
    if len(not_lower) > 0:
        if len(not_upper) == 0:
            if gStyle:
                style("Parameters should be declared in lower-case")
        else:
            if len(not_lower) == 1 and len(all_lower) > 5:
                bad = not_lower[0]
                warning("Parameters have inconsistent case (all lower-case except '%s')" % bad)
            elif len(not_upper) == 1 and len(all_upper) > 5:
                bad = not_upper[0]
                warning("Parameters have inconsistent case (all upper-case except '%s')" % bad)
            else:
                warning("Parameters have inconsistent case (prefer all lower-case)")

    for var in gVariables.values():
        if var.name.startswith("FUNCTION::"):
            var.assign = 0 # suppress printing in summary
        else:
            if var.funcNameAndArg[0] != "":
                # was assigned as a function output variable
                if var.used:
                    fname = var.funcNameAndArg[0]
                    argno = var.funcNameAndArg[1]
                    if fname in gUserFunctions:
                        func = gUserFunctions[fname]
                        if not func.outarg_used[argno]:
                            func.outarg_used[argno] = True
                            if gDebug:
                                print("    Function '%s' output argument '%s' is used because variable '%s' is used"
                                      % (fname,func.args[argno],var.name))
                    else: # pragma: no cover
                        fatal("Programming error: function arguments")
                else:
                    var.used = True # don't report as unused
            if var.assign == -1:
                gFileName.append(var.declare[0])
                gLineNo.append(var.declare[1])
                if var.oppt:
                    error("Operating-point variable '%s' was never assigned a value" % var.name)
                elif var.used:
                    warning("Variable '%s' was never assigned a value" % var.name)
                else:
                    warning("Variable '%s' was never set and never used" % var.name)
                gLineNo.pop()
                gFileName.pop()
            elif not var.used and var.assign > 0:
                warning("Variable '%s' was never used" % var.name)
            elif var.used and var.name == "$vt":
                warning("$vt not recommended; value differs slightly between simulators")
            if var.oppt and var.name.lower() in param_names_lower:
                parname = ""
                for par in gParameters.values():
                    if par.name.lower() == var.name.lower():
                        parname = par.name
                        break
                warning("Operating-point variable '%s' differs only in case from parameter '%s'" % (var.name, parname))
    for func in gUserFunctions.values():
        if func.used:
            if len(func.outputs) > 0:
                for i in range(len(func.args)):
                    if not func.outarg_used[i]:
                        warning("Function '%s' output argument '%s' is never used"
                                % (func.name, func.args[i]) )
        elif gVerbose:
            notice("User-defined function '%s' is never used" % func.name)
# end of check checkParmsAndVars


def getIdentifierType( name ):
    type = ""
    scope = getCurrentScope()
    vname = scope + name
    if vname in gVariables:
        type = gVariables[vname].type
    elif name in gVariables:
        type = gVariables[name].type
    elif name in gParameters:
        type = gParameters[name].type
    return type


def checkDependencies( deplist, unset_message, bias_dep_in, is_condition, do_warn ):
    [conds, bias_dep] = getCurrentConditions(Expression("NOTHING"), bias_dep_in)
    biases = []
    if is_condition:
        # condition for an if, for, while, case statement
        # don't use bias_dep from getCurrentConditions
        bias_dep = bias_dep_in
    for dep in deplist:
        if dep.find("(") > 0:
            biases.append(dep)
            if bias_dep < 2:
                bias_dep = 2
            is_set = True
            found = True
            scope = ""
        elif dep == "ddx":
            bias_dep = 4
            is_set = True
            found = True
            scope = ""
        else:
            is_set = False
            found  = False
            scope = getCurrentScope()
        while not found:
            vn = scope + dep
            if vn in gVariables:
                found = True
                gVariables[vn].used = True
                if bias_dep < gVariables[vn].bias_dep:
                    bias_dep = gVariables[vn].bias_dep
                for bias in gVariables[vn].biases:
                    if not bias in biases:
                        biases.append(bias)
                assign = gVariables[vn].assign
                if assign >= 0:
                    is_set = True
                elif assign <= -2:
                    is_set = True
                    if assign == -2 and not unset_message.startswith("Task") and do_warn:
                        warning("%s %s, which should not appear in a compact model" % (unset_message, dep))
                if len(gVariables[vn].conditions) > 0:
                    if not conditionCovered(conds, gVariables[vn].conditions):
                        if do_warn and gDebug:
                            if conds == "":
                                print("Hidden state: '%s' used (always)" % vn)
                            else:
                                used_conds = reformatConditions(conds)
                                print("Hidden state: '%s' used if %s" % (vn, used_conds))
                            assigned_conds = reformatConditions(gVariables[vn].conditions)
                            print("    assigned under conditions: %s" % assigned_conds)
                        if not vn in gHiddenState:
                            gHiddenState[vn] = vn
                break
            if scope != "":
                scopes = scope.split(".")
                scope = ""
                if len(scopes) > 2:
                    for sc in scopes[0:-2]:
                        scope += sc + "."
            else:
                break
        if not found:
            if dep in gParameters:
                if gParameters[dep].is_alias:
                    par = gParameters[dep].defv.e1
                    if par in gParameters:
                        error("Reference to aliasparam '%s' should use parameter '%s'" % (dep,par))
                    else:
                        error("Reference to aliasparam '%s'" % dep)
                gParameters[dep].used = True
                is_set = True
            elif dep in gPortnames or dep in gNodenames:
                if bias_dep == 0:
                    bias_dep = 1
                if do_warn:
                    error("Reference to node '%s' without access function" % dep)
                is_set = True
            elif dep in gBranches:
                if bias_dep == 0:
                    bias_dep = 1
                if do_warn:
                    error("Reference to branch '%s' without access function" % dep)
                is_set = True
        if not is_set and do_warn:
            if dep in gVAMSkeywords:
                error("%s Verilog-AMS keyword '%s'" % (unset_message, dep))
            # $xposition, $yposition, $angle, $hflip, $vflip
            elif dep in gVAMShiersysparm:
                warning("%s hierarchical system parameter '%s'" % (unset_message, dep))
            else:
                error("%s %s, which has not been set" % (unset_message, dep))
    return [bias_dep, biases]
# end of checkDependencies


# split on "&&", but watch out for parentheses
def splitConditions( str ):
    retval = []
    num_parens = 0
    cond = ""
    i = 0
    while  i < len(str):
        ch = str[i]
        if ch == '(':
            num_parens += 1
        elif ch == ')':
            num_parens -= 1
        if num_parens == 0 and i+1 < len(str) and ch == '&' and str[i+1] == '&':
            retval.append(cond)
            cond = ""
            i += 2
        else:
            cond += ch
            i += 1
    if cond != "":
        retval.append(cond)
    if num_parens != 0: # pragma: no cover
        fatal("Programming error in splitConditions")
    return retval


# variable used under condition 'test'
# see if it was set under the same or less restrictive conditions
def conditionCovered( test, condlist ):
    retval = False
    if test in condlist:
        retval = True
    else:
        for cond in condlist:
            newparts = splitConditions(test)
            oldparts = splitConditions(cond)
            both = []
            newonly = []
            oldonly = []
            for op in oldparts:
                if op in newparts:
                    both.append(op)
                else:
                    oldonly.append(op)
            for op in newparts:
                if op in oldparts:
                    if not op in both: # pragma: no cover
                        fatal("Programming error in conditionCovered")
                else:
                    newonly.append(op)
            if len(oldonly) == 0:
                #print("cond      %s\nsubset of %s" % (test,cond))
                retval = True
                break
    return retval


def mergeConditions( oldlist, new ):
    retlist = []
    add_new = True
    did_merge = False
    for old in oldlist:
        if new == "NOT(" + old + ")":
            # A and !(A)
            add_new = False
            retlist = []
            break
        elif new == old:
            # already have this condition
            add_new = False
            retlist.append(old)
        else:
            oldparts = splitConditions(old)
            newparts = splitConditions(new)
            both = []
            oldonly = []
            newonly = []
            for op in oldparts:
                if op in newparts:
                    both.append(op)
                else:
                    oldonly.append(op)
            for op in newparts:
                if op in oldparts:
                    if not op in both: # pragma: no cover
                        fatal("Programming error in mergeConditions")
                else:
                    newonly.append(op)
            replace = False
            retval = ""
            if len(newonly) == 0:
                # dropped some number of old conditions
                replace = True
                retval = new
                add_new = False
            elif len(oldonly) == 0:
                # new condition is more restrictive
                add_new = False
            elif len(oldonly) == 1 and len(newonly) == 1 and \
                    (oldonly[0] == "NOT(" + newonly[0] + ")" or \
                     newonly[0] == "NOT(" + oldonly[0] + ")"):
                # A..&&B and A..&&!(B) -> A..
                # A..&&!(B) and A..&&B -> A..
                for op in both:
                    if retval != "":
                        retval += "&&"
                    retval += op
                replace = True
                add_new = False
                did_merge = True
            elif len(newonly) == 1:
                olds = oldonly[0]
                for ol in oldonly[1:]:
                    olds += "&&" + ol
                if newonly[0] == "NOT(" + olds + ")":
                    # A&&B and !(A&&B)
                    # A&&B&&C and !(A&&B&&C)
                    for op in both:
                        if retval != "":
                            retval += "&&"
                        retval += op
                    replace = True
                    add_new = False
                    did_merge = True
            if replace:
                retlist.append(retval)
            else:
                retlist.append(old)
    if add_new:
        retlist.append(new)
    elif did_merge and len(retlist) > 1:
        retlist = mergeConditions(retlist[:-1],retlist[-1])
    if "" in retlist:
        retlist = []
    return retlist


# replace NOT(x) with !(x)
def reformatConditions( conds ):
    if isinstance(conds, list):
        retlist = []
        for cond in conds:
            cond = reformatConditions(cond)
            retlist.append(cond)
        return retlist
    else:
        cond = conds
        idx = cond.find("NOT(")
        while idx >= 0:
            if (idx > 0):
                cond = cond[0:idx] + "!(" + cond[idx+4:]
            else:
                cond = "!(" + cond[idx+4:]
            idx = cond.find("NOT(")
        return cond


def markVariableAsSet( varname, bias_dep_in, cond_in, biases_in, for_var, set_by_func, fname, argnum ):
    global gLineNo

    found = False
    scope = getCurrentScope()
    while not found:
        vn = scope + varname
        if vn in gVariables:
            found = True
            [conds, bias_dep] = getCurrentConditions(cond_in, bias_dep_in)
            biases = biases_in
            if gVariables[vn].type == "integer":
                if bias_dep > 1:
                    bias_dep = 1
                    biases = []
            if gVariables[vn].assign < 0:
                # first assignment
                gVariables[vn].assign = gLineNo[-1]
                if conds != "":
                    gVariables[vn].conditions = [conds]
                gVariables[vn].bias_dep = bias_dep
                gVariables[vn].biases = biases
            elif len(gVariables[vn].conditions) > 0:
                if conds == "":
                    # no conditions for this assignment
                    gVariables[vn].conditions = []
                    gVariables[vn].bias_dep = bias_dep
                    gVariables[vn].biases   = biases
                else:
                    oldconds = gVariables[vn].conditions
                    gVariables[vn].conditions = mergeConditions(oldconds, conds)
            if conds == "" or gVariables[vn].bias_dep < bias_dep:
                gVariables[vn].bias_dep = bias_dep
            if conds == "": # replace previous biases
                gVariables[vn].biases   = biases
            else: # append these biases, possibly overkill
                for bias in biases:
                    if not bias in gVariables[vn].biases:
                        gVariables[vn].biases.append(bias)
            if set_by_func:
                gVariables[vn].funcNameAndArg = [fname, argnum]
            break
        if scope != "":
            scopes = scope.split(".")
            scope = ""
            if len(scopes) > 2:
                for sc in scopes[0:-2]:
                    scope += sc + "."
        else:
            break

    if not found:
        if for_var:
            if varname in gParameters:
                error("Parameter '%s' cannot be used as a for loop variable" % varname)
            else:
                error("Identifier '%s' is not valid as a for loop variable" % varname)
        elif set_by_func:
            if varname in gParameters:
                error("Parameter '%s' cannot be used as an output argument of function" % varname)
            else:
                error("Identifier '%s' is not valid as an output argument of function" % varname)
        else:
            error("Assignment to '%s' which is not a variable" % varname)
    return found
# end of markVariableAsSet


# try to parse term as something like L{pname} * Inv_L
def checkBinningTerm(t_in, pname):
    term = deepcopy(t_in) # term may be modified below
    fact = ""
    pfx  = ""
    sfx  = ""
    oper = "*"
    if term.type == "*":
        # consider Inv_L * Inv_W * P{param}
        if term.e1.type == "*" and term.e2.type == "NAME":
            if term.e1.e1.type == "NAME" and term.e1.e2.type == "NAME":
                if term.e2.e1 in gParameters:
                    term.e1.type = "NAME"
                    term.e1.e1 = term.e1.e1.e1 + "*" + term.e1.e2.e1
                elif term.e1.e1.e1 in gParameters:
                    sparm = term.e1.e1.e1
                    term.e1.type = "NAME"
                    term.e1.e1 = term.e1.e2.e1 + "*" + term.e2.e1
                    term.e2.e1 = sparm
                elif term.e1.e2.e1 in gParameters:
                    sparm = term.e1.e2.e1
                    term.e1.type = "NAME"
                    term.e1.e1 = term.e1.e1.e1 + "*" + term.e2.e1
                    term.e2.e1 = sparm
        elif term.e2.type == "*" and term.e1.type == "NAME":
            if term.e2.e1.type == "NAME" and term.e2.e2.type == "NAME":
                if term.e1.e1 in gParameters:
                    term.e2.type = "NAME"
                    term.e2.e1 = term.e2.e1.e1 + "*" + term.e2.e2.e1
                elif term.e2.e1.e1 in gParameters:
                    sparm = term.e2.e1.e1
                    term.e2.type = "NAME"
                    term.e2.e1 = term.e2.e2.e1 + "*" + term.e1.e1
                    term.e1.e1 = sparm
                elif term.e2.e2.e1 in gParameters:
                    sparm = term.e2.e2.e1
                    term.e2.type = "NAME"
                    term.e2.e1 = term.e2.e1.e1 + "*" + term.e1.e1
                    term.e1.e1 = sparm
    if term.type == "/":
        # consider 1.0e-6 * LXL / L
        if term.e1.type == "*" and term.e2.type == "NAME" and (term.e2.e1 == "L" or term.e2.e1 == "W"):
            if term.e1.e1.type == "NUMBER" and term.e1.e2.type == "NAME" and term.e1.e2.e1 in gParameters:
                # convert to 1.0e-6 / L * LXL
                term.type = "*"
                term.e1.type = "/"
                lparm = term.e2 # L or W
                term.e2 = term.e1.e2
                term.e1.e2 = lparm
        # consider PXL * 1.0e-6 / (L * NFIN)
        if term.e1.type == "*" and term.e2.type == "*":
            if term.e1.e1.type == "NAME" and term.e1.e1.e1 in gParameters and term.e1.e2.type == "NUMBER":
                if term.e2.e1.type == "NAME" and term.e2.e2.type == "NAME":
                    # convert to PXL * (1.0e-6 / (L * NFIN))
                    newfact = Expression("/")
                    newfact.e1 = term.e1.e2
                    newfact.e2 = term.e2
                    term.type = "*"
                    term.e1 = term.e1.e1
                    term.e2 = newfact
    if (term.type == "*" or term.type == "/") and (term.e1.type == "NAME" or term.e2.type == "NAME"):
        sname = ""
        oper = term.type
        if term.e1.type == "NAME" and term.e1.e1 in gParameters:
            sname = term.e1.e1
            fact = term.e2.asString()
        elif term.e2.type == "NAME" and term.e2.e1 in gParameters:
            sname = term.e2.e1
            fact = term.e1.asString()
        elif term.e2.type == "NUMBER" and term.e2.e1 == 0.0:
            # deliberately zero binning term
            pfx = "0.0"
        if sname != "":
            if sname in gParameters:
                spar = gParameters[sname]
                if spar.defv:
                    if spar.defv.type == "NUMBER":
                        if spar.defv.e1 != 0:
                            warning("Non-zero default value (%e) for binning parameter '%s'" % (spar.defv.e1, sname))
                    elif spar.defv.type == "NAME" and spar.defv.e1 in gParameters:
                        dpar = gParameters[spar.defv.e1]
                        if sname[0] != dpar.name[0] or (sname[0:2] == "P2" and dpar.name[0:2] != "P2"):
                            # WPAR default set from PPAR
                            warning("Unexpected default '%s' for binning parameter '%s'" % (dpar.name, sname))
                        while dpar.defv.type == "NAME" and dpar.defv.e1 in gParameters:
                            dpar = gParameters[dpar.defv.e1]
                        if dpar.defv.type != "NUMBER" or dpar.defv.e1 != 0:
                            warning("Non-zero default value (%s) for binning parameter '%s'" % (spar.defv.e1, sname))
                    else:
                        warning("Unexpected default value (type %s) for binning parameter '%s'" % (sname, spar.defv.type))
            if len(sname) > len(pname) and pname == sname[len(sname)-len(pname):]:
                # LPAR, WPAR, PPAR, P2PAR, etc.
                pfx = sname[0:len(sname)-len(pname)]
            elif (pname[-1].lower() == 'o' or pname[-1] == '0') \
                and pname[:-1] == sname[0:len(pname)-1]:
                # PARO and PARL, PARWL, etc.
                sfx = sname[len(pname)-1:]
            elif pname[0:2].lower() == 'po' and pname[2:] == sname[len(sname)-len(pname)+2:]:
                # POPAR and LPAR, WPAR, etc.
                pfx = sname[0:2]
    return [fact, oper, pfx, sfx]
# end of checkBinningTerm

# keep track of binning equation patterns
gBinningPatterns = [[],[]]

def checkBinning(vname, t0, t1, t2, t3, t4, t5, line):
    global gBinningPatterns, gFileName, gLineNo
    if t0.type == "NAME" and t0.e1 in gParameters:
        pname = t0.e1
        pname_l = pname.lower()
        plen = len(pname)
        vname_l = vname.lower()
        vlen = len(vname)
        vpfx = ""
        vsfx = ""
        pfx0 = ""
        sfx0 = ""
        if vlen > plen:
            if vname_l[0:plen] == pname_l:
                # VFB_i = VFB + ...
                vsfx = vname[plen:]
            elif vname_l[vlen-plen:] == pname_l:
                # pParam_VFB = VFB + ...
                vpfx = vname[0:vlen-plen]
            elif (pname_l[-1] == 'o' or pname_l[-1] == '0') \
                and pname_l[:-1] == vname_l[0:plen-1]:
                # VFB_i = VFBO + ...
                vsfx = vname[plen-1:]
                sfx0 = pname[-1]
        elif vlen < plen and (pname_l[-1] == 'o' or pname_l[-1] == '0') \
            and pname_l[:-1] == vname_l:
            # KUOWE = KUOWEO + ...
            vpfx = "NONE"
            sfx0 = pname[-1]
        elif pname_l[0:2] == 'po' and pname_l[2:] == vname_l[0:plen-2]:
            vsfx = vname[plen-2:]
            pfx0 = pname[0:2]
            if vsfx == "":
                vsfx = "NONE"
        [fact1, op1, pfx1, sfx1] = checkBinningTerm(t1, pname)
        [fact2, op2, pfx2, sfx2] = checkBinningTerm(t2, pname)
        [fact3, op3, pfx3, sfx3] = checkBinningTerm(t3, pname)
        if t4:
            [fact4, op4, pfx4, sfx4] = checkBinningTerm(t4, pname)
        else:
            [fact4, op4, pfx4, sfx4] = ["", "", "", ""]
        if t5:
            [fact5, op5, pfx5, sfx5] = checkBinningTerm(t5, pname)
        else:
            [fact5, op5, pfx5, sfx5] = ["", "", "", ""]
        warned = False
        if fact1 != "" and fact2 != "" and fact3 != "":
            if pfx1+sfx1 != "" and pfx2+sfx2 != "" and pfx3+sfx3 != "":
                if gBinning and vpfx == "" and vsfx == "":
                    warning("Did not find parameter name '%s' in binning variable name '%s'" % (pname, vname))
                    warned = True
                if t4 and pfx4+sfx4 == "":
                    warning("Binning equation to set '%s' involves '%s'" % (vname, t4.asString()))
                    warned = True
                if t5 and pfx5+sfx5 == "":
                    warning("Binning equation to set '%s' involves '%s'" % (vname, t5.asString()))
                    warned = True
            elif pfx1+sfx1 == "" and (pfx2+sfx2 != "" or pfx3+sfx3 != ""):
                warning("Binning equation to set '%s' involves '%s'" % (vname, t1.asString()))
                warned = True
            elif pfx2+sfx2 == "" and (pfx1+sfx1 != "" or pfx3+sfx3 != ""):
                warning("Binning equation to set '%s' involves '%s'" % (vname, t2.asString()))
                warned = True
            elif pfx3+sfx3 == "" and (pfx1+sfx1 != "" or pfx2+sfx2 != ""):
                warning("Binning equation to set '%s' involves '%s'" % (vname, t3.asString()))
                warned = True
        new_pat = [vpfx+vsfx, pfx0+sfx0, fact1+op1, pfx1+sfx1, fact2+op2, pfx2+sfx2, fact3+op3, pfx3+sfx3, fact4+op4, pfx4+sfx4, fact5+op5, pfx5+sfx5]
        if gBinning and (gVerbose or gDebug):
            if fact3 == "":
                t3 = "0"
            else:
                t3 = pfx3 + "[par]" + sfx3 + " " + op3 + " " + fact3
            if fact4 == "":
                t4 = "0"
            else:
                t4 = pfx4 + "[par]" + sfx4 + " " + op4 + " " + fact4
            if fact5 == "":
                t5 = "0"
            else:
                t5 = pfx5 + "[par]" + sfx5 + " " + op5 + " " + fact5
        do_print = False
        if len(gBinningPatterns[0]) == 0:
            if (vpfx != "" or vsfx != "") and (pfx1+sfx1 != "" or pfx2+sfx2 != "" or pfx3+sfx2 != ""):
                gBinningPatterns[0] = [gFileName[-1], gLineNo[-1]] + new_pat
                if gBinning and (gVerbose or gDebug):
                    do_print = True
        elif gBinning and not warned:
            matched = False
            for old_pat in gBinningPatterns:
                if old_pat[2:] == new_pat:
                    matched = True
            if not matched:
                warning("Binning equation does not match that on line %d of %s" % (gBinningPatterns[0][1], gBinningPatterns[0][0]))
                gBinningPatterns.append( [gFileName[-1], gLineNo[-1]] + new_pat )
                if gVerbose or gDebug:
                    do_print = True
        if do_print:
            print("Binning pattern: %s[par]%s = %s[par]%s + %s[par]%s %s %s + %s[par]%s %s %s + %s + %s + %s" \
                  % (vpfx, vsfx, pfx0, sfx0, pfx1, sfx1, op1, fact1, pfx2, sfx2, op2, fact2, t3, t4, t5))
            if gDebug:
                print("    derived from %s" % line)
# end of checkBinning


# for if and while
gConditionForNextStmt = Expression("NOTHING")
gConditionForLastStmt = Expression("NOTHING")
gCondBiasDForNextStmt = 0
gCondBiasDForLastStmt = 0

def parseOther( line ):
    global gScopeList, gFileName, gLineNo, gDebug
    global gBlocknames
    global gConditions, gCondBiasDep
    global gConditionForNextStmt, gConditionForLastStmt
    global gCondBiasDForNextStmt, gCondBiasDForLastStmt
    global gStatementInCurrentBlock
    global gAnalogInitial

    this_cond = gConditionForNextStmt
    gConditionForNextStmt = Expression("NOTHING")
    last_cond = gConditionForLastStmt
    gConditionForLastStmt = Expression("NOTHING")
    this_c_bd = gCondBiasDForNextStmt
    gCondBiasDForNextStmt = 0
    last_c_bd = gCondBiasDForLastStmt
    gCondBiasDForLastStmt = 0
    is_analog_initial = False

    parser = Parser(line)
    val = parser.lex()
    keywd = ""

    if parser.isIdentifier():
        keywd = parser.getString()
        if keywd == "analog":
            gStatementInCurrentBlock = True
            val = parser.lex()
            keywd = parser.getString()
            if keywd == "initial":
                val = parser.lex()
                warning("Restrictions in 'analog initial' not checked")
                is_analog_initial = True

    if parser.isNumber() or (keywd == "default" and parser.peekChar() == ":"):
        if len(gScopeList) > 0 and gScopeList[-1].startswith("CASE::"):
            last_val = 0
            if parser.isNumber():
                # case 0: begin ...
                scope = gScopeList[-1]
                this_cond = Expression("==")
                this_cond.e1 = Expression("NAME")
                this_cond.e1.e1 = scope[6:]
                this_cond.e2 = Expression("NUMBER")
                last_val = parser.getNumber()
                this_cond.e2.number = last_val
            else:
                # default: (assume no hidden state if "default" is present)
                this_cond = Expression("NOTHING")
            val = parser.lex()
            while val == ord(','):
                val = parser.lex()
                if parser.isNumber():
                    old_cond = this_cond
                    new_cond = Expression("==")
                    new_cond.e1 = Expression("NAME")
                    new_cond.e1.e1 = scope[6:]
                    new_cond.e2 = Expression("NUMBER")
                    last_val = parser.getNumber()
                    new_cond.e2.number = last_val
                    this_cond = Expression("||")
                    this_cond.e1 = old_cond
                    this_cond.e2 = new_cond
                    val = parser.lex()
                else:
                    error("Expected number after ',' in case value list")
            if val == ord(':'):
                val = parser.lex()
            else:
                error("Expected ':' after case value %d" % last_val)


    # if, else, for, while, repeat, begin, end
    if parser.isIdentifier():
        keywd = parser.getString()

        if keywd in ["end", "else", "if", "for", "repeat", "while", "begin", "case"]:
            gStatementInCurrentBlock = True

        if keywd == "end":
            #warning("pop scope '%s'" % gScopeList[-1])
            bad_scope = False
            if len(gScopeList) > 0:
                gScopeList.pop()
            if len(gConditions) > 0:
                last_cond = gConditions[-1]
                gConditions.pop()
            if len(gCondBiasDep) > 0:
                last_c_bd = gCondBiasDep[-1]
                gCondBiasDep.pop()
            if bad_scope:
                error("Found 'end' without corresponding 'begin'")
            if len(gAnalogInitial) > 0:
                gAnalogInitial.pop()
            if len(gAnalogInitial) > 0:
                is_analog_initial = gAnalogInitial[-1]
            else:
                is_analog_initial = False
            val = parser.lex()
            if parser.isIdentifier():
                keywd = parser.getString()
            elif val == 0:
                # possible else on next line
                this_cond = last_cond
                last_cond = 0
                this_c_bd = last_c_bd
                last_c_bd = 0

        if keywd == "else":
            if last_cond.type != "NOTHING":
                if last_cond.type == "&&" and last_cond.e1.type == "NOT":
                    this_cond = Expression("&&")
                    this_cond.e1 = last_cond.e1
                    this_cond.e2 = Expression("NOT")
                    this_cond.e2.e1 = last_cond.e2
                else:
                    this_cond = Expression("NOT")
                    this_cond.e1 = last_cond
                this_c_bd = last_c_bd
            else: # pragma: no cover
                error("Programming error: 'else' condition")
            val = parser.lex()
            if parser.isIdentifier():
                keywd = parser.getString()
            elif val == 0:
                # statement on next line
                gConditionForNextStmt = this_cond
                gCondBiasDForNextStmt = this_c_bd

        if keywd == "if":
            new_cond = parser.getExpression(True)
            deps = new_cond.getDependencies(False, False)
            [bias_dep, biases] = checkDependencies(deps, "If condition depends on", 0, True, True)
            if this_cond.type == "NOTHING":
                this_cond = new_cond
            else:
                old_cond = this_cond
                this_cond = Expression("&&")
                this_cond.e1 = old_cond
                this_cond.e2 = new_cond
            if bias_dep > 1 and (new_cond.type == "==" or new_cond.type == "!="):
                bias_dep = 3
                warning("Bias-dependent '%s': if %s" % (new_cond.type, new_cond.asString()))
            if bias_dep == 2:
                bias_dep = 1
            this_c_bd = bias_dep
            val = parser.lex()
            if parser.isIdentifier():
                keywd = parser.getString()
                # expect begin, but could be if (cond) var = value;
                if keywd == "if":
                    # if (A) if (B) ... becomes if (A&&B)
                    # TODO: doesn't handle subsequent else
                    new_cond = parser.getExpression(True)
                    deps = new_cond.getDependencies(False, False)
                    [bias_dep, biases] = checkDependencies(deps, "If condition depends on", 0, True, True)
                    old_cond = this_cond
                    this_cond = Expression("&&")
                    this_cond.e1 = old_cond
                    this_cond.e2 = new_cond
                    if bias_dep > 1 and (new_cond.type == "==" or new_cond.type == "!="):
                        bias_dep = 3
                        warning("Bias-dependent '%s': if %s" % (new_cond.type, new_cond.asString()))
                    if bias_dep == 2:
                        bias_dep = 1
                    this_c_bd = bias_dep
                    val = parser.lex()
            if parser.isIdentifier():
                keywd = parser.getString()
            elif val == ord(';'):
                # null statement
                val = parser.lex()
            elif val == 0:
                # statement on next line
                gConditionForNextStmt = this_cond
                gCondBiasDForNextStmt = this_c_bd

        if keywd == "for":
            # for (init; cond; incr)
            varname = ""
            val = parser.lex()
            if val == ord('('):
                val = parser.lex()
            else:
                error("Missing '(' after 'for'")
            if val == ord(';'):
                error("Missing initializer in 'for' loop")
            else:
                if parser.isIdentifier():
                    varname = parser.getString()
                    val = parser.lex()
                else:
                    error("Missing variable name in 'for' loop initialization")
                    if val != ord('='):
                        val = parser.lex()
                bad_init = False
                if val == ord('='):
                    init_expr = parser.getExpression()
                    if init_expr.type == "NOTHING":
                        error("Missing initializer value in 'for' loop")
                    elif varname != "":
                        markVariableAsSet(varname, False, this_cond, [], True, False, "", 0)
                    if init_expr.type == "NUMBER":
                        pass
                    elif init_expr.type == "NAME":
                        deps = [init_expr.e1]
                        checkDependencies(deps, "For loop initialization depends on", 0, False, True)
                    else:
                        # deps = init_expr.getDependencies(False, False)
                        bad_init = True
                        if gDebug:
                            notice("For loop init type %s" % init_expr.type)
                else:
                    bad_init = True
                if bad_init:
                    error("Unexpected initialization of 'for' loop")
                if parser.peekChar() == ';':
                    val = parser.lex()
            if val != ord(';'):
                error("Missing ';' in for loop")
            this_cond = parser.getExpression(True)
            if this_cond.type != "NOTHING":
                deps = this_cond.getDependencies(False, False)
                [bias_dep, biases] = checkDependencies(deps, "For loop condition depends on", 0, True, True)
                val = parser.lex()
            else:
                error("Missing condition in 'for' loop")
                val = parser.lex()
            if val == ord(';'):
                val = parser.lex()
            else:
                error("Missing ';' in for loop")
            good_incr = True
            if parser.isIdentifier():
                incr_var = parser.getString()
                if varname != "" and varname != incr_var:
                    error("For loop initializes '%s' but increments '%s'" % (varname, incr_var))
            elif val == ord(')'):
                error("Missing increment statement in 'for' loop")
            else:
                good_incr = False
            while parser.peekChar().isspace():
                parser.getChar()
            oper = parser.peekOper()
            if oper == "=":
                oper = parser.getOper()
                incr_rhs = parser.getExpression()
                if incr_rhs.type != "+" and incr_rhs.type != "-":
                    good_incr = False
            elif oper == "++":
                oper = parser.getOper()
                # error printed by peekOper()
            elif oper == "--":
                oper = parser.getOper()
                # error printed by peekOper()
            elif oper == "+=" or oper == "-=":
                oper = parser.getOper()
                # error printed by peekOper()
                val = parser.lex()
            elif oper != "":
                oper = parser.getOper()
                good_incr = False
                val = parser.lex()
            if not good_incr:
                error("Unexpected increment statement for 'for' loop")
            if oper != "":
                val = parser.lex()
            if val == ord(')'):
                val = parser.lex()
            else:
                error("Missing ')' in for loop")
            if parser.isIdentifier():
                keywd = parser.getString()
                if keywd != "begin":
                    error("expected 'begin' after 'for'")
            elif val == 0:
                # statement on next line
                gConditionForNextStmt = this_cond
                gCondBiasDForNextStmt = this_c_bd

        if keywd == "repeat":
            this_cond = parser.getExpression()
            if this_cond.type != "NUMBER" or not this_cond.is_int:
                error("Expression for 'repeat' should be an integer")
            val = parser.lex()
            if parser.isIdentifier():
                keywd = parser.getString()
                if keywd != "begin":
                    error("expected 'begin' after 'repeat'")
            elif val == 0:
                # statement on next line
                gConditionForNextStmt = this_cond
                gCondBiasDForNextStmt = this_c_bd

        if keywd == "while":
            this_cond = parser.getExpression(True)
            deps = this_cond.getDependencies(False, False)
            [bias_dep, biases] = checkDependencies(deps, "While condition depends on", 0, True, True)
            this_c_bd = bias_dep
            val = parser.lex()
            if parser.isIdentifier():
                keywd = parser.getString()
                if keywd != "begin":
                    error("expected 'begin' after 'while'")
            elif val == 0:
                # statement on next line
                gConditionForNextStmt = this_cond
                gCondBiasDForNextStmt = this_c_bd

        if keywd == "begin":
            bname = ""
            val = parser.lex()
            if val == ord(':'):
                val = parser.lex()
                if parser.isIdentifier():
                    bname = parser.getString()
                    escaped = parser.isEscaped()
                    scope = getCurrentScope()
                    scoped_name = scope + bname
                    valid = checkIdentifierCollisions(bname, scoped_name, escaped, "Block name")
                    if valid:
                        gBlocknames[bname] = 1
                    val = parser.lex()
                    # named block: can declare variables
                    gStatementInCurrentBlock = False
                else:
                    error("Missing block name after ':'")
            elif parser.isIdentifier():
                keywd = parser.getString()
                if keywd == "if":
                    new_cond = parser.getExpression(True)
                    deps = new_cond.getDependencies(False, False)
                    [bias_dep, biases] = checkDependencies(deps, "If condition depends on", 0, True, True)
                    if this_cond.type == "NOTHING":
                        this_cond = new_cond
                    else:
                        old_cond = this_cond
                        this_cond = Expression("&&")
                        this_cond.e1 = old_cond
                        this_cond.e2 = new_cond
                    if bias_dep > 1 and (new_cond.type == "==" or new_cond.type == "!="):
                        bias_dep = 3
                        warning("Bias-dependent '%s': if %s" % (new_cond.type, new_cond.asString()))
                    if bias_dep == 2:
                        bias_dep = 1
                    this_c_bd = bias_dep
                    val = parser.lex()
            elif val != 0:
                error("Unexpected token %d after 'begin'" % val)
            #warning("push scope '%s'" % bname)
            gScopeList.append(bname)
            gConditions.append(this_cond)
            gCondBiasDep.append(this_c_bd)
            this_cond = Expression("NOTHING")
            this_c_bd = 0
            gAnalogInitial.append(is_analog_initial)

        if keywd == "case":
            var = parser.getExpression()
            if var.type == "NAME":
                varname = var.e1
                [bias_dep, biases] = checkDependencies([varname], "Case statement depends on", 0, True, True)
                if bias_dep:
                    error("Case statement variable '%s' is bias-dependent" % varname)
                elif getIdentifierType(varname) != "integer":
                    error("Case statement variable '%s' should be integer type" % varname)
            else:
                varname = "NOTHING"
                error("Expected variable for 'case' statement")
            gScopeList.append("CASE::" + varname)
            val = parser.lex()

        elif keywd == "endcase":
            if len(gScopeList) > 0 and gScopeList[-1].startswith("CASE::"):
                gScopeList.pop()
            else:
                error("Found 'endcase' without corresponding 'case'")
            val = parser.lex()

    if parser.isIdentifier():
        keywd = parser.getString()

        if keywd in gDisciplines:
            checkDeclarationContext("Discipline")
            bus_range = []
            val = parser.lex()
            if val == ord('['):
                bus_range = parser.getBusRange()
                val = parser.token
            while parser.isIdentifier():
                pname = parser.getString()
                escaped = parser.isEscaped()
                if pname in gPortnames:
                    node = gPortnames[pname]
                    if node.discipline != "":
                        error("Duplicate specification of discipline for port (terminal) '%s'" % pname)
                    node.discipline = keywd
                    if bus_range != []:
                        if node.is_bus:
                            if node.msb != bus_range[0] or node.lsb != bus_range[1]:
                                error("Discipline declaration for port '%s' has different bus range than port declaration" % pname)
                        else:
                            error("Discipline declaration for scalar port '%s' has bus range [%d:%d]" % (pname, bus_range[0], bus_range[1]))
                    elif node.is_bus:
                        error("Discipline declaration for '%s' should also have bus range [%d:%d]" % (pname, node.msb, node.lsb))
                else:
                    # internal node
                    if pname in gNodenames and gNodenames[pname].type == "ground":
                        node = gNodenames[pname]
                        if node.discipline == "":
                            node.discipline = keywd
                        else:
                            error("Duplicate discipline declaration for ground '%s'" % pname)
                    else:
                        valid = checkIdentifierCollisions(pname, pname, escaped, "Node name")
                        if valid:
                            node = Port(pname, gFileName[-1], gLineNo[-1])
                            node.type = "internal"
                            node.discipline = keywd
                            if bus_range != []:
                                node.is_bus = True
                                node.msb = bus_range[0]
                                node.lsb = bus_range[1]
                            gNodenames[pname] = node
                val = parser.lex()
                if val == ord(','):
                    val = parser.lex()
                elif val == ord('='):
                    error("Net initializers not supported")
                    defv = parser.getExpression()
                    val = parser.lex()
                elif parser.isIdentifier():
                    error("Missing ',' in list of nodes")
            if val == ord(';'):
                val = parser.lex()
            else:
                error("Missing ';' after list of nodes")

        elif keywd == "ground":
            checkDeclarationContext("Ground")
            val = parser.lex()
            while parser.isIdentifier():
                pname = parser.getString()
                escaped = parser.isEscaped()
                if pname in gNodenames:
                    node = gNodenames[pname]
                    if node.type == "internal":
                        node.type = "ground"
                    elif node.type == "ground":
                        error("Duplicate 'ground' declaration for '%s'" % pname)
                    else:
                        error("Cannot specify 'ground' for port '%s'" % pname)
                else:
                    valid = checkIdentifierCollisions(pname, pname, escaped, "Node name")
                    if valid:
                        node = Port(pname, gFileName[-1], gLineNo[-1])
                        node.type = "ground"
                        node.discipline = ""
                        gNodenames[pname] = node
                val = parser.lex()
                if val == ord(','):
                    val = parser.lex()
                elif val == ord('='):
                    error("Net initializers not supported for grounds")
                    defv = parser.getExpression()
                    val = parser.lex()
                elif parser.isIdentifier():
                    error("Missing ',' in list of nodes")
            if val == ord(';'):
                val = parser.lex()
            else:
                error("Missing ';' after list of nodes")

        elif keywd == "branch":
            checkDeclarationContext("Branch")
            bname = ""
            node1 = ""
            node2 = ""
            val = parser.lex()
            if val == ord('('):
                val = parser.lex()
            else:
                error("Missing '(' in branch declaration")
            if parser.isIdentifier():
                node1 = parser.getString()
            else:
                error("Expected node name in branch declaration")
            if val != ord(')'):
                val = parser.lex()
                if val == ord(','):
                    val = parser.lex()
                    if parser.isIdentifier():
                        node2 = parser.getString()
                    else:
                        error("Expected node name in branch declaration")
            if val != ord(')'):
                val = parser.lex()
            if val == ord(')'):
                val = parser.lex()
            else:
                error("Missing ')' in branch declaration")
            while parser.isIdentifier():
                bname = parser.getString()
                escaped = parser.isEscaped()
                valid = checkIdentifierCollisions(bname, bname, escaped, "Branch")
                if valid:
                    branch = Branch(bname)
                    branch.node1 = node1
                    branch.node2 = node2
                    gBranches[bname] = branch
                    disc1 = ""
                    if node1 in gPortnames:
                        disc1 = gPortnames[node1].discipline
                    elif node1 in gNodenames:
                        disc1 = gNodenames[node1].discipline
                    elif node1 != "":
                        error("Unknown node '%s' in branch declaration" % node1)
                    if node2 != "":
                        disc2 = ""
                        if node2 in gPortnames:
                            disc2 = gPortnames[node2].discipline
                        elif node2 in gNodenames:
                            disc2 = gNodenames[node2].discipline
                        else:
                            error("Unknown node '%s' in branch declaration" % node1)
                        if disc1 != "" and disc2 != "" and disc1 != disc2:
                            error("Conflicting disciplines of nodes in branch '%s'" % bname)
                    branch.discipline = disc1
                val = parser.lex()
                if val == ord(','):
                    val = parser.lex()
            if bname == "":
                error("Missing name for branch declaration")
            if val == ord(';'):
                val = parser.lex()
            else:
                error("Missing ';' after branch declaration")

        elif keywd in ["cmos", "rcmos", "bufif0", "bufif1", "notif0", "notif1", \
                       "nmos", "pmos", "rnmos", "rpmos", \
                       "and", "nand", "or", "nor", "xor", "xnor", "buf", "not", \
                       "tranif0", "tranif1", "rtranif0", "rtranif1", \
                       "tran", "rtran", "pulldown", "pullup"]: # pragma: no cover
            error("Gate instantiation '%s' not expected in compact model" % keywd)

        elif keywd in ["always", "initial"]: # pragma: no cover
            error("Digital '%s' blocks not supported in Verilog-A" % keywd)
            rest = parser.getRestOfLine()

        elif keywd in ["task", "endtask"]: # pragma: no cover
            error("Tasks not supported in compact models")

        elif keywd in ["supply0", "supply1", "tri", "triand", "trior", "tri0", "tri1", \
                       "uwire", "wire", "wand", "wor", "trireg", "wreal", \
                       "realtime", "reg", "time", "specparam"]: # pragma: no cover
            error("Unsupported declaration: '%s'" % keywd)

        else:
            if parser.peekChar() == '[':
                parser.checkBusIndex(keywd)
            val = parser.lex()
            if val == ord('='):
                all_zero_assigns = True
                while val == ord('='):
                    # variable assignment
                    gStatementInCurrentBlock = True
                    varname = keywd
                    keywd = ""
                    gConditions.append(this_cond)
                    rhs = parser.getExpression()
                    if rhs.type != "NUMBER" or rhs.number != 0:
                        all_zero_assigns = False
                    gConditions.pop()
                    # ensure RHS values have all been set
                    deps = rhs.getDependencies(True, False)
                    [bias_dep, biases] = checkDependencies(deps, "Variable %s depends on" % varname, this_c_bd, False, True)
                    found = markVariableAsSet(varname, bias_dep, this_cond, biases, False, False, "", 0)
                    if rhs.type == "+" and (rhs.e1.type == "+" or rhs.e2.type == "+"):
                        # possible binning equation P_i = P + LP*Inv_L + WP*Inv_W + PP*Inv_P
                        t3 = rhs.e2
                        if rhs.e1.type == "+" and (t3.type == "*" or t3.type == "/"):
                            t2 = rhs.e1.e2
                            if rhs.e1.e1.type == "+" and (t2.type == "*" or t2.type == "/"):
                                t1 = rhs.e1.e1.e2
                                t0 = rhs.e1.e1.e1
                                if t0.type == "NAME" and (t1.type == "*" or t1.type == "/"):
                                    checkBinning(varname, t0, t1, t2, t3, 0, 0, line)
                                elif t0.type == "+" and (t1.type == "*" or t1.type == "/"):
                                    # P_i = P + Inv_L*LP + Inv_NFIN*NP + Inv_LNFIN*PP + Inv_W*WP + Inv_WL*WLP;
                                    t5 = t3
                                    t4 = t2
                                    t3 = t1
                                    t2 = t0.e2
                                    if t0.e1.type == "+" and (t2.type == "*" or t2.type == "/"):
                                        t1 = t0.e1.e2
                                        t0 = t0.e1.e1
                                        if t0.type == "NAME" and (t1.type == "*" or t1.type == "/"):
                                            checkBinning(varname, t0, t1, t2, t3, t4, t5, line)
                    val = parser.lex()
                    if val == ord(';'):
                        val = parser.lex()
                    else:
                        error("Missing ';' after assignment")
                    if parser.isIdentifier():
                        keywd = parser.getString()
                        if keywd == "else":
                            # if (A) X=0; else X=1;
                            rest = keywd + " " + parser.getRestOfLine()
                            keywd = ""
                            gConditionForLastStmt = this_cond
                            gCondBiasDForLastStmt = this_c_bd
                            parseOther(rest)
                            this_cond = gConditionForLastStmt
                            this_c_bd = gCondBiasDForLastStmt
                        elif keywd == "if":
                            parser.ungetChars("if")
                            if gStyle:
                                style("Multiple statements on a single line")
                        elif gStyle and not all_zero_assigns:
                            style("Multiple assignments on a single line")
                        if parser.peekChar() == '[':
                            parser.checkBusIndex(keywd)
                        if keywd != "if":
                            val = parser.lex()
                val = parser.lex()
                if val == 0 and keywd != "":
                    if keywd == "end":
                        gConditionForLastStmt = this_cond
                        gCondBiasDForLastStmt = this_c_bd
                        parseOther(keywd)
                        this_cond = gConditionForLastStmt
                        this_c_bd = gCondBiasDForLastStmt
                    else: # pragma: no cover
                        fatal("Unexpected '%s' after assignment" % keywd)

            elif val == ord('('):
                if keywd in gAccessFuncs:
                    gStatementInCurrentBlock = True
                    acc = keywd
                    if is_analog_initial:
                        error("Branch access in analog initial block")
                    nname1 = ""
                    nname2 = ""
                    val = parser.lex()
                    is_port_flow = False
                    if val == ord('<'):
                        is_port_flow = True
                        val = parser.lex()
                    if parser.isIdentifier():
                        nname1 = parser.getString()
                        val = parser.lex()
                        if val == ord('['):
                            parser.checkBusIndex(nname1)
                            val = parser.lex()
                    else:
                        error("Missing node/branch name in potential or flow access")
                    if is_port_flow:
                        if val == ord('>'):
                            val = parser.lex()
                        else:
                            error("Missing '>' in port flow access")
                    elif val == ord(','):
                        val = parser.lex()
                        if parser.isIdentifier():
                            nname2 = parser.getString()
                            val = parser.lex()
                            if val == ord('['):
                                parser.checkBusIndex(nname2)
                                val = parser.lex()
                        else:
                            error("Missing node/branch name after ','")
                    if val == ord(')'):
                        val = parser.lex()
                    else:
                        error("Missing ')' in potential or flow access")
                    dname = ""
                    if nname1 != "":
                        if nname1 in gBranches:
                            dname = gBranches[nname1].discipline
                            if nname2 != "":
                                if nname2 in gBranches:
                                    error("Invalid potential or flow access %s(branch,branch)" % acc)
                                else:
                                    error("Invalid potential or flow access %s(branch,node)" % acc)
                        else:
                            if nname1 in gPortnames:
                                dname = gPortnames[nname1].discipline
                            elif nname1 in gNodenames:
                                dname = gNodenames[nname1].discipline
                            else:
                                error("Identifier '%s' is not a port, node, or branch in potential or flow access" % nname1)
                            if nname2 != "":
                                dname2 = ""
                                if nname2 in gBranches:
                                    dname2 = gBranches[nname2].discipline
                                    error("Invalid potential or flow access %s(node,branch)" % acc)
                                elif nname2 in gPortnames:
                                    dname2 = gPortnames[nname2].discipline
                                elif nname2 in gNodenames:
                                    dname2 = gNodenames[nname2].discipline
                                else:
                                    error("Identifier '%s' is not a port, node, or branch in potential or flow access" % nname2)
                                if dname2 != "" and dname2 != dname:
                                    error("Nodes '%s' and '%s' belong to different disciplines" % (nname1, nname2))
                    if dname in gDisciplines:
                        disc = gDisciplines[dname]
                        fname = disc.flow
                        pname = disc.potential
                        if fname in gNatures and acc == gNatures[fname].access:
                            registerContrib(acc, "flow", nname1, nname2, this_cond, this_c_bd)
                        elif pname in gNatures and acc == gNatures[pname].access:
                            registerContrib(acc, "potential", nname1, nname2, this_cond, this_c_bd)
                        else:
                            error("Incorrect access function '%s' for %s '%s'" % (acc, dname, nname1))
                    if val == ord('<') and parser.peekChar() == '+':
                        parser.getChar()
                        rhs = parser.getExpression()
                        deps = rhs.getDependencies(False, False)
                        [bias_dep, biases] = checkDependencies(deps, "Branch contribution depends on", this_c_bd, False, True)
                        # check dependencies again, but ignore bias_dep in _noise calls
                        deps = rhs.getDependencies(False, True)
                        [bias_dep, biases] = checkDependencies(deps, "Branch contribution depends on", this_c_bd, False, False)
                        if bias_dep == 4:
                            error("Branch contribution depends on quantity obtained from ddx (requires second derivatives)")
                        elif bias_dep >= 3:
                            error("Branch contribution depends on quantity with bad derivative (from bias-dependent == or !=)")
                    val = parser.lex()
                    if val == ord(';'):
                        val = parser.lex()
                        if gStyle and parser.isIdentifier():
                            keywd = parser.getString()
                            if keywd in gAccessFuncs:
                                style("Multiple contributions on a single line")
                    else:
                        error("Missing ';' after contribution")

                elif keywd in ["$strobe", "$display", "$monitor", "$write", "$debug", \
                               "$fatal", "$error", "$warning", "$finish", "$stop"]:
                    gStatementInCurrentBlock = True
                    if keywd in ["$finish", "$stop"]:
                        notice("$error() preferred over %s()" % keywd)
                    elif keywd == "$debug":
                        warning("%s() may degrade performance" % keywd)

                    args = parser.parseArgList()
                    deps = []
                    for arg in args:
                        deps += arg.getDependencies(False, False)
                    checkDependencies(deps, "Task %s references" % keywd, this_c_bd, False, True)
                    val = parser.lex()
                    if val == ord(')'):
                        val = parser.lex()
                    else:
                        error("Missing ')' in call of task %s" % keywd)
                    if keywd in ["$finish", "$stop"]:
                        if len(args) == 1:
                            if args[0].type != "NUMBER":
                                error("Argument to %s must be a number" % keywd)
                        elif len(args) > 1:
                            error("Task %s takes at most one argument" % keywd)
                    elif keywd == "$fatal":
                        if len(args) >= 1:
                            if args[0].type != "NUMBER":
                                error("Argument to %s must be a number" % keywd)
                    else:
                        if len(args) == 0:
                            error("Expected at least 1 argument for %s" % keywd)
                        elif args[0].type != "STRING":
                            error("First argument to %s must be a string" % keywd)
                    if val == ord(';'):
                        val = parser.lex()
                    else:
                        error("Missing ';' after task %s" % keywd)
                else:
                    error("Unexpected '%s'" % keywd)
                    rest = parser.getRestOfLine()
                    val = 0
            else:
                reported = False
                if val in [ord('+'), ord('-'), ord('*'), ord('/')]:
                    next = parser.lex()
                    if next in [ord('+'), ord('-'), ord('*'), ord('/'), ord('=')]:
                        oper = chr(val) + chr(next)
                        if oper in ["++", "--", "+=", "-=", "*=", "/="]:
                            error("Operator '%s' not valid in Verilog-A" % oper)
                            reported = True
                rest = parser.getRestOfLine()
                if not reported:
                    error("Unexpected '%s'" % line)
                    val = 0

    elif val == ord('@'):
        error("Events (@) should not be used in compact models")
        #  @ (initial_step or initial_step("static","pss") or initial_step("static","pdisto")) begin
        val = parser.lex()
        args = []
        if val == ord('('):
            val = parser.lex()
            while parser.isIdentifier():
                arg = Expression("NAME")
                arg.e1 = parser.getString()
                args.append(arg)
                val = parser.lex()
                if val == ord('('):
                    arg.type = "FUNCCALL"
                    if arg.e1 in ["initial_step", "final_step"]:
                        val = parser.lex()
                        while parser.isString():
                            sub = Expression("STRING")
                            sub.e1 = parser.getString()
                            arg.args.append(sub)
                            val = parser.lex()
                            if val == ord(','):
                                val = parser.lex()
                    else: # cross, above
                        while val != ord(')'):
                            sub = parser.getExpression()
                            arg.args.append(sub)
                            val = parser.lex()
                            if val == ord(','):
                                val = parser.lex()
                    if val == ord(')'):
                        val = parser.lex()
                    else:
                        error("Missing ')' in event specification")
                    if parser.isIdentifier():
                        if parser.getString() == "or":
                            val = parser.lex()
                            if not parser.isIdentifier():
                                error("Missing event after 'or'")
            if val == ord(')'):
                val = parser.lex()
        if parser.isIdentifier():
            keywd = parser.getString()
            if keywd == "begin":
                bname = ""
                val = parser.lex()
                if val == ord(':'):
                    val = parser.lex()
                    if parser.isIdentifier():
                        bname = parser.getString()
                        gStatementInCurrentBlock = False
                        val = parser.lex()
                    else:
                        error("Missing block name after ':'")
                gScopeList.append(bname)
                event = Expression("EVENT")
                event.args = args
                gConditions.append(event)
                gCondBiasDep.append(0)

    elif val == ord('`'):
        rest = parser.getRestOfLine()
        error("Unhandled macro `%s" % rest)
    elif val == ord('#'):
        rest = parser.getRestOfLine()
        error("Delays (#) should not be used in compact models")
    elif val == 0:
        pass # end of line
    elif parser.isNumber():
        num = parser.getNumber()
        error("Unexpected number %d" % num)
    else:
        error("Unexpected %s" % format_char(val))

    gConditionForLastStmt = this_cond
    gCondBiasDForLastStmt = this_c_bd

    if val != 0:
        rest = ""
        if parser.isIdentifier():
            rest = parser.getString()
        rest += parser.getRestOfLine()
        parseOther(rest)
# end of parseOther


def getCurrentScope():
    global gScopeList
    scope = ""
    if len(gScopeList) > 1:
        scope = gScopeList[1]
        for sc in gScopeList[2:]:
            if sc != "":
                if scope != "":
                    scope += "."
                scope += sc
    if scope != "":
        scope += "."
    return scope

def getCurrentConditions( cond_in, c_bd_in ):
    global gConditions, gCondBiasDep
    conds = ""
    bias_dep = c_bd_in
    if cond_in.type != "NOTHING":
        conds = cond_in.asString()
    for cond in gConditions:
        if cond.type != "NOTHING":
            if conds != "":
                conds += "&&"
            conds += cond.asString()
    for c_bd in gCondBiasDep:
        if c_bd > bias_dep:
            bias_dep = c_bd
    return [conds, bias_dep]


# create disciplines, else it's pointless to continue
def addBasicDisciplines():
    warning("Creating basic disciplines and natures to continue parsing")
    #
    volt = Nature("Voltage", True)
    volt.access = "V"
    gNatures[volt.name] = volt
    curr = Nature("Current", True)
    curr.access = "I"
    gNatures[curr.name] = curr
    elec = Discipline("electrical")
    elec.potential = "Voltage"
    elec.flow = "Current"
    gDisciplines[elec.name] = elec
    #
    pwr = Nature("Power", True)
    pwr.access = "Pwr"
    gNatures[pwr.name] = pwr
    temp = Nature("Temperature", True)
    temp.access = "Temp"
    gNatures[temp.name] = temp
    thermal = Discipline("thermal")
    thermal.potential = "Temperature"
    thermal.flow = "Power"
    gDisciplines[thermal.name] = thermal
    #
    gAccessFuncs["V"] = 1
    gAccessFuncs["I"] = 1
    gAccessFuncs["Pwr"] = 1
    gAccessFuncs["Temp"] = 1


def handleCompilerDirectives( line ):
    global gIfDefStatus, gIncDir, gMissingConstantsFile

    line = line.strip()
    if line == "":
        return line

    # compiler directives
    did_compiler_directive = False

    if line.startswith("`ifdef") or line.startswith("`ifndef"):
        if line.startswith("`ifdef"):
            start = 6
            true_status = "TRUE"
            false_status = "FALSE"
        else:
            start = 7
            false_status = "TRUE"
            true_status = "FALSE"
        if not line[start].isspace():
            error("Missing space after %s" % line[0:start])
        while line[start].isspace():
            start += 1
        stop = start + 1
        while stop < len(line) and not line[stop].isspace():
            stop += 1
        token = line[start:stop]
        if stop < len(line):
            error("Extra characters after `ifdef")
        if token in gMacros:
            gIfDefStatus.append(true_status)
        else:
            gIfDefStatus.append(false_status)
        did_compiler_directive = True
    elif line.startswith("`else"):
        if len(gIfDefStatus) > 0:
            if gIfDefStatus[-1] == "TRUE":
                gIfDefStatus[-1] = "FALSE"
            elif gIfDefStatus[-1] == "FALSE":
                gIfDefStatus[-1] = "TRUE"
            else: # pragma: no cover
                fatal("Invalid `ifdef status")
        else:
            error("Found `else without `ifdef")
        did_compiler_directive = True
    elif line.startswith("`endif"):
        if len(gIfDefStatus) > 0:
            gIfDefStatus.pop()
        else:
            error("Unmatched `endif")
        did_compiler_directive = True

    # check ifdef status
    ifdef_status = True
    for stat in gIfDefStatus:
        if stat == "FALSE":
            ifdef_status = False

    if ifdef_status:
        if line.startswith("`define"):
            parseMacro( line )
            did_compiler_directive = True
        elif line.startswith("`undef"):
            parseUndef( line )
            did_compiler_directive = True

        elif line.startswith("`include"):
            rest = line[8:].strip()
            while rest != "":
                if rest[0] == '"':
                    fname = rest[1:]
                    qt = fname.find("\"")
                    if qt > 0:
                        rest  = fname[qt+1:].strip()
                        fname = fname[0:qt]
                    elif qt == 0:
                        error("Empty filename for `include")
                        rest = fname[1:].strip()
                        fname = ""
                    else:
                        error("Missing end-quote for `include")
                        rest = ""
                else:
                    error("Filename should be enclosed in quotes")
                    fname = rest
                    qt = fname.find(" ")
                    if qt > 0:
                        rest  = fname[qt+1:].strip()
                        fname = fname[0:qt]
                    else:
                        rest = ""
                if fname != "":
                    fpath = os.path.dirname(gFileName[-1])
                    f_fpn = os.path.join(fpath, fname)
                    found = False
                    if os.path.exists(f_fpn) and os.path.isfile(f_fpn):
                        found = True
                        parseFile( f_fpn )
                    else:
                        for fpath in gIncDir:
                            f_fpn = os.path.join(fpath, fname)
                            if os.path.exists(f_fpn) and os.path.isfile(f_fpn):
                                found = True
                                parseFile( f_fpn )
                                break
                    if not found:
                        error("File not found: %s" % line)
                        if fname == "discipline.h" or fname == "disciplines.vams":
                            addBasicDisciplines()
                        elif fname == "constants.h" or fname == "constants.vams":
                            gMissingConstantsFile = fname
                if len(rest) > 0:
                    if rest.startswith("`include"):
                        error("Multiple `include directives on one line")
                        rest = rest[8:].strip()
                    else:
                        error("Unexpected characters after `include: %s" % rest)
                        rest = ""
            did_compiler_directive = True

        elif line.startswith("`begin_keywords")      or line.startswith("`end_keywords") \
          or line.startswith("`celldefine")          or line.startswith("`endcelldefine") \
          or line.startswith("`default_discipline")  or line.startswith("`default_nettype") \
          or line.startswith("`default_transition") \
          or line.startswith("`nounconnected_drive") or line.startswith("`unconnected_drive") \
          or line.startswith("`pragma")              or line.startswith("`resetall") \
          or line.startswith("`timescale"):
            error("Unsupported compiler directive: %s" % line)
            did_compiler_directive = True

    if did_compiler_directive or not ifdef_status:
        line = "" # done with this line

    if line.find("`") >= 0:
        line = expandMacro( line )

    return line
# end of handleCompilerDirectives


def preProcessNatures( lines ):
    # preprocess so idt_natures are defined
    for line in lines:
        gLineNo[-1] += 1
        if line.startswith("nature") and len(gScopeList) == 0:
            parseNatureDecl(line, False)


def parseLines( lines ):
    global gScopeList, gLineNo
    global gStatementInCurrentBlock
    prev_attribs = []

    j = 0
    while j < len(lines):
        gLineNo[-1] += 1
        line = lines[j]
        j += 1

        line = handleCompilerDirectives(line)
        line = line.strip()
        if line == "":
            continue

        # macros may have multiple lines, want to parse them separately
        parts = line.split("\n")
        i = 0
        while i < len(parts):
            part = parts[i]
            i += 1
            part = part.strip()
            if part == "":
                continue
            retarray = getAttributes(part)
            part = retarray[0]
            attribs = prev_attribs + retarray[1:]
            if part == "":
                if i < len(parts):
                    part = parts[i]
                    i += 1
                    retarray = getAttributes(part)
                    part = retarray[0]
                    attribs.append(retarray[1:])
                else:
                    prev_attribs = attribs
                    continue
            elif part[-1] == ';':
                prev_attribs = []

            # natures and disciplines
            #
            if part.startswith("nature"):
                if len(gScopeList) == 0:
                    parseNatureDecl(part, True)
                else: # pragma: no cover
                    error("Natures should be defined at top level")
            elif len(gScopeList) == 1 and gScopeList[0].startswith("NATURE::"):
                parseNatureLine(part)
            elif part.startswith("discipline"):
                if len(gScopeList) == 0:
                    parseDisciplineDecl(part)
                else: # pragma: no cover
                    error("Disciplines should be defined at top level")
            elif len(gScopeList) == 1 and gScopeList[0].startswith("DISCIPLINE::"):
                parseDisciplineLine(part)
            elif part.startswith("endnature") or part.startswith("enddiscipline"):
                verifyEndScope(part[3:])

            # module
            #
            elif part.startswith("module") or part.startswith("macromodule"):
                if part.startswith("macromodule"):
                    error("Macromodels not supported")
                if len(gScopeList) == 0:
                    gScopeList.append("MODULE")
                    parseModuleDecl(part)
                    gStatementInCurrentBlock = False
                elif len(gScopeList) == 1: # pragma: no cover
                    if gScopeList[-1] == "MODULE":
                        fatal("Nested module declaration")
                    else:
                        fatal("Module declared in unexpected context")
                else: # pragma: no cover
                    scope = getCurrentScope()
                    fatal("Unexpected 'module' in scope %s" % scope)
            elif part.startswith("endmodule"):
                verifyEndScope("module")
                checkPorts()
                checkParmsAndVars()

            # functions
            #
            elif part.startswith("analog function") or part.startswith("function"):
                parseFunction(part)
                checkDeclarationContext("Function")
                gStatementInCurrentBlock = False
            elif part.startswith("endfunction"):
                gStatementInCurrentBlock = False
                if len(gScopeList) > 0 and gScopeList[-1].startswith("FUNCTION::"):
                    if gCurrentFunction:
                        fname = gCurrentFunction.name
                        fscope = gScopeList[-1]
                        if fscope == "FUNCTION::" + fname:
                            fscope += "."
                        else: # pragma: no cover
                            fatal("Programming error: endfunction %s but scope %s" % (fname, fscope))
                        func_assign = False
                        for var in gVariables.values():
                            if var.name.startswith(fscope):
                                vn = var.name[len(fscope):]
                                if vn == fname:
                                    if var.assign > 0:
                                        func_assign = True
                                elif var.assign < 0:
                                    gFileName.append(var.declare[0])
                                    gLineNo.append(var.declare[1])
                                    if var.used:
                                        warning("In function %s, variable '%s' was never assigned a value" % (fname,vn))
                                    else:
                                        warning("In function %s, variable '%s' was never set and never used" % (fname,vn))
                                    gLineNo.pop()
                                    gFileName.pop()
                                elif not var.used and not var.assign == 0:
                                    warning("In function %s, variable '%s' was never used" % (fname,vn))
                        if not func_assign:
                            if len(gCurrentFunction.outputs) > 0:
                                warning("Function '%s' is not assigned a return value" % gCurrentFunction.name)
                            else:
                                error("Function '%s' is not assigned a return value" % gCurrentFunction.name)
                    else: # pragma: no cover
                        fatal("Programming error: no current function")
                    gScopeList.pop()
                    if len(gScopeList) == 1 and gScopeList[0].startswith("FUNCTION::"):
                        # dummy scope created for function declared out of module
                        gScopeList.pop()
                else:
                    error("Found 'endfunction' without corresponding 'function'")

            # parameters, variables, ports
            #
            elif part.startswith("parameter") or part.startswith("aliasparam") or part.startswith("localparam") or \
                    (part.startswith("real") and part[4].isspace()) or (part.startswith("integer") and part[7].isspace()) or \
                    (part.startswith("genvar") and part[6].isspace()):
                if part[-1] != ';' and i >= len(parts):
                    tmpline = part
                    tmpj = j
                    while tmpline[-1] != ';' and tmpj < len(lines):
                        nextln = lines[tmpj]
                        nextln = handleCompilerDirectives(nextln)
                        if nextln != "":
                            nparts = nextln.split()
                            if nparts[0] in gVAMSkeywords and not nparts[0] in gMathFunctions:
                                break
                            tmpline += " " + nextln
                            tmpline = tmpline.strip()
                        tmpj += 1
                    if tmpline[-1] == ';' and tmpj > j:
                        part = tmpline
                        j = tmpj
                        gLineNo[-1] = j
                if part.startswith("parameter") or part.startswith("aliasparam") or part.startswith("localparam"):
                    parseParamDecl(part, attribs)
                else:
                    parseVariableDecl(part, attribs)
            elif part.startswith("inout") or part.startswith("input") or part.startswith("output"):
                parsePortDirection(part)

            # statements (and other content)
            #
            else:
                while part[-1] != ';' and i < len(parts):
                    # may want to join the next line
                    # in the case of multi-line assignments
                    # or if conditions split across lines
                    num_parens = 0
                    if part.find("if") >= 0:
                        for ch in part:
                            if ch == '(':
                                num_parens += 1
                            elif ch == ')':
                                num_parens -= 1
                    if num_parens > 0 or \
                           (part.find("if") < 0 and part.find("else") < 0 and \
                            part.find("begin") < 0 and part.find("end") < 0 and \
                            part.find("case") < 0):
                        part += parts[i]
                        part = part.strip()
                        i += 1
                    else:
                        break
                if i >= len(parts) and part[-1] != ';' and part != "end":
                    tmpj = j
                    tmpline = part
                    while tmpj < len(lines):
                        num_parens = 0
                        for ch in tmpline:
                            if ch == '(':
                                num_parens += 1
                            elif ch == ')':
                                num_parens -= 1
                        if num_parens > 0:
                            nextln = lines[tmpj]
                            nextln = handleCompilerDirectives(nextln)
                            tmpline += " " + nextln
                            tmpj += 1
                        elif num_parens == 0:
                            if tmpj > j:
                                part = tmpline
                                j = tmpj
                                gLineNo[-1] = j
                            break
                        elif num_parens < 0: # pragma: no cover
                            fatal("Unexpected ')'")

                    if tmpline[-1] != ';':
                        parser = Parser(tmpline)
                        val = parser.lex()
                        if parser.isIdentifier():
                            do_it = False
                            keywd = parser.getString()
                            val = parser.lex()
                            if val == ord('='):
                                do_it = True # assignment -> need ;
                            elif val == ord('(') and keywd in gAccessFuncs:
                                do_it = True # contribution -> need ;
                            if do_it:
                                while tmpline[-1] != ';' and tmpj < len(lines):
                                    nextln = lines[tmpj]
                                    nextln = handleCompilerDirectives(nextln)
                                    if nextln != "":
                                        nparts = nextln.split()
                                        if nparts[0] in gVAMSkeywords and not nparts[0] in gMathFunctions:
                                            break
                                        tmpline += " " + nextln
                                        tmpline = tmpline.strip()
                                    tmpj += 1
                                if tmpline[-1] == ';' and tmpj > j:
                                    part = tmpline
                                    j = tmpj
                                    gLineNo[-1] = j
                parseOther(part)
# end of parseLines


# deal with comments: // and /* */ (watch out for quotes)
def removeComments( line, in_comment ):
    retstr = ""
    in_quote = False
    i = 0
    while i < len(line):
        if in_comment and not in_quote and i+1 < len(line) and line[i:i+2] == "*/":
            in_comment = False
            i += 2
        if not in_comment and not in_quote and i+1 < len(line) and line[i:i+2] == "/*":
            in_comment = True
            i += 2
        if not in_comment and not in_quote and i+1 < len(line) and line[i:i+2] == "//":
            return [retstr, in_comment]
        ch = line[i]
        if ch == '"':
            in_quote = not in_quote
        if not in_comment:
            retstr += ch
        i += 1
    return [retstr, in_comment]


def getIndent( line ):
    num_spaces = 0
    for ch in line:
        if ch == '\t': # TAB discouraged
            return -1
        elif ch == ' ':
            num_spaces += 1
        else:
            break
    return num_spaces


def makeIndent( num_spaces ):
    indent = ""
    for i in range(num_spaces):
        indent += " "
    return indent


def fixSpaces(line, keywd):
    ret = ""
    i = 0
    klen = len(keywd)
    while i < len(line):
        if line[i:i+klen] == keywd:
            ret += line[i:i+klen]
            i += klen
            break
        else:
            ret += line[i]
            i += 1
    ret += ' '
    while i < len(line) and line[i].isspace():
        i += 1
    if line[i] == '(':
        num_parens = 0
        while i < len(line):
            ch = line[i]
            if ch == '(':
                num_parens += 1
            elif ch == ')':
                num_parens -= 1
            ret += line[i]
            i += 1
            if num_parens == 0:
                break
        if i < len(line):
            ret += ' '
            while i < len(line) and line[i].isspace():
                i += 1
    ret += line[i:]
    return ret


def parseFile( filename ):
    global gFileName, gLineNo, gScopeList, gVAMScompdir
    print("Reading %s" % filename)

    check_style = gStyle
    if filename.find("disciplines.vams") >= 0 or filename.find("discipline.h") >= 0 \
            or filename.find("constants.vams") >= 0 or filename.find("constants.h") >= 0:
        # turn off style-checking for standard header files
        check_style = False

    # read the lines from the file
    try:
        if sys.version_info >= (3,0):
            with open(filename, newline='') as IF:
                linesOfCode=IF.readlines()
        else:
            with open(filename) as IF:
                linesOfCode=IF.readlines()

    except:
        try:
            with open(filename, encoding='windows-1252') as IF:
                linesOfCode=IF.readlines()
            warning("Non-standard encoding in file %s (windows-1252)" % filename)
        except:
            fatal("Failed to open file: %s" % filename)
    gFileName.append(filename)
    gLineNo.append(0)

    # variables to track indentation style
    fixedLines = []
    did_fix = False
    required_indent = 0
    optional_indent = 0
    optional_indent_reason = ""
    single_line_indent = 0
    previous_single_indent = 0
    indent_module = -1
    ifdef_depth = 0
    indents_in_ifdefs = []
    indent_ifdef = -1
    first_ifdef_indent = -1
    indent_inside_ifdef = -1
    first_inside_ifdef = -1
    indent_multiline_define = -1
    indent_define = -1
    first_define_indent = -1
    special_define_indent = 0
    indent_this_define = -1
    first_multiline_define = -1
    in_multiline_define = False
    in_case_block = False
    continued_line = ""
    unclosed_parens = 0
    dos_format = -1

    # handle comments and line-continuation characters
    linesToParse = []
    lineToParse = ""
    in_comment = False

    for line in linesOfCode:
        gLineNo[-1] += 1
        if line.endswith("\r\n"):
            # dos format; change to unix
            line = line[:-2] + "\n"
            if dos_format == -1:
                dos_format = 1
                style("MS-DOS file format (\\r\\n)")
            elif dos_format == 0:
                style("Carriage return (\\r or ^M) before newline")
        elif dos_format == 1:
            style("Inconsistent line termination (\\n vs \\r\\n)")
        elif dos_format == -1:
            dos_format = 0
        orig_line = line
        if check_style:
            if line.find("\t") >= 0:
                style("Use of TAB characters discouraged")
            if len(line) > 1 and line[-2] == ' ':
                style("Space character at end of line")
        if in_comment:
            if line.find("*/") >= 0:
                retval = removeComments(line, in_comment)
                line = retval[0]
                in_comment = retval[1]
            else:
                line = ""
            if gFixIndent:
                fixed_line = orig_line.rstrip() + "\n"
                fixedLines.append(fixed_line)
                if fixed_line != orig_line:
                    did_fix = True
        else:
            if line.find("//") >= 0 or line.find("/*") >= 0:
                retval = removeComments(line, in_comment)
                line = retval[0]
                in_comment = retval[1]
            if check_style:
                # getIndent from orig_line to check indent of comments
                indent = getIndent(orig_line)
                fixed_line = orig_line.rstrip() + "\n"
                stripped = line.strip()
                if stripped == "":
                    if optional_indent > 0 and orig_line.strip() != "":
                        if indent >= required_indent + optional_indent:
                            required_indent += optional_indent
                            if ifdef_depth > 0 and len(indents_in_ifdefs) == ifdef_depth:
                                depths = indents_in_ifdefs[ifdef_depth-1]
                                depths[0] += optional_indent
                            optional_indent = 0
                            optional_indent_reason = ""
                            if indent_inside_ifdef == -1:
                                indent_inside_ifdef = 1
                                first_inside_ifdef = gLineNo[-1]
                            elif indent_inside_ifdef != 1:
                                style("Inconsistent indentation inside `ifdef (compare with line %d)" % first_inside_ifdef )
                    stripped = orig_line.strip()
                    if stripped != "":
                        fixed_line = makeIndent(required_indent+single_line_indent) + stripped + "\n"
                    else:
                        fixed_line = "\n"
                else:
                    if orig_line.endswith("\\\n"):
                        #orig_line = orig_line[:-2] + "\n"
                        if stripped.endswith("\\"):
                            stripped = stripped[:-1].strip()
                    parser = Parser(line)
                    val = parser.lex()
                    bad_indent = False
                    keywd = ""
                    is_vams_compdir = False
                    if parser.isIdentifier():
                        keywd = parser.getString()
                        if keywd in ["end", "endcase", "enddiscipline", "endfunction", "endmodule", "endnature"]:
                            if required_indent >= gSpcPerInd:
                                required_indent -= gSpcPerInd
                            # else extra end, or incorrect indent on a previous line
                            if keywd == "endcase":
                                in_case_block = False
                            elif keywd == "end":
                                rest = parser.peekRestOfLine().strip()
                                if rest.startswith("else if"):
                                    rest = rest[7:]
                                    if len(rest) < 2 or rest[0] != ' ' or rest[1] != '(':
                                        style("Prefer single space after 'if'")
                                        if gFixIndent:
                                            fixed_line = fixSpaces(fixed_line, "else if")
                        elif keywd in ["if", "for", "while", "repeat", "case"]:
                            rest = parser.peekRestOfLine()
                            if len(rest) < 2 or rest[0] != ' ' or rest[1] != '(':
                                style("Prefer single space after '%s'" % keywd)
                                if gFixIndent:
                                    fixed_line = fixSpaces(fixed_line, keywd)
                        elif keywd == "else":
                            rest = parser.peekRestOfLine().strip()
                            if rest.startswith("if"):
                                rest = rest[2:]
                                if len(rest) < 2 or rest[0] != ' ' or rest[1] != '(':
                                    style("Prefer single space after 'if'")
                                    if gFixIndent:
                                        fixed_line = fixSpaces(fixed_line, "if")
                    elif val == ord('@'):
                        keywd = "@"
                    this_indent = required_indent
                    if val == ord('`') and this_indent >= gSpcPerInd:
                        if (indent_ifdef != 0 and stripped.startswith("`else")) or \
                               (indent_ifdef == -1 and stripped.startswith("`endif")):
                            this_indent -= gSpcPerInd
                    elif in_multiline_define:
                        this_indent += special_define_indent
                    previous_single_indent = single_line_indent
                    if single_line_indent > 0:
                        if keywd != "begin":
                            this_indent += single_line_indent
                            if continued_line != "" and stripped.find("begin") > 0:
                                required_indent += single_line_indent
                                single_line_indent = 0
                        if continued_line == "":
                            single_line_indent = 0
                    if indent != this_indent:
                        bad_indent = True
                        #if indent >= this_indent + optional_indent:
                        if indent > this_indent:
                            if indent == this_indent + optional_indent:
                                bad_indent = False
                            required_indent += optional_indent
                            this_indent += optional_indent
                            if optional_indent_reason == "module":
                                indent_module = 1
                            elif optional_indent_reason == "define":
                                if indent_multiline_define == -1:
                                    indent_multiline_define = 1
                                    first_multiline_define = gLineNo[-1]
                                elif indent_multiline_define != 1:
                                    style("Inconsistent indentation of multi-line `define (compare with line %d)" % first_multiline_define )
                                    redo_indent = True
                                indent_this_define = 1
                            elif optional_indent_reason == "ifdef":
                                if ifdef_depth > 0 and len(indents_in_ifdefs) == ifdef_depth:
                                    depths = indents_in_ifdefs[ifdef_depth-1]
                                    depths[0] += optional_indent
                                if indent_inside_ifdef == -1:
                                    indent_inside_ifdef = 1
                                    first_inside_ifdef = gLineNo[-1]
                                elif indent_inside_ifdef != 1:
                                    style("Inconsistent indentation inside `ifdef (compare with line %d)" % first_inside_ifdef )
                            optional_indent = 0
                            optional_indent_reason = ""
                        #elif indent == -1:
                        #    # TAB handled above
                        #    bad_indent = False
                    elif indent == this_indent and optional_indent > 0:
                        if optional_indent_reason == "module":
                            indent_module = 0
                        elif optional_indent_reason == "define":
                            if indent_multiline_define == -1:
                                indent_multiline_define = 0
                                first_multiline_define = gLineNo[-1]
                            elif indent_multiline_define != 0:
                                style("Inconsistent indentation of multi-line `define (compare with line %d)" % first_multiline_define )
                                redo_indent = True
                            indent_this_define = 0
                        elif optional_indent_reason == "ifdef":
                            if val == ord('`') and orig_line.startswith("`else"):
                                pass # only comments since `ifdef
                            elif indent_inside_ifdef == -1:
                                indent_inside_ifdef = 0
                                first_inside_ifdef = gLineNo[-1]
                            elif indent_inside_ifdef > 0:
                                style("Inconsistent indentation inside `ifdef (compare with line %d)" % first_inside_ifdef )
                        optional_indent = 0
                        optional_indent_reason = ""
     
                    may_be_assign_or_contrib = 0
                    if keywd != "":
                        if keywd == "module" or keywd == "macromodule":
                            optional_indent = gSpcPerInd
                            optional_indent_reason = "module"
                        elif keywd in ["analog", "if", "else", "for", "repeat", "while", "@"]:
                            single_line_indent = previous_single_indent + gSpcPerInd
                            num_parens = 0
                            while val != 0:
                                val = parser.lex()
                                if parser.isIdentifier():
                                    str = parser.getString()
                                    if str == "begin" or str == "function":
                                        optional_indent = 0
                                        required_indent += gSpcPerInd
                                        single_line_indent = 0
                                    elif str == "end":
                                        required_indent -= gSpcPerInd
                                elif val == ord('('):
                                    num_parens += 1
                                elif val == ord(')'):
                                    num_parens -= 1
                                elif val == ord(';') and num_parens == 0:
                                    single_line_indent = 0
                        elif keywd in ["begin", "case", "discipline", "function", "nature"]:
                            optional_indent = 0
                            required_indent += gSpcPerInd
                            if keywd == "case":
                                in_case_block = True
                        elif keywd == "end" or (in_case_block and keywd == "default"):
                            # required_indent decreased above for "end"
                            while val != 0:
                                val = parser.lex()
                                if parser.isIdentifier():
                                    if parser.getString() == "else":
                                        single_line_indent = previous_single_indent + gSpcPerInd
                                    elif parser.getString() == "begin":
                                        optional_indent = 0
                                        required_indent += gSpcPerInd
                                        single_line_indent = 0
                        else:
                            val = parser.lex()
                            if val == ord('='):
                                may_be_assign_or_contrib = 1
                            elif val == ord('('):
                                rest = parser.peekRestOfLine()
                                if rest.find("<+"):
                                    may_be_assign_or_contrib = 2
                    elif in_case_block and parser.isNumber():
                        while val != 0:
                            val = parser.lex()
                            if parser.isIdentifier():
                                if parser.getString() == "begin":
                                    optional_indent = 0
                                    required_indent += gSpcPerInd
                    elif val == ord('`'):
                        val = parser.lex()
                        if parser.isIdentifier():
                            keywd = parser.getString()
                            if keywd in gVAMScompdir:
                                bad_indent = False
                                is_vams_compdir = True
                                if gFixIndent and indent > this_indent:
                                    fixed_line = makeIndent(this_indent) + fixed_line.strip() + "\n"
                            if keywd == "define":
                                if required_indent > 0:
                                    if indent == 0:
                                        if indent_define == -1:
                                            indent_define = 0
                                            first_define_indent = gLineNo[-1]
                                        elif indent_define != 0:
                                            style("Inconsistent indentation of `define (compare with line %d)" % first_define_indent )
                                    else:
                                        if indent_define == -1:
                                            indent_define = 1
                                            first_define_indent = gLineNo[-1]
                                        elif indent_define != 1:
                                            style("Inconsistent indentation of `define (compare with line %d)" % first_define_indent )
                                if indent != 0 and indent != required_indent:
                                    if required_indent != 0:
                                        style("Incorrect indent: %d, should be %d (or 0)" % (indent, required_indent), "Incorrect indent")
                                    else:
                                        style("Incorrect indent: %d, should be 0" % indent, "Incorrect indent")
                                if line[-2] == '\\':
                                    optional_indent = gSpcPerInd
                                    optional_indent_reason = "define"
                                    in_multiline_define = True
                                    if indent == 0 and indent != required_indent:
                                        special_define_indent = indent - required_indent
                                    else:
                                        special_define_indent = 0
                            elif keywd == "ifdef" or keywd == "ifndef":
                                # optional indent inside ifdef
                                if indent_inside_ifdef == 1:
                                    required_indent += gSpcPerInd
                                else:
                                    optional_indent += gSpcPerInd
                                    optional_indent_reason = "ifdef"
                                ifdef_depth += 1
                                indents_in_ifdefs.append([required_indent, -1])
                                if indent_ifdef == -1:
                                    # first encounter
                                    if this_indent > 0:
                                        if indent == 0:
                                            indent_ifdef = 0
                                            first_ifdef_indent = gLineNo[-1]
                                            this_indent = 0
                                        else:
                                            indent_ifdef = gSpcPerInd
                                            first_ifdef_indent = gLineNo[-1]
                                            if indent != this_indent:
                                                bad_indent = True
                                    # if this_indent == 0, can't decide yet
                                elif (indent_ifdef == 0 and indent != 0) or \
                                     (indent_ifdef > 0 and indent != this_indent):
                                    style("Inconsistent indentation for `ifdef (compare with line %d)" % first_ifdef_indent, "Incorrect indent")
                                if indent_ifdef == 0:
                                    this_indent = 0
                                fixed_line = makeIndent(this_indent) + fixed_line.strip() + "\n"
                            elif keywd == "else":
                                if indent_inside_ifdef == -1:
                                    # only comments in the `ifdef block
                                    optional_indent += gSpcPerInd
                                    optional_indent_reason = "ifdef"
                                if ifdef_depth > 0 and len(indents_in_ifdefs) == ifdef_depth:
                                    depths = indents_in_ifdefs[ifdef_depth-1]
                                    depths[1] = required_indent
                                    if depths[0] != required_indent:
                                        required_indent = depths[0]
                                if (indent_ifdef == 0 and indent != 0) or \
                                        (indent_ifdef > 0 and indent != this_indent):
                                    bad_indent = True
                                    if indent_ifdef == 0:
                                        this_indent = 0
                            elif keywd == "endif":
                                if ifdef_depth > 0:
                                    if ifdef_depth > 0 and len(indents_in_ifdefs) == ifdef_depth:
                                        depths = indents_in_ifdefs[ifdef_depth-1]
                                        if depths[1] == -1: # no `else
                                            if depths[0] != required_indent:
                                                style("Conditional indent inside `ifdef, but no `else")
                                        elif depths[1] > required_indent:
                                            style("Extra indent in `ifdef compared with `else block")
                                        elif depths[1] < required_indent:
                                            style("Extra indent in `else compared with `ifdef block")
                                        indents_in_ifdefs.pop()
                                    ifdef_depth -= 1
                                    if indent_inside_ifdef > 0:
                                        if this_indent >= gSpcPerInd:
                                            this_indent -= gSpcPerInd
                                        if required_indent >= gSpcPerInd:
                                            required_indent -= gSpcPerInd
                                    elif indent_inside_ifdef == -1:
                                        optional_indent_reason = ""
                                    if optional_indent >= gSpcPerInd:
                                        optional_indent -= gSpcPerInd
                                # else error reported later
                                if (indent_ifdef == 0 and indent != 0) or \
                                        (indent_ifdef > 0 and indent != required_indent):
                                    bad_indent = True
                                    if indent_ifdef == 0:
                                        this_indent = 0
                    if in_multiline_define:
                        if line[-2] != '\\':
                            in_multiline_define = False
                            if indent_this_define > 0:
                                if required_indent >= gSpcPerInd:
                                    required_indent -= gSpcPerInd
                    if bad_indent and continued_line != "" and indent > this_indent:
                        bad_indent = False
                    if bad_indent and indent != -1:
                        style("Incorrect indent: %d, should be %d" % (indent, this_indent), "Incorrect indent")
                    if bad_indent or (indent_this_define == 0 and not is_vams_compdir):
                        if indent_this_define == 0:
                            this_indent += gSpcPerInd
                        fixed_line = makeIndent(this_indent) + fixed_line.strip() + "\n"
                    # try to determine if this line continues
                    if continued_line == "no_semicolon":
                        if stripped != "" and stripped[-1] == ';':
                            continued_line = ""
                    elif continued_line == "parentheses":
                        for ch in line:
                            if ch == '(':
                                unclosed_parens += 1
                            elif ch == ')':
                                unclosed_parens -= 1
                        if unclosed_parens <= 0:
                            unclosed_parens = 0
                            continued_line = ""
                    else:
                        continued_line = ""
                        if stripped != "" and stripped[-1] != ';':
                            if may_be_assign_or_contrib:
                                continued_line = "no_semicolon"
                            else:
                                num_parens = 0
                                for ch in line:
                                    if ch == '(':
                                        num_parens += 1
                                    elif ch == ')':
                                        num_parens -= 1
                                if num_parens > 0:
                                    continued_line = "parentheses"
                                    unclosed_parens = num_parens
                if gFixIndent:
                    fixedLines.append(fixed_line)
                    if fixed_line != orig_line:
                        did_fix = True
            # end of check_style

        lineToParse += line
        if line.endswith("\\\n"):
            lineToParse = lineToParse[:-2] + "\n"
            linesToParse.append("") # preserve line numbering
        else:
            linesToParse.append(lineToParse)
            lineToParse = ""

    # write file, if there are any changes
    if gFixIndent:
        fixfile = filename + ".fixed"
        if did_fix:
            with open(fixfile, "w") as OF:
                OF.writelines(fixedLines)
            print("Wrote %s with fixed indentation" % fixfile)
        else:
            if os.path.exists(fixfile) and os.path.isfile(fixfile):
                print("%s did not need fixes to indentation; removing %s" % (filename, fixfile))
                os.remove(fixfile)
            elif gVerbose:
                print("%s did not need fixes to indentation" % filename)

    # empty the arrays
    fixedLines = []
    linesOfCode = []

    # process the lines
    gLineNo[-1] = 0
    preProcessNatures(linesToParse)
    gLineNo[-1] = 0
    parseLines(linesToParse)

    # done with this file
    gFileName.pop()
    gLineNo.pop()
    if len(gFileName) == 0 and len(gScopeList) > 0:
        error("Missing 'endmodule'")
# end of parseFile




################################################################################
# main routine

class infoAction(argparse.Action):
    def __init__(self, option_strings, dest, nargs=None, help=None):
        argparse.Action.__init__(self, option_strings=option_strings, dest=dest, nargs=nargs, help=help)
        return
    def __call__(self, parser, namespace, values, option_string=None):
        print("""
This program runs checks on Verilog-A models to detect common problems, including:

- hidden state (variables used before assigned)
- unused parameters or variables
- division by zero for parameters (1/parm where parm's range allows 0)
- integer division (1/2 = 0)
- incorrect ddx() usage
- ports without direction and/or discipline
- incorrect access functions for discipline
- unnamed noise sources
- poor coding style
""")
        sys.exit()

class versionAction(argparse.Action):
    def __init__(self, option_strings, dest, nargs=None, help=None):
        argparse.Action.__init__(self, option_strings=option_strings, dest=dest, nargs=nargs, help=help)
        return
    def __call__(self, parser, namespace, values, option_string=None):
        print("""
VAMPyRE version 1.5
""")
        sys.exit()


#
#   Process command line arguments
#

parser = argparse.ArgumentParser(description='Parse and run basic checks on Verilog-A models')
parser.add_argument('main_file',             help='name of top-level Verilog-A file')
parser.add_argument('-a', '--all',           help='equivalent to --max_num 0', action='store_true')
parser.add_argument('-b', '--binning',       help='analyze binning equations', action='store_true')
parser.add_argument('-d', '--debug',         help='turn on debug mode', action='store_true')
parser.add_argument('-f', '--fix_indent',    help='fix indentation problems', action='store_true')
parser.add_argument('-i', '--info',          help='print detailed usage information', action=infoAction, nargs=0)
parser.add_argument('-I', '--inc_dir',       help='add search path for include', action='append', default=[])
parser.add_argument(      '--max_num',       help='maximum number of each type of error/warning message', type=int, default=5)
parser.add_argument('-n', '--no_style',      help='turn off style checking', action='store_false')
parser.add_argument('-s', '--indent_spaces', help='number of spaces to indent', type=int, default=4)
parser.add_argument('-v', '--verbose',       help='turn on verbose mode', action='store_true')
parser.add_argument(      '--version',       help='print version number', action=versionAction, nargs=0)
args       = parser.parse_args()
gIncDir    = args.inc_dir
gDebug     = args.debug
gBinning   = args.binning
gFixIndent = args.fix_indent
if args.all:
    gMaxNum = 0
else:
    gMaxNum = args.max_num
gSpcPerInd = args.indent_spaces
gStyle     = args.no_style
gVerbose   = args.verbose
if gFixIndent and not gStyle:
    fatal("Cannot fix indentation without checking style")
parseFile(args.main_file)

# Summary
print("\nModule: %s" % gModuleName)
if len(gParameters):
    InstParms = {}
    ModelParms = {}
    for (k,p) in gParameters.items():
        if p.inst:
            InstParms[k] = p
        else:
            ModelParms[k] = p
    if len(InstParms):
        if gVerbose or len(InstParms) < 20:
            print("\nInstance parameters (%d):" % len(InstParms))
            print_list(InstParms.keys(), "    ", 70)
        else:
            print("\nInstance parameters (%d)" % len(InstParms))
        if gVerbose or len(ModelParms) < 20:
            print("\nModel parameters (%d):" % len(ModelParms))
            print_list(ModelParms.keys(), "    ", 70)
        else:
            print("\nModel parameters (%d)" % len(ModelParms))
    else:
        if gVerbose or len(gParameters) < 40:
            print("\nParameters (%d):" % len(gParameters))
            print_list(gParameters.keys(), "    ", 70)
        else:
            print("\nParameters (%d)" % len(gParameters))

if len(gVariables):
    OpPtVars = {}
    vars_to_pop = []
    for (k,v) in gVariables.items():
        if v.assign == 0 or v.assign < -1:
            vars_to_pop.append(v.name)
        elif v.oppt:
            OpPtVars[k] = v
    for var in vars_to_pop:
        gVariables.pop(var)
    if len(OpPtVars):
        if gVerbose or len(OpPtVars) < 50:
            print("\nOperating-point variables (%d):" % len(OpPtVars))
            print_list(OpPtVars.keys(), "    ", 70)
        else:
            print("\nOperating-point variables (%d)" % len(OpPtVars))
    else:
        print("\nNo operating-point variables")
    if gVerbose:
        print("Variables (%d):" % len(gVariables))
        print_list(gVariables.keys(), "    ", 70)

if len(gHiddenState):
    print("Possible hidden-state variables:")
    print_list(gHiddenState.keys(), "    ", 70)
