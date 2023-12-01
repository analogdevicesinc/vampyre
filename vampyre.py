#!/usr/bin/env python
""" VAMPyRE
    Verilog-A Model Pythonic Rule Enforcer
    version 1.9.3, 01-Dec-2023

    intended for checking for issues like:
     1) hidden state (variables used before assigned)
        or conditionally-assigned operating point values
     2) unused parameters, variables, ports, or nodes
     3) division by zero for parameters (1/parm where parm's range allows 0)
        or domain errors for ln() and sqrt()
     4) integer division (1/2 = 0)
     5) incorrect ddx() usage
     6) ports without direction and/or discipline
     7) incorrect access functions for discipline
     8) unnamed noise sources
     9) use of features not appropriate for compact models
    10) various issues of poor coding style
    11) compliance with CMC Verilog-A Code Standards:
       - use of lowercase identifiers
       - proper use of tref and dtemp
       - proper use of the multiplicity attribute
       - proper implementation of mult_i, mult_q, mult_fn
    12) misuse of limexp
    13) various problems with binning equations and units
    14) nonlinear ddt() expressions and C(V)*ddt(V) formulations
    15) superfluous assignments

    Copyright (c) 2023 Analog Devices, Inc.

    Licensed under Educational Community License, Version 2.0 (the "License");
    you may not use this file except in compliance with the License. You may
    obtain a copy of the license at http://opensource.org/licenses/ECL-2.0

    Unless required by applicable law or agreed to in writing, software
    distributed under the License is distributed on an "AS IS" BASIS,
    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
    See the License for the specific language governing permissions and
    limitations under the License.
"""


################################################################################
# setup

import os
import sys
import argparse
from copy import deepcopy

################################################################################
# global variables

gVersionNumber = "1.9.3"
gModuleName    = ""
gNatures       = {}
gDisciplines   = {}
gParameters    = {}
gVariables     = {}
gHiddenState   = {}
gUserFunctions = {}
gAccessFuncs   = {}
gPortnames     = {}  # terminals
gPortlist      = ""
gNodenames     = {}  # internal nodes
gBranches      = {}
gMultScaled    = {}
gContribs      = {}  # contributions
gBlocknames    = {}
gMacros        = {}
gIfDefStatus   = []
gScopeList     = []
gLoopTypes     = []
gConditions    = []  # if-conditions in play
gCondBiasDep   = []  # whether conditions are bias-dependent:
                     # 0 (no), 1 (from if cond), 2 (yes), 3 (deriv error), 4 (ddx)
gAnalogBlock   = []  # analog or analog initial
gMissingConstantsFile = ""
gStatementInCurrentBlock = False
gLastDisplayTask  = []
gLastLineWasEvent = False
gCurrentFunc      = 0
gFileName         = []
gLineNo           = []
gSubline          = 0
gIncDir           = []
gDebug            = False
gDefines          = []
gFixIndent        = False
gMaxNum           = 5
gSpcPerInd        = 4
gStyle            = False
gVerbose          = False
gBinning          = False
gBinningPatterns  = [[], []]
gBinningFactors   = {}
gNoiseTypes       = []
gSuperfluous      = False
gCompactModel     = True
gCheckMultFactors = 0
gUsesDDT          = False
gPrDefVals        = False
gPreProcess       = False

# dictionaries to track how many of each type of message
gErrorMsgDict   = {}
gWarningMsgDict = {}
gNoticeMsgDict  = {}
gStyleMsgDict   = {}

# tokens
# (0-255 are ascii and latin-1)
gTokenUnused     = 256
gTokenNumber     = 257
gTokenIdentifier = 258
gTokenString     = 259

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
gMathFunctions = {"abs": 1, "limexp": 1, "log": 1, "$log10": 1, "min": 2, "max": 2}
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
gUnitsMultiply = ["A", "Amp", "Amps", "F", "Farads", "C", "coul", "Coulomb", "Coulombs",
                  "W", "Watt", "Watts", "A/V", "S", "mho", "A*V", "A*V/K", "A/K", "C/K",
                  "J/K", "s*W/K", "A/V^2", "A^2/Hz"]
gUnitsDivide   = ["Ohm", "Ohms", "V/A", "K/W", "V^2/Hz"]

# MULT_* factors
gMultFactors = ["mult_i", "mult_q", "mult_fn", "MULT_I", "MULT_Q", "MULT_FN"]

# for use in asString()
gExprTypeNoExtraParen = ["==", "!=", "<", ">", "<=", ">=", "&&", "||", "!", "NOT", "NAME", "NUMBER", "FUNCCALL"]

################################################################################
# classes

class Nature:
    """ Verilog-AMS Nature """
    def __init__(self, name, defn):
        self.name = name
        self.defined = defn  # False during pre-processing
        self.units = ""
        self.access = ""
        self.idt_nature = ""


class Discipline:
    """ Verilog-AMS Discipline """
    def __init__(self, name):
        self.name = name
        self.potential = ""
        self.flow = ""
        self.domain = ""


class Parameter:
    """ Verilog-AMS Parameter """
    def __init__(self, name, type, file, line):
        self.name = name
        self.type = type         # real, integer, string
        self.declare = [file, line]
        self.defv = 0
        self.units = 0
        self.range = []          # array parameters
        self.value_range = "nb"  # nb, cc, co, oc, oo, cz, oz, sw
        self.max = float("inf")
        self.min = float("-inf")
        self.exclude = []
        self.inst = False        # instance parameter (or both)
        self.model = False       # model parameter (or both)
        self.scale = ""          # linear or quadratic
        self.used = 0            # number of times used
        self.is_alias = False    # True for aliasparam (should not be used)

    def getValueRangeStr(self):
        ret = ""
        if self.value_range == "nb":
            ret = "(-inf:inf)"
        elif self.value_range in ["cc", "co", "oc", "oo"]:
            if self.value_range[0] == "c":
                lend = "["
            else:
                lend = "("
            if self.value_range[1] == "c":
                rend = "]"
            else:
                rend = ")"
            ret = lend + self.min.asString() + ":" + self.max.asString() + rend
        #elif self.value_range == "cz":
        #    ret = "[0:inf)"
        #elif self.value_range == "oz":
        #    ret = "(0:inf)"
        #elif self.value_range == "sw":
        #    ret = "[0:1]"
        return ret


class Variable:
    """ Verilog-AMS variable """
    def __init__(self, name, type, oppt, file, line):
        self.name = name
        self.type = type              # real, integer, string
        self.range = []
        self.oppt = oppt              # operating-point variable
        self.units = ""               # units (only for oppt)
        self.declare = [file, line]
        self.assign = -1              # line where value assigned
        self.prev_assign = []         # details about previous assignment
        self.non_zero_assign = False  # true if assigned a non-zero value
        self.conditions = []          # conditions when assigned
        self.accesses = []            # details about assignments and uses
        self.used = False             # whether used
        self.bias_dep = 0             # bias-dependent: 0=no, 1=from if cond, 2=yes, 3=deriv error, 4=ddx
        self.biases = []              # node voltages (for ddx check)
        self.zeros = []               # parameters that, when zero, make the variable zero
        self.factors = []             # factors (for MULT_* checks)
        self.ddt_use = 0              # involves ddt: 0=no, 1=linear, 2=nonlinear
        self.func_name_arg = ["", 0]  # tracing whether function output arg is used


class Access:
    """ Variable access info """
    def __init__(self, type, file, line, subl, cond, zini=0):
        self.type = type              # "SET" or "USE"
        self.file = file
        self.line = line
        self.subl = subl
        self.cond = cond              # conditions in effect
        self.loop = False             # if used in loop
        self.zini = zini              # initialization to zero
        self.used = False


class Port:
    """ Verilog-AMS port """
    def __init__(self, name, file, line):
        self.name = name
        self.type = ""
        self.direction = ""   # inout, input, output
        self.discipline = ""  # electrical, thermal, ...
        self.declare = [file, line]
        self.is_bus = False
        self.msb = 0
        self.lsb = 0
        self.mult = False     # whether contributions to this node need mult scaling
        self.used = False


class Branch:
    """ Verilog-AMS branch """
    def __init__(self, name):
        self.name = name
        self.node1 = ""
        self.node2 = ""
        self.discipline = ""  # electrical, thermal, ...
        self.conds = []       # conditions under which contributed to
        self.lhs_flow = 0     # flow contrib: 0=no, 1=yes, 2=bias-dep condition
        self.lhs_pot  = 0     # potential contrib: 0=no, 1=yes, 2=bias-dep condition
        self.mult = False     # whether contributions to this branch need mult scaling
        self.used = False


class Contribution:
    """ Verilog-AMS contribution """
    def __init__(self, acc):
        self.acc = acc
        self.func = ""
        self.node1 = ""
        self.node2 = ""
        self.branch = ""
        self.discipline = ""


class Macro:
    """ Verilog-AMS macro """
    def __init__(self, name, text):
        self.name = name
        self.text = text
        self.args = []
        self.used = False      # whether used


class Function:
    """ Verilog-AMS function """
    def __init__(self, name, type):
        self.name = name
        self.type = type       # real, integer, string
        self.args = []
        self.inputs = []
        self.outputs = []
        self.used = False      # whether used
        self.outarg_used = []  # whether output args are used in any call to the function


class Expression:
    """ Verilog-AMS expression """
    def __init__(self, type):
        self.type = type       # NAME, STRING, NUMBER, FUNCCALL, or operator
        self.e1 = 0
        self.e2 = 0
        self.e3 = 0
        self.number = 0
        self.is_int = False    # is the expression an integer?
        self.args = []

    def getDependencies(self, assign_context, branch_contrib):
        deps = []
        if self.type == "NAME":
            deps = [self.e1]
        elif self.type in ["NUMBER", "STRING"]:
            deps = []
        elif self.type == "ARRAY":
            deps = []
            for arg in self.args:
                deps += arg.getDependencies(assign_context, branch_contrib)
        elif self.type == "FUNCCALL":
            if self.e1 in gAccessFuncs:
                if len(gAnalogBlock) > 0 and gAnalogBlock[-1] == "initial":
                    error("Branch access in analog initial block")
                nargs = len(self.args)
                if nargs > 0:
                    nname1 = self.args[0].e1
                    nname2 = ""
                    if nname1 in gPortnames:
                        dname = gPortnames[nname1].discipline
                        gPortnames[nname1].used = True
                    elif nname1 in gNodenames:
                        dname = gNodenames[nname1].discipline
                        gNodenames[nname1].used = True
                    elif nname1 in gBranches:
                        dname = gBranches[nname1].discipline
                        gBranches[nname1].used = True
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
                            error("Incorrect access function '%s' for %s '%s'"
                                  % (self.e1, dname, nname1))
                    elif dname == "" and not branch_contrib:
                        # if branch_contrib, error was reported already
                        error("Cannot determine discipline of '%s'" % nname1)
                    if nargs > 2 or (nargs > 1 and (nname1 in gBranches or self.args[1].e1 in gBranches)):
                        error("Invalid potential or flow access")
                    elif nargs == 2:
                        nname2 = self.args[1].e1
                        dname2 = ""
                        if nname2 in gPortnames:
                            dname2 = gPortnames[nname2].discipline
                            gPortnames[nname2].used = True
                        elif nname2 in gNodenames:
                            dname2 = gNodenames[nname2].discipline
                            gNodenames[nname2].used = True
                        if dname2 != "" and dname2 != dname and not branch_contrib:
                            # if branch_contrib, error was reported already
                            error("Nodes '%s' and '%s' belong to different disciplines"
                                  % (nname1, nname2))
                    if nname1 != "":
                        if nname1 in gBranches:
                            dep = self.e1 + "(" + nname1 + ")"
                            deps.append(dep)
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
            elif self.e1 in ["$port_connected", "$param_given"]:
                if len(self.args) == 1:
                    nname = self.args[0].e1
                    if self.e1 == "$port_connected":
                        if nname in gPortnames:
                            gPortnames[nname].used = True
                        else:
                            error("Invalid argument '%s' to %s, must be a port" % (nname, self.e1))
                    elif self.e1 == "$param_given":
                        if nname in gParameters:
                            gParameters[nname].used += 1
                        else:
                            error("Invalid argument '%s' to %s, must be a parameter" % (nname, self.e1))
                else:
                    error("Incorrect number of arguments to %s" % self.e1)
            elif self.e1 == "limexp":
                if len(self.args) == 1:
                    arg = self.args[0]
                    deps = arg.getDependencies(assign_context, branch_contrib)
                    [bias_dep, biases, ddt] = checkDependencies(deps, "Call to function '%s' depends on"
                                                                % self.e1, 0, False, True, False)
                    if not bias_dep:
                        warning("Call to %s(%s), but argument is not bias-dependent"
                                % (self.e1, arg.asString()))
                # else: warning issued elsewhere
            elif self.e1 == "$limit":
                if len(self.args) == 0:
                    error("Require at least one argument to %s" % self.e1)
                else:
                    access = self.args[0]
                    deps += access.getDependencies(assign_context, branch_contrib)
                    if access.type != "FUNCCALL" or access.e1 not in gAccessFuncs:
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
                        if fname not in ["\"pnjlim\"", "\"fetlim\""]:
                            notice("Non-standard limiting function %s" % fname)
                    else:
                        warning("Second argument to $limit should be a string or identifier")
                    for arg in self.args[2:]:
                        argdeps = arg.getDependencies(assign_context, branch_contrib)
                        [bias_dep, biases, ddt] = checkDependencies(argdeps,
                                                                    "Call to function '%s' depends on"
                                                                    % self.e1, 0, False, True, False)
                        deps += argdeps
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
                        [bias_dep, biases, ddt] = checkDependencies(argdeps, "Call to ddx() depends on",
                                                                    0, False, True, False)
                        if ddt:
                            error("Call to ddx() where argument %s involves ddt()" % arg.asString())
                        if accstr not in biases:
                            warning("Call to ddx(), but %s does not depend on %s"
                                    % (arg.asString(), acc.asString()))
                            if gDebug:
                                print("    Depends on: %s" % biases)
                # else error printed during parsing
            elif assign_context and self.e1 in ["white_noise", "flicker_noise", "noise_table", "noise_table_log"]:
                warning("Small-signal noise source '%s' in assignment" % self.e1)
            elif branch_contrib and self.e1 in ["white_noise", "flicker_noise", "noise_table", "noise_table_log",
                                                "laplace_zp", "laplace_zd", "laplace_np", "laplace_nd",
                                                "zi_zp", "zi_zd", "zi_np", "zi_nd"]:
                # ignore dependencies in noise functions and laplace and z operators
                pass
            else:
                out_arg_pos = []
                if (assign_context or branch_contrib) and self.e1 == "ddt":
                    deps = ["ddt"]
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
                        pass  # null argument
                    else:
                        deps += arg.getDependencies(assign_context, branch_contrib)
                if assign_context and len(out_arg_pos) > 0:
                    [bias_dep, biases, ddt] = checkDependencies(deps, "Call to function '%s' depends on"
                                                                      % self.e1, 0, False, False, True)
                    for i in range(len(self.args)):
                        if i in out_arg_pos:
                            arg = self.args[i]
                            markVariableAsSet(arg.e1, None, [], bias_dep, Expression("NOTHING"),
                                              biases, 0, False, True, self.e1, i)

        elif self.type in ["!", "~"]:
            deps = self.e1.getDependencies(assign_context, branch_contrib)
        elif self.type in ["+", "-"]:
            dep1 = self.e1.getDependencies(assign_context, branch_contrib)
            if self.e2:
                dep2 = self.e2.getDependencies(assign_context, branch_contrib)
            else:  # unary +/-
                dep2 = []
            deps = dep1 + dep2
        elif self.type in ["*", "/", "%", "**",
                           "==", "!=", "<", ">", "<=", ">=",
                           "&&", "||", "&", "|", "^"]:
            dep1 = self.e1.getDependencies(assign_context, branch_contrib)
            dep2 = self.e2.getDependencies(assign_context, branch_contrib)
            deps = dep1 + dep2
            if self.type == "*" and "ddt" in deps:
                # check for C(V) * ddt(V)
                ddt_biases = []
                if self.e1.type == "FUNCCALL" and self.e1.e1 == "ddt":
                    ddt_biases = checkDependencies(dep1, "", 0, False, False, False)[1]
                    cdep = dep2
                elif self.e2.type == "FUNCCALL" and self.e2.e1 == "ddt":
                    ddt_biases = checkDependencies(dep2, "", 0, False, False, False)[1]
                    cdep = dep1
                if ddt_biases:
                    biases = checkDependencies(cdep, "", 0, False, False, False)[1]
                    bad_biases = []
                    for bias in biases:
                        if bias in ddt_biases:
                            bad_biases.append(bias)
                    if bad_biases:
                        bias = simplifyBiases(bad_biases)
                        error("Non-charge-conserving f(%s) * ddt(%s)" % (bias, bias))
        elif self.type in ["<<", ">>", "<<<", ">>>", "===", "!==",
                           "^~", "~^", "~&", "~|"]:
            # not expected in Verilog-A
            dep1 = self.e1.getDependencies(assign_context, branch_contrib)
            dep2 = self.e2.getDependencies(assign_context, branch_contrib)
            deps = dep1 + dep2
        elif self.type in ["?:"]:
            dep1 = self.e1.getDependencies(assign_context, branch_contrib)
            dep2 = self.e2.getDependencies(assign_context, branch_contrib)
            dep3 = self.e3.getDependencies(assign_context, branch_contrib)
            deps = dep1 + dep2 + dep3
        else:  # pragma: no cover
            fatal("Unhandled expression type '%s' in getDependencies" % self.type)
        return deps

    def ddtCheck(self):
        ret = 0
        if self.type == "NAME":
            if self.e1 in gParameters:
                ret = 0
            else:
                vn = findVariableInScope(self.e1)
                if vn in gVariables:
                    ret = gVariables[vn].ddt_use
                else:  # pragma: no cover
                    fatal("Unexpected identifier '%s' in ddtCheck" % self.e1)
        elif self.type in ["NUMBER", "STRING"]:
            ret = 0
        elif self.type == "FUNCCALL":
            if self.e1 in gAccessFuncs:
                ret = 0
            elif self.e1 == "ddt":
                ret = 1 + self.args[0].ddtCheck()
            else:
                for arg in self.args:
                    ddt1 = arg.ddtCheck()
                    if ddt1 > 0:
                        ret = 2  # nonlinear
        elif self.type in ["+", "-"]:
            ddt1 = self.e1.ddtCheck()
            if self.e2:
                ddt2 = self.e2.ddtCheck()
            else:
                ddt2 = 0
            if ddt1 > ddt2:
                ret = ddt1
            else:
                ret = ddt2
        elif self.type in ["*", "/", "%", "**"]:
            ddt1 = self.e1.ddtCheck()
            ddt2 = self.e2.ddtCheck()
            if ddt1 == 0 and ddt2 == 0:
                ret = 0
            elif ddt1 == 1 and ddt2 == 0 and self.type in ["*", "/"]:
                ret = 1
            elif ddt1 == 0 and ddt2 == 1 and self.type == "*":
                ret = 1
            else:
                ret = 2  # assume nonlinear
        elif self.type in ["?:"]:
            ddt1 = self.e1.ddtCheck()
            ddt2 = self.e2.ddtCheck()
            ddt3 = self.e3.ddtCheck()
            if ddt1 == 0:
                if ddt2 > ddt3:
                    ret = ddt2
                else:
                    ret = ddt3
            else:
                error("ddt() dependence in condition for ?: operator")
        elif self.type in ["!", "~"]:
            error("ddt() dependence in operand for '%s' operator" % self.type)
        elif self.type in ["==", "!=", "<", ">", "<=", ">=",
                           "&&", "||", "&", "|", "^",
                           "<<", ">>", "<<<", ">>>", "===", "!==",
                           "^~", "~^", "~&", "~|"]:
            error("ddt() dependence in operands for '%s' operator" % self.type)
        else:  # pragma: no cover
            fatal("Unhandled expression type '%s' in ddtCheck" % self.type)
        return ret

    def asString(self):
        ret = ""
        if self.type in ["NAME", "STRING"]:
            ret = self.e1
        elif self.type == "NUMBER":
            if isinstance(self.e1, str):  # original parsed string
                ret = self.e1
            elif self.is_int:
                ret = str(int(self.number))
            else:
                ret = str(self.number)
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
        elif self.type in ["!", "~"]:
            ret = self.type + "(" + self.e1.asString() + ")"
        elif self.type == "NOT":
            ret = "NOT(" + self.e1.asString() + ")"
        elif self.type in ["+", "-"]:
            if self.e2:
                ret = "(" + self.e1.asString() + self.type + self.e2.asString() + ")"
            else:  # unary +/-
                ret = self.type + "(" + self.e1.asString() + ")"
        elif self.type == "&&":
            if self.e1.type in gExprTypeNoExtraParen and self.e2.type in gExprTypeNoExtraParen:
                ret = self.e1.asString() + "&&" + self.e2.asString()
            else:
                ret = "(" + self.e1.asString() + self.type + self.e2.asString() + ")"
        elif self.type in ["*", "/", "%", "**",
                           "==", "!=", "<", ">", "<=", ">=",
                           "||", "&", "|", "^"]:
            ret = "(" + self.e1.asString() + self.type + self.e2.asString() + ")"
        elif self.type in ["<<", ">>", "<<<", ">>>", "===", "!==",
                           "^~", "~^", "~&", "~|"]:
            # not expected in Verilog-A
            ret = "(" + self.e1.asString() + self.type + self.e2.asString() + ")"
        elif self.type in ["?:"]:
            ret = "(" + self.e1.asString() + ")?(" + self.e2.asString() \
                   + "):(" + self.e3.asString() + ")"
        elif self.type == "ARRAY":
            ret = "'{"
            first = True
            for arg in self.args:
                if first:
                    first = False
                else:
                    ret += ","
                ret += arg.asString()
            ret += "}"
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
        #elif self.type == "PORT_FLOW":
        #    return "<" + self.e1 + ">"
        #elif self.type == "NOTHING":
        #    return "TODO"
        else:  # pragma: no cover
            fatal("Unhandled expression type '%s' in asString" % self.type)
        return ret


class Parser:
    """ Parser for Verilog-AMS source code """
    def __init__(self, line):
        self.line = line
        self.chpt = 0
        self.end  = len(line)
        self.token = 0
        self.number = 0
        self.ternary_cond = Expression("NOTHING")
        self.is_int = False
        self.escaped = False
        self.string = ""
        self.parse_units = False

    def isNumber(self):
        return self.token == gTokenNumber

    def isIdentifier(self):
        return self.token == gTokenIdentifier

    def isString(self):
        return self.token == gTokenString

    def isEscaped(self):
        return self.escaped

    def getNumber(self):
        return self.number

    def getString(self):
        return self.string

    def getRestOfLine(self):
        retstr = self.line[self.chpt:]
        self.chpt = self.end
        return retstr

    def peekRestOfLine(self):
        return self.line[self.chpt:]

    def eatSpace(self):
        while self.peekChar().isspace():
            self.getChar()

    def peekToken(self):
        ret = ""
        pcp = self.chpt
        while pcp < self.end and self.line[pcp].isspace():
            pcp += 1
        if pcp < self.end:
            ret = self.line[pcp]
        return ret

    def peekChar(self):
        ret = ""
        if self.chpt < self.end:
            ret = self.line[self.chpt]
        return ret

    def getChar(self):
        ret = ""
        if self.chpt < self.end:
            ch = self.line[self.chpt]
            self.chpt += 1
            ret = ch
        return ret

    def ungetChar(self, ch):
        if self.chpt > 0:
            self.chpt -= 1
            if self.line[self.chpt] != ch:  # pragma: no cover
                fatal("Parse failure in ungetChar")

    def ungetChars(self, chars):
        for i in range(len(chars), 0, -1):
            ch = chars[i-1]
            self.ungetChar(ch)

    def lexNumber(self):
        value = 0
        string = ""
        sci_not = ""
        is_int = True
        while self.peekChar().isdigit():
            string += self.getChar()
        if self.peekChar() == '.':
            string += self.getChar()
            is_int = False
            while self.peekChar().isdigit():
                string += self.getChar()
        if self.peekChar() == '.':  # pragma: no cover
            error("Invalid number")
        ch = self.peekChar()
        if ch in ['e', 'E']:
            sci_not = self.getChar()
            ch = self.peekChar()
            last_ch = ""
            if ch in ['+', '-']:
                sci_not += self.getChar()
            while self.peekChar().isdigit():
                last_ch = self.getChar()
                sci_not += last_ch
            if last_ch.isdigit():
                string += sci_not
                is_int = False
            else:
                self.ungetChars(sci_not)
                sci_not = ""
        value = float(string)

        # if not sci not, look for SI prefixes
        if sci_not == "":
            ch = self.peekChar()
            if ch == 'T':
                value *= 1e12
                string += self.getChar()
            elif ch == 'G':
                value *= 1e9
                string += self.getChar()
            elif ch == 'M':
                value *= 1e6
                string += self.getChar()
            elif ch in ['K', 'k']:
                value *= 1e3
                string += self.getChar()
            elif ch == 'm':
                value *= 1e-3
                string += self.getChar()
                is_int = False
            elif ch == 'u':
                value *= 1e-6
                string += self.getChar()
                is_int = False
            elif ch == 'n':
                value *= 1e-9
                string += self.getChar()
                is_int = False
            elif ch == 'p':
                value *= 1e-12
                string += self.getChar()
                is_int = False
            elif ch == 'f':
                value *= 1e-15
                string += self.getChar()
                is_int = False
            elif ch == 'a':
                value *= 1e-18
                string += self.getChar()
                is_int = False

        return [value, is_int, string]
    # end of lexNumber

    def lexName(self):
        string = ""
        while self.peekChar() == '$':
            string += self.getChar()
        while self.peekChar().isalnum() or self.peekChar() == '_':
            string += self.getChar()
        return string

    def lexString(self, qstr):
        string = self.getChar()
        ch = ""
        escaped = False
        while ch != qstr or escaped:
            ch = self.getChar()
            if ch == "":
                error("Missing end-quote %s" % qstr)
                if string[-1] == ';':
                    self.ungetChar(';')
                    string = string[0:-1]
                break
            if escaped:
                string += ch
                ch = " "
                escaped = False
            elif ch == "\\":
                escaped = True
            else:
                string += ch
        return string

    def lex(self):
        self.eatSpace()
        ch = self.peekChar()
        if ch == "":
            self.token = 0
        elif ch.isdigit():
            number = self.lexNumber()
            self.number = number[0]
            self.is_int = number[1]
            self.string = number[2]
            self.token  = gTokenNumber
        elif ch.isalpha() or ch == '_' or ch == '$':
            self.string = self.lexName()
            self.token = gTokenIdentifier
        elif ch == '\\':
            # escaped identifier
            self.string = ""
            ch = self.getChar()
            ch = self.getChar()
            while ch != "" and not ch.isspace():
                self.string += ch
                ch = self.getChar()
            if self.string != "":
                self.token = gTokenIdentifier
                self.escaped = True
            else:
                ch = '\\'
                self.token = ord(ch)
        elif ch == '"':
            self.string = self.lexString(ch)
            self.token = gTokenString
        else:
            ch = self.getChar()
            if ord(ch) > 127 and not self.parse_units:
                error("Encountered non-ASCII character '%s' (decimal %d)" % (ch, ord(ch)))
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
            if self.chpt < self.end-1:
                nextch = self.line[self.chpt+1]
            else:
                nextch = ""
            if ch == "+" and nextch == "+":
                oper += nextch
            elif ch == "-" and nextch == "-":
                oper += nextch
            elif ch == "*" and nextch == "*":
                oper += nextch
            elif ch in "+-*/" and nextch == "=":  # pragma: no cover
                oper += nextch
            elif ch == "&" and nextch == "&":
                oper += nextch
            elif ch == "|" and nextch == "|":
                oper += nextch
            elif ch == "=" and nextch == "=":
                oper += nextch
            elif ch == "!" and nextch == "=":
                oper += nextch
            elif ch == "<" and nextch == "<":
                oper += nextch
            elif ch == "<" and nextch == "=":
                oper += nextch
            elif ch == ">" and nextch == ">":
                oper += nextch
            elif ch == ">" and nextch == "=":
                oper += nextch
            elif ch == "~" and nextch == "^":
                oper += nextch
            elif ch == "^" and nextch == "~":
                oper += nextch
            elif ch == "~" and nextch == "&":
                oper += nextch
            elif ch == "~" and nextch == "|":
                oper += nextch
            if len(oper) == 2 and self.chpt < self.end-2:
                nextnext = self.line[self.chpt+2]
                if oper == "<<" and nextnext == "<":
                    oper = "<<<"
                elif oper == ">>" and nextnext == ">":
                    oper = "<<<"
                elif oper == "==" and nextnext == "=":
                    oper = "==="
                elif oper == "!=" and nextnext == "=":
                    oper = "!=="
            # error printed by getOper() if it's used
            #if oper in ["++", "+=", "--", "-=", "*=", "/=", "<<", ">>", "<<<", ">>>", \
            #            "===", "!==", "^~", "~^", "~&", "~|"]:
            #    error("Operator '%s' not valid in Verilog-A" % oper)
        return oper

    def getOper(self):
        oper = ""
        ch = self.peekChar()
        if ch in "+-*/%><!&|=~^?:":
            oper = self.getChar()
            nextch = self.peekChar()
            if ch == "+" and nextch == "+":
                oper += self.getChar()
            elif ch == "-" and nextch == "-":
                oper += self.getChar()
            elif ch == "*" and nextch == "*":
                oper += self.getChar()
            elif ch in "+-*/" and nextch == "=":  # pragma: no cover
                oper += self.getChar()
            elif ch == "&" and nextch == "&":
                oper += self.getChar()
            elif ch == "|" and nextch == "|":
                oper += self.getChar()
            elif ch == "=" and nextch == "=":
                oper += self.getChar()
            elif ch == "!" and nextch == "=":
                oper += self.getChar()
            elif ch == "<" and nextch == "<":
                oper += self.getChar()
            elif ch == "<" and nextch == "=":
                oper += self.getChar()
            elif ch == ">" and nextch == ">":
                oper += self.getChar()
            elif ch == ">" and nextch == "=":
                oper += self.getChar()
            elif ch == "~" and nextch == "^":
                oper += self.getChar()
            elif ch == "^" and nextch == "~":
                oper += self.getChar()
            elif ch == "~" and nextch == "&":
                oper += self.getChar()
            elif ch == "~" and nextch == "|":
                oper += self.getChar()
            if len(oper) == 2 and self.chpt < self.end:
                nextnext = self.peekChar()
                if oper == "<<" and nextnext == "<":
                    oper += self.getChar()
                elif oper == ">>" and nextnext == ">":
                    oper += self.getChar()
                elif oper == "==" and nextnext == "=":
                    oper += self.getChar()
                elif oper == "!=" and nextnext == "=":
                    oper += self.getChar()
            if oper in ["++", "+=", "--", "-=", "*=", "/=", "<<", ">>", "<<<", ">>>",
                        "===", "!==", "^~", "~^", "~&", "~|"]:
                if not self.parse_units:
                    error("Operator '%s' not valid in Verilog-A" % oper)
        self.eatSpace()
        return oper

    # expression parsing - watch for order of operations

    def parseArgList(self, is_cond=False, allow_repl=False):
        args = []
        expr = self.getExpression(is_cond)
        if expr.type == "NOTHING" and self.peekChar() == ')':
            return []
        self.eatSpace()
        ch = self.peekChar()
        if ch == '{' and allow_repl:
            self.getChar()
            sub_args = self.parseArgList()
            self.eatSpace()
            if self.peekChar() == '}':
                self.getChar()
                self.eatSpace()
            else:
                error("Missing '}'")
            if not expr.is_int:
                error("Non-integer multiplier '%s' in replication" % expr.asString())
            repl = 1
            if expr.type == "NUMBER":
                repl = int(expr.number)
            elif expr.type == "NAME" and expr.e1 in gParameters:
                gParameters[expr.e1].used += 1
                defv = gParameters[expr.e1].defv
                if defv.type == "NUMBER":
                    repl = int(defv.number)
                else:
                    warning("Cannot determine value of multiplier '%s' in replication" % expr.e1)
            else:
                error("Invalid multiplier '%s' in replication" % expr.asString())
            args = repl * sub_args
        else:
            args.append(expr)
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
            while self.peekChar() == '.':
                # hierarchical identifier
                expr.e1 += self.getChar()
                val = self.lex()
                if self.isIdentifier():
                    expr.e1 += self.getString()
                else:  # pragma: no cover
                    error("Invalid hierarchical identifier %s" % expr.e1)
        elif self.isNumber():
            expr = Expression("NUMBER")
            expr.number = self.getNumber()
            expr.is_int = self.is_int
            expr.e1     = self.getString()
        elif self.isString():
            expr = Expression("STRING")
            expr.e1 = self.getString()
        elif val == ord('('):
            expr = self.getExpression(is_cond)
            self.eatSpace()
            if self.parse_units and self.peekChar().isalpha():
                expr1 = expr
                expr2 = self.getExpression()
                expr = Expression("*")
                expr.e1 = expr1
                expr.e2 = expr2
                self.eatSpace()
            if self.peekChar() == ')':
                self.getChar()
                self.eatSpace()
            elif not self.parse_units:
                error("Missing ')'")
        elif val == ord('\''):
            if self.peekChar() == '{':
                self.getChar()
                expr = Expression("ARRAY")
                expr.args = self.parseArgList(is_cond, True)
                self.eatSpace()
                if self.peekChar() == '}':
                    self.getChar()
                    self.eatSpace()
                elif not self.parse_units:
                    error("Missing '}'")
            else:  # pragma: no cover
                self.ungetChar(chr(val))
                expr = Expression("NOTHING")
        elif val == ord('{'):
            expr = Expression("ARRAY")
            expr.args = self.parseArgList(is_cond, True)
            self.eatSpace()
            if self.peekChar() == '}':
                self.getChar()
                self.eatSpace()
            elif not self.parse_units:
                error("Missing '}'")
        elif val == ord('`'):
            # undefined macro
            self.lex()  # drop `
            expr = Expression("NAME")
            expr.e1 = "`" + self.getString()

        else:
            if val:
                self.ungetChar(chr(val))
            expr = Expression("NOTHING")
        self.eatSpace()
        return expr

    def parsePostfix(self, is_cond):
        expr = self.parsePrimary(is_cond)
        if expr.type == "NAME":
            ch = self.peekChar()
            if ch == '(':
                self.getChar()
                fname = expr.e1
                self.eatSpace()
                ch = self.peekChar()
                if fname in gAccessFuncs and ch == '<':
                    self.getChar()
                    pname = ""
                    ch = self.peekChar()
                    while ch not in [0, '>', ')']:
                        pname += self.getChar()
                        ch = self.peekChar()
                    if ch == '>':
                        self.getChar()
                        if pname not in gPortnames:
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
                        if gCompactModel:
                            warning("Function '%s' should not be used in a compact model" % fname)
                    elif nargs != len(args):
                        error("Incorrect number of arguments (%d) for function %s (expect %d)"
                              % (len(args), fname, nargs))
                    if fname == "log":
                        warning("Found base-10 log(); use ln() for natural log")
                    check_arg = False
                    allow_zero = False
                    if fname in ["ln", "$ln", "log", "$log10", "sqrt", "$sqrt"] and len(args) == 1:
                        check_arg = True
                    elif fname in ["pow", "$pow"] and len(args) == 2:
                        expon = args[1]
                        if expon.type == "NUMBER" and expon.number < 0:
                            check_arg = True
                        elif expon.type == "NUMBER" and expon.number > 0 and not expon.is_int:
                            check_arg = True
                            allow_zero = True
                    if check_arg:
                        arg = args[0]
                        if fname in ["sqrt", "$sqrt"]:
                            allow_zero = True
                        zeros = getParamZeros(arg, True)
                        if len(zeros) > 0:
                            for pname in zeros:
                                arg_str = ""
                                if pname[0] == '-':
                                    pname = pname[1:]
                                    if not checkParamRangeNeg(pname, allow_zero):
                                        if not checkConditionsRequireNegative(pname, allow_zero, self.ternary_cond):
                                            arg_str = arg.asString()
                                else:
                                    if not checkParamRangePos(pname, allow_zero):
                                        if not checkConditionsRequirePositive(pname, allow_zero, self.ternary_cond):
                                            arg_str = arg.asString()
                                if arg_str != "":
                                    if allow_zero:
                                        bad_val = "negative"
                                    else:
                                        bad_val = "non-positive"
                                    if fname in ["pow", "$pow"]:
                                        warning("Range of parameter '%s' allows %s argument for function '%s(arg,%g)'"
                                                % (pname, bad_val, fname, args[1].number))
                                    else:
                                        if pname == arg_str:
                                            warning("Parameter range allows %s argument for function '%s(%s)'"
                                                    % (bad_val, fname, pname))
                                        else:
                                            warning("Range of parameter '%s' allows %s argument for function '%s()'"
                                                    % (pname, bad_val, fname))
                        elif arg.type == "NUMBER":
                            if arg.number < 0:
                                warning("Negative argument to function '%s(%s)'" % (fname, arg.asString()))
                            elif arg.number == 0 and not allow_zero:
                                warning("Zero argument to function '%s(%s)'" % (fname, arg.asString()))
                        elif checkFuncArgNoConditions(arg):
                            warning("Parameter range(s) allow negative argument to function '%s(%s)'"
                                    % (fname, arg.asString()))
                elif fname in gUserFunctions:
                    if gUserFunctions[fname].type == "integer":
                        expr.is_int = True
                elif fname == "white_noise":
                    if len(args) == 1:
                        error("Missing name for white_noise()")
                    elif len(args) > 2 or len(args) == 0:
                        error("Incorrect number of arguments (%d) for function %s (expect 1 or 2)"
                              % (len(args), fname))
                    if "white_noise" not in gNoiseTypes:
                        gNoiseTypes.append("white_noise")
                elif fname == "flicker_noise":
                    if len(args) == 2:
                        error("Missing name for flicker_noise()")
                    elif len(args) > 3 or len(args) < 2:
                        error("Incorrect number of arguments (%d) for function %s (expect 2 or 3)"
                              % (len(args), fname))
                    if "flicker_noise" not in gNoiseTypes:
                        gNoiseTypes.append("flicker_noise")
                elif fname == "ddt":
                    if len(args) != 1:
                        error("Incorrect number of arguments (%d) for function %s (expect 1)"
                              % (len(args), fname))
                elif fname == "ddx":
                    if len(args) == 2:
                        acc = args[1]
                        if acc.type != "FUNCCALL" or acc.e1 not in gAccessFuncs:
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
                                defv = args[1]
                                if defv.type != "NUMBER" or defv.number != 0:
                                    warning("$simparam(\"gmin\", %s) should use 0 for default value" % defv.asString())
                    else:
                        error("Incorrect number of arguments (%d) for function %s (expect 1 or 2)"
                              % (len(args), fname))

                self.eatSpace()
                if self.peekChar() == ')':
                    self.getChar()
                    self.eatSpace()
                else:
                    error("Missing ')' in function call")
            elif ch == '[':
                # bus select
                self.checkBusIndex(expr.e1)
            else:
                if ch in ['+', '-']:
                    oper = self.peekOper()
                    if oper in ["++", "--"]:
                        # error printed by getOper(); just discard it
                        self.getOper()
                vname = expr.e1
                found = False
                vn = findVariableInScope(vname)
                if vn in gVariables:
                    found = True
                    if gVariables[vn].type == "integer":
                        expr.is_int = True
                if not found and vname in gParameters:
                    if gParameters[vname].type == "integer":
                        expr.is_int = True
        return expr

    def parsePrefix(self, is_cond):
        oper = self.peekOper()
        if oper == "-":               # Unary minus
            oper = self.getOper()
            do_warn = ""
            if self.peekOper() == "-":
                do_warn = "-"
            elif self.peekOper() == "+":
                do_warn = "-+"
            arg = self.parsePrefix(is_cond)
            if arg.type == "NUMBER":
                expr = Expression("NUMBER")
                expr.number = -arg.number
                expr.is_int = arg.is_int
                if isinstance(arg.e1, str) and do_warn == "":
                    expr.e1 = oper + arg.e1
            else:
                expr = Expression(oper)
                expr.e1 = arg
                expr.is_int = arg.is_int
            if do_warn != "":
                warning("Adjacent unary operators: '%s%s'" % (do_warn, arg.asString()))
        elif oper == "+":             # Unary plus
            oper = self.getOper()
            do_warn = ""
            if self.peekOper() == "+":
                do_warn = "++"
            elif self.peekOper() == "-":
                do_warn = "+"
            expr = self.parsePrefix(is_cond)
            if do_warn != "":
                warning("Adjacent unary operators: '%s%s'" % (do_warn, expr.asString()))
        elif oper == "!":             # Logical not
            oper = self.getOper()
            do_warn = ""
            if self.peekOper() == "!":
                do_warn = "!"
            arg = self.parsePrefix(is_cond)
            expr = Expression(oper)
            expr.e1 = arg
            if do_warn != "":
                warning("Adjacent unary operators: '%s%s'" % (do_warn, arg.asString()))
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
        if self.peekOper() == "**" or (self.parse_units and self.peekOper() == "^"):
            oper = self.getOper()
            lhs = expr
            rhs = self.parsePower(is_cond)
            expr = Expression(oper)
            expr.e1 = lhs
            expr.e2 = rhs
            expr.is_int = lhs.is_int and rhs.is_int
        return expr

    def checkExprZero(self, expr):
        zero_factors = []
        if expr.type in ["+", "-"] and expr.e2:
            pname = ""
            if expr.e1.type == "NAME" and expr.e1.e1 in gParameters:
                pname = expr.e1.e1
                other = expr.e2
            elif expr.e2.type == "NAME" and expr.e2.e1 in gParameters:
                pname = expr.e2.e1
                other = expr.e1
            if pname != "":
                if other.type == "NUMBER":
                    if expr.type == "+":
                        excl = -other.number
                    else:
                        excl = other.number
                    if not checkParameterRangeExclude(pname, excl):
                        if not checkConditionsExclude(pname, excl, self.ternary_cond):
                            msg = "expression 1/" + expr.asString()
                            if gVerbose:
                                msg += "\n    range of " + pname + " is "
                                msg += gParameters[pname].getValueRangeStr()
                            zero_factors.append(msg)
                elif other.type == "NAME" and other.e1 in gParameters:
                    pname2 = other.e1
                    if expr.type == "+":
                        if (checkParamRangePos(pname, False) and checkParamRangePos(pname2, True)) or \
                           (checkParamRangePos(pname, True)  and checkParamRangePos(pname2, False)) or \
                           (checkParamRangeNeg(pname, False) and checkParamRangeNeg(pname2, True)) or \
                           (checkParamRangeNeg(pname, True)  and checkParamRangeNeg(pname2, False)):
                            bad = False  # both positive or both negative
                        else:
                            bad = False
                            if not checkParameterRangeExclude(pname, 0) and \
                               not checkParameterRangeExclude(pname2, 0) and \
                               not checkConditionsExclude(pname, 0, self.ternary_cond) and \
                               not checkConditionsExclude(pname2, 0, self.ternary_cond):
                                bad = True
                            else:
                                can_check = True
                                for cond in gConditions:
                                    if cond.type != "NOTHING":
                                        can_check = False
                                if can_check:
                                    par1 = gParameters[pname]
                                    par2 = gParameters[pname2]
                                    if par1.value_range == "nb" or par2.value_range == "nb":
                                        bad = True
                                    elif par1.value_range in ["cc", "co", "oc", "oo"] and \
                                         par2.value_range in ["cc", "co", "oc", "oo"]:
                                        [min1, got_min1] = getRangeValue(par1.min)
                                        [max1, got_max1] = getRangeValue(par1.max)
                                        [min2, got_min2] = getRangeValue(par2.min)
                                        [max2, got_max2] = getRangeValue(par2.max)
                                        if got_max1 and max1 > 0 and got_min2 and min2 < 0:
                                            bad = True
                                        elif got_max2 and max2 > 0 and got_min1 and min1 < 0:
                                            bad = True
                            if bad:
                                msg = "expression 1/" + expr.asString()
                                if gVerbose:
                                    msg += "\n    range of " + pname + " is "
                                    msg += gParameters[pname].getValueRangeStr()
                                    if pname2 != pname:
                                        msg += "\n    range of " + pname2 + " is "
                                        msg += gParameters[pname2].getValueRangeStr()
                                zero_factors.append(msg)
                    else:
                        if (checkParamRangePos(pname, False) and checkParamRangeNeg(pname2, True)) or \
                           (checkParamRangePos(pname, True)  and checkParamRangeNeg(pname2, False)) or \
                           (checkParamRangeNeg(pname, False) and checkParamRangePos(pname2, True)) or \
                           (checkParamRangeNeg(pname, True)  and checkParamRangePos(pname2, False)):
                            bad = False  # (pos)-(neg) or (neg)-(pos)
                        else:
                            msg = "expression 1/" + expr.asString()
                            if gVerbose:
                                msg += "\n    range of " + pname + " is "
                                msg += gParameters[pname].getValueRangeStr()
                                if pname2 != pname:
                                    msg += "\n    range of " + pname2 + " is "
                                    msg += gParameters[pname2].getValueRangeStr()
                            zero_factors.append(msg)
        else:
            if expr.type in ["+", "-"] and not expr.e2:
                zeros = getParamZeros(expr.e1, False)
            else:
                zeros = getParamZeros(expr, False)
            zero_pnames = []
            for pname in zeros:
                if pname not in zero_pnames and not checkParameterRangeExclude(pname, 0):
                    if not checkConditionsExclude(pname, 0, self.ternary_cond):
                        zero_pnames.append(pname)
                        msg = "parameter " + pname
                        if gVerbose:
                            msg += "\n    range of " + pname + " is "
                            msg += gParameters[pname].getValueRangeStr()
                        zero_factors.append(msg)
        for pname in zero_factors:
            warning("Possible division by zero for %s" % pname)

    def parseMultiplicative(self, is_cond):
        expr = self.parsePower(is_cond)
        while self.peekOper() in ["*", "/", "%"]:
            oper = self.getOper()
            lhs = expr
            rhs = self.parsePower(is_cond)
            expr = Expression(oper)
            expr.e1 = lhs
            expr.e2 = rhs
            expr.is_int = lhs.is_int and rhs.is_int
            if oper == "/":
                if expr.is_int and not self.parse_units:
                    warning("Integer divide")
                self.checkExprZero(expr.e2)
        return expr

    def parseAdditive(self, is_cond):
        expr = self.parseMultiplicative(is_cond)
        while self.peekOper() in ["+", "-"]:
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
        if self.peekOper() == '?':
            self.ternary_cond = cond
            self.getOper()
            second = self.parseTernary(is_cond)
            has_colon = False
            if self.peekOper() == ':':
                self.getOper()
                has_colon = True
                self.ternary_cond = Expression("!")
                self.ternary_cond.e1 = cond
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
            self.ternary_cond = Expression("NOTHING")
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
            if name in gParameters:
                var = gParameters[name]
            else:
                vn = findVariableInScope(name)
                if vn in gVariables:
                    var = gVariables[vn]
                else:
                    var = False
            if not port and not var:
                error("Identifier '%s' is not a vector port or array, cannot index" % name)
            sel = self.getExpression()
            if sel.type == "NUMBER":
                index = sel.number
                if port:
                    if port.is_bus:
                        if (index > port.msb and index > port.lsb) or (index < port.msb and index < port.lsb):
                            error("Invalid index %d, port '%s' range is [%d:%d]" % (index, name, port.msb, port.lsb))
                    else:
                        error("Port '%s' is not a bus, cannot index" % name)
                elif var:
                    if len(var.range) == 2:
                        msb = var.range[0]
                        lsb = var.range[1]
                        if (index > msb and index > lsb) or (index < msb and index < lsb):
                            error("Invalid index %d, variable '%s' range is [%d:%d]" % (index, name, msb, lsb))
                    elif var.range == []:
                        error("Identifier '%s' is not an array, cannot index" % name)
                    else:  # pragma: no cover
                        error("Unexpected range of variable '%s%s'" % (name, var.range))
            elif gVerbose:
                notice("Not validating vector index '%s[%s]'" % (name, sel.asString()))
            if self.peekChar() == ']':
                self.getChar()
                self.eatSpace()
            else:
                error("missing ']' after bus index")
        return True

    def getParmRange(self, pname, r_or_ex):
        val = self.lex()
        lend = ""
        if val == ord('['):
            lend = "c"
        elif val == ord('('):
            lend = "o"
        elif val == ord(';'):
            error("Missing range after 'from'")
            return [val, float("-inf"), float("inf"), "nb"]
        else:
            error("Invalid %s in %s" % (formatChar(val), r_or_ex))
        pmin = self.getExpression()
        if pmin.type == "NOTHING":
            error("Missing lower bound in %s" % r_or_ex)
        elif pmin.type == "NAME" and pmin.e1 in gParameters:
            gParameters[pmin.e1].used += 1
        val = self.lex()
        if val == ord(';'):
            error("Incomplete %s" % r_or_ex)
            return [val, pmin, float("inf"), "nb"]
        if val != ord(':'):
            error("Invalid %s for range separator" % formatChar(val))
        pmax = self.getExpression()
        if pmax.type == "NOTHING":
            error("Missing upper bound in %s" % r_or_ex)
        elif pmax.type == "NAME" and pmax.e1 in gParameters:
            gParameters[pmax.e1].used += 1
        val = self.lex()
        rend = ""
        if val == ord(']'):
            rend = "c"
        elif val == ord(')'):
            rend = "o"
        elif val == ord(';'):
            error("Missing ')' or ']' in %s" % r_or_ex)
            return [val, pmin, pmax, "nb"]
        else:
            error("Invalid %s in %s" % (formatChar(val), r_or_ex))
        pvalrng = lend + rend
        if rend == "c" and pmax.type == "NAME" and pmax.e1 == "inf":
            if lend == "c" and pmin.type == "-" and pmin.e1.type == "NAME" and pmin.e1.e1 == "inf":
                error("Invalid %s [-inf:inf] for parameter '%s'; should be (-inf:inf)"
                      % (r_or_ex, pname), "Invalid range [-inf:inf]")
            else:
                if lend == "c":
                    lend = "["
                elif lend == "o":
                    lend = "("
                lend = lend + pmin.asString()
                error("Invalid %s %s:inf] for parameter '%s'; should be %s:inf)"
                      % (r_or_ex, lend, pname, lend), "Invalid range inf]")
        elif lend == "c" and pmin.type == "-" and pmin.e1.type == "NAME" and pmin.e1.e1 == "inf":
            if rend == "c":
                rend = "]"
            elif rend == "o":
                rend = ")"
            rend = pmax.asString() + rend
            error("Invalid %s [-inf:%s for parameter '%s'; should be (-inf:%s"
                  % (r_or_ex, rend, pname, rend), "Invalid range [-inf:inf]")
        val = self.lex()
        return [val, pmin, pmax, pvalrng]

    def getBusRange(self):
        bus_range = []
        msb = 0
        lsb = 0
        valid = True
        if valid:
            val = self.getExpression()
            [msb, valid] = getSimpleValue(val, "msb", "[")
        if valid:
            val = self.lex()
            if val != ord(':'):
                error("Expected ':' in bus range")
                valid = False
        if valid:
            val = self.getExpression()
            [lsb, valid] = getSimpleValue(val, "lsb", ":")
        if valid:
            val = self.lex()
            if val != ord(']'):  # pragma: no cover
                valid = False
        if valid:
            val = self.lex()
            bus_range = [msb, lsb]
        else:
            while val != ord(']'):
                val = self.lex()
                if val == 0:  # pragma: no cover
                    fatal("Parse error looking for ']'")
        return bus_range

# end of class Parser


################################################################################
# functions

def fatal( message ):  # pragma: no cover
    if len(gLineNo) > 0:
        print("FATAL: File %s, line %d: %s" % (gFileName[-1], gLineNo[-1], message))
    else:
        print("FATAL: %s" % message)
    sys.exit(1)


def error( message, type=None ):
    if type is None:
        type = message
    count = 0
    if type in gErrorMsgDict:
        count = gErrorMsgDict[type]
    count += 1
    gErrorMsgDict[type] = count
    if (count <= gMaxNum or gMaxNum == 0) and not gPreProcess:
        if len(gLineNo) > 0:
            print("ERROR in file %s, line %d: %s" % (gFileName[-1], gLineNo[-1], message))
        else:  # pragma: no cover
            print("ERROR: %s" % message)
        if count == gMaxNum:
            print("    Further errors of this type will be suppressed")


def warning( message, type=None ):
    if type is None:
        type = message
    count = 0
    if type in gWarningMsgDict:
        count = gWarningMsgDict[type]
    count += 1
    gWarningMsgDict[type] = count
    if (count <= gMaxNum or gMaxNum == 0) and not gPreProcess:
        if len(gLineNo) > 0:
            print("WARNING in file %s, line %d: %s" % (gFileName[-1], gLineNo[-1], message))
        else:  # pragma: no cover
            print("WARNING: %s" % message)
        if count == gMaxNum:
            print("    Further warnings of this type will be suppressed")


def style( message, type=None ):
    if type is None:
        type = message
    count = 0
    if type in gStyleMsgDict:
        count = gStyleMsgDict[type]
    count += 1
    gStyleMsgDict[type] = count
    if (count <= gMaxNum or gMaxNum == 0) and not gPreProcess:
        if len(gLineNo) > 0:
            print("STYLE in file %s, line %d: %s" % (gFileName[-1], gLineNo[-1], message))
        else:  # pragma: no cover
            print("STYLE: %s" % message)
        if count == gMaxNum:
            print("    Further style comments of this type will be suppressed")


def notice( message, type=None ):
    if type is None:
        type = message
    count = 0
    if type in gNoticeMsgDict:
        count = gNoticeMsgDict[type]
    count += 1
    gNoticeMsgDict[type] = count
    if (count <= gMaxNum or gMaxNum == 0) and not gPreProcess:
        if len(gLineNo) > 0:
            print("NOTICE in file %s, line %d: %s" % (gFileName[-1], gLineNo[-1], message))
        else:  # pragma: no cover
            print("NOTICE: %s" % message)
        if count == gMaxNum:
            print("    Further notices of this type will be suppressed")


def printList( keys, start, maxlen ):
    outstr = start
    for item in sorted(keys):
        if outstr == start:
            outstr += item
        elif len(outstr) + len(item) + 2 > maxlen:
            outstr += ","
            print(outstr)
            outstr = start + item
        else:
            outstr += ", " + item
    print(outstr)


def formatChar( chrnum ):  # pragma: no cover
    if chrnum > 127:
        retstr = "non-ASCII character"
    elif chrnum < 32 or chrnum == 127:
        retstr = "non-printable character"
    elif chrnum in [34, 39, 96]:
        # various quotation marks: " ' `
        retstr = ("character %s" % chr(chrnum))
    else:
        retstr = ("character '%s'" % chr(chrnum))
    return retstr


# get value of constant expression (for bus range)
def getSimpleValue( expr, m_or_l, char ):
    value = 0
    valid = True
    if expr.type == "NOTHING":
        error("Expected constant expression for %s after '%s'" % (m_or_l, char))
        valid = False
    elif expr.type == "NUMBER":
        value = expr.number
    elif expr.type == "NAME":
        name = expr.e1
        if name in gParameters:
            gParameters[name].used += 1
            defv = gParameters[name].defv
            if defv.type == "NUMBER":
                value = defv.number
            else:
                warning("Cannot determine range %s from '%s'" % (m_or_l, name))
        else:
            error("Expected constant expression for %s after '%s', got '%s'" % (m_or_l, char, name))
    elif expr.type in ["+", "-", "*", "/"]:
        if valid:
            [e1, valid] = getSimpleValue(expr.e1, m_or_l, char)
        if valid:
            if expr.e2:
                [e2, valid] = getSimpleValue(expr.e2, m_or_l, char)
            elif expr.type in ["+", "-"]:  # unary +/-
                e2 = e1
                e1 = 0
            else:  # pragma: no cover
                error("Missing second operand for %s" % expr.type)
                valid = False
        if valid:
            if expr.type == "+":
                value = e1 + e2
            elif expr.type == "-":
                value = e1 - e2
            elif expr.type == "*":
                value = e1 * e2
            elif expr.type == "/" and e2 != 0:
                value = e1 / e2
    else:
        warning("Cannot determine range %s from '%s'" % (m_or_l, expr.asString()))
    return [value, valid]


# get a numerical value for range (min/max)
def getRangeValue( expr ):
    value = 0
    valid = False
    if expr.type == "NUMBER":
        value = expr.number
        valid = True
    elif expr.type == "NAME" and expr.e1 == "inf":
        value = 1e300
        valid = True
    elif expr.type == "-" and expr.e1.type == "NAME" and expr.e1.e1 == "inf":
        value = -1e300
        valid = True
    return [value, valid]


# get a list of all parameters that, when zero, make expr zero
# (for checking division by zero)
def getParamZeros( expr, signed ):
    ret = []
    if expr.type == "NAME":
        if expr.e1 in gParameters:
            ret.append(expr.e1)
    elif expr.type == "*":
        ret += getParamZeros(expr.e1, signed)
        ret += getParamZeros(expr.e2, signed)
    elif expr.type == "/":
        ret += getParamZeros(expr.e1, signed)
        # denom checked when parsing this subexpression
    elif expr.type == "-" and not expr.e2:
        fac1 = getParamZeros(expr.e1, signed)
        for fac in fac1:
            if signed:
                ret.append("-" + fac)
            else:
                ret.append(fac)
    elif expr.type in ["+", "-"]:
        fac1 = getParamZeros(expr.e1, signed)
        fac2 = getParamZeros(expr.e2, signed)
        for fac in fac1:
            if fac in fac2:
                ret.append(fac)
    elif expr.type == "FUNCCALL" and expr.e1 in ["sqrt", "$sqrt", "pow", "$pow",
            "sin", "$sin", "asin", "$asin", "sinh", "$sinh", "asinh", "$asinh",
            "tan", "$tan", "atan", "$atan", "tanh", "$tanh", "atanh", "$atanh"]:
        ret += getParamZeros(expr.args[0], signed)
    return ret


# check if parameter range excludes the value 'val'
def checkParameterRangeExclude( pname, val ):
    excludes = False
    parm = gParameters[pname]
    if val in parm.exclude:
        excludes = True
    if parm.value_range == "cc":
        if parm.min.type == "NUMBER" and parm.min.number > val:
            excludes = True
        elif parm.max.type == "NUMBER" and parm.max.number < val:
            excludes = True
    elif parm.value_range == "co":
        if parm.min.type == "NUMBER" and parm.min.number > val:
            excludes = True
        elif parm.max.type == "NUMBER" and parm.max.number <= val:
            excludes = True
    elif parm.value_range == "oc":
        if parm.min.type == "NUMBER" and parm.min.number >= val:
            excludes = True
        elif parm.max.type == "NUMBER" and parm.max.number < val:
            excludes = True
    elif parm.value_range == "oo":
        if parm.min.type == "NUMBER" and parm.min.number >= val:
            excludes = True
        elif parm.max.type == "NUMBER" and parm.max.number <= val:
            excludes = True
    # else "nb"
    if not excludes and val == 0 and parm.value_range in ["cc", "co", "oc", "oo"]:
        if parm.min.type == "NAME" and parm.min.e1 in gParameters:
            if parm.value_range in ["cc", "co"]:
                allow_zero = False
            else:
                allow_zero = True
            if checkParamRangePos(parm.min.e1, allow_zero):
                excludes = True
        if parm.max.type == "NAME" and parm.max.e1 in gParameters:
            if parm.value_range in ["cc", "oc"]:
                allow_zero = False
            else:
                allow_zero = True
            if checkParamRangeNeg(parm.max.e1, allow_zero):
                excludes = True
    return excludes


def checkParamRangePos( pname, allow_zero ):
    parm = gParameters[pname]
    if parm.value_range != "nb":
        if parm.min.type == "NUMBER":
            a = parm.min.number
            if parm.value_range in ["cc", "co"]:
                # from [a:inf), a>0
                if a > 0 or (allow_zero and a == 0):
                    return True
            elif parm.value_range in ["oc", "oo"]:
                # from (a:inf), a>=0
                if a >= 0:
                    return True
        elif parm.min.type == "NAME" and parm.min.e1 in gParameters:
            if parm.value_range in ["oc", "oo"]:
                allow_zero = True
            if checkParamRangePos(parm.min.e1, allow_zero):
                return True
    return False


def checkParamRangeNeg( pname, allow_zero ):
    parm = gParameters[pname]
    if parm.value_range != "nb":
        if parm.max.type == "NUMBER":
            b = parm.max.number
            if parm.value_range in ["cc", "oc"]:
                # from (-inf:b], b<0
                if b < 0 or (allow_zero and b == 0):
                    return True
            elif parm.value_range in ["co", "oo"]:
                # from (-inf:b), b<=0
                if b <= 0:
                    return True
        elif parm.max.type == "NAME" and parm.max.e1 in gParameters:
            if allow_zero or parm.value_range in ["co", "oo"]:
                allow_zero = True
            if checkParamRangeNeg(parm.max.e1, allow_zero):
                return True
    return False


# not comprehensive
def exprIsPositive( expr ):
    ret = False
    if expr.type == "NUMBER":
        if expr.number > 0:
            ret = True
    elif expr.type in ["+", "*", "/"]:
        if exprIsPositive(expr.e1) and exprIsPositive(expr.e2):
            ret = True
    elif expr.type == "NAME":
        if expr.e1 in gParameters and checkParamRangePos(expr.e1, False):
            ret = True
    return ret


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
        do_abs = 0
        if exval == 0:
            # special case for excluding 0, if inequality involves
            # a variable whose factors include the parameter
            if e1.type == "NAME" and e1.e1 not in gParameters:
                vname = findVariableInScope(e1.e1)
                if vname in gVariables:
                    zeros = gVariables[vname].zeros
                    if pname in zeros:
                        pname = vname
            elif e2.type == "NAME" and e2.e1 not in gParameters:
                vname = findVariableInScope(e2.e1)
                if vname in gVariables:
                    zeros = gVariables[vname].zeros
                    if pname in zeros:
                        pname = vname
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
        elif e1.type == "FUNCCALL" and e1.e1 == "abs" and e2.type == "NUMBER" and e2.number > 0:
            do_abs = 1
            arg = e1.args[0]
        elif e2.type == "FUNCCALL" and e2.e1 == "abs" and e1.type == "NUMBER" and e1.number > 0:
            do_abs = 2
            arg = e2.args[0]
        if do_abs:
            if do_abs == 1:
                oplist1 = [">", ">="]
                oplist2 = ["<", "<="]
            else:
                oplist2 = [">", ">="]
                oplist1 = ["<", "<="]
            if arg.type in ["+", "-"] and not arg.e2:
                # abs(unary +/-)
                arg = arg.e1
            if arg.type in ["+", "-"] and arg.e1.type == "NAME" and arg.e1.e1 == pname:
                # abs(pname - exval) > val
                if arg.e2 and arg.e2.type == "NUMBER":
                    if (arg.type == "-" and arg.e2.number == exval) or (arg.type == "+" and arg.e2.number == -exval):
                        if not negated and oper in oplist1:
                            excludes = True
                        elif negated and oper in oplist2:
                            excludes = True
            elif arg.type == "NAME" and arg.e1 == pname:
                if exval == 0:
                    if not negated and oper in oplist1:
                        excludes = True
                    elif negated and oper in oplist2:
                        excludes = True
        if val:
            if val_eq:
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
                if val.type == "NUMBER" and val.number > exval:
                    excludes = True
                elif val.type == "NAME" and val.e1 in gParameters and exval == 0:
                    if checkParamRangePos(val.e1, False):
                        excludes = True
            elif val_lt:
                if val.type == "NUMBER" and val.number < exval:
                    excludes = True
                elif val.type == "NAME" and val.e1 in gParameters and exval == 0:
                    if checkParamRangeNeg(val.e1, False):
                        excludes = True
            if not excludes and val.type != "NUMBER" and exval == 0:
                # pname > val or pname >= val, with val > 0
                if (val_gt or val_ge) and exprIsPositive(val):
                    excludes = True
                #elif val_gt and exprIsPositiveOrZero(val):
                #    excludes = True
                #elif (val_lt or val_le) and exprIsNegative(val):
                #    excludes = True
    elif oper == "&&" and not negated:
        c1 = checkConditionsExcludeOper(e1.type, e1.e1, e1.e2, negated, pname, exval, False)
        c2 = checkConditionsExcludeOper(e2.type, e2.e1, e2.e2, negated, pname, exval, False)
        excludes = c1 or c2
    elif oper == "||" and negated and recurse:  # else of (c1 || c2)
        c1 = checkConditionsExcludeOper(e1.type, e1.e1, e1.e2, negated, pname, exval, False)
        c2 = checkConditionsExcludeOper(e2.type, e2.e1, e2.e2, negated, pname, exval, False)
        excludes = c1 or c2
    elif oper in ["!", "NOT"]:
        excludes = checkConditionsExcludeOper(e1.type, e1.e1, e1.e2, not negated, pname, exval, True)
    return excludes


# check if conditions in effect exclude the value 'exval'
def checkConditionsExclude( pname, exval, extra ):
    excludes = False
    for cond in gConditions:
        if cond.type != "NOTHING":
            excludes |= checkConditionsExcludeOper(cond.type, cond.e1, cond.e2, False, pname, exval, True)
    if extra.type != "NOTHING":
        # ternary condition
        excludes |= checkConditionsExcludeOper(extra.type, extra.e1, extra.e2, False, pname, exval, True)
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
            if (oper == ">" and not negated) or (oper == "<=" and negated):
                val_ge = True
            elif (oper == ">=" and not negated) or (oper == "<" and negated):
                val_gt = True
        elif e2.type == "NAME" and pname == e2.e1:
            # val < PAR, val <= PAR
            val = e1
            if (oper == "<" and not negated) or (oper == ">=" and negated):
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
    elif oper == "||" and negated and recurse:  # else of (c1 || c2)
        c1 = checkConditionsRequirePosOper(e1.type, e1.e1, e1.e2, negated, pname, allow_zero, False)
        c2 = checkConditionsRequirePosOper(e2.type, e2.e1, e2.e2, negated, pname, allow_zero, False)
        requires = c1 or c2
    elif oper in ["!", "NOT"]:
        requires = checkConditionsRequirePosOper(e1.type, e1.e1, e1.e2, not negated, pname, allow_zero, True)
    return requires


# check if conditions in effect require positive parameter value
def checkConditionsRequirePositive( pname, allow_zero, extra ):
    requires = False
    for cond in gConditions:
        if cond.type != "NOTHING":
            requires |= checkConditionsRequirePosOper(cond.type, cond.e1, cond.e2, False, pname, allow_zero, True)
    if extra.type != "NOTHING":
        # ternary condition
        requires |= checkConditionsRequirePosOper(extra.type, extra.e1, extra.e2, False, pname, allow_zero, True)
    return requires


# recursive!
def checkConditionsRequireNegOper( oper, e1, e2, negated, pname, allow_zero, recurse ):
    requires = False
    if oper in [">=", ">", "<=", "<"]:
        val_lt = False
        val_le = False
        val = 0
        if e1.type == "NAME" and pname == e1.e1:
            # PAR < val, PAR <= val
            val = e2
            if (oper == "<" and not negated) or (oper == ">=" and negated):
                val_le = True
            elif (oper == "<=" and not negated) or (oper == ">" and negated):
                val_lt = True
        elif e2.type == "NAME" and pname == e2.e1:
            # val > PAR, val >= PAR
            val = e1
            if (oper == ">" and not negated) or (oper == "<=" and negated):
                val_le = True
            elif (oper == ">=" and not negated) or (oper == "<" and negated):
                val_lt = True
        if val:
            if val_le:
                if val.type == "NUMBER" and val.number <= 0:
                    requires = True
                elif val.type == "NAME" and val.e1 in gParameters:
                    if checkParamRangeNeg(val.e1, True):
                        requires = True
            elif val_lt:
                if val.type == "NUMBER" and (val.number < 0 or (allow_zero and val.number == 0)):
                    requires = True
                elif val.type == "NAME" and val.e1 in gParameters:
                    if checkParamRangeNeg(val.e1, allow_zero):
                        requires = True

    elif oper == "&&" and not negated:
        c1 = checkConditionsRequireNegOper(e1.type, e1.e1, e1.e2, negated, pname, allow_zero, False)
        c2 = checkConditionsRequireNegOper(e2.type, e2.e1, e2.e2, negated, pname, allow_zero, False)
        requires = c1 or c2
    elif oper == "||" and negated and recurse:  # else of (c1 || c2)
        c1 = checkConditionsRequireNegOper(e1.type, e1.e1, e1.e2, negated, pname, allow_zero, False)
        c2 = checkConditionsRequireNegOper(e2.type, e2.e1, e2.e2, negated, pname, allow_zero, False)
        requires = c1 or c2
    elif oper in ["!", "NOT"]:
        requires = checkConditionsRequireNegOper(e1.type, e1.e1, e1.e2, not negated, pname, allow_zero, True)
    return requires


# check if conditions in effect require negative parameter value
def checkConditionsRequireNegative( pname, allow_zero, extra ):
    requires = False
    for cond in gConditions:
        if cond.type != "NOTHING":
            requires |= checkConditionsRequireNegOper(cond.type, cond.e1, cond.e2, False, pname, allow_zero, True)
    if extra.type != "NOTHING":
        # ternary condition
        requires |= checkConditionsRequireNegOper(extra.type, extra.e1, extra.e2, False, pname, allow_zero, True)
    return requires


# check if parameter ranges allow negative argument
# (abort checking if conditions are in play, too hard to check in general)
def checkFuncArgNoConditions( arg ):
    allows = False
    can_check = True
    for cond in gConditions:
        if cond.type != "NOTHING":
            can_check = False
    if can_check:
        if arg.type in ["+", "-", "*", "/"] and arg.e2:
            e1_can_be_neg = False
            e1_must_be_neg = False
            if arg.e1.type == "NAME" and arg.e1.e1 in gParameters:
                if not checkParamRangePos(arg.e1.e1, True):
                    e1_can_be_neg = True
                if checkParamRangeNeg(arg.e1.e1, False):
                    e1_must_be_neg = True
            elif arg.e1.type == "NUMBER":
                if arg.e1.number < 0:
                    e1_can_be_neg = True
                    e1_must_be_neg = True
            elif arg.e1.type == "-" and not arg.e1.e2:
                if arg.e1.e1.type == "NAME" and arg.e1.e1.e1 in gParameters:
                    if not checkParamRangeNeg(arg.e1.e1.e1, True):
                        e1_can_be_neg = True
                    if checkParamRangePos(arg.e1.e1.e1, False):
                        e1_must_be_neg = True
            else:
                can_check = False
            e2_can_be_neg = False
            e2_must_be_neg = False
            e2_can_be_pos = False
            if arg.e2.type == "NAME" and arg.e2.e1 in gParameters:
                if not checkParamRangePos(arg.e2.e1, True):
                    e2_can_be_neg = True
                if checkParamRangeNeg(arg.e2.e1, False):
                    e2_must_be_neg = True
                elif not checkParamRangeNeg(arg.e2.e1, True):
                    e2_can_be_pos = True
            elif arg.e2.type == "NUMBER":
                if arg.e2.number < 0:
                    e2_can_be_neg = True
                    e2_must_be_neg = True
                elif arg.e2.number > 0:
                    e2_can_be_pos = True
            else:
                can_check = False
            if can_check:
                if arg.type == "+" and (e1_can_be_neg or e2_can_be_neg):
                    allows = True
                elif arg.type == "-" and (e1_can_be_neg or e2_can_be_pos):
                    allows = True
                elif arg.type in ["*", "/"]:
                    if (e1_can_be_neg or e2_can_be_neg) and not (e1_must_be_neg and e2_must_be_neg):
                        allows = True
    return allows


def checkIdentifierCollisions( basename, name, escaped, type ):
    valid = True
    if basename in gVAMSkeywords and not escaped:
        error("%s '%s' collides with Verilog-AMS keyword" % (type, name))
        valid = False
    elif name == gModuleName:
        error("%s '%s' collides with module name" % (type, name))
        valid = False
    elif name in gNatures and gNatures[name].defined:
        if type == "Nature":
            error("Duplicate declaration of nature '%s'" % name)
        else:
            error("%s '%s' collides with previously-declared nature" % (type, name))
        valid = False
    elif name in gDisciplines:
        if type == "Discipline":
            error("Duplicate declaration of discipline '%s'" % name)
        else:
            error("%s '%s' collides with previously-declared discipline" % (type, name))
        valid = False
    elif name in gParameters:
        if type == "Parameter":
            error("Duplicate declaration of parameter '%s'" % name)
        else:
            error("%s '%s' collides with previously-declared parameter" % (type, name))
        valid = False
    elif name in gVariables:
        if type == "Variable":
            error("Duplicate declaration of variable '%s'" % name)
        else:
            error("%s '%s' collides with previously-declared variable" % (type, name))
        valid = False
    elif name in gNodenames:
        if type == "Node name":
            error("Duplicate declaration of internal node '%s'" % name)
        else:
            error("%s '%s' collides with previously-declared node" % (type, name))
        valid = False
    elif name in gBranches:
        if type == "Branch":
            error("Duplicate declaration of branch '%s'" % name)
        else:
            error("%s '%s' collides with previously-declared branch" % (type, name))
        valid = False
    elif name in gUserFunctions:
        if type == "Function":
            error("Duplicate user-defined function '%s'" % name)
        else:
            error("%s '%s' collides with previously-defined function" % (type, name))
        valid = False
    elif name in gBlocknames:
        if type == "Block name":
            error("Duplicate named block '%s'" % name)
        else:  # pragma: no cover
            error("%s '%s' collides with previously-found named block" % (type, name))
        valid = False
    elif name in gPortnames and type != "Port name":
        error("%s '%s' collides with previously-declared port (terminal)" % (type, name))
        valid = False
    elif type == "Port name" and name not in gPortnames:
        if len(gScopeList) == 1:
            error("Identifier '%s' is not a port (terminal), cannot set direction" % name)
        valid = False
    return valid


def getAttributes( line, in_attrib ):
    retval = []
    do_append = True
    had_comma = False
    got_attrib = False
    if in_attrib == "":
        start  = line.find("(*")
        start2 = start + 2
    else:
        start  = 0
        start2 = 0
    stop = line.find("*)")
    if stop == 0:
        line = line[stop+2:]
        line = line.strip()
        if len(line) > 0 and line[0] == ';':
            warning("Attribute should prefix the item to which it applies")
        start = line.find("(*")
        start2 = start + 2
        stop = line.find("*)")
        in_attrib = ""
    if start >= 0 and stop < 0:
        # attribute continues on next line
        retval = [" "]
    else:
        retval = [""]
    while start >= 0:
        semi = line.find(";")
        if semi >= 0 and semi < start:
            do_append = False
            error("Attribute found after ';', please add line break")
        stop = line.find("*)")
        if stop > 0:
            quot = line.find("\"")
            if quot >= 0 and quot < stop:
                stop = -1
                in_quote = False
                for i in range(len(line)-1):
                    ch = line[i]
                    if ch == '"':
                        in_quote = not in_quote
                    elif not in_quote and ch == '*' and line[i+1] == ')':
                        stop = i
                        break
        if stop >= 0:
            i = stop+2
            while i < len(line) and line[i].isspace():
                i += 1
            if i < len(line) and line[i] == ';':
                warning("Attribute should prefix the item to which it applies")
        if stop > start or stop == -1:
            if stop == -1:
                attrib = line[start2:]
                line = line[:start]
            else:
                attrib = line[start2:stop]
                line = line[:start] + line[stop+2:]
            line = line.strip()
            attrib = attrib.strip()
            if attrib != "":
                if in_attrib == " ":
                    if attrib[0] != ',':
                        error("Missing ',' in attribute list")
                elif attrib[0] == ',':
                    if in_attrib in [',', ""]:
                        error("Extra comma before attributes")
                        attrib = attrib[1:]
                nested = attrib.find("(*")
                if nested > 0:
                    quot = attrib.find("\"")
                    if quot >= 0 and quot < nested:
                        nested = -1
                        in_quote = False
                        for i in range(len(attrib)-1):
                            ch = attrib[i]
                            if ch == '"':
                                in_quote = not in_quote
                            elif not in_quote and ch == '(' and attrib[i+1] == '*':  # pragma: no cover
                                nested = i
                                break
                if nested >= 0:
                    countParentheses(attrib, 0, 1)
                    attrib = attrib[:nested]
                    if line != "":
                        nested = line.find("*)")
                        if nested >= 0:
                            line = line[nested+2:]
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
                if name == "":
                    error("Missing name for attribute")
                if do_append:
                    retval.append(name)
                while i < len(attrib) and attrib[i].isspace():
                    i += 1
                got_attrib = True
                get_value = False
                value = ""
                if i < len(attrib):
                    if attrib[i] == '=':
                        i += 1
                        get_value = True
                    elif attrib[i] == '"':
                        error("Missing '=' before value for attribute '%s'" % name)
                        get_value = True
                if get_value:
                    # (* name = value *)
                    parser = Parser(attrib[i:])
                    expr = parser.getExpression()
                    if expr.type != "NOTHING":
                        value = expr.asString()
                    i = len(attrib)-len(parser.getRestOfLine())
                    if value == "":
                        error("Missing value for attribute '%s'" % name)
                # else:
                #     (* name *)
                if do_append:
                    retval.append(value)
                had_comma = False
                if i < len(attrib):
                    if attrib[i] == ',':
                        had_comma = True
                        i += 1
                        while i < len(attrib) and attrib[i].isspace():
                            i += 1
                            if attrib[i] == ',':
                                error("Extra comma in attributes")
                                i += 1
                    elif attrib[i] == ';' and name == "inherited_mfactor":
                        error("Archaic attribute syntax: %s" % attrib)
                        i += 1
                    elif attrib[i].isalnum() or attrib[i] == '_':
                        if value == "":
                            error("Missing '=' or ',' in attribute")
                        else:
                            error("Missing ',' in attribute list")
                    else:
                        error("Unexpected character '%s' in attribute" % attrib[i])
                        i += 1
        else:  # pragma: no cover
            fatal("Malformed attribute")
        start = line.find("(*")
        start2 = start + 2
    if had_comma:
        if retval[0] == " ":
            retval[0] = ","
        else:
            error("Extra comma at end of attributes")
    elif retval[0] == " " and not got_attrib:
        # (* by itself on a line, pretend we got a comma
        retval[0] = ","
    retval.insert(0, line)
    return retval
# end of getAttributes


def parseMacro( line ):
    if line.startswith("`define"):
        i = len("`define")
        max_i = len(line)
        if i < max_i:
            if not line[i].isspace():
                error("Missing space after `define")
            while i < max_i and line[i].isspace():
                i += 1
            name = ""
            if i < max_i and not (line[i].isalpha() or line[i] == '_'):
                error("Macro name must start with a letter or underscore (_)")
            while i < max_i and (line[i].isalnum() or line[i] == '_'):
                name += line[i]
                i += 1
            if len(name) > 7 and name[0:7] == "__VAMS_":
                error("Macro name may not start with __VAMS_")
            args = []
            if i < max_i:
                if line[i] == '(':
                    while i < max_i and line[i] != ')':
                        i += 1
                        while i < max_i and line[i].isspace():
                            i += 1
                        arg = ""
                        if i < max_i and not (line[i].isalpha() or line[i] == '_'):
                            error("Macro formal arguments must start with a letter or underscore (_)")
                        while i < max_i and (line[i].isalnum() or line[i] == '_'):
                            arg += line[i]
                            i += 1
                        args.append(arg)
                        while i < max_i and line[i].isspace():
                            i += 1
                        if i < max_i and line[i] != ',' and line[i] != ')':
                            error("Macro formal arguments should be separated by a comma (,)")
                    if i < max_i and line[i] == ')':
                        i += 1
                    else:
                        error("Missing ')' for macro definition")
            if i < max_i:
                value = " " + line[i:].strip()
            else:
                value = ""
            if name in gVAMScompdir:
                error("Invalid macro name '%s'" % name)
            else:
                if name in gMacros:
                    warning("Redefining macro '%s'" % name)
                mac = Macro(name, value)
                mac.args = args
                gMacros[name] = mac
                filename = gFileName[-1]
                if filename.find("disciplines.vams") >= 0 or filename.find("discipline.h") >= 0 \
                        or filename.find("disciplines.h") >= 0 \
                        or filename.find("constants.vams") >= 0 or filename.find("constants.h") >= 0:
                    # don't complain about unused macros from standard header files
                    mac.used = True
        else:
            error("Missing macro name after `define")
# end of parseMacro


def parseUndef( line ):
    if line.startswith("`undef"):
        token = "`undef"
        if line.startswith("`undefine"):
            token = "`undefine"
            error("Invalid compiler directive `undefine (use `undef)")
        i = len(token)
        max_i = len(line)
        if i < max_i:
            if not line[i].isspace():
                error("Missing space after %s" % token)
            while i < max_i and line[i].isspace():
                i += 1
            name = ""
            while i < max_i and (line[i].isalnum() or line[i] == '_'):
                name += line[i]
                i += 1
            if i < max_i:
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
    # check if actual argument contains match for later formal argument
    if check_collision:
        collide = False
        for i in range(len(actuals)):
            act = actuals[i]
            for j in range(i+1, len(args)):
                arg = args[j]
                if act.find(arg) >= 0:
                    collide = True
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
                                string = text[0:j]
                                string += actuals[i]
                                string += text[j+len(arg):]
                                text = string
                                j += len(actuals[i]) - 1
                else:
                    word_bnd = True
            j += 1
    return text


def findFirstUnquotedTic( line, start ):
    pos = -1
    i = start
    max_i = len(line)
    in_quote = False
    escaped = False
    while i < max_i:
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
            gMacros[name].used = True
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
                else:  # pragma: no cover
                    fatal("Macro `%s called with %d arguments, expected %d" % (name, len(actuals), len(args)))
            elif len(args) > 0:
                error("Macro `%s called without arguments, expected %d" % (name, len(args)))
            if len(text) > 0 and text[-1] == ';' and len(line) > i and line[i] == ';':
                warning("Expansion of macro '%s' results in repeated ';'" % name)
            string = line[0:pos]
            string += text
            string += line[i:]
            line = string
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
            if t2 < 0:  # pragma: no cover
                fatal("Found `%s in macro text without corresponding `endif" % name)
            else:
                i = t2 + 6
            rest = line[te:t2]
            if rest.find("`ifdef") > 0 or rest.find("`ifndef") > 0:  # pragma: no cover
                fatal("Found nested `ifdef in macro text, not supported")
            elif rest.find("`elsif") > 0:  # pragma: no cover
                fatal("Found `elsif in macro text, not supported")
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
            string = line[0:pos]
            string += ifdef_text
            string += line[i:]
            line = string
        elif name in ["else", "elsif", "endif"]:  # pragma: no cover
            fatal("Found `%s in macro text without corresponding `ifdef" % name)
        elif name in gVAMScompdir:
            error("Can't handle compiler directive `%s in macro text" % name)
            start = i
        elif name == "elseif":
            # error printed in handleCompilerDirectives
            start = i
        else:
            if gMissingConstantsFile != "" and name[0:2] in ["M_", "P_"]:
                error("Undefined macro `%s; please check %s" % (name, gMissingConstantsFile))
            else:
                error("Undefined macro `%s" % name)
            start = i
        pos = findFirstUnquotedTic(line, start)
    return line
# end of expandMacro


def verifyEndScope( what ):
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
    else:  # pragma: no cover
        scope = getCurrentScope()
        error("Unexpected '%s' while still in scope %s" % (endname, scope))


def parseNatureDecl( line, defined ):
    parser = Parser(line)
    parser.lex()
    if not parser.isIdentifier() or parser.getString() != "nature":  # pragma: no cover
        fatal("Parse error looking for 'nature'")

    name = ""
    escaped = False
    parser.lex()
    if parser.isIdentifier():
        name = parser.getString()
        escaped = parser.isEscaped()
    else:  # pragma: no cover
        fatal("Parse error looking for nature name")

    if name != "":
        if defined:
            valid = checkIdentifierCollisions(name, name, escaped, "Nature")
            if valid:
                nature = Nature(name, defined)
                gNatures[name] = nature
            gScopeList.append("NATURE::" + name)
        elif name not in gNatures:
            nature = Nature(name, defined)
            gNatures[name] = nature


def parseNatureLine( line ):
    name = (gScopeList[0])[len("NATURE::"):]
    parser = Parser(line)
    parser.lex()
    attrib = ""
    if parser.isIdentifier():
        attrib = parser.getString()
    else:  # pragma: no cover
        error("Parse error in nature '%s'" % name)

    if attrib == "end":
        error("Found 'end' instead of 'endnature'")
        attrib = "endnature"

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
                    if idt not in gNatures:
                        error("In nature '%s', idt_nature references an unknown nature '%s'" % (name, idt))
                else:
                    error("Expected identifier for %s in nature '%s'" % (attrib, name))

        val = parser.lex()
        if name != "" and val != ord(';'):
            error("Missing ';' at end of item in nature '%s'" % name)
# end of parseNatureLine


def parseDisciplineDecl( line ):
    parser = Parser(line)
    parser.lex()
    if not parser.isIdentifier() or parser.getString() != "discipline":  # pragma: no cover
        fatal("Parse error looking for 'discipline'")

    name = ""
    escaped = False
    parser.lex()
    if parser.isIdentifier():
        name = parser.getString()
        escaped = parser.isEscaped()
    else:  # pragma: no cover
        fatal("Parse error looking for discipline name")

    valid = checkIdentifierCollisions(name, name, escaped, "Discipline")
    if valid:
        disc = Discipline(name)
        gDisciplines[name] = disc
    gScopeList.append("DISCIPLINE::"+name)


def parseDisciplineLine( line ):
    disc = (gScopeList[0])[len("DISCIPLINE::"):]
    parser = Parser(line)
    val = parser.lex()
    attrib = ""
    if parser.isIdentifier():
        attrib = parser.getString()
    else:
        error("Parse error in discipline '%s'" % disc)

    if attrib == "end":
        error("Found 'end' instead of 'enddiscipline'")
        attrib = "enddiscipline"

    if attrib == "enddiscipline":
        if gVerbose:
            if disc in gDisciplines and gDisciplines[disc].domain == "":
                if gDisciplines[disc].potential == "":
                    warning("No potential nature specified for discipline '%s'" % disc)
                if gDisciplines[disc].flow == "":
                    warning("No flow nature specified for discipline '%s'" % disc)
        verifyEndScope("discipline")
    else:
        if attrib in ["potential", "flow"]:
            val = parser.lex()
            nature = ""
            if parser.isIdentifier() or parser.isString():
                nature = parser.getString()
                if nature not in gNatures:
                    error("Undefined nature '%s' for discipline '%s'" % (nature, disc))
            elif val == ord(';'):
                error("Missing nature for '%s' in discipline '%s'" % (attrib, disc))
            else:  # pragma: no cover
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
                        if units != "" and units not in gUnitsMultiply:
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


# reset all macros
def resetMacros():
    global gMacros
    gMacros      = {}
    # predefined macros
    mac1 = Macro("__VAMS_ENABLE__", "")
    mac1.used = True
    gMacros[mac1.name] = mac1
    mac2 = Macro("__VAMS_COMPACT_MODELING__", "")
    mac2.used = True
    gMacros[mac2.name] = mac2

    # -D command-line defines
    for tok in gDefines:
        parts = tok.split("=")
        if len(parts) > 1:
            mac = Macro(parts[0], parts[1])
        else:
            mac = Macro(parts[0], "")
        mac.used = True
        gMacros[mac.name] = mac


# clear all module-scope items for new module
def initializeModule():
    global gParameters, gVariables, gHiddenState, gPortnames, gPortlist, gNodenames, gBranches, \
           gMultScaled, gUsesDDT, gContribs, gBlocknames, gBinningPatterns, gBinningFactors, gNoiseTypes
    gParameters      = {}
    gVariables       = {}
    gHiddenState     = {}
    gPortnames       = {}
    gPortlist        = ""
    gNodenames       = {}
    gBranches        = {}
    gMultScaled      = {}
    gUsesDDT         = False
    gContribs        = {}
    gBlocknames      = {}
    gBinningPatterns = [[], []]
    gBinningFactors  = {}
    gNoiseTypes      = []
    resetMacros()

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


# print parameter defaults
def printParamDefaults( plist ):
    for p in plist:
        defv = gParameters[p].defv
        if defv.type == "NUMBER":
            defv = defv.asString()
        elif defv.type == "NAME":
            pname = defv.e1
            if pname in gParameters:
                defv = gParameters[pname].defv.asString()
            else:  # pragma: no cover
                # should have been caught during parsing
                defv = "UNKNOWN"
        else:
            defv = defv.asString()
        print("    %s = %s" % (p, defv))


# print summary
def printModuleSummary():
    if gModuleName:
        print("\nModule: %s" % gModuleName)
    else:  # pragma: no cover
        print("\nNo module declaration found")
    if gVerbose:
        print("\nPorts: %s" % gPortlist)
        if gNodenames:
            if len(gNodenames) < 10:
                printList(gNodenames.keys(), "Internal nodes: ", 70)
            else:
                print("Internal nodes:")
                printList(gNodenames.keys(), "    ", 70)
        if len(gMultScaled):
            start = "Ports/nodes that get mult scaling: "
            printList(gMultScaled.keys(), start, 70)
    if len(gParameters):
        inst_parms = {}
        model_parms = {}
        lin_scale = []
        quad_scale = []
        for (k, p) in gParameters.items():
            if not p.is_alias:
                if p.inst:
                    inst_parms[k] = p
                if p.model:  # can be both model and inst
                    model_parms[k] = p
            if p.scale == "\"linear\"":
                lin_scale.append(p.name)
            elif p.scale == "\"quadratic\"":
                quad_scale.append(p.name)
        if len(inst_parms):
            if gVerbose or len(inst_parms) < 20:
                print("\nInstance parameters (%d):" % len(inst_parms))
                printList(inst_parms.keys(), "    ", 70)
            else:
                print("\nInstance parameters (%d)" % len(inst_parms))
            if len(lin_scale) > 0:
                print("  Linear scaling:")
                printList(lin_scale, "    ", 70)
            if len(quad_scale) > 0:
                print("  Quadratic scaling:")
                printList(quad_scale, "    ", 70)
            if gPrDefVals:
                printParamDefaults(inst_parms)
            if len(model_parms) == 0:
                if gVerbose:
                    print("\nNo model parameters")
            elif gVerbose or len(model_parms) < 20:
                print("\nModel parameters (%d):" % len(model_parms))
                printList(model_parms.keys(), "    ", 70)
            else:
                print("\nModel parameters (%d)" % len(model_parms))
            if gPrDefVals:
                printParamDefaults(model_parms)
        else:
            if gVerbose or len(gParameters) < 40:
                print("\nParameters (%d):" % len(gParameters))
                printList(gParameters.keys(), "    ", 70)
            else:
                print("\nParameters (%d)" % len(gParameters))
            if gPrDefVals:
                print("Default values:")
                printParamDefaults(gParameters)

    if len(gVariables):
        op_pt_vars = {}
        vars_to_pop = []
        for (k, v) in gVariables.items():
            if v.assign == 0 or v.assign < -1:
                vars_to_pop.append(v.name)
            elif v.oppt:
                op_pt_vars[k] = v
        for var in vars_to_pop:
            gVariables.pop(var)
        if len(op_pt_vars):
            if gVerbose or len(op_pt_vars) < 50:
                print("\nOperating-point variables (%d):" % len(op_pt_vars))
                printList(op_pt_vars.keys(), "    ", 70)
            else:
                print("\nOperating-point variables (%d)" % len(op_pt_vars))
        elif len(gVariables):
            print("\nNo operating-point variables")
        if gVerbose:
            print("\nVariables (%d):" % len(gVariables))
            printList(gVariables.keys(), "    ", 70)

    if gCompactModel and len(gHiddenState):
        print("Possible hidden-state variables:")
        printList(gHiddenState.keys(), "    ", 70)
# end of printModuleSummary


def parseModuleDecl( line, attribs ):
    global gModuleName, gPortlist

    parser = Parser(line)
    val = parser.lex()
    mtype = ""
    if parser.isIdentifier():
        mtype = parser.getString()
    if mtype not in ["module",  "macromodule"]:  # pragma: no cover
        fatal("Parse error looking for 'module'")

    mname = ""
    val = parser.lex()
    if parser.isIdentifier():
        mname = parser.getString()
    else:  # pragma: no cover
        fatal("Parse error looking for module name")

    if mname in gVAMSkeywords:
        error("Module name '%s' collides with Verilog-AMS keyword" % mname)
    else:
        if gModuleName != "":
            if gModuleName == mname:
                warning("File contains duplicate definition of module '%s'" % gModuleName)
            else:
                warning("File contains multiple modules: %s and %s" % (gModuleName, mname))
                printModuleSummary()
                initializeModule()
                print("")
        gModuleName = mname

    if attribs:
        for i in range(0, len(attribs)-1, 2):
            if attribs[i] == "instance_parameter_list":
                warning("Module declared with non-standard attribute (* %s *)" % attribs[i])
            elif gVerbose and attribs[i] in ["compact_module", "compact_model"]:
                notice("Module declared as (* %s *)" % attribs[i])

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
                if gPortlist != "":
                    gPortlist += ", "
                gPortlist += pname
            val = parser.lex()
            if val == ord(','):
                val = parser.lex()
            elif val != ord(')') and val != ord(';') and not parser.isIdentifier():
                error("Unexpected %s in port list" % formatChar(val))
                val = parser.lex()
        if val == ord(')'):
            val = parser.lex()
        else:
            error("Missing ')' after list of ports in module declaration")
    if val != ord(';'):
        error("Missing ';' at end of module declaration")
    elif parser.peekChar() == ';':
        error("Extra ';' at end of module declaration")
# end of parseModuleDecl


def parseFunction( line ):
    global gCurrentFunc
    analog = False

    parser = Parser(line)
    val = parser.lex()
    if parser.isIdentifier() and parser.getString() == "analog":
        analog = True
        val = parser.lex()

    if not parser.isIdentifier() or parser.getString() != "function":  # pragma: no cover
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

    if ftype not in ["real", "integer", "string"]:
        error("Invalid type '%s' for function '%s'" % (ftype, fname))

    gScopeList.append("FUNCTION::" + fname)
    gCurrentFunc = Function(fname, ftype)
    valid = checkIdentifierCollisions(fname, fname, escaped, "Function")
    if valid:
        gUserFunctions[fname] = gCurrentFunc
        # implicit declaration of return variable
        scope = getCurrentScope()
        vname = scope + fname
        retval = Variable(vname, ftype, False, gFileName[-1], gLineNo[-1])
        gVariables[vname] = retval
# end of parseFunction


def parseParamDecl( line, attribs ):
    global gCheckMultFactors
    in_module = True

    parm_or_alias = ""
    parser = Parser(line)
    val = parser.lex()
    if parser.isIdentifier():
        parm_or_alias = parser.getString()
    if parm_or_alias not in ["parameter", "aliasparam", "localparam"]:  # pragma: no cover
        fatal("Parse error looking for 'parameter'")

    inst = False
    model = True
    units = 0
    scale = ""
    if attribs:
        for i in range(0, len(attribs)-1, 2):
            if attribs[i] == "type":
                val = attribs[i+1].lower()
                if val in ["instance", "\"instance\""]:
                    inst = True
                    model = False
                elif val in ["both", "\"both\""]:
                    inst = True
            elif attribs[i] == "units":
                units = attribs[i+1].strip('"').strip()
            elif attribs[i] == "scale":
                val = attribs[i+1]
                if val in ["\"linear\"", "\"quadratic\""]:
                    if scale == "":
                        scale = val
                    else:
                        error("Scaling specified twice: %s and %s" % (val, scale))
                elif val in ["linear", "quadratic"]:
                    error("Attribute value should be in quotes")
                else:
                    error("Invalid attribute: %s=%s" % (attribs[i], val))

    ptype = ""
    if parm_or_alias == "aliasparam":
        ptype = "alias"
    else:
        val = parser.lex()
        if parser.isIdentifier():
            ptype = parser.getString()
        else:  # pragma: no cover
            fatal("Parse error looking for parameter type")

    val = ord(',')
    while val == ord(','):
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
        parm.model = model
        parm.units = units
        if parm.inst:
            parm.scale = scale
        elif scale != "":
            error("May not specify %s-scaling for model parameter '%s'" % (scale, pname))

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
            if ptype not in ["real", "integer", "string"]:
                error("Invalid type '%s' for parameter '%s'" % (ptype, pname))

        if pname != "":
            valid = checkIdentifierCollisions(pname, pname, escaped, "Parameter")
        else:
            valid = False
        # insert into global map
        if valid and in_module:
            gParameters[pname] = parm
            if pname in gMultFactors:
                gCheckMultFactors = 1

        if val == ord('['):
            parm.range = parser.getBusRange()
            val = parser.token

        if val != ord('='):
            error("Incomplete %s specification" % parm_or_alias)
            return

        parm.defv = parser.getExpression()
        if parm.defv.type == "ARRAY" or len(parm.range) > 0:
            if parm.defv.type != "ARRAY":
                error("Array parameter '%s' needs array of default values" % pname)
            elif len(parm.range) == 0:
                error("Non-array parameter '%s' has array of default values" % pname)
            else:
                nump = parm.range[1] - parm.range[0]
                if nump > 0:
                    nump += 1
                else:
                    nump = 1 - nump
                if nump != len(parm.defv.args):
                    error("Array parameter '%s' has different dimension than array of default values" % pname)

        if parm_or_alias == "aliasparam":
            if parm.defv.type == "NAME" and parm.defv.e1 in gParameters:
                parm.type = gParameters[parm.defv.e1].type
                gParameters[parm.defv.e1].used += 1
            else:
                error("Default value for aliasparam '%s' must be a parameter" % pname)
            val = parser.lex()
            parm.used += 1
            parm.is_alias = True

        else:
            if parm.defv.type == "NAME" and parm.defv.e1 in gParameters:
                dparm = gParameters[parm.defv.e1]
                if dparm.units and units and dparm.units != units:
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
                if (dparm.units or units) and dparm.units != units:
                    if units:
                        u1 = "units '" + units + "'"
                    else:
                        u1 = "no units"
                    if dparm.units:
                        u2 = "units '" + dparm.units + "'"
                    else:
                        u2 = "no units"
                    notice("Parameter '%s' with %s has default '%s' with %s" % (pname, u1, dparm.name, u2))
            deps = parm.defv.getDependencies(False, False)
            for dep in deps:
                if dep == pname:
                    error("Parameter '%s' default value may not depend on itself" % pname)
                elif dep in gParameters:
                    gParameters[dep].used += 1
                    if parm.model and not gParameters[dep].model:
                        if parm.defv.type == "NAME" and parm.defv.e1 == dep:
                            error("Default value of model parameter '%s' is an instance parameter '%s'"
                                  % (pname, dep))
                        else:
                            error("Default value of model parameter '%s' depends on instance parameter '%s'"
                                  % (pname, dep))
                else:
                    if dep in gVariables:
                        error("Parameter '%s' default value may not depend on variable '%s'" % (pname, dep))
                    elif dep in gPortnames:
                        error("Parameter '%s' default value may not depend on port '%s'" % (pname, dep))
                    elif dep in gNodenames:
                        error("Parameter '%s' default value may not depend on node '%s'" % (pname, dep))
                    elif dep in gBranches:
                        error("Parameter '%s' default value may not depend on branch '%s'" % (pname, dep))
                    elif dep == "inf":
                        error("Infinity (inf) may only be used in ranges; not a valid default value for parameter '%s'"
                              % pname)
                    else:
                        error("Parameter '%s' default value depends on unknown identifier '%s'" % (pname, dep))

            val = parser.lex()
            parse_error = False
            while parser.isIdentifier():
                string = parser.getString()
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
                        if parser.isNumber():
                            parse_error = True
                            error("Invalid number in range for string parameter")
                            break
                    else:
                        parse_error = True
                        error("Invalid %s in range" % formatChar(val))
                        break
                    if val == ord('}'):
                        val = parser.lex()
                    else:
                        parse_error = True
                        error("Invalid %s in range" % formatChar(val))
                        break
                elif string == "from":
                    [val, parm.min, parm.max, parm.value_range] = parser.getParmRange(parm.name, "range")
                elif string == "exclude":
                    parser.eatSpace()
                    val = parser.peekChar()
                    if val in ['[', '(']:
                        # exclude range, parsed but ignored
                        val = parser.getParmRange(parm.name, "exclude")[0]
                    else:
                        exclude = parser.getExpression()
                        if exclude.type == "NUMBER":
                            parm.exclude.append(exclude.number)
                        elif exclude.type == "NAME":
                            if exclude.e1 in gParameters:
                                gParameters[exclude.e1].used += 1
                            else:
                                warning("Parameter range 'exclude %s' is not valid" % exclude.e1)
                        elif exclude.type == "NOTHING":
                            error("Missing value for 'exclude'")
                        else:
                            warning("Unsupported exclusion of type %s" % exclude.type)
                        val = parser.lex()
                else:  # pragma: no cover
                    fatal("Unexpected '%s' in parameter declaration" % string)

        if parm.defv.type == "NUMBER":
            num = parm.defv.number
            if num in parm.exclude:
                warning("Default value of %s=%d is excluded from range" % (parm.name, num))
            elif parm.value_range != "nb":
                if parm.min.type == "NUMBER" and (parm.min.number > num \
                        or (parm.value_range in ["oc", "oo"] and parm.min.number == num)):
                    warning("Default value of %s=%d is too small" % (parm.name, num))
                elif parm.max.type == "NUMBER" and (parm.max.number < num \
                        or (parm.value_range in ["co", "oo"] and parm.max.number == num)):
                    warning("Default value of %s=%d is too large" % (parm.name, num))

    if val != ord(';'):
        if not parse_error:
            error("Missing ';' at end of %s declaration" % parm_or_alias)
    elif parser.peekChar() == ';':
        error("Extra ';' at end of %s declaration" % parm_or_alias)
# end of parseParamDecl


def parseVariableDecl( line, attribs ):
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
                if vtype not in ["real", "integer", "genvar", "string"]:
                    err_msg = ("Unexpected variable type '%s'" % vtype)
                val = parser.lex()

        if parser.isIdentifier():
            string = parser.getString()
            escaped = parser.isEscaped()
            scope = getCurrentScope()
            vname = scope + string
            valid = checkIdentifierCollisions(string, vname, escaped, "Variable")
            # insert into global map
            if valid and in_module:
                var = Variable(vname, vtype, oppt, gFileName[-1], gLineNo[-1])
                if oppt:
                    var.used = True
                    var.units = units
                    if units in gUnitsMultiply:
                        if multi != "\"multiply\"":
                            warning("Operating-point variable %s should have multiplicity=\"multiply\""
                                    % var.name, "multiplicity")
                    elif units in gUnitsDivide:
                        if multi != "\"divide\"":
                            warning("Operating-point variable %s should have multiplicity=\"divide\""
                                    % var.name, "multiplicity")
                    elif multi != "":
                        if units == "":
                            unitstr = "no units"
                        else:
                            unitstr = "units '" + units + "'"
                        warning("Unexpected multiplicity=%s for operating-point variable '%s' with %s"
                                % (multi, var.name, unitstr), "multiplicity")
                gVariables[vname] = var
                if len(gScopeList) > 0 and gScopeList[-1].startswith("FUNCTION::"):
                    if gCurrentFunc:
                        if string in gCurrentFunc.inputs:
                            # hack, not really where it was assigned a value
                            var.assign = gLineNo[-1]
                        elif string in gCurrentFunc.outputs:
                            # marks the variable within the function as used
                            # (separate question whether the output is used outside the function)
                            var.used = True
            val = parser.lex()
            if val == ord('['):
                gVariables[vname].range = parser.getBusRange()
                val = parser.token
            if val == ord('='):
                # real var = 1;
                error("Variable initializers not supported")
                parser.getExpression()
                val = parser.lex()
        else:
            if vtype == "":  # pragma: no cover
                fatal("Parse error looking for variable type")
            elif vtype in ["real", "integer", "genvar", "string"]:
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
            while val == ord(';'):
                error("Extra ';' at end of variable declaration")
                val = parser.lex()
        elif parser.isIdentifier():
            string = parser.getString()
            if string in ["inout", "input", "output", "real", "integer", "genvar", "string"]:
                error("Missing ';' after list of variables")
                vtype = string
                val = parser.lex()
            else:
                error("Missing ',' in list of variables")
        elif val != 0:
            error("Invalid %s in variable declaration" % formatChar(val))
            val = parser.lex()
    if line[-1] != ';':
        error("Missing ';' at end of variable declaration")
# end of parseVariableDecl


def parsePortDirection( line ):
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
                if ptype in ["input", "output"]:
                    if len(gScopeList) == 1 and gScopeList[-1] == "MODULE::":
                        if gCompactModel:
                            warning("Port direction '%s' should be 'inout'" % ptype)
                elif ptype == "inout":
                    if len(gScopeList) > 0 and gScopeList[-1].startswith("FUNCTION::"):
                        warning("'inout' function arguments are not recommended")
                elif ptype in ["real", "integer", "string"]:
                    # input x; real x; OK for functions
                    if gStyle and (len(gScopeList) == 0 or not gScopeList[-1].startswith("FUNCTION::")):
                        style("Variable declaration on same line as port declaration")
                elif ptype in gDisciplines:
                    # inout a; electrical a;
                    pass
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
                if gCurrentFunc:
                    if pname in gCurrentFunc.args:
                        if ptype in ["real", "integer", "string"]:
                            # input x; real x;
                            scope = getCurrentScope()
                            vname = scope + pname
                            valid = checkIdentifierCollisions(pname, vname, escaped, "Variable")
                            # insert into global map
                            if valid and in_module:
                                var = Variable(vname, ptype, False, gFileName[-1], gLineNo[-1])
                                gVariables[vname] = var
                                if pname in gCurrentFunc.inputs:
                                    # hack, not really where it was assigned a value
                                    var.assign = gLineNo[-1]
                                elif pname in gCurrentFunc.outputs:
                                    # marks the variable within the function as used
                                    # (separate question whether the output is used outside the function)
                                    var.used = True
                        else:
                            error("Function '%s' already had an argument '%s'"
                                  % (gCurrentFunc.name, pname))
                    else:
                        gCurrentFunc.args.append(pname)
                        if ptype in ["input", "inout"]:
                            gCurrentFunc.inputs.append(pname)
                            gCurrentFunc.outarg_used.append(True)
                        else:
                            gCurrentFunc.outputs.append(pname)
                            gCurrentFunc.outarg_used.append(False)
                else:  # pragma: no cover
                    fatal("Programming error %s" % gScopeList[-1])
            elif ptype in ["real", "integer", "string"]:
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
                    if ptype in ["inout", "input", "output"]:
                        if port.direction != "":
                            error("Port (terminal) '%s' already had direction '%s'" % (pname, port.direction))
                        else:
                            port.direction = ptype
                        if len(bus_range) == 2:
                            port.is_bus = True
                            port.msb = bus_range[0]
                            port.lsb = bus_range[1]
                    elif ptype in gDisciplines:
                        if port.discipline != "":
                            error("Duplicate specification of discipline for port (terminal) '%s'" % pname)
                        port.discipline = ptype
                        if ptype == "electrical":
                            port.mult = True
                        if bus_range:
                            if port.is_bus:
                                if port.msb != bus_range[0] or port.lsb != bus_range[1]:
                                    error("Port '%s' has different bus ranges in port and discipline declarations"
                                          % pname)
                            elif port.direction != "":
                                error("Discipline declaration for scalar port '%s' has bus range [%d:%d]"
                                      % (pname, bus_range[0], bus_range[1]))
                            # else have not seen direction declaration
                        elif port.is_bus:
                            error("Discipline declaration for port '%s' should also have bus range [%d:%d]"
                                  % (pname, port.msb, port.lsb))
                    else:
                        error("Unexpected identifier '%s' in port declaration" % ptype)
            val = parser.lex()
        elif parser.isNumber():
            num = parser.getNumber()
            if ptype in ["real", "integer", "string"]:
                error("Unexpected number %g in variable declaration" % num)
            else:
                error("Unexpected number %g in port declaration" % num)
            val = parser.lex()
        else:
            if ptype == "":  # pragma: no cover
                fatal("Parse error looking for port direction")
            elif ptype in ["inout", "input", "output", "real", "integer", "string"]:
                error("Missing identifier after '%s' " % ptype)
            else:
                error("Missing port direction before '%s' " % ptype)

        if val == ord(','):
            val = parser.lex()
        elif val == ord(';'):
            ptype = ""
            bus_range = []
            val = parser.lex()
            while val == ord(';'):
                error("Extra ';' at end of port direction declaration")
                val = parser.lex()
        elif parser.isIdentifier():
            string = parser.getString()
            if string in ["inout", "input", "output", "real", "integer", "string"]:
                error("Missing ';' after list of ports")
                ptype = string
                val = parser.lex()
            else:
                error("Missing ',' in list of ports")
        elif val != 0:
            error("Invalid %s in port declaration" % formatChar(val))
            val = parser.lex()
    if line[-1] != ';':
        error("Missing ';' at end of port direction declaration")
# end of parsePortDirection


def checkDeclarationContext( decl_type ):
    if gStatementInCurrentBlock:
        if len(gScopeList) > 1 and gScopeList[-1] == "":
            error("%s declarations only allowed in named blocks" % decl_type)
        else:
            error("%s declarations must preceed all statements" % decl_type)


def classifyContribs():
    global gCheckMultFactors
    for contrib in gContribs.values():
        acc = contrib.acc
        func = ""
        node1 = ""
        node2 = ""
        bran = ""
        st = ""
        for ch in acc:
            if ch == '(':
                if func or node1 or node2:
                    # if (cond) I(a,b) <+ expr -> restart parsing
                    node1 = ""
                    node2 = ""
                func = st
                st = ""
            elif ch == ',':
                node1 = st
                st = ""
            elif ch == ')':
                if node1 == "":
                    node1 = st
                else:
                    node2 = st
                st = ""
            elif ch.isspace():
                if func == "":
                    st = ""
            elif ch in ['<','>']:
                pass
            else:
                st += ch
        if node1 != "" and node1.find('[') > 0:
            parts = node1.split('[')
            node1 = parts[0]
        if node2 != "" and node2.find('[') > 0:
            parts = node2.split('[')
            node2 = parts[0]
        disc = ""
        if node1 in gPortnames:
            disc = gPortnames[node1].discipline
        elif node1 in gNodenames:
            disc = gNodenames[node1].discipline
        elif node1 in gBranches:
            bran = node1
            br = gBranches[bran]
            disc = br.discipline
            node1 = br.node1
            node2 = br.node2
        contrib.func = func
        contrib.node1 = node1
        contrib.node2 = node2
        contrib.branch = bran
        contrib.discipline = disc
    more = True
    while more:
        more = False
        for contrib in gContribs.values():
            if contrib.discipline == "electrical":
                node1 = contrib.node1
                node2 = contrib.node2
                scale = False
                scale = ""
                if node1 in gPortnames:
                    scale = True
                    scale = node1
                    port = gPortnames[node1]
                    if not port.mult:
                        port.mult = True
                        gMultScaled[node1] = 1
                        more = True
                elif node1 in gNodenames and gNodenames[node1].mult:
                    scale = True
                    scale = node1
                if node2 in gPortnames:
                    scale = True
                    scale = node2
                    port = gPortnames[node2]
                    if not port.mult:
                        port.mult = True
                        gMultScaled[node2] = 1
                        more = True
                elif node2 in gNodenames and gNodenames[node2].mult:
                    scale = True
                    scale = node2
                if scale:
                    if node1 in gNodenames:
                        node = gNodenames[node1]
                        if not node.mult:
                            node.mult = True
                            gMultScaled[node1] = 1
                            more = True
                    if node2 in gNodenames:
                        node = gNodenames[node2]
                        if not node.mult:
                            node.mult = True
                            gMultScaled[node2] = 1
                            more = True
                    if contrib.branch in gBranches:
                        gBranches[contrib.branch].mult = True
    gCheckMultFactors = 2
# end of classifyContribs


def registerContrib( type, node1, node2, cond_in, bias_dep_in ):
    [conds, bias_dep] = getCurrentConditions(cond_in, bias_dep_in)
    bname = ""
    bprint = node1
    if node2 == "":
        if node1 in gBranches:
            bname = node1
        else:
            bname = "UNNAMED_" + node1 + "_GND"
    else:
        bprint = node1 + "," + node2
        bname = "UNNAMED_" + node2 + "_" + node1
        if not bname in gBranches:
            # only one unnamed branch between two nets
            bname = "UNNAMED_" + node1 + "_" + node2
    if bname in gBranches:
        branch = gBranches[bname]
    else:
        branch = Branch(bname)
        branch.node1 = node1
        branch.node2 = node2
        gBranches[bname] = branch
    if branch.conds:
        if conds:
            if conds not in branch.conds:
                branch.conds = mergeConditions(branch.conds, conds)
        else:
            branch.conds = []
    elif conds and branch.lhs_flow == 0 and branch.lhs_pot == 0:
        # first contribution
        branch.conds.append(conds)
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
    if branch.lhs_flow and branch.lhs_pot and conds == "":
        warning("Switch branch (%s) collapsed by %s contribution" % (bprint, type))


def getFactors( expr ):
    ret = []
    if expr.type == "NAME":
        if expr.e1 in gParameters:
            ret.append(expr.e1)
        else:
            vname = findVariableInScope(expr.e1)
            if vname in gVariables:
                ret += gVariables[vname].factors
    elif expr.type == "FUNCCALL":
        ret.append(expr.e1)
        if expr.e1 in ["ddt", "ddx", "white_noise", "flicker_noise"]:
            ret += getFactors(expr.args[0])
    elif expr.type == "*":
        ret += getFactors(expr.e1)
        ret += getFactors(expr.e2)
    elif expr.type == "/":
        ret += getFactors(expr.e1)
    elif expr.type == "-" and not expr.e2:
        ret += getFactors(expr.e1)
    elif expr.type in ["+", "-"]:
        fac1 = getFactors(expr.e1)
        fac2 = getFactors(expr.e2)
        for mname in gMultFactors:
            if (mname in fac1 and mname not in fac2) or \
               (mname in fac2 and mname not in fac1):
                error("Combining factors with different %s dependence" % mname)
        facs = [fac for fac in fac1 if fac in fac2]
        ret += facs
    return ret


def countMultFactors( mname, factors ):
    count = 0
    errs = []
    for fac in factors:
        if fac == mname:
            count += 1
        elif fac in gMultFactors and fac not in errs:
            error("Contribution has unexpected factor of %s (expected %s)" % (fac, mname))
            errs.append(fac)
    if count > 1:
        error("Contribution has extra factor of %s" % mname)


def checkMultFactors( rhs ):
    parts = []
    more = [rhs]

    # split all additive terms
    while len(parts) != len(more):
        parts = more
        more = []
        for part in parts:
            if part.type in ["+", "-"]:
                more.append(part.e1)
                if part.e2:
                    more.append(part.e2)
            else:
                more.append(part)
    parts = more

    # classify
    for part in parts:
        factors = getFactors(part)
        if "flicker_noise" in factors:
            if "mult_fn" in factors:
                countMultFactors("mult_fn", factors)
            elif "MULT_FN" in factors:
                countMultFactors("MULT_FN", factors)
            else:
                if "MULT_I" in factors:
                    error("Flicker noise contribution uses MULT_I instead of MULT_FN")
                elif "mult_i" in factors:
                    error("Flicker noise contribution uses mult_i instead of mult_fn")
                elif "MULT_FN" in gParameters:
                    error("Flicker noise contribution missing MULT_FN")
                else:
                    error("Flicker noise contribution missing mult_fn")
        elif "ddt" in factors:
            if "mult_q" in factors:
                countMultFactors("mult_q", factors)
            elif "MULT_Q" in factors:
                countMultFactors("MULT_Q", factors)
            else:
                if "MULT_I" in factors:
                    error("Charge contribution uses MULT_I instead of MULT_Q")
                elif "mult_i" in factors:
                    error("Charge contribution uses mult_i instead of mult_q")
                elif "MULT_Q" in gParameters:
                    error("Charge contribution missing MULT_Q")
                else:
                    error("Charge contribution missing mult_q")
        else:
            if "mult_i" in factors:
                countMultFactors("mult_i", factors)
            elif "MULT_I" in factors:
                countMultFactors("MULT_I", factors)
            elif "MULT_I" in gParameters:
                error("Current contribution missing MULT_I")
            else:
                error("Current contribution missing mult_i")


def checkPorts():
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
            elif gStyle and gCompactModel:
                style("Port names should be lowercase ('%s' not '%s')" % (port.name.lower(), port.name))
        gLineNo.pop()
        gFileName.pop()


def handleEqualityConditions(expr, parms, excls):
    pname = ""
    value = 0
    handled = True
    if expr.type == "FUNCCALL" and expr.e1 == "NOT" and len(expr.args) == 1:
        expr = expr.args[0]
        if expr.type == "==":
            expr.type = "!="
        elif expr.type == "!=":
            expr.type = "=="
        else:
            handled = False
    elif expr.type == "NAME":
        if expr.e1 in gParameters:
            pname = expr.e1
            expr.e1 = Expression("NAME")
            expr.e1.e1 = pname
            expr.e2 = Expression("NUMBER")
            expr.e2.value = 0
            expr.type = "!="
        else:
            handled = False
    if handled:
        if expr.e1.type == "NAME" and expr.e1.e1 in gParameters:
            pname = expr.e1.e1
            value = expr.e2
        elif expr.e2.type == "NAME" and expr.e2.e1 in gParameters:
            pname = expr.e2.e1
            value = expr.e1
    if pname:
        if value.type == "NUMBER":
            val = value.number
            if expr.type == "==":
                if pname in parms:
                    vals = parms[pname]
                    if val not in vals:
                        vals.append(val)
                else:
                    vals = [val]
                parms[pname] = vals
            else:
                if pname in excls:
                    vals = excls[pname]
                    if val not in vals:
                        vals.append(val)
                else:
                    vals = [val]
                excls[pname] = vals
        else:
            handled = False
    else:
        handled = False
    return [handled, parms, excls]


def handleAndConditions(expr, parms, excls):
    if expr.e1.type in ["==", "!=", "FUNCCALL", "NAME"]:
        [ha1, parms, excls] = handleEqualityConditions(expr.e1, parms, excls)
    else:
        ha1 = False
    if expr.e2.type in ["==", "!=", "FUNCCALL", "NAME"]:
        [ha2, parms, excls] = handleEqualityConditions(expr.e2, parms, excls)
    else:
        ha2 = False
    handled = False
    if ha1 and ha2 and len(parms) == 1 and len(excls) == 1:
        # par==1 and par!=0 is the same as par==1
        popkey = False
        for (k1, v1) in parms.items():
            for (k2, v2) in excls.items():
                if k1 == k2 and len(v2) == 1 and v2[0] not in v1:
                    handled = True
                    popkey = k2
        if popkey:
            excls.pop(popkey)
    return [handled, parms, excls]


def handleOrConditions(expr, parms, excls):
    if expr.e1.type == "||":
        [ha1, parms, excls] = handleOrConditions(expr.e1, parms, excls)
    elif expr.e1.type in ["==", "!=", "FUNCCALL", "NAME"]:
        [ha1, parms, excls] = handleEqualityConditions(expr.e1, parms, excls)
    else:
        ha1 = False
    if expr.e2.type == "||":
        [ha2, parms, excls] = handleOrConditions(expr.e2, parms, excls)
    elif expr.e2.type in ["==", "!=", "FUNCCALL", "NAME"]:
        [ha2, parms, excls] = handleEqualityConditions(expr.e2, parms, excls)
    else:
        ha2 = False
    handled = ha1 and ha2
    return [handled, parms, excls]


def checkInternalNodes():
    for node in gNodenames.values():
        gFileName.append(node.declare[0])
        gLineNo.append(node.declare[1])
        conds = []
        found_contrib = 0
        bad_branch = False
        for branch in gBranches.values():
            if node.name in [branch.node1, branch.node2] and (branch.lhs_pot or branch.lhs_flow):
                if branch.conds:
                    found_contrib = 1
                    for co in branch.conds:
                        conds = mergeConditions(conds, co, True)
                        if not conds:
                            found_contrib = 2
                            break
                    if branch.lhs_pot == 2 and branch.lhs_flow == 0:
                        # bias-dependent conditional voltage contribution
                        found_contrib = 3
                        bad_branch = branch
                else:
                    found_contrib = 2
                if found_contrib == 2:
                    break
        if found_contrib == 1:
            do_report = True
            parms = {}
            excls = {}
            for co in conds:
                parser = Parser(co)
                expr = parser.getExpression()
                handled = True
                if expr.type in ["==", "!=", "FUNCCALL", "NAME"]:
                    [handled, parms, excls] = handleEqualityConditions(expr, parms, excls)
                elif expr.type == "||":
                    [handled, parms, excls] = handleOrConditions(expr, parms, excls)
                elif expr.type == "&&":
                    [handled, parms, excls] = handleAndConditions(expr, parms, excls)
                elif expr.type in [">", ">=", "<", "<="]:
                    pname = ""
                    if expr.e1.type == "NAME" and expr.e1.e1 in gParameters:
                        pname = expr.e1.e1
                    elif expr.e2.type == "NAME" and expr.e2.e1 in gParameters:
                        pname = expr.e2.e1
                    if pname:
                        parms[pname] = []
                    else:
                        handled = False
                else: # pragma: no cover
                    handled = False
                if not handled:
                    do_report = False
            if do_report:
                for (pname, vals) in parms.items():
                    parm = gParameters[pname]
                    if parm.type == "integer" and parm.min.type == "NUMBER" and parm.max.type == "NUMBER":
                        start = int(parm.min.number)
                        stop = int(parm.max.number)
                        if parm.value_range in ["oc", "oo"]:
                            start += 1
                        if parm.value_range in ["co", "oo"]:
                            stop -= 1
                        if parm.min.number >= -10 and parm.max.number <= 10:
                            missing = []
                            for v in range(start,stop+1):
                                if not v in vals and v not in parm.exclude:
                                    missing.append(v)
                                    if pname in excls:
                                        excl = excls[pname]
                                        if v not in excl:
                                            missing.pop()
                            if not missing:
                                pass
                            elif len(missing) == 1:
                                warning("No contributions to node '%s' when parameter %s=%d"
                                        % (node.name, pname, missing[0]))
                            elif missing:
                                warning("No contributions to node '%s' when parameter %s=%s"
                                        % (node.name, pname, missing))
                        else:
                            if len(vals) < len(range(start,stop+1)) - len(parm.exclude):
                                warning("Contributions to node '%s' only for some values of parameter '%s'"
                                        % (node.name, pname))
                    elif parm.type == "real" and len(conds) == 1:
                        warning("Contributions to node '%s' only when %s" % (node.name, conds[0]))
                for (pname, vals) in excls.items():
                    if len(vals) == 1 and pname not in parms:
                        warning("No contributions to node '%s' when parameter %s=%d" % (node.name, pname, vals[0]))
            elif gVerbose:
                notice("Node '%s' has contributions under conditions %s" % (node.name, conds))
        elif found_contrib == 3:
            bprint = bad_branch.node1
            if bad_branch.node2:
                bprint += "," + bad_branch.node2
            warning("Switch branch (%s) with bias-dependent condition" % bprint)
        elif found_contrib == 0:
            error("No contributions found for node '%s'" % node.name)
        gLineNo.pop()
        gFileName.pop()
# end of checkInternalNodes


def checkMultDeclaration( parm ):
    if parm.model or not parm.inst or parm.type != "real":
        error("Parameter %s should be a real instance parameter" % parm.name)
    if parm.name == "mult_fn":
        if parm.defv.type != "NAME" or parm.defv.e1 != "mult_i":
            error("Parameter mult_fn default value is %s, should be mult_i" % parm.defv.asString())
    elif parm.name == "MULT_FN":
        if parm.defv.type != "NAME" or parm.defv.e1 != "MULT_I":
            error("Parameter MULT_FN default value is %s, should be MULT_I" % parm.defv.asString())
    else:
        if parm.defv.type != "NUMBER" or parm.defv.number != 1.0:
            error("Parameter %s default value is %s, should be 1.0" % (parm.name, parm.defv.asString()))
    if parm.value_range != "co" or parm.min.type != "NUMBER" or parm.min.number != 0 \
            or parm.max.type != "NAME" or parm.max.e1 != "inf":
        error("Parameter %s range should be [0:inf)" % parm.name)


def getPrimaryParameter( name, alt ):
    ret = False
    if name in gParameters:
        ret  = gParameters[name]
    if alt in gParameters:
        if ret and not ret.is_alias and not gParameters[alt].is_alias:
            error("Both '%s' and '%s' declared" % (name, alt))
        else:
            ret = gParameters[alt]
    return ret


def checkParmsAndVars():
    temp = gVariables["$temperature"]
    has_tref = False
    has_dtemp = False
    all_lower = []
    all_upper = []
    not_lower = []
    not_upper = []
    some_units = False
    some_inst_parms = False
    param_names_lower = []
    for branch in gBranches.values():
        n1 = branch.node1
        n2 = branch.node2
        if n1 in gPortnames:
            gPortnames[n1].used = True
        elif n1 in gNodenames:
            gNodenames[n1].used = True
        if n2 in gPortnames:
            gPortnames[n2].used = True
        elif n2 in gNodenames:
            gNodenames[n2].used = True
    for port in gPortnames.values():
        if not port.used:
            gFileName.append(port.declare[0])
            gLineNo.append(port.declare[1])
            warning("Port '%s' was never used" % port.name)
            gLineNo.pop()
            gFileName.pop()
    for node in gNodenames.values():
        if not node.used:
            gFileName.append(node.declare[0])
            gLineNo.append(node.declare[1])
            warning("Node '%s' was never used" % node.name)
            gLineNo.pop()
            gFileName.pop()
    for par in gParameters.values():
        if par.inst:
            some_inst_parms = True
            break
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
                            warning("Reference temperature parameter '%s' (aliased as %s) should not be %s"
                                    % (par2.name, par.name, "an instance parameter"))
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
                            warning("Temperature offset parameter '%s' (aliased as %s) should be %s"
                                    % (par2.name, par.name, "an instance parameter"))
            elif some_inst_parms and not par.inst:
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
        if par.units:
            some_units = True
        gLineNo.pop()
        gFileName.pop()
    if temp.used > 0:
        if not has_tref:
            warning("Module uses $temperature but does not have 'tref'")
        if not has_dtemp:
            warning("Module uses $temperature but does not have 'dtemp'")
    if len(not_lower) > 0:
        if len(not_upper) == 0:
            if gStyle and gCompactModel:
                style("Parameters should be declared in lowercase")
        else:
            if len(not_lower) == 1 and len(all_lower) > 5:
                bad = not_lower[0]
                warning("Parameters have inconsistent case (all lowercase except '%s')" % bad)
            elif len(not_upper) == 1 and len(all_upper) > 5:
                bad = not_upper[0]
                warning("Parameters have inconsistent case (all uppercase except '%s')" % bad)
            else:
                warning("Parameters have inconsistent case (prefer all lowercase)")
    if not some_units and len(gParameters) > 0 and gCompactModel:
        warning("No units specified for any parameters")
    mult_i  = getPrimaryParameter("mult_i",  "MULT_I")
    mult_q  = getPrimaryParameter("mult_q",  "MULT_Q")
    mult_fn = getPrimaryParameter("mult_fn", "MULT_FN")
    if mult_i or mult_q or mult_fn:
        if not mult_i or (not mult_q and gUsesDDT) or (not mult_fn and "flicker_noise" in gNoiseTypes):
            if len(not_lower) == 0:
                error("Parameters mult_i, mult_q, mult_fn must be implemented together.")
            else:
                error("Parameters MULT_I, MULT_Q, MULT_FN must be implemented together.")
        if mult_i:
            checkMultDeclaration(mult_i)
        if mult_q:
            checkMultDeclaration(mult_q)
        if mult_fn:
            checkMultDeclaration(mult_fn)
    elif gVerbose and gCompactModel:
        if len(not_lower) == 0:
            notice("This compact model does not implement mult_i, mult_q, mult_fn")
        else:
            notice("This compact model does not implement MULT_I, MULT_Q, MULT_FN")

    for var in gVariables.values():
        if var.func_name_arg[0] != "":
            # was assigned as a function output variable
            if var.used:
                fname = var.func_name_arg[0]
                argno = var.func_name_arg[1]
                if fname in gUserFunctions:
                    func = gUserFunctions[fname]
                    if not func.outarg_used[argno]:
                        func.outarg_used[argno] = True
                        if gDebug:
                            print("    Function '%s' output argument '%s' is used because variable '%s' is used"
                                  % (fname, func.args[argno], var.name))
                else:  # pragma: no cover
                    fatal("Programming error: function arguments")
            else:
                var.used = True  # don't report as unused
        if var.name.startswith("FUNCTION::"):
            var.assign = 0  # suppress printing in summary
        else:
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
            elif var.oppt and var.conditions != []:
                gFileName.append(var.declare[0])
                gLineNo.append(var.declare[1])
                error("Operating-point variable '%s' is only conditionally assigned a value" % var.name)
                gLineNo.pop()
                gFileName.pop()
            if gSuperfluous and var.used and var.accesses and not var.range and var.name not in gHiddenState:
                if var.oppt:
                    # oppt variable used at end
                    acc = Access("USE", var.declare[0], 999999, 0, [])
                    var.accesses.append(acc)
                for i in range(len(var.accesses)):
                    acc_i = var.accesses[i]
                    if acc_i.type == "USE":
                        have_unused = False
                        for j in range(i):
                            acc_j = var.accesses[j]
                            if not acc_j.used and acc_j.type == "SET":
                                have_unused = True
                                break
                        if have_unused:
                            set_conditions = []
                            for j in range(i-1,-1,-1):
                                acc_j = var.accesses[j]
                                done = False
                                if acc_j.type == "SET":
                                    for cond in acc_j.cond:
                                        if cond in acc_i.cond:
                                            # used under same conditions as set
                                            acc_j.used = True
                                            done = True
                                            break
                                        if not cond:
                                            # set with no conditions
                                            acc_j.used = True
                                            done = True
                                            break
                                        acc_j.used = True
                                        if set_conditions:
                                            set_conditions = mergeConditions(set_conditions, cond)
                                            if set_conditions == [] and (len(cond) < 3 or cond[0:3] != "NOT"):
                                                done = True
                                                break
                                        else:
                                            set_conditions = acc_j.cond
                                if done:
                                    break
                            if acc_i.loop:
                                set_in_loop = False
                                for j in range(i):
                                    acc_j = var.accesses[j]
                                    if acc_j.loop == acc_i.loop and acc_j.type == "SET":
                                        acc_j.used = True
                                        set_in_loop = True
                                if not set_in_loop:
                                    # look for last assignment in the loop
                                    for j in range(len(var.accesses)-1,i,-1):
                                        acc_j = var.accesses[j]
                                        if acc_j.loop == acc_i.loop and acc_j.type == "SET" and \
                                                (acc_j.file != acc_i.file or \
                                                 acc_j.line != acc_i.line or acc_j.subl != acc_i.subl):
                                            acc_j.used = True
                                            break
                zini = False
                init = True
                for acc in var.accesses:
                    if acc.type == "SET":
                        if acc.zini == 1:
                            if zini:
                                if gVerbose or isinstance(zini, str):
                                    gFileName.append(acc.file)
                                    gLineNo.append(acc.line)
                                    notice("Repeated initialization of %s=0" % var.name, "RepeatedInit")
                                else:
                                    gFileName.append(zini.file)
                                    gLineNo.append(zini.line)
                                    if zini.file == acc.file:
                                        notice("Repeated initialization of %s=0, also line %d" % (var.name, acc.line))
                                    else:
                                        notice("Repeated initialization of %s=0" % var.name)
                                        gFileName[-1] = acc.file
                                        gLineNo[-1] = acc.line
                                        notice("Repeated initialization of %s=0" % var.name)
                                gFileName.pop()
                                gLineNo.pop()
                                zini = "DONE"
                            else:
                                zini = acc
                        if not acc.used:
                            gFileName.append(acc.file)
                            gLineNo.append(acc.line)
                            if acc.zini != 1:
                                if acc.zini == 2 and init:
                                    what = "initialization of"
                                else:
                                    what = "assignment to"
                                    init = False
                                if gVerbose and acc.subl > 0:
                                    notice("Superfluous %s %s (macro line %d)" % (what, var.name, acc.subl))
                                else:
                                    notice("Superfluous %s %s" % (what, var.name))
                            elif gVerbose:
                                if acc.subl > 0:
                                    notice("Unnecessary initialization of %s=0 (macro line %d)" % (var.name, acc.subl))
                                else:
                                    notice("Unnecessary initialization of %s=0" % var.name)
                            gFileName.pop()
                            gLineNo.pop()
                        else:
                            init = False
                    else:
                        zini = False
            if var.oppt and var.name.lower() in param_names_lower:
                parname = ""
                for par in gParameters.values():
                    if par.name.lower() == var.name.lower():
                        parname = par.name
                        break
                warning("Operating-point variable '%s' differs only in case from parameter '%s'" % (var.name, parname))
            if var.oppt and (mult_i or mult_q):
                if var.units in ["A", "A/V", "mho"]:
                    if not mult_i or not mult_i.name in var.factors:
                        if mult_i:
                            mname = mult_i.name
                        elif len(not_upper) == 0:  # pragma: no cover
                            mname = "MULT_I"
                        else:
                            mname = "mult_i"
                        warning("Operating-point variable '%s' with units '%s' not scaled by %s"
                                % (var.name, var.units, mname))
                elif var.units in ["C", "F"]:
                    if not mult_q or not mult_q.name in var.factors:
                        if mult_q:
                            mname = mult_q.name
                        elif len(not_upper) == 0:
                            mname = "MULT_Q"
                        else:
                            mname = "mult_q"
                        warning("Operating-point variable '%s' with units '%s' not scaled by %s"
                                % (var.name, var.units, mname))

    for func in gUserFunctions.values():
        if func.used:
            if len(func.outputs) > 0:
                for i in range(len(func.args)):
                    if not func.outarg_used[i]:
                        warning("Function '%s' output argument '%s' is never used"
                                % (func.name, func.args[i]) )
        elif gVerbose:
            notice("User-defined function '%s' is never used" % func.name)
    for mac in gMacros.values():
        if gVerbose and not mac.used:
            notice("Macro '%s' is never used" % mac.name)
# end of checkParmsAndVars


def findVariableInScope( name ):
    scope = getCurrentScope()
    while 1:
        vname = scope + name
        if vname in gVariables:
            return vname
        if scope != "":
            scopes = scope.split(".")
            scope = ""
            if len(scopes) > 2:
                for sc in scopes[0:-2]:
                    scope += sc + "."
        else:
            break
    return ""


def getIdentifierType( name ):
    type = ""
    vname = findVariableInScope(name)
    if vname in gVariables:
        type = gVariables[vname].type
    elif name in gParameters:
        type = gParameters[name].type
    return type


def checkDependencies( deplist, unset_message, bias_dep_in, is_condition, no_ddt, do_warn ):
    global gUsesDDT
    [conds, bias_dep] = getCurrentConditions(Expression("NOTHING"), bias_dep_in)
    biases = []
    ddt_call = 0
    if is_condition:
        # condition for an if, for, while, case statement
        # don't use bias_dep from getCurrentConditions
        bias_dep = bias_dep_in
    for dep in deplist:
        if dep.find("(") > 0:
            biases.append(dep)
            bias_dep = max(bias_dep, 2)
            is_set = True
            found = True
        elif dep == "ddx":
            bias_dep = 4
            is_set = True
            found = True
        elif dep == "ddt":
            is_set = True
            ddt_call = 1
            gUsesDDT = True
            found = True
        else:
            is_set = False
            found  = False
        vn = findVariableInScope(dep)
        if vn in gVariables:
            found = True
            gVariables[vn].used = True
            if do_warn:
                accesses = gVariables[vn].accesses
                if len(accesses) and accesses[-1].type == "USE" and accesses[-1].cond == [conds]:
                    pass
                else:
                    access = Access("USE", gFileName[-1], gLineNo[-1], gSubline, [conds])
                    for lt in gLoopTypes:
                        if lt[0] in ["for", "repeat", "while"]:
                            access.loop = lt
                            break
                    gVariables[vn].accesses.append(access)
            if bias_dep < gVariables[vn].bias_dep:
                bias_dep = gVariables[vn].bias_dep
            for bias in gVariables[vn].biases:
                if bias not in biases:
                    biases.append(bias)
            if ddt_call < gVariables[vn].ddt_use:
                ddt_call = gVariables[vn].ddt_use
            assign = gVariables[vn].assign
            if assign >= 0:
                is_set = True
            elif assign <= -2:
                is_set = True
                if assign == -2 and not unset_message.startswith("Task") and do_warn and gCompactModel:
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
                    if vn not in gHiddenState:
                        gHiddenState[vn] = vn
        if not found:
            if dep in gParameters:
                if gParameters[dep].is_alias:
                    par = gParameters[dep].defv.e1
                    if par in gParameters:
                        error("Reference to aliasparam '%s' should use parameter '%s'" % (dep, par))
                    else:
                        error("Reference to aliasparam '%s'" % dep)
                if do_warn:
                    gParameters[dep].used += 1
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
            elif dep.find(".") > 0:
                parts = dep.split(".")
                dname = ""
                if parts[0] in gPortnames:
                    dname = gPortnames[parts[0]].discipline
                elif parts[0] in gNodenames:
                    dname = gNodenames[parts[0]].discipline
                if dname != "" and dname in gDisciplines:
                    #disc = gDisciplines[dname]
                    if len(parts) == 3 and parts[1] in ["potential", "flow"] and parts[2] == "abstol":
                        is_set = True
        if not is_set and do_warn:
            if dep in gVAMSkeywords:
                error("%s Verilog-AMS keyword '%s'" % (unset_message, dep))
            # $xposition, $yposition, $angle, $hflip, $vflip
            elif dep in gVAMShiersysparm:
                warning("%s hierarchical system parameter '%s'" % (unset_message, dep))
            else:
                error("%s %s, which has not been set" % (unset_message, dep))
    if ddt_call and no_ddt:
        error("%s argument involving ddt()" % unset_message)
    return [bias_dep, biases, ddt_call]
# end of checkDependencies


def registerLoopVariableUse( loop_cond, loop_info ):
    deps = loop_cond.getDependencies(False, False)
    conds = getCurrentConditions(Expression("NOTHING"), 0)[0]
    done = []
    for dep in deps:
        vn = findVariableInScope(dep)
        if vn in gVariables:
            access = Access("USE", gFileName[-1], gLineNo[-1], gSubline, [conds])
            gVariables[vn].accesses.append(access)
            done.append(vn)
    if loop_info:
        for (vn, var) in gVariables.items():
            if vn not in done:
                for i in range(len(var.accesses)):
                    acc = var.accesses[i]
                    if acc.loop:
                        if acc.type == "USE":
                            if acc.file == loop_info[1] and acc.line > loop_info[2] \
                                and (i+1 == len(var.accesses) or var.accesses[i+1].type != "SET" \
                                     or acc.line != var.accesses[i+1].line):
                                access = Access("USE", gFileName[-1], gLineNo[-1], gSubline, [conds])
                                var.accesses.append(access)
                        elif acc.type == "SET" and acc.file == loop_info[1] and acc.line > loop_info[2]:
                            break


# split on "&&", but watch out for parentheses
def splitConditions( string ):
    retval = []
    num_parens = 0
    cond = ""
    i = 0
    while i < len(string):
        ch = string[i]
        if ch == '(':
            num_parens += 1
        elif ch == ')':
            num_parens -= 1
        if num_parens == 0 and i+1 < len(string) and ch == '&' and string[i+1] == '&':
            retval.append(cond)
            cond = ""
            i += 2
        else:
            cond += ch
            i += 1
    if cond != "":
        retval.append(cond)
    if num_parens != 0:  # pragma: no cover
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
                    if op not in both:  # pragma: no cover
                        fatal("Programming error in conditionCovered")
                else:
                    newonly.append(op)
            if len(oldonly) == 0:
                retval = True
                break
    return retval


def mergeConditions( oldlist, new, do_or=False ):
    retlist = []
    add_new = new
    did_merge = False
    if new in oldlist:
        # already have this condition
        return oldlist
    newparts = splitConditions(new)
    if len(newparts) == len(oldlist):
        all_found = True
        for old in oldlist:
            not_old = "NOT(" + old + ")"
            if not_old not in newparts:
                all_found = False
                break
        if all_found:
            return []
    for old in oldlist:
        if new == "NOT(" + old + ")":
            # A and !(A)
            add_new = False
            retlist = []
            break
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
                if op not in both:  # pragma: no cover
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
                (oldonly[0] == "NOT(" + newonly[0] + ")" or
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
            elif do_or and oldonly[0] == "NOT(" + newonly[0] + ")":
                # !A&&B.. and A -> B&&.. and A
                retval = ""
                for op in oldonly[1:]:
                    if retval != "":
                        retval += "&&"
                    retval += op
                for op in both:
                    if retval != "":
                        retval += "&&"
                    retval += op
                replace = True
        elif do_or and len(oldparts) == 1 and len(newparts) == 2:
            if newparts[0] == "NOT(" + oldparts[0] + ")":
                # A and !A&&B -> new condition is B
                add_new = newparts[1]
        if replace:
            retlist.append(retval)
        else:
            retlist.append(old)
    if add_new:
        retlist.append(add_new)
    elif did_merge and len(retlist) > 1:
        retlist = mergeConditions(retlist[:-1], retlist[-1])
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
    cond = conds
    idx = cond.find("NOT(")
    while idx >= 0:
        if idx > 0:
            cond = cond[0:idx] + "!(" + cond[idx+4:]
        else:
            cond = "!(" + cond[idx+4:]
        idx = cond.find("NOT(")
    return cond


def checkCondForBiasDepEquality( cond ):
    ret = False
    if cond.type in ["==", "!="]:
        is_integer_test = 0
        if cond.e1.type == "NUMBER" and cond.e1.is_int:
            is_integer_test = 1
        elif cond.e1.type == "NAME":
            if getIdentifierType(cond.e1.e1) == "integer":
                is_integer_test = 1
        elif cond.e1.type == "FUNCCALL":
            fname = cond.e1.e1
            if fname in gUserFunctions and gUserFunctions[fname].type == "integer":
                is_integer_test = 1
        if cond.e2.type == "NUMBER" and cond.e2.is_int:
            is_integer_test += 2
        elif cond.e2.type == "NAME":
            if getIdentifierType(cond.e2.e1) == "integer":
                is_integer_test += 2
        elif cond.e2.type == "FUNCCALL":
            fname = cond.e2.e1
            if fname in gUserFunctions and gUserFunctions[fname].type == "integer":
                is_integer_test += 2
        if is_integer_test != 3:
            deps = cond.getDependencies(False, False)
            bias_dep = checkDependencies(deps, "", 0, True, True, False)[0]
            if bias_dep > 1:
                ret = True
    elif cond.type == "&&":
        if checkCondForBiasDepEquality(cond.e1) or checkCondForBiasDepEquality(cond.e2):
            ret = True
    return ret


def simplifyBiases(bias_list):
    ret = bias_list
    if len(bias_list) == 3:
        pfx0 = bias_list[0].find("(")
        if pfx0 > 0:
            bname = bias_list[0][pfx0+1:-1]
            if bname in gBranches:
                branch = gBranches[bname]
                pfx1 = bias_list[1].find("(")
                pfx2 = bias_list[2].find("(")
                if pfx1 == pfx2 and bias_list[1][0:pfx1] == bias_list[2][0:pfx2]:
                    nname1 = bias_list[1][pfx1+1:-1]
                    nname2 = bias_list[2][pfx2+1:-1]
                    if nname1 == branch.node1 and nname2 == branch.node2:
                        ret = bias_list[0]
    elif len(bias_list) == 2:
        pfx0 = bias_list[0].find("(")
        pfx1 = bias_list[1].find("(")
        if pfx0 == pfx1 and bias_list[0][0:pfx0] == bias_list[1][0:pfx1]:
            ret = bias_list[0][0:-1] + "," + bias_list[1][pfx1+1:]
    return ret


def markVariableAsSet( varname, rhs, deps, bias_dep_in, cond_in, bias_in, ddt_in, for_var, set_by_func, fname, argnum ):
    found = False
    vn = findVariableInScope(varname)
    if vn in gVariables:
        var = gVariables[vn]
        found = True
        [conds, bias_dep] = getCurrentConditions(cond_in, bias_dep_in)
        biases = bias_in
        if var.type == "integer":
            if bias_dep > 1:
                bias_dep = 1
                biases = []
            if rhs and rhs.type == "NUMBER" and rhs.is_int == 0:
                warning("Integer variable '%s' assigned a real value %s" % (varname, rhs.asString()))
        elif var.type == "real" and rhs:
            if rhs.type in ["<", ">", "<=", ">=", "==", "!=", "&&", "||"]:
                warning("Real variable '%s' assigned a boolean (result of %s)" % (varname, rhs.type))
        first_assign = False
        if var.assign < 0:
            first_assign = True
            var.assign = gLineNo[-1]
            if conds != "":
                var.conditions = [conds]
            var.bias_dep = bias_dep
            var.biases = biases
        elif len(var.conditions) > 0:
            if conds in ["", "NOT(0)"]:
                # no conditions for this assignment
                var.conditions = []
                var.bias_dep = bias_dep
                var.biases   = biases
            else:
                oldconds = var.conditions
                var.conditions = mergeConditions(oldconds, conds)
        zero_init = 0
        if rhs and rhs.type == "NUMBER" and rhs.number == 0:
            if conds:
                zero_init = 2
            else:
                zero_init = 1
        access = Access("SET", gFileName[-1], gLineNo[-1], gSubline, [conds], zero_init)
        for lt in gLoopTypes:
            if lt[0] in ["for", "repeat", "while"]:
                access.loop = lt
                break
        var.accesses.append(access)
        if conds == "" or var.bias_dep < bias_dep:
            var.bias_dep = bias_dep
        elif bias_dep == 0 and var.bias_dep > 0:
            if var.prev_assign != [] and var.prev_assign[1] == 0 and \
                    conds == "NOT(" + var.prev_assign[0] + ")":
                var.bias_dep = 0
        if conds == "":  # replace previous biases
            var.biases   = biases
        else:  # append these biases, possibly overkill
            for bias in biases:
                if bias not in var.biases:
                    var.biases.append(bias)
        if rhs and rhs.type != "NONE":
            zeros = getParamZeros(rhs, False)
            factors = getFactors(rhs)
            if len(deps) > 0:
                for mname in gMultFactors:
                    if mname in deps and mname not in factors:
                        error("Variable '%s' has an invalid dependence on %s" % (varname, mname))
            if not first_assign and conds != "" and (len(var.factors) > 0 or len(factors) > 0):
                if var.non_zero_assign and (rhs.type != "NUMBER" or rhs.number != 0):
                    for mname in gMultFactors:
                        if (mname in factors and mname not in var.factors) or \
                           (mname not in factors and mname in var.factors):
                            error("Variable '%s' depends on %s differently in if/else blocks" % (varname, mname))
            if conds == "" or first_assign:
                var.zeros = zeros
                var.factors = factors
            else:
                if len(var.zeros) > 0 and len(zeros) > 0:
                    var.zeros = [zer for zer in var.zeros if zer in zeros]
                else:
                    var.zeros = []
                if len(var.factors) > 0 and len(factors) > 0:
                    var.factors = [fac for fac in var.factors if fac in factors]
                else:
                    var.factors = []
            if rhs.type != "NUMBER" or rhs.number != 0:
                var.non_zero_assign = True
        else:
            var.zeros = []
            var.factors = []
        if set_by_func:
            var.func_name_arg = [fname, argnum]
        if conds == "" or var.ddt_use < ddt_in:
            var.ddt_use = ddt_in
        elif ddt_in == 0 and var.ddt_use > 0:
            if var.prev_assign != [] and var.prev_assign[2] == 0 and \
                    conds == "NOT(" + var.prev_assign[0] + ")":
                var.ddt_use = 0
        var.prev_assign = [conds, bias_dep, ddt_in]

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


# recursive parsing of units expression
def getUnitFactors( expr ):
    num = []
    den = []
    mct = 0
    if expr.type == "NAME":
        num = [expr.e1]
        if expr.e1 == "m":
            mct += 1
    elif expr.type == "NUMBER":
        if gDebug and expr.number != 1.0:
            notice("Unexpected number %g in units expression" % expr.number)
    elif expr.type in ["*", "-"]:
        # "-" shows up in "ohm-um^WR"
        fac1 = getUnitFactors(expr.e1)
        fac2 = getUnitFactors(expr.e2)
        num = fac1[0] + fac2[0]
        den = fac1[1] + fac2[1]
        mct += fac1[2] + fac2[2]
    elif expr.type == "/":
        fac1 = getUnitFactors(expr.e1)
        fac2 = getUnitFactors(expr.e2)
        num = fac1[0] + fac2[1]
        den = fac1[1] + fac2[0]
        mct += fac1[2] - fac2[2]
    elif expr.type == "^":
        fac = getUnitFactors(expr.e1)
        if expr.e2.is_int:
            cnt = int(expr.e2.number)
            if cnt > 0:
                num = cnt * fac[0]
                den = cnt * fac[1]
            elif cnt < 0:
                num = (-cnt) * fac[1]
                den = (-cnt) * fac[0]
            mct += cnt * fac[2]
        else:
            # V^0.5 or m^(LLN+1)
            num = [expr.asString()]
            if expr.e2.type == "NUMBER":
                mct += expr.e2.number * fac[2]
    elif expr.type == "FUNCCALL" and len(expr.args) == 1:
        # ohm(um)^WR
        fname = expr.e1
        num.append(fname)
        if fname == "m":
            mct += 1
        fac = getUnitFactors(expr.args[0])
        num += fac[0]
        den += fac[1]
        mct += fac[2]
    return [num, den, mct]


# parse units string
def parseUnitsString( units ):
    num = []
    den = []
    mct = 0
    if units and units != "-" and units[0] != "?":
        if units[0] == '/':
            units = "1" + units
        parser = Parser(units)
        parser.parse_units = True
        expr = parser.getExpression()
        [num, den, mct] = getUnitFactors(expr)
    return [num, den, mct]


def compareUnits( num1, den1, cnt1, num2, den2, cnt2, add ):
    [numa, dena, cnta] = parseUnitsString(add)
    for fac in numa:
        if fac in num2:
            idx = num2.index(fac)
            num2[idx] = ""
        elif fac in den1:
            idx = den1.index(fac)
            den1[idx] = ""
        else:
            num1.append(fac)
    for fac in dena:
        if fac in den2:
            idx = den2.index(fac)
            den2[idx] = ""
        elif fac in num1:
            idx = num1.index(fac)
            num1[idx] = ""
        else:
            den1.append(fac)
    for i in range(len(num1)):
        fac = num1[i]
        if fac in num2:
            idx = num2.index(fac)
            num2[idx] = ""
            num1[i] = ""
    for i in range(len(den1)):
        fac = den1[i]
        if fac in den2:
            idx = den2.index(fac)
            den2[idx] = ""
            den1[i] = ""
    bad_num = []
    for fac in (num1+den2):
        if fac:
            bad_num.append(fac)
    bad_den = []
    for fac in (den1+num2):
        if fac:
            bad_den.append(fac)
    if bad_num or bad_den:
        if cnt1 + cnta == cnt2:
            # powers of 'm' match; see if there are any other units
            for fac in (bad_num+bad_den):
                if fac != "m" and \
                        (len(fac) < 2 or fac[0:2] != "m^") and \
                        (len(fac) < 3 or fac[0:3] != "(m^"):
                    return False
        else:
            return False
    return True


# check units of binning parameter against main parameter
def checkBinUnits( spar, prefix, name, units ):
    if isinstance(units, str) or isinstance(spar.units, str):
        [nump, denp, cntp] = parseUnitsString(units)
        [nums, dens, cnts] = parseUnitsString(spar.units)
        bad_units = False
        expect = "?"
        prefix = prefix.upper()
        if gBinningFactors and prefix in gBinningFactors:
            expect = gBinningFactors[prefix]

            if units in ("", "-"):
                if spar.units != expect:
                    bad_units = expect
            elif expect == "":
                if units != spar.units:
                    bad_units = units
            elif not compareUnits(nump, denp, cntp, nums, dens, cnts, expect):
                if units[0] == '/':
                    bad_units = expect + units
                else:
                    bad_units = units + "*" + expect

            if bad_units:
                if units:
                    if bad_units == units:
                        should = "(should be '%s' to match units for '%s')" % (bad_units, name)
                    else:
                        should = "(should be '%s' since units for '%s' are '%s')" % (bad_units, name, units)
                else:
                    should = "(should be '%s' since '%s' has no units)" % (bad_units, name)
                gFileName.append(spar.declare[0])
                gLineNo.append(spar.declare[1])
                if spar.units and spar.units != "-":
                    notice("Unexpected units '%s' for binning parameter '%s' %s"
                           % (spar.units, spar.name, should), "Unexpected units")
                else:
                    notice("Missing units for binning parameter '%s' %s"
                           % (spar.name, should), "Unexpected units")
                gFileName.pop()
                gLineNo.pop()
        else:
            if spar.units == units:
                expect = ""
            elif spar.units in (units+"*m", "m*"+units):
                expect = "m"
            elif spar.units == "m^2" and units == "m":
                expect = "m"
            elif spar.units in (units+"*m^2", "m^2*"+units, "(m^2)*"+units):
                expect = "m^2"
            elif spar.units == units + "/m":
                expect = "1/m"
            elif spar.units == units + "/m^2":
                expect = "1/m^2"
            elif len(units) > 2 and units[0:2] == "1/":
                if spar.units == "m" + units[1:]:
                    expect = "m"
                elif spar.units == "m^2" + units[1:]:
                    expect = "m^2"
            else:
                if cnts - cntp == 1:
                    expect = "m"
                elif cnts - cntp == 2:
                    expect = "m^2"

            if expect != "?":
                gBinningFactors[prefix] = expect
                if gDebug:
                    notice("Binning factor for %s-parameters set to '%s'" % (prefix, expect))
            elif units and (not spar.units or spar.units == "-"):
                gFileName.append(spar.declare[0])
                gLineNo.append(spar.declare[1])
                should = "(note that '%s' has units %s)" % (name, units)
                notice("Missing units for binning parameter '%s' %s" % (spar.name, should), "Unexpected units")
                gFileName.pop()
                gLineNo.pop()
# end of checkBinUnits


# try to parse term as something like L{pname} * Inv_L
def checkBinningTerm( t_in, pname, punit ):
    term = deepcopy(t_in)  # term may be modified below
    fact = ""
    pfx  = ""
    sfx  = ""
    oper = "*"
    if term.type == "*":
        # consider Inv_L * Inv_W * P{param}, group "Inv_L*Inv_W" into a single factor
        if term.e1.type == "*" and term.e2.type == "NAME":
            if term.e1.e1.type == "NAME" and term.e1.e2.type == "NAME":
                if term.e2.e1 in gParameters:
                    # Inv_L * Inv_W * P{param}
                    term.e1.type = "NAME"
                    term.e1.e1 = term.e1.e1.e1 + "*" + term.e1.e2.e1
                elif term.e1.e1.e1 in gParameters:
                    # P{param} * Inv_L * Inv_W
                    sparm = term.e1.e1.e1
                    term.e1.type = "NAME"
                    term.e1.e1 = term.e1.e2.e1 + "*" + term.e2.e1
                    term.e2.e1 = sparm
                elif term.e1.e2.e1 in gParameters:
                    # Inv_L * P{param} * Inv_W
                    sparm = term.e1.e2.e1
                    term.e1.type = "NAME"
                    term.e1.e1 = term.e1.e1.e1 + "*" + term.e2.e1
                    term.e2.e1 = sparm
        elif term.e2.type == "*" and term.e1.type == "NAME":
            # consider P{param} * (Inv_L * Inv_W) or permutations thereof
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
        # consider 1.0e-6 * LXL / L or WXW * 1.0e-6 / WGAA
        if term.e1.type == "*" and term.e2.type == "NAME" and term.e2.e1 in gParameters and gParameters[term.e2.e1].inst:
            if term.e1.e1.type == "NUMBER" and term.e1.e2.type == "NAME" and term.e1.e2.e1 in gParameters:
                # convert to 1.0e-6 / L * LXL
                term.type = "*"
                term.e1.type = "/"
                lparm = term.e2  # L or W
                term.e2 = term.e1.e2
                term.e1.e2 = lparm
            elif term.e1.e2.type == "NUMBER" and term.e1.e1.type == "NAME" and term.e1.e1.e1 in gParameters:
                # convert to 1.0e-6 / WGAA * WXW
                term.type = "*"
                term.e1.type = "/"
                lparm = term.e2  # L or W
                term.e2 = term.e1.e1
                term.e1.e1 = term.e1.e2
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
    if (term.type in ["*", "/"]) and (term.e1.type == "NAME" or term.e2.type == "NAME"):
        sname = ""
        oper = term.type
        if term.e1.type == "NAME" and term.e1.e1 in gParameters:
            sname = term.e1.e1
            fact = term.e2.asString()
        elif oper == "*" and term.e2.type == "NAME" and term.e2.e1 in gParameters:
            sname = term.e2.e1
            fact = term.e1.asString()
        elif (term.e1.type == "NUMBER" and term.e1.number == 0.0) or \
             (term.e2.type == "NUMBER" and term.e2.number == 0.0):
            # deliberately zero binning term
            pfx = "0.0"
        if sname != "":
            if sname in gParameters:
                spar = gParameters[sname]
                gFileName.append(spar.declare[0])
                gLineNo.append(spar.declare[1])
                if spar.defv:
                    if spar.defv.type == "NUMBER":
                        if spar.defv.number != 0:
                            warning("Non-zero default value (%e) for binning parameter '%s'"
                                    % (spar.defv.number, sname))
                    elif spar.defv.type == "NAME" and spar.defv.e1 in gParameters:
                        dpar = gParameters[spar.defv.e1]
                        if sname[0] != dpar.name[0] or (sname[0:2] == "P2" and dpar.name[0:2] != "P2"):
                            # WPAR default set from PPAR
                            warning("Unexpected default '%s' for binning parameter '%s'"
                                    % (dpar.name, sname))
                        while dpar.defv.type == "NAME" and dpar.defv.e1 in gParameters:
                            if dpar == gParameters[dpar.defv.e1]:
                                # default is itself - infinite loop, error reported elsewhere
                                dpar = False
                                break
                            dpar = gParameters[dpar.defv.e1]
                        if dpar and (dpar.defv.type != "NUMBER" or dpar.defv.number != 0):
                            warning("Non-zero default value (%s) for binning parameter '%s'"
                                    % (spar.defv.e1, sname))
                    else:
                        warning("Unexpected default value (type %s) for binning parameter '%s'"
                                % (spar.defv.type, sname))
                gFileName.pop()
                gLineNo.pop()
            if len(sname) > len(pname) and pname == sname[len(sname)-len(pname):]:
                # LPAR, WPAR, PPAR, P2PAR, etc.
                pfx = sname[0:len(sname)-len(pname)]
                checkBinUnits(spar, pfx, pname, punit)
            elif (pname[-1].lower() == 'o' or pname[-1] == '0') and \
                    pname[:-1] == sname[0:len(pname)-1]:
                # PARO and PARL, PARWL, etc.
                sfx = sname[len(pname)-1:]
                checkBinUnits(spar, sfx, pname, punit)
            elif pname[0:2].lower() == 'po' and pname[2:] == sname[len(sname)-len(pname)+2:]:
                # POPAR and LPAR, WPAR, etc.
                pfx = sname[0:2]
    return [fact, oper, pfx, sfx]
# end of checkBinningTerm



def checkBinning( vname, t_0, t_1, t_2, t_3, t_4, t_5, line ):
    if t_0.type == "NAME" and t_0.e1 in gParameters:
        pname = t_0.e1
        pname_l = pname.lower()
        plen = len(pname)
        punit = gParameters[pname].units
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
            elif (pname_l[-1] == 'o' or pname_l[-1] == '0') and \
                    pname_l[:-1] == vname_l[0:plen-1]:
                # VFB_i = VFBO + ...
                vsfx = vname[plen-1:]
                sfx0 = pname[-1]
        elif vlen < plen and (pname_l[-1] == 'o' or pname_l[-1] == '0') and \
                pname_l[:-1] == vname_l:
            # KUOWE = KUOWEO + ...
            vpfx = "NONE"
            sfx0 = pname[-1]
        elif pname_l[0:2] == 'po' and pname_l[2:] == vname_l[0:plen-2]:
            vsfx = vname[plen-2:]
            pfx0 = pname[0:2]
            if vsfx == "":
                vsfx = "NONE"
        [fact1, op1, pfx1, sfx1] = checkBinningTerm(t_1, pname, punit)
        [fact2, op2, pfx2, sfx2] = checkBinningTerm(t_2, pname, punit)
        [fact3, op3, pfx3, sfx3] = checkBinningTerm(t_3, pname, punit)
        if t_4:
            [fact4, op4, pfx4, sfx4] = checkBinningTerm(t_4, pname, punit)
        else:
            [fact4, op4, pfx4, sfx4] = ["", "", "", ""]
        if t_5:
            [fact5, op5, pfx5, sfx5] = checkBinningTerm(t_5, pname, punit)
        else:
            [fact5, op5, pfx5, sfx5] = ["", "", "", ""]
        warned = False
        if fact1 != "" and fact2 != "" and fact3 != "":
            if pfx1+sfx1 != "" and pfx2+sfx2 != "" and pfx3+sfx3 != "":
                if gBinning and vpfx == "" and vsfx == "":
                    warning("Did not find parameter name '%s' in binning variable name '%s'" % (pname, vname))
                    warned = True
                if t_4 and pfx4+sfx4 == "":
                    warning("Binning equation to set '%s' involves '%s'" % (vname, t_4.asString()))
                    warned = True
                if t_5 and pfx5+sfx5 == "":
                    warning("Binning equation to set '%s' involves '%s'" % (vname, t_5.asString()))
                    warned = True
            elif pfx1+sfx1 == "" and (pfx2+sfx2 != "" or pfx3+sfx3 != ""):
                warning("Binning equation to set '%s' involves '%s'" % (vname, t_1.asString()))
                warned = True
            elif pfx2+sfx2 == "" and (pfx1+sfx1 != "" or pfx3+sfx3 != ""):
                warning("Binning equation to set '%s' involves '%s'" % (vname, t_2.asString()))
                warned = True
            elif pfx3+sfx3 == "" and (pfx1+sfx1 != "" or pfx2+sfx2 != ""):
                warning("Binning equation to set '%s' involves '%s'" % (vname, t_3.asString()))
                warned = True
        new_pat = [vpfx+vsfx, pfx0+sfx0, fact1+op1, pfx1+sfx1, fact2+op2, pfx2+sfx2, fact3+op3, pfx3+sfx3,
                   fact4+op4, pfx4+sfx4, fact5+op5, pfx5+sfx5]
        if gBinning and (gVerbose or gDebug):
            t_3 = "0"
            if fact3:
                t_3 = pfx3 + "[par]" + sfx3 + " " + op3 + " " + fact3
            t_4 = "0"
            if fact4:
                t_4 = pfx4 + "[par]" + sfx4 + " " + op4 + " " + fact4
            t_5 = "0"
            if fact5:
                t_5 = pfx5 + "[par]" + sfx5 + " " + op5 + " " + fact5
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
                warning("Binning equation does not match that on line %d of %s"
                        % (gBinningPatterns[0][1], gBinningPatterns[0][0]))
                gBinningPatterns.append( [gFileName[-1], gLineNo[-1]] + new_pat )
                if gVerbose or gDebug:
                    do_print = True
        if do_print:
            print("Binning pattern: %s[par]%s = %s[par]%s + %s[par]%s %s %s + %s[par]%s %s %s + %s + %s + %s"
                  % (vpfx, vsfx, pfx0, sfx0, pfx1, sfx1, op1, fact1, pfx2, sfx2, op2, fact2, t_3, t_4, t_5))
            if gDebug:
                print("    derived from %s" % line)
# end of checkBinning


# for if, while, for
gLoopTypeForNextStmt = []
gLoopTypeForLastStmt = []
gConditionForNextStmt = []
gConditionForLastStmt = []
gCondBiasDForNextStmt = []
gCondBiasDForLastStmt = []


def parseOther( line ):
    global gStatementInCurrentBlock
    global gLastDisplayTask, gLastLineWasEvent

    if gLoopTypeForNextStmt:
        this_loop = gLoopTypeForNextStmt.pop()
    else:
        this_loop = ""
    last_loop = gLoopTypeForLastStmt
    if gConditionForNextStmt:
        this_cond = gConditionForNextStmt.pop()
    else:
        this_cond = Expression("NOTHING")
    if gConditionForLastStmt:
        last_cond = gConditionForLastStmt.pop()
    else:
        last_cond = Expression("NOTHING")
    if gCondBiasDForNextStmt:
        this_c_bd = gCondBiasDForNextStmt.pop()
    else:
        this_c_bd = 0
    if gCondBiasDForLastStmt:
        last_c_bd = gCondBiasDForLastStmt.pop()
    else:
        last_c_bd = 0
    analog_block = "None"

    parser = Parser(line)
    val = parser.lex()
    if val == ord(';') and (gLastLineWasEvent or this_cond.type != "NOTHING"):
        # null statement
        val = parser.lex()
    keywd = ""

    got_paren = False
    if val == ord('(') and len(gScopeList) > 0 and gScopeList[-1].startswith("CASE::"):
        got_paren = True
        val = parser.lex()

    if parser.isIdentifier():
        keywd = parser.getString()
        if keywd == "analog":
            gStatementInCurrentBlock = True
            analog_block = "analog"
            if gCheckMultFactors == 1:
                classifyContribs()
            val = parser.lex()
            keywd = parser.getString()
            if keywd == "initial":
                val = parser.lex()
                warning("Restrictions in 'analog initial' not checked")
                analog_block = "initial"

    if parser.isNumber() or (keywd == "default" and parser.peekToken() == ":"):
        if len(gScopeList) > 0 and gScopeList[-1].startswith("CASE::"):
            case_val = 0
            if parser.isNumber():
                # case 0: begin ...
                scope = gScopeList[-1]
                this_cond = Expression("==")
                this_cond.e1 = Expression("NAME")
                this_cond.e1.e1 = scope[6:]
                this_cond.e2 = Expression("NUMBER")
                case_val = parser.getNumber()
                this_cond.e2.number = case_val
            else:
                # default: (assume no hidden state if "default" is present)
                this_cond = Expression("NOT")
                this_cond.e1 = Expression("NUMBER")
                this_cond.e1.number = 0
            val = parser.lex()
            while val == ord(','):
                val = parser.lex()
                if parser.isNumber():
                    old_cond = this_cond
                    new_cond = Expression("==")
                    new_cond.e1 = Expression("NAME")
                    new_cond.e1.e1 = scope[6:]
                    new_cond.e2 = Expression("NUMBER")
                    case_val = parser.getNumber()
                    new_cond.e2.number = case_val
                    this_cond = Expression("||")
                    this_cond.e1 = old_cond
                    this_cond.e2 = new_cond
                    val = parser.lex()
                else:
                    error("Expected number after ',' in case value list")
            if got_paren:
                if val == ord(')'):
                    val = parser.lex()
                else:
                    error("Missing ')' after case value %d" % case_val)
            if val == ord(':'):
                val = parser.lex()
            else:
                error("Expected ':' after case value %d" % case_val)

    # if, else, for, while, repeat, begin, end
    if parser.isIdentifier():
        keywd = parser.getString()

        if keywd in ["end", "else", "if", "for", "repeat", "while", "begin", "case", "generate"]:
            gStatementInCurrentBlock = True

        if keywd == "end":
            if len(gScopeList) > 0:
                scope = gScopeList.pop()
                col = scope.find("::")
                if col > 0:
                    error("Found 'end' instead of 'end%s'" % scope[0:col].lower())
            else:
                error("Found 'end' without corresponding 'begin'")
            if len(gLoopTypes) > 0:
                last_loop = gLoopTypes.pop()
            if len(gConditions) > 0:
                last_cond = gConditions.pop()
            if len(gCondBiasDep) > 0:
                last_c_bd = gCondBiasDep.pop()
            if len(gAnalogBlock) > 0:
                gAnalogBlock.pop()
            if len(gAnalogBlock) > 0:
                analog_block = gAnalogBlock[-1]
            else:
                analog_block = "None"
            val = parser.lex()
            if last_loop[0] in ["for", "repeat", "while"] and last_cond.type != "NOTHING":
                registerLoopVariableUse(last_cond, last_loop)
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
            else:
                # this "else" corresponds to an "if" that's always true, eg if(1)
                this_cond = Expression("NUMBER")
                this_cond.number = 0
            val = parser.lex()
            if parser.isIdentifier():
                keywd = parser.getString()
            elif val == 0:
                # statement on next line
                gConditionForNextStmt.append(this_cond)
                gCondBiasDForNextStmt.append(this_c_bd)

        if keywd == "if":
            new_cond = parser.getExpression(True)
            if new_cond.type == "NUMBER" and new_cond.number != 0:
                # if (1)
                new_cond = Expression("NOTHING")
                deps = []
            elif new_cond.type == "NOTHING":
                error("Missing condition for 'if'")
                deps = []
            else:
                deps = new_cond.getDependencies(False, False)
            [bias_dep, biases, ddt] = checkDependencies(deps, "If condition depends on", 0, "if", True, True)
            if this_cond.type == "NOTHING":
                this_cond = new_cond
            elif new_cond.type == "NOTHING":
                pass
            else:
                old_cond = this_cond
                this_cond = Expression("&&")
                this_cond.e1 = old_cond
                this_cond.e2 = new_cond
            if bias_dep > 1 and checkCondForBiasDepEquality(new_cond):
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
                    rest = "if " + parser.getRestOfLine()
                    gConditionForNextStmt.append(this_cond)
                    gCondBiasDForNextStmt.append(this_c_bd)
                    hold_cond = this_cond
                    hold_c_bd = this_c_bd
                    parseOther(rest)
                    this_cond = hold_cond
                    this_c_bd = hold_c_bd
                    val = parser.lex()
                    # if (A) if (B) ... becomes if (A&&B)
                    # TODO: doesn't handle subsequent else
            if parser.isIdentifier():
                keywd = parser.getString()
            elif val == ord(';'):
                # null statement
                val = parser.lex()
            elif val == 0:
                # statement on next line
                gConditionForNextStmt.append(this_cond)
                gCondBiasDForNextStmt.append(this_c_bd)

        if keywd == "for":
            # for (init; cond; incr)
            this_loop = keywd
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
                        markVariableAsSet(varname, None, [], False, this_cond, [], 0, True, False, "", 0)
                    if init_expr.type == "NUMBER":
                        pass
                    elif init_expr.type == "NAME":
                        deps = [init_expr.e1]
                        checkDependencies(deps, "For loop initialization depends on", 0, False, True, True)
                    else:
                        # deps = init_expr.getDependencies(False, False)
                        bad_init = True
                        if gDebug:  # pragma: no cover
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
                [bias_dep, biases, ddt] = checkDependencies(deps, "For loop condition depends on", 0, "for", True, True)
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
                if varname not in ["", incr_var]:
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
                if incr_rhs.type not in ["+", "-"]:
                    good_incr = False
            elif oper == "++":
                oper = parser.getOper()
                # error printed by peekOper()
            elif oper == "--":
                oper = parser.getOper()
                # error printed by peekOper()
            elif oper in ["+=", "-="]:
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
                gLoopTypeForNextStmt.append(this_loop)
                this_loop = ""
                gConditionForNextStmt.append(this_cond)
                gCondBiasDForNextStmt.append(this_c_bd)

        if keywd == "repeat":
            this_loop = keywd
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
                gLoopTypeForNextStmt.append(this_loop)
                this_loop = ""
                gConditionForNextStmt.append(this_cond)
                gCondBiasDForNextStmt.append(this_c_bd)

        if keywd == "while":
            this_loop = keywd
            this_cond = parser.getExpression(True)
            deps = this_cond.getDependencies(False, False)
            [bias_dep, biases, ddt] = checkDependencies(deps, "While condition depends on", 0, "while", True, True)
            this_c_bd = bias_dep
            val = parser.lex()
            if parser.isIdentifier():
                keywd = parser.getString()
                if keywd != "begin":
                    error("expected 'begin' after 'while'")
            elif val == 0:
                # statement on next line
                gLoopTypeForNextStmt.append(this_loop)
                this_loop = ""
                gConditionForNextStmt.append(this_cond)
                gCondBiasDForNextStmt.append(this_c_bd)

        if keywd == "generate":
            # generate i (start, stop [, incr])
            val = parser.lex()
            if parser.isIdentifier():
                warning("Archaic 'generate' statement")
                name = parser.getString()
                scope = getCurrentScope()
                vname = scope + name
                if vname in gVariables:
                    notice("generate statement redeclares variable '%s'" % name)
                    markVariableAsSet(vname, None, [], False, this_cond, [], 0, False, False, "", 0)
                else:
                    var = Variable(vname, "integer", False, gFileName[-1], gLineNo[-1])
                    var.assign = 0  # don't report in summary
                    gVariables[vname] = var
                val = parser.lex()
            else:
                error("Missing variable name for generate")
            valid = True
            if val == ord('('):
                start = parser.getExpression()
                if start.type == "NOTHING":
                    error("Missing (null) start expression in generate statement")
                val = parser.lex()
                if val != ord(','):
                    valid = False
                if val != ord(')'):
                    stop = parser.getExpression()
                    if stop.type == "NOTHING":
                        error("Missing (null) stop expression in generate statement")
                    val = parser.lex()
                if val == ord(','):
                    incr = parser.getExpression()
                    if incr.type == "NOTHING":
                        error("Missing (null) increment in generate statement")
                    val = parser.lex()
                if val == ord(')'):
                    val = parser.lex()
                    if parser.isIdentifier():
                        keywd = parser.getString()
                else:  # pragma: no cover
                    valid = False
            else:
                valid = False
            if not valid:
                error("Expected (start, stop) for generate")

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
                    [bias_dep, biases, ddt] = checkDependencies(deps, "If condition depends on", 0, "if", True, True)
                    if this_cond.type == "NOTHING":
                        this_cond = new_cond
                    else:
                        old_cond = this_cond
                        this_cond = Expression("&&")
                        this_cond.e1 = old_cond
                        this_cond.e2 = new_cond
                    if bias_dep > 1 and checkCondForBiasDepEquality(new_cond):
                        bias_dep = 3
                        warning("Bias-dependent '%s': if %s" % (new_cond.type, new_cond.asString()))
                    if bias_dep == 2:
                        bias_dep = 1
                    this_c_bd = bias_dep
                    val = parser.lex()
            elif val != 0:
                error("Unexpected token %d after 'begin'" % val)
            gScopeList.append(bname)
            gLoopTypes.append([this_loop, gFileName[-1], gLineNo[-1]])
            gConditions.append(this_cond)
            gCondBiasDep.append(this_c_bd)
            this_loop = ""
            this_cond = Expression("NOTHING")
            this_c_bd = 0
            gAnalogBlock.append(analog_block)

        if keywd == "case":
            var = parser.getExpression()
            if var.type == "NAME":
                varname = var.e1
                [bias_dep, biases, ddt] = checkDependencies([varname], "Case statement depends on", 0, "case", True, True)
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
                    if bus_range:
                        if node.is_bus:
                            if node.msb != bus_range[0] or node.lsb != bus_range[1]:
                                error("Port '%s' has different bus ranges in port and discipline declarations" % pname)
                        elif node.direction != "":
                            error("Discipline declaration for scalar port '%s' has bus range [%d:%d]"
                                  % (pname, bus_range[0], bus_range[1]))
                        # else have not seen direction declaration
                    elif node.is_bus:
                        error("Discipline declaration for port '%s' should also have bus range [%d:%d]"
                              % (pname, node.msb, node.lsb))
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
                            if bus_range:
                                node.is_bus = True
                                node.msb = bus_range[0]
                                node.lsb = bus_range[1]
                            gNodenames[pname] = node
                val = parser.lex()
                if val == ord(','):
                    val = parser.lex()
                elif val == ord('='):
                    error("Net initializers not supported")
                    parser.getExpression()
                    val = parser.lex()
                elif parser.isIdentifier():
                    error("Missing ',' in list of nodes")
            if val == ord(';'):
                val = parser.lex()
                while val == ord(';'):
                    error("Extra ';' after list of nodes")
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
                    else:  # pragma: no cover
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
                    parser.getExpression()
                    val = parser.lex()
                elif parser.isIdentifier():
                    error("Missing ',' in list of nodes")
            if val == ord(';'):
                val = parser.lex()
                while val == ord(';'):
                    error("Extra ';' after list of nodes")
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
                while val == ord(';'):
                    error("Extra ';' after branch declaration")
                    val = parser.lex()
            else:
                error("Missing ';' after branch declaration")

        elif keywd in ["cmos", "rcmos", "bufif0", "bufif1", "notif0", "notif1",
                       "nmos", "pmos", "rnmos", "rpmos",
                       "and", "nand", "or", "nor", "xor", "xnor", "buf", "not",
                       "tranif0", "tranif1", "rtranif0", "rtranif1",
                       "tran", "rtran", "pulldown", "pullup"]:  # pragma: no cover
            if gCompactModel:
                error("Gate instantiation '%s' not expected in compact model" % keywd)

        elif keywd in ["resistor", "capacitor", "inductor", "tline"]:  # pragma: no cover
            if analog_block == "None":
                warning("Spice instantiation '%s' not recommended" % keywd)
            else:
                error("Spice instantiation '%s' inside analog block" % keywd)
            val = 0

        elif keywd in ["always", "initial"]:  # pragma: no cover
            error("Digital '%s' blocks not supported in Verilog-A" % keywd)
            val = parser.lex()

        elif keywd in ["task", "endtask"]:  # pragma: no cover
            if gCompactModel:
                error("Tasks not supported in compact models")

        elif keywd in ["supply0", "supply1", "tri", "triand", "trior", "tri0", "tri1",
                       "uwire", "wire", "wand", "wor", "trireg", "wreal",
                       "realtime", "reg", "time", "specparam"]:  # pragma: no cover
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
                    [bias_dep, biases, ddt] = checkDependencies(deps, "Variable %s depends on"
                                                                % varname, this_c_bd, False, False, True)
                    if ddt:
                        ddt = rhs.ddtCheck()
                        if ddt > 1:
                            error("Assignment to '%s' is nonlinear in ddt()" % varname)
                    markVariableAsSet(varname, rhs, deps, bias_dep, this_cond, biases, ddt, False, False, "", 0)
                    if not bias_dep and rhs.type == "+" and (rhs.e1.type == "+" or rhs.e2.type == "+"):
                        # possible binning equation P_i = P + LP*Inv_L + WP*Inv_W + PP*Inv_P
                        t3 = rhs.e2
                        if rhs.e1.type == "+" and (t3.type in ["*", "/"]):
                            t2 = rhs.e1.e2
                            if rhs.e1.e1.type == "+" and (t2.type in ["*", "/"]):
                                t1 = rhs.e1.e1.e2
                                t0 = rhs.e1.e1.e1
                                if t0.type == "NAME" and (t1.type in ["*", "/"]):
                                    checkBinning(varname, t0, t1, t2, t3, 0, 0, line)
                                elif t0.type == "+" and (t1.type in ["*", "/"]):
                                    if t0.e1.type == "NAME":
                                        # 4 binning parameters
                                        t4 = t3
                                        t3 = t2
                                        t2 = t1
                                        t1 = t0.e2
                                        t0 = t0.e1
                                        checkBinning(varname, t0, t1, t2, t3, t4, 0, line)
                                    else:
                                        # 5 binning parameters
                                        t5 = t3
                                        t4 = t2
                                        t3 = t1
                                        t2 = t0.e2
                                        if t0.e1.type == "+" and (t2.type in ["*", "/"]):
                                            t1 = t0.e1.e2
                                            t0 = t0.e1.e1
                                            if t0.type == "NAME" and (t1.type in ["*", "/"]):
                                                checkBinning(varname, t0, t1, t2, t3, t4, t5, line)
                    val = parser.lex()
                    if val == ord(';'):
                        val = parser.lex()
                        if val == ord(';'):
                            error("Extra ';' after assignment")
                            val = parser.lex()
                    else:
                        error("Missing ';' after assignment")
                    if parser.isIdentifier():
                        keywd = parser.getString()
                        if keywd == "else":
                            # if (A) X=0; else X=1;
                            rest = keywd + " " + parser.getRestOfLine()
                            keywd = ""
                            gConditionForLastStmt.append(this_cond)
                            gCondBiasDForLastStmt.append(this_c_bd)
                            parseOther(rest)
                            this_cond = gConditionForLastStmt.pop()
                            this_c_bd = gCondBiasDForLastStmt.pop()
                        elif keywd == "if":
                            parser.ungetChars("if")
                            if gStyle:
                                style("Multiple statements on a single line")
                        elif keywd == "end" and parser.peekRestOfLine().strip() == "":
                            pass
                        elif gStyle and not all_zero_assigns:
                            style("Multiple assignments on a single line")
                        if parser.peekChar() == '[':
                            parser.checkBusIndex(keywd)
                        if keywd != "if":
                            val = parser.lex()
                val = parser.lex()
                if val == 0 and keywd != "":
                    if keywd == "end":
                        gConditionForLastStmt.append(this_cond)
                        gCondBiasDForLastStmt.append(this_c_bd)
                        parseOther(keywd)
                        this_cond = gConditionForLastStmt.pop()
                        this_c_bd = gCondBiasDForLastStmt.pop()
                    else:  # pragma: no cover
                        fatal("Unexpected '%s' after assignment" % keywd)

            elif val == ord('(') or (val == ord(';') and keywd in ["$fatal", "$error", "$warning", "$finish", "$stop"]):
                if keywd in gAccessFuncs:
                    gStatementInCurrentBlock = True
                    acc = keywd
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
                            gBranches[nname1].used = True
                            if nname2 != "":
                                if nname2 in gBranches:
                                    error("Invalid potential or flow access %s(branch,branch)" % acc)
                                else:
                                    error("Invalid potential or flow access %s(branch,node)" % acc)
                        else:
                            if nname1 in gPortnames:
                                dname = gPortnames[nname1].discipline
                                gPortnames[nname1].used = True
                            elif nname1 in gNodenames:
                                dname = gNodenames[nname1].discipline
                                gNodenames[nname1].used = True
                            else:
                                error("Identifier '%s' is not a port, node, or branch in potential or flow access"
                                      % nname1)
                            if nname2 != "":
                                dname2 = ""
                                if nname2 in gBranches:
                                    dname2 = gBranches[nname2].discipline
                                    error("Invalid potential or flow access %s(node,branch)" % acc)
                                elif nname2 in gPortnames:
                                    dname2 = gPortnames[nname2].discipline
                                    gPortnames[nname2].used = True
                                elif nname2 in gNodenames:
                                    dname2 = gNodenames[nname2].discipline
                                    gNodenames[nname2].used = True
                                else:
                                    error("Identifier '%s' is not a port, node, or branch in potential or flow access"
                                          % nname2)
                                if dname2 not in ["", dname]:
                                    error("Nodes '%s' and '%s' belong to different disciplines" % (nname1, nname2))
                    if dname in gDisciplines:
                        disc = gDisciplines[dname]
                        fname = disc.flow
                        pname = disc.potential
                        if fname in gNatures and acc == gNatures[fname].access:
                            registerContrib("flow", nname1, nname2, this_cond, this_c_bd)
                        elif pname in gNatures and acc == gNatures[pname].access:
                            registerContrib("potential", nname1, nname2, this_cond, this_c_bd)
                        else:
                            error("Incorrect access function '%s' for %s '%s'" % (acc, dname, nname1))
                    if val == ord('<') and parser.peekChar() == '+':
                        parser.getChar()
                        rhs = parser.getExpression()
                        deps = rhs.getDependencies(False, False)
                        acc_ref = acc + "(" + nname1
                        if nname2 != "":
                            acc_ref = acc_ref + "," + nname2
                        acc_ref = acc_ref + ")"
                        if acc_ref in deps:
                            error("Branch contribution to '%s' depends on '%s'" % (acc_ref, acc_ref))
                        [bias_dep, biases, ddt] = checkDependencies(deps, "Branch contribution depends on",
                                                                    this_c_bd, False, False, True)
                        if ddt:
                            ddt = rhs.ddtCheck()
                            if ddt > 1:
                                error("Contribution to '%s' is nonlinear in ddt()" % acc_ref)
                        if gCheckMultFactors:
                            if acc == "I" and \
                                    ((nname1 in gPortnames and gPortnames[nname1].mult) or \
                                     (nname1 in gNodenames and gNodenames[nname1].mult) or \
                                     (nname1 in gBranches and gBranches[nname1].mult)) :
                                checkMultFactors(rhs)
                            else:
                                for mult in gMultFactors:
                                    if mult in deps:
                                        error("Unexpected '%s' in branch contribution" % mult)
                        # check dependencies again, but ignore bias_dep in _noise calls
                        deps = rhs.getDependencies(False, True)
                        [bias_dep, biases, ddt] = checkDependencies(deps, "Branch contribution depends on",
                                                                    this_c_bd, False, False, False)
                        if bias_dep == 4:
                            error("Branch contribution depends on quantity obtained from ddx (requires second derivatives)")
                        elif bias_dep >= 3:
                            error("Branch contribution depends on quantity with bad derivative (from bias-dependent == or !=)")
                    elif val == ord(':'):
                        indirect = parser.getExpression()
                        deps = indirect.getDependencies(False, False)
                        [bias_dep, biases, ddt] = checkDependencies(deps, "Indirect branch contribution depends on",
                                                                    this_c_bd, False, False, True)
                        if ddt:  # pragma: no cover
                            ddt = indirect.ddtCheck()
                            if ddt > 1:
                                error("Indirect branch contribution is nonlinear in ddt()")
                        invalid = False
                        if indirect.type == "==" and indirect.e1.type == "FUNCCALL":
                            fname = indirect.e1.e1
                            if not (fname in ["ddt", "idt"] or fname in gAccessFuncs):
                                invalid = True
                        else:
                            invalid = True
                        if invalid:
                            error("Invalid indirect branch contribution")
                    val = parser.lex()
                    if val == ord(';'):
                        val = parser.lex()
                        if gStyle and parser.isIdentifier():
                            keywd = parser.getString()
                            if keywd in gAccessFuncs:
                                style("Multiple contributions on a single line")
                        elif val == ord(';'):
                            error("Extra ';' after contribution")
                            val = parser.lex()
                    else:
                        error("Missing ';' after contribution")

                elif keywd in ["$strobe", "$display", "$write", "$monitor", "$debug", "$bound_step", "$discontinuity",
                               "$fstrobe", "$fdisplay", "$fwrite", "$fmonitor",
                               "$fatal", "$error", "$warning", "$info", "$finish", "$stop"]:
                    gStatementInCurrentBlock = True
                    if val == ord('('):
                        args = parser.parseArgList()
                    else:
                        args = []
                    if keywd in ["$finish", "$stop"]:
                        arg = ""
                        if len(gLastDisplayTask) > 2 and gFileName[-1] == gLastDisplayTask[1] \
                                                     and gLineNo[-1] >= gLastDisplayTask[2]:
                            if len(args) > 0:
                                arg = "(" + args[0].asString() + ")"
                            notice("$error(msg) preferred over %s(msg); %s%s;" % (gLastDisplayTask[0], keywd, arg))
                        else:
                            if val == ord('('):
                                arg = ()
                            notice("$error() preferred over %s%s" % (keywd, arg))
                    elif keywd in ["$debug", "$monitor"]:
                        warning("%s() may degrade performance" % keywd)
                    elif keywd in ["$strobe", "$display", "$write"]:
                        gLastDisplayTask = [keywd, gFileName[-1], gLineNo[-1]]

                    if val == ord('('):
                        deps = []
                        for arg in args:
                            deps += arg.getDependencies(False, False)
                        bias_dep = checkDependencies(deps, "Task %s references" % keywd, this_c_bd, False, False, True)[0]
                        if bias_dep and keywd in ["$strobe", "$display", "$write", "$warning", "$info",
                                                  "$fstrobe", "$fdisplay", "$fwrite", "$fmonitor"]:
                            warning("bias-dependent %s() may degrade performance" % keywd)
                        val = parser.lex()
                        if val == ord(')'):
                            val = parser.lex()
                        else:
                            error("Missing ')' in call of task %s" % keywd)
                    if keywd in ["$finish", "$stop", "$discontinuity"]:
                        if len(args) == 1:
                            if args[0].type != "NUMBER":
                                error("Argument to %s must be a number" % keywd)
                        elif len(args) > 1:
                            error("Task %s takes at most one argument" % keywd)
                    elif keywd == "$fatal":
                        if len(args) >= 1:
                            if args[0].type != "NUMBER":
                                error("Argument to %s must be a number" % keywd)
                    elif keywd == "$bound_step":
                        if gCompactModel:
                            warning("Task '%s' should not be used in a compact model" % keywd)
                    else:
                        strpos = 0
                        posstr = "First"
                        if keywd in ["$fstrobe", "$fdisplay", "$fwrite", "$fmonitor"]:
                            badarg = False
                            if args[0].type == "NAME":
                                vn = findVariableInScope(args[0].e1)
                                if vn in gVariables:
                                    var = gVariables[vn]
                                    if var.type != "integer":
                                        badarg = True
                                else:
                                    badarg = True
                            else:
                                badarg = True
                            if badarg:
                                error("First argument to %s must be a file descriptor (integer variable)" % keywd)
                            strpos = 1
                            posstr = "Second"
                        if len(args) == 0:
                            error("Expected at least 1 argument for %s" % keywd)
                        elif len(args) <= strpos:
                            error("Expected at least %d arguments for %s" % (strpos+1,keywd))
                        elif args[strpos].type != "STRING":
                            error("%s argument to %s must be a string" % (posstr, keywd))
                    if val == ord(';'):
                        val = parser.lex()
                        if val == ord(';'):
                            error("Extra ';' after task %s" % keywd)
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
                    nextch = parser.lex()
                    if nextch in [ord('+'), ord('-'), ord('*'), ord('/'), ord('=')]:
                        oper = chr(val) + chr(nextch)
                        if oper in ["++", "--", "+=", "-=", "*=", "/="]:
                            error("Operator '%s' not valid in Verilog-A" % oper)
                            reported = True
                rest = parser.getRestOfLine()
                if not reported:
                    error("Unexpected '%s'" % line)
                    val = 0

    elif val == ord('@'):
        if gCompactModel:
            error("Events (@) should not be used in compact models")
        gLastLineWasEvent = True
        #  @ (initial_step or initial_step("static","pss") or initial_step("static","pdisto")) begin
        val = parser.lex()
        args = []
        keywd = ""
        if val == ord('('):
            val = parser.lex()
            while parser.isIdentifier() and keywd == "":
                arg = Expression("NAME")
                arg.e1 = parser.getString()
                args.append(arg)
                val = parser.lex()
                if val == ord('('):
                    arg.type = "FUNCCALL"
                    arg.args = parser.parseArgList()
                    for farg in arg.args:
                        deps = farg.getDependencies(False, False)
                        checkDependencies(deps, "%s event references" % arg.e1, this_c_bd, False, False, True)
                    if parser.peekChar() == ')':
                        parser.getChar()
                        val = parser.lex()
                    else:
                        error("Missing ')' in event")
                if parser.isIdentifier():
                    keywd = parser.getString()
                    if keywd == "or":
                        keywd = ""
                        val = parser.lex()
                        if not parser.isIdentifier():
                            error("Missing event after 'or'")
            if val == ord(')'):
                parser.lex()
            else:
                error("Missing ')' in event")
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
        if rest.startswith("elseif"):
            pass # error printed in handleCompilerDirectives
        else:
            error("Unhandled macro `%s" % rest)
    elif val == ord('#'):
        rest = parser.getRestOfLine()
        if gCompactModel:
            error("Delays (#) should not be used in compact models")
    elif val == 0:
        pass  # end of line
    elif parser.isNumber():
        num = parser.getNumber()
        error("Unexpected number %d" % num)
    elif val == ord(';') and len(gScopeList) > 0 and gScopeList[-1].startswith("CASE::"):
        pass  # null case statement
    else:
        error("Unexpected %s" % formatChar(val))

    gLoopTypeForLastStmt.append(this_loop)
    if this_loop in ["for", "repeat", "while"] and this_cond.type != "NOTHING":
        registerLoopVariableUse(this_cond, False)
    gConditionForLastStmt.append(this_cond)
    gCondBiasDForLastStmt.append(this_c_bd)

    if val != 0:
        rest = ""
        if parser.isIdentifier():
            rest = parser.getString()
        rest += parser.getRestOfLine()
        parseOther(rest)
# end of parseOther


def getCurrentScope():
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
    conds = ""
    bias_dep = c_bd_in
    if cond_in.type != "NOTHING":
        conds = cond_in.asString()
    for cond in gConditionForNextStmt:
        if cond.type != "NOTHING":
            if conds != "":
                conds += "&&"
            conds += cond.asString()
    for cond in gConditions:
        if cond.type != "NOTHING":
            if conds != "":
                conds += "&&"
            conds += cond.asString()
    for c_bd in gCondBiasDForNextStmt:
        if c_bd > bias_dep:
            bias_dep = c_bd
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
    global gMissingConstantsFile

    indent = 0
    if gPreProcess:
        indent = getIndent(line)
    line = line.strip()
    if line == "":
        return line

    # compiler directives
    did_compiler_directive = False

    if line.startswith("`ifdef") or line.startswith("`ifndef") or line.startswith("`elsif"):
        if line.startswith("`ifndef"):
            start = 7
            false_status = "TRUE"
            true_status = "FALSE"
        else:
            start = 6
            true_status = "TRUE"
            false_status = "FALSE"
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
        if line.startswith("`elsif"):
            if len(gIfDefStatus) > 0:
                if gIfDefStatus[-1] == "FALSE":
                    if token in gMacros:
                        gMacros[token].used = True
                        gIfDefStatus[-1] = "TRUE:elsif"
                    else:
                        gIfDefStatus[-1] = "FALSE"
                else: # TRUE, TRUE:elsif, FALSE:elsif means was TRUE at some point
                    gIfDefStatus[-1] = "FALSE:elsif"
            else:
                error("Found `elsif without `ifdef")
        else:
            if token in gMacros:
                gMacros[token].used = True
                gIfDefStatus.append(true_status)
            else:
                gIfDefStatus.append(false_status)
        did_compiler_directive = True
    elif line.startswith("`elseif"):
        if "elseif" not in gMacros:
            error("Undefined macro `elseif (was `elsif intended?)")
    elif line.startswith("`else"):
        if len(gIfDefStatus) > 0:
            if gIfDefStatus[-1] == "FALSE":
                gIfDefStatus[-1] = "TRUE"
            else: # TRUE, TRUE:elsif, FALSE:elsif means was TRUE at some point
                gIfDefStatus[-1] = "FALSE"
        else:
            error("Found `else without `ifdef")
        did_compiler_directive = True
    elif line.startswith("`endif"):
        start = 6
        if len(line) > start and not line[start].isspace():
            error("Missing space after %s" % line[0:start])
        if len(gIfDefStatus) > 0:
            gIfDefStatus.pop()
        else:
            error("Unmatched `endif")
        did_compiler_directive = True

    # check ifdef status
    ifdef_status = True
    for stat in gIfDefStatus:
        if stat in ["FALSE", "FALSE:elsif"]:
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
                        lines = parseFile( f_fpn )
                        if gPreProcess:
                            line = "//" + line + "\n"
                            for li in lines:
                                line = line + li + "\n"
                    else:
                        for fpath in gIncDir:
                            f_fpn = os.path.join(fpath, fname)
                            if os.path.exists(f_fpn) and os.path.isfile(f_fpn):
                                found = True
                                parseFile( f_fpn )
                                break
                    if not found and fname in ["discipline.h", "disciplines.h", "disciplines.vams",
                                               "constants.h", "constants.vams"]:
                        # check directory where this script resides
                        fpath = sys.path[0]
                        f_fpn = os.path.join(fpath, fname)
                        if os.path.exists(f_fpn) and os.path.isfile(f_fpn):
                            found = True
                            parseFile( f_fpn )
                    if not found and not gPreProcess:
                        if fname in ["discipline.h", "disciplines.h", "disciplines.vams"]:
                            loc = "https://accellera.org/images/downloads/standards/v-ams/disciplines_2-4.vams"
                            error("File not found: %s\n  Download from %s" % (line, loc))
                            addBasicDisciplines()
                        elif fname in ["constants.h", "constants.vams"]:
                            loc = "https://accellera.org/images/downloads/standards/v-ams/constants_2-4.vams"
                            error("File not found: %s\n  Download from %s" % (line, loc))
                            gMissingConstantsFile = fname
                        else:
                            error("File not found: %s" % line)
                if len(rest) > 0:
                    if rest.startswith("`include"):
                        error("Multiple `include directives on one line")
                        rest = rest[8:].strip()
                    else:
                        error("Unexpected characters after `include: %s" % rest)
                        rest = ""
            did_compiler_directive = not gPreProcess

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
        line = ""  # done with this line
    elif gPreProcess:
        if indent < 0: # tab
            indent = gSpcPerInd
        line = (indent * " ") + line

    if line.find("`") >= 0:
        line = expandMacro( line )

    return line
# end of handleCompilerDirectives


def preProcessNaturesAndContribs( lines ):
    # preprocess so idt_natures are defined
    for line in lines:
        gLineNo[-1] += 1
        if line.startswith("nature") and len(gScopeList) == 0:
            parseNatureDecl(line, False)
        elif line.find("<+") > 0:
            parts = line.split("<+")
            contrib = Contribution(parts[0].strip())
            key = len(gContribs)
            gContribs[key] = contrib


def preProcessLines( lines ):
    processed = []
    for line in lines:
        gLineNo[-1] += 1
        line = handleCompilerDirectives(line)
        if line == "":
            continue
        parts = line.split("\n")
        for part in parts:
            processed.append(part)
    return processed


def parseLines( lines ):
    global gStatementInCurrentBlock
    global gSubline
    prev_line = ""
    prev_attribs = []
    in_attribute = ""

    j = 0
    while j < len(lines):
        gLineNo[-1] += 1
        gSubline = 0
        line = lines[j]
        j += 1

        line = handleCompilerDirectives(line)
        num_sublines = len(line.split("\n"))
        line = line.strip()
        if line == "":
            continue

        # macros may have multiple lines, want to parse them separately
        parts = line.split("\n")
        i = 0
        while i < len(parts):
            part = parts[i]
            i += 1
            if num_sublines > 1:
                gSubline += 1
            part = part.strip()
            if part == "":
                continue
            retarray = getAttributes(part, in_attribute)
            part = retarray[0]
            in_attribute = retarray[1]
            attribs = prev_attribs + retarray[2:]
            while part == "":
                if i < len(parts):
                    part = parts[i]
                    i += 1
                    retarray = getAttributes(part, in_attribute)
                    part = retarray[0]
                    in_attribute = retarray[1]
                    attribs += retarray[2:]
                else:
                    prev_attribs = attribs
                    break
            if part == "" or in_attribute != "":
                prev_line = prev_line + part
                continue
            if part[-1] == ';':
                prev_attribs = []
            part = prev_line + part
            prev_line = ""

            # natures and disciplines
            #
            if part.startswith("nature"):
                if len(gScopeList) == 0:
                    parseNatureDecl(part, True)
                else:  # pragma: no cover
                    error("Natures should be defined at top level")
            elif len(gScopeList) == 1 and gScopeList[0].startswith("NATURE::"):
                parseNatureLine(part)
            elif part.startswith("discipline"):
                if len(gScopeList) == 0:
                    parseDisciplineDecl(part)
                else:  # pragma: no cover
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
                    gScopeList.append("MODULE::")
                    parseModuleDecl(part, attribs)
                    gStatementInCurrentBlock = False
                elif len(gScopeList) == 1:  # pragma: no cover
                    if gScopeList[-1] == "MODULE::":
                        fatal("Nested module declaration")
                    else:
                        fatal("Module declared in unexpected context")
                else:  # pragma: no cover
                    scope = getCurrentScope()
                    fatal("Unexpected 'module' in scope %s" % scope)
            elif part.startswith("endmodule"):
                verifyEndScope("module")
                checkPorts()
                checkParmsAndVars()
                checkInternalNodes()

            # functions
            #
            elif part.startswith("analog function") or part.startswith("function"):
                parseFunction(part)
                checkDeclarationContext("Function")
                gStatementInCurrentBlock = False
            elif part.startswith("endfunction"):
                gStatementInCurrentBlock = False
                if len(gScopeList) > 0 and gScopeList[-1].startswith("FUNCTION::"):
                    if gCurrentFunc:
                        fname = gCurrentFunc.name
                        fscope = gScopeList[-1]
                        if fscope == "FUNCTION::" + fname:
                            fscope += "."
                        else:  # pragma: no cover
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
                                        warning("In function %s, variable '%s' was never assigned a value"
                                                % (fname, vn))
                                    else:
                                        warning("In function %s, variable '%s' was never set and never used"
                                                % (fname, vn))
                                    gLineNo.pop()
                                    gFileName.pop()
                                elif not var.used and not var.assign == 0:
                                    warning("In function %s, variable '%s' was never used" % (fname, vn))
                        if not func_assign:
                            if len(gCurrentFunc.outputs) > 0:
                                warning("Function '%s' is not assigned a return value" % gCurrentFunc.name)
                            else:
                                error("Function '%s' is not assigned a return value" % gCurrentFunc.name)
                    else:  # pragma: no cover
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
                     part.startswith("inout") or part.startswith("input") or part.startswith("output") or \
                    (part.startswith("real") and part[4].isspace()) or \
                    (part.startswith("integer") and part[7].isspace()) or \
                    (part.startswith("genvar") and part[6].isspace()) or \
                    (part.startswith("string") and part[6].isspace()):
                is_port_dir = False
                if part.startswith("inout") or part.startswith("input") or part.startswith("output"):
                    is_port_dir = True
                if part[-1] != ';' and i >= len(parts):
                    tmpline = part
                    tmpj = j
                    while tmpline[-1] != ';' and tmpj < len(lines):
                        nextln = lines[tmpj]
                        nextln = handleCompilerDirectives(nextln)
                        if nextln != "":
                            nparts = nextln.split()
                            if nparts[0] in gVAMSkeywords and nparts[0] not in gMathFunctions:
                                break
                            if is_port_dir and (nparts[0] in ["inout", "input", "output"] or nparts[0] in gDisciplines):
                                # missing ; so don't join this part
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
                elif part.startswith("inout") or part.startswith("input") or part.startswith("output"):
                    parsePortDirection(part)
                else:
                    parseVariableDecl(part, attribs)

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
                           (part.find("if") < 0 and part.find("else") < 0 and
                            part.find("begin") < 0 and part.find("end") < 0 and
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
                        nextln = lines[tmpj].strip()
                        if num_parens > 0 and nextln == "end":
                            break
                        if num_parens > 0 or tmpline[-1] in "+-*/%><!&|=~^?" or \
                                (nextln != "" and nextln[0] in "+-*/%><!&|=~^?"):
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
                        elif num_parens < 0:  # pragma: no cover
                            fatal("Unexpected ')'")

                    if tmpline[-1] != ';':
                        parser = Parser(tmpline)
                        val = parser.lex()
                        if parser.isIdentifier():
                            do_it = False
                            keywd = parser.getString()
                            val = parser.lex()
                            if val == ord('='):
                                do_it = True  # assignment -> need ;
                            elif val == ord('(') and keywd in gAccessFuncs:
                                do_it = True  # contribution -> need ;
                            if do_it:
                                while tmpline[-1] != ';' and tmpj < len(lines):
                                    nextln = lines[tmpj]
                                    nextln = handleCompilerDirectives(nextln)
                                    if nextln != "":
                                        nparts = nextln.split()
                                        if nparts[0] in gVAMSkeywords and nparts[0] not in gMathFunctions:
                                            break
                                        tmpline += " " + nextln
                                        tmpline = tmpline.strip()
                                    tmpj += 1
                                if tmpline[-1] == ';' and tmpj > j:
                                    part = tmpline
                                    j = tmpj
                                    gLineNo[-1] = j
                retarray = getAttributes(part, in_attribute)
                part = retarray[0]
                in_attribute = retarray[1]
                attribs += retarray[2:]
                parseOther(part)
# end of parseLines


# deal with comments: // and /* */ (watch out for quotes)
def removeComments( line, in_comment ):
    retstr = ""
    in_quote = False
    i = 0
    last_ch = ""
    while i < len(line):
        if in_comment and not in_quote and i+1 < len(line) and line[i:i+2] == "*/":
            in_comment = False
            i += 2
            if last_ch == ' ' and line[i] == ' ':
                i += 1
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
            last_ch = ch
        i += 1
    return [retstr, in_comment]


def getIndent( line ):
    num_spaces = 0
    for ch in line:
        if ch == '\t':  # TAB discouraged
            return -1
        if ch == ' ':
            num_spaces += 1
        else:
            break
    return num_spaces


def checkSingleSpaceAfterParen( line, keywd ):
    ret = True
    if line.find(")") >= 0 and line.find(keywd) > 0:
        ret = False
        i = 0
        klen = len(keywd)
        while i < len(line):
            if line[i] == ')' and len(line) > i+1+klen:
                if line[i+1:i+2+klen] == " " + keywd:
                    ret = True
                    break
            i += 1
    return ret


def fixSpaces( line, keywd ):
    ret = ""
    i = 0
    klen = len(keywd)
    while i < len(line):
        if line[i:i+klen] == keywd:
            ret += line[i:i+klen]
            i += klen
            break
        ret += line[i]
        i += 1
    ret += ' ' # single space after if
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
            while i < len(line) and line[i].isspace():
                i += 1
            if i < len(line):
                ret += ' '
    if i < len(line):
        ret += line[i:]
    else:
        ret += "\n"
    return ret


def countParentheses( line, num_parens, num_attrib ):
    in_quote = False
    for i in range(len(line)):
        ch = line[i]
        if ch == '"':
            in_quote = not in_quote
        if not in_quote:
            if ch == '(':
                if i+1 < len(line) and line[i+1] == '*':
                    if num_attrib > 0:
                        error("Nested attribute")
                    num_attrib += 1
                else:
                    num_parens += 1
            elif ch == ')':
                if i > 0 and line[i-1] == '*':
                    if num_attrib <= 0:
                        error("Unexpected '*)'")
                    num_attrib -= 1
                else:
                    num_parens -= 1
    return [num_parens, num_attrib]


def parseFile( filename ):
    if not gPreProcess:
        print("Reading %s" % filename)

    if gStyle:
        if filename.find("disciplines.vams") >= 0 or filename.find("discipline.h") >= 0 \
                or filename.find("constants.vams") >= 0 or filename.find("constants.h") >= 0:
            # turn off style-checking for standard header files
            check_style = False
        else:
            check_style = True
    else:
        check_style = False

    # read the lines from the file
    try:
        if sys.version_info >= (3, 0):
            with open(filename, newline='') as in_file:
                lines_of_code = in_file.readlines()
        else:  # pragma: no cover
            with open(filename) as in_file:
                lines_of_code = in_file.readlines()

    except:
        try:
            with open(filename, encoding='windows-1252') as in_file:
                lines_of_code = in_file.readlines()
            if not gPreProcess:
                warning("Non-standard encoding in file %s (windows-1252)" % filename)
        except:
            fatal("Failed to open file: %s" % filename)
    gFileName.append(filename)
    gLineNo.append(0)
    initial_ifdef = len(gIfDefStatus)

    # variables to track indentation style
    fixed_lines = []
    did_fix = False
    required_indent = 0
    optional_indent = 0
    optional_indent_reason = ""
    single_line_indent = 0
    single_line_keywd = [0, ""]
    previous_single_indent = 0
    previous_single_keywd = [0, ""]
    indent_module = -1
    first_module_indent = 0
    ifdef_depth = 0
    indents_in_ifdefs = []
    indent_ifdef = -1
    first_ifdef_indent = -1
    indent_inside_ifdef = -1
    first_inside_ifdef = -1
    indent_multiline_define = -1
    indent_define = -1
    indent_before_define = 0
    first_define_indent = -1
    special_define_indent = 0
    indent_this_define = -1
    first_multiline_define = -1
    in_multiline_define = False
    in_case_block = False
    continued_line = ""
    unclosed_parens = 0
    paren_continued_lines = 0
    dos_format = -1

    # handle comments and line-continuation characters
    lines_to_parse = []
    line_to_parse = ""
    in_comment = False

    for line in lines_of_code:
        gLineNo[-1] += 1
        if line.endswith("\r\n"):
            # dos format; change to unix
            line = line[:-2] + "\n"
            if dos_format == -1:
                dos_format = 1
                style("MS-DOS file format (\\r\\n)")
            elif dos_format == 0 and not gPreProcess:
                style("Carriage return (\\r or ^M) before newline")
        elif dos_format == 1 and not gPreProcess:
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
                [line, in_comment] = removeComments(line, in_comment)
            else:
                line = ""
            if gFixIndent:
                fixed_line = orig_line.rstrip() + "\n"
                fixed_lines.append(fixed_line)
                if fixed_line != orig_line:
                    did_fix = True
        else:
            if continued_line == "parentheses":
                paren_continued_lines += 1
            else:
                paren_continued_lines = 0
            if line.find("//") >= 0 or line.find("/*") >= 0:
                [line, in_comment] = removeComments(line, in_comment)
            if check_style:
                # getIndent from orig_line to check indent of comments
                indent = getIndent(orig_line)
                fixed_line = orig_line.rstrip() + "\n"
                stripped = line.strip()
                if stripped == "":
                    if in_multiline_define:
                        in_multiline_define = False
                        required_indent = indent_before_define
                        indent_this_define = -1
                    if optional_indent > 0 and orig_line.strip() != "":
                        if indent >= required_indent + optional_indent:
                            if indent_inside_ifdef == 0:
                                style("Inconsistent indentation inside `ifdef (compare with line %d)"
                                      % first_inside_ifdef )
                            else:
                                required_indent += optional_indent
                                if ifdef_depth > 0 and len(indents_in_ifdefs) == ifdef_depth:
                                    depths = indents_in_ifdefs[ifdef_depth-1]
                                    depths[0] += optional_indent
                                optional_indent = 0
                                optional_indent_reason = ""
                                if indent_inside_ifdef == -1:
                                    indent_inside_ifdef = 1
                                    first_inside_ifdef = gLineNo[-1]
                    stripped = orig_line.strip()
                    if stripped != "":
                        this_indent = required_indent + single_line_indent
                        if stripped[0:2] == "//":
                            if indent_define == 0 and stripped.find("`define") > 0:
                                this_indent = 0
                            if indent_ifdef == 0:
                                if stripped.find("`ifdef") > 0 or stripped.find("`ifndef") > 0 \
                                        or stripped.find("`else") > 0 or stripped.find("`elsif") > 0 \
                                        or stripped.find("`endif") > 0:
                                    this_indent = 0
                        fixed_line = (this_indent * " ") + stripped + "\n"
                    else:
                        fixed_line = "\n"
                else:
                    if orig_line.endswith("\\\n"):
                        #orig_line = orig_line[:-2] + "\n"
                        if stripped.endswith("\\"):
                            stripped = stripped[:-1].strip()
                        if len(orig_line) > 2 and orig_line[-3] != ' ':
                            style("Prefer space before \\")
                            if gFixIndent:
                                fixed_line = fixed_line[:-2] + " \\\n"
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
                                    need_fix = False
                                    if len(rest) < 2 or rest[0] != ' ' or rest[1] != '(':
                                        style("Prefer single space after 'if'")
                                        need_fix = gFixIndent
                                    if not checkSingleSpaceAfterParen( rest, "begin" ):
                                        style("Prefer single space between ')' and 'begin'")
                                        need_fix = gFixIndent
                                    if need_fix:
                                        fixed_line = fixSpaces(fixed_line, "else if")
                        elif keywd in ["if", "for", "while", "repeat", "case"]:
                            rest = parser.peekRestOfLine()
                            need_fix = False
                            if len(rest) < 2 or rest[0] != ' ' or rest[1] != '(':
                                style("Prefer single space after '%s'" % keywd)
                                need_fix = gFixIndent
                            if not checkSingleSpaceAfterParen( rest, "begin" ):
                                style("Prefer single space between ')' and 'begin'")
                                need_fix = gFixIndent
                            if need_fix:
                                fixed_line = fixSpaces(fixed_line, keywd)
                        elif keywd == "else":
                            rest = parser.peekRestOfLine().strip()
                            if rest.startswith("if"):
                                rest = rest[2:]
                                need_fix = False
                                if len(rest) < 2 or rest[0] != ' ' or rest[1] != '(':
                                    style("Prefer single space after 'if'")
                                    need_fix = gFixIndent
                                if not checkSingleSpaceAfterParen( rest, "begin" ):
                                    style("Prefer single space between ')' and 'begin'")
                                    need_fix = gFixIndent
                                if need_fix:
                                    fixed_line = fixSpaces(fixed_line, "if")
                    elif val == ord('@'):
                        keywd = "@"
                    this_indent = required_indent
                    if val == ord('`') and this_indent >= gSpcPerInd:
                        if (indent_ifdef != 0 and stripped.startswith("`else")) or \
                               (indent_ifdef == -1 and stripped.startswith("`endif")):
                            this_indent -= gSpcPerInd
                    if in_multiline_define:
                        this_indent += special_define_indent
                    previous_single_indent = single_line_indent
                    previous_single_keywd = single_line_keywd
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
                        if indent > this_indent:
                            if optional_indent > 0 and optional_indent_reason == "ifdef" and indent_inside_ifdef == 0:
                                style("Inconsistent indentation inside `ifdef (compare with line %d)"
                                      % first_inside_ifdef )
                            else:
                                if indent == this_indent + optional_indent:
                                    bad_indent = False
                                required_indent += optional_indent
                                this_indent += optional_indent
                                if optional_indent_reason == "module":
                                    indent_module = 1
                                    first_module_indent = gLineNo[-1]
                                elif optional_indent_reason == "define":
                                    if indent_multiline_define == -1:
                                        indent_multiline_define = 1
                                        first_multiline_define = gLineNo[-1]
                                    elif indent_multiline_define != 1:
                                        style("Inconsistent indentation of multi-line `define (compare with line %d)"
                                              % first_multiline_define )
                                    indent_this_define = 1
                                elif optional_indent_reason == "ifdef":
                                    if ifdef_depth > 0 and len(indents_in_ifdefs) == ifdef_depth:
                                        depths = indents_in_ifdefs[ifdef_depth-1]
                                        depths[0] += optional_indent
                                    if indent_inside_ifdef == -1:
                                        indent_inside_ifdef = 1
                                        first_inside_ifdef = gLineNo[-1]
                                optional_indent = 0
                                optional_indent_reason = ""
                    elif indent == this_indent and optional_indent > 0:
                        if optional_indent_reason == "module":
                            if indent_module == 1:
                                if ifdef_depth > 0 and len(indents_in_ifdefs) == ifdef_depth:
                                    depths = indents_in_ifdefs[ifdef_depth-1]
                                    if depths[1] == -1:
                                        depths[0] += optional_indent
                                    else:
                                        required_indent += optional_indent
                                if val != ord('`'):
                                    style("Inconsistent indentation inside module block (compare with line %d)"
                                          % first_module_indent )
                                    if gFixIndent:
                                        fixed_line = (required_indent * " ") + fixed_line.strip() + "\n"
                            else:
                                indent_module = 0
                        elif optional_indent_reason == "define":
                            if indent_multiline_define == -1:
                                indent_multiline_define = 0
                                first_multiline_define = gLineNo[-1]
                            elif indent_multiline_define != 0:
                                style("Inconsistent indentation of multi-line `define (compare with line %d)"
                                      % first_multiline_define )
                            indent_this_define = 0
                        elif optional_indent_reason == "ifdef":
                            if val == ord('`') and orig_line.startswith("`else"):
                                pass  # only comments since `ifdef
                            elif indent_inside_ifdef == -1:
                                indent_inside_ifdef = 0
                                first_inside_ifdef = gLineNo[-1]
                        optional_indent = 0
                        optional_indent_reason = ""

                    may_be_assign_or_contrib = 0
                    if val == ord('`'):
                        val = parser.lex()
                        if parser.isIdentifier():
                            token = parser.getString()
                            if token in gVAMScompdir:
                                bad_indent = False
                                is_vams_compdir = True
                                if gFixIndent and indent > this_indent:
                                    fixed_line = (this_indent * " ") + fixed_line.strip() + "\n"
                            if token not in ["ifdef", "ifndef", "else", "elsif", "endif"]:
                                if continued_line == "parentheses":
                                    warning("Unclosed parenthesis before this point (opened on line %d)"
                                            % (gLineNo[-1]-paren_continued_lines))
                                    unclosed_parens = 0
                                    continued_line = ""
                            if token == "define":
                                indent_before_define = required_indent
                                if required_indent > 0:
                                    if indent == 0:
                                        if indent_define == -1:
                                            indent_define = 0
                                            first_define_indent = gLineNo[-1]
                                        elif indent_define != 0:
                                            style("Inconsistent indentation of `define (compare with line %d)"
                                                  % first_define_indent )
                                            if gFixIndent:
                                                fixed_line = (this_indent * " ") + fixed_line.strip() + "\n"
                                    else:
                                        if indent_define == -1:
                                            indent_define = 1
                                            first_define_indent = gLineNo[-1]
                                        elif indent_define != 1:
                                            style("Inconsistent indentation of `define (compare with line %d)"
                                                  % first_define_indent )
                                            if gFixIndent:
                                                fixed_line = fixed_line.strip() + "\n"
                                if indent not in [0, required_indent]:
                                    if required_indent != 0:
                                        style("Incorrect indent: %d, should be %d (or 0)"
                                              % (indent, required_indent), "Incorrect indent")
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
                            elif token in ["ifdef", "ifndef"]:
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
                                    style("Inconsistent indentation for `ifdef (compare with line %d)"
                                          % first_ifdef_indent, "Incorrect indent")
                                if indent_ifdef == 0:
                                    this_indent = 0
                                fixed_line = (this_indent * " ") + fixed_line.strip() + "\n"
                            elif token in ["else", "elsif"]:
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
                            elif token == "endif":
                                if ifdef_depth > 0:
                                    if ifdef_depth > 0 and len(indents_in_ifdefs) == ifdef_depth:
                                        depths = indents_in_ifdefs[ifdef_depth-1]
                                        if depths[1] == -1:  # no `else
                                            pass  # spurious warnings
                                            #if depths[0] != required_indent:
                                            #    style("Conditional indent inside `ifdef, but no `else")
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
                            elif not is_vams_compdir:
                                # `MODEL begin:
                                val = parser.lex()
                                if parser.isIdentifier():
                                    keywd = parser.getString()
                                else:
                                    val = 0
                    if keywd != "":
                        if keywd in ["module", "macromodule"]:
                            optional_indent = gSpcPerInd
                            optional_indent_reason = "module"
                        elif keywd in ["analog", "if", "else", "for", "repeat", "while", "generate", "@"]:
                            single_line_indent = previous_single_indent + gSpcPerInd
                            single_line_keywd = [gLineNo[-1], keywd]
                            num_parens = 0
                            while val != 0:
                                val = parser.lex()
                                if parser.isIdentifier():
                                    string = parser.getString()
                                    if string in ["begin", "function"]:
                                        optional_indent = 0
                                        required_indent += gSpcPerInd
                                        single_line_indent = 0
                                        if previous_single_indent > 0:
                                            gLineNo.append(previous_single_keywd[0])
                                            style("Missing begin/end for '%s'" % previous_single_keywd[1])
                                            gLineNo.pop()
                                    elif string == "end":
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
                            num_parens = 0
                            while val != 0:
                                val = parser.lex()
                                if parser.isIdentifier():
                                    if parser.getString() == "else":
                                        single_line_indent = previous_single_indent + gSpcPerInd
                                    elif parser.getString() == "begin":
                                        optional_indent = 0
                                        required_indent += gSpcPerInd
                                        single_line_indent = 0
                                    elif string == "end":
                                        required_indent -= gSpcPerInd
                                elif val == ord('('):
                                    num_parens += 1
                                elif val == ord(')'):
                                    num_parens -= 1
                                elif val == ord(';') and num_parens == 0:
                                    single_line_indent = 0
                        else:
                            val = parser.lex()
                            if val == ord('['):
                                rest = parser.peekRestOfLine()
                                if rest.find("]"):
                                    while val != ord(']'):
                                        val = parser.lex()
                                    val = parser.lex()
                            if val == ord('='):
                                may_be_assign_or_contrib = 1
                            elif val == ord('('):
                                rest = parser.peekRestOfLine()
                                if rest.find("<+"):
                                    may_be_assign_or_contrib = 2
                    elif in_case_block:
                        if val == ord('('):
                            val = parser.lex()
                        if parser.isNumber():
                            while val != 0:
                                val = parser.lex()
                                if parser.isIdentifier():
                                    if parser.getString() == "begin":
                                        optional_indent = 0
                                        required_indent += gSpcPerInd
                    if bad_indent and continued_line != "" and indent > this_indent:
                        bad_indent = False
                    elif indent_this_define == 0 and not is_vams_compdir:
                        bad_indent = True
                    if bad_indent and indent != -1 and indent != this_indent:
                        style("Incorrect indent: %d, should be %d" % (indent, this_indent), "Incorrect indent")
                    if bad_indent:
                        if indent_this_define == 0:
                            this_indent += gSpcPerInd
                        fixed_line = (this_indent * " ") + fixed_line.strip() + "\n"
                    if in_multiline_define and line[-2] != '\\':
                        in_multiline_define = False
                        required_indent = indent_before_define
                        indent_this_define = -1
                    # try to determine if this line continues
                    if continued_line == "no_semicolon":
                        if stripped != "" and stripped[-1] == ';':
                            continued_line = ""
                    elif continued_line in ["parentheses", "attributes"]:
                        [unclosed_parens, num_attrib] = countParentheses(line, unclosed_parens, 1)
                        if continued_line == "parentheses":
                            if unclosed_parens <= 0:
                                unclosed_parens = 0
                                continued_line = ""
                            elif stripped == "end" or (stripped != "" and stripped[-1] == ';'):
                                warning("Unclosed parenthesis before this point (opened on line %d)"
                                        % (gLineNo[-1]-paren_continued_lines))
                                unclosed_parens = 0
                                continued_line = ""
                        else:
                            if num_attrib <= 0:
                                continued_line = ""
                    else:
                        continued_line = ""
                        if stripped != "" and stripped[-1] != ';':
                            if may_be_assign_or_contrib:
                                continued_line = "no_semicolon"
                            else:
                                [num_parens, num_attrib] = countParentheses(line, 0, 0)
                                if num_parens > 0:
                                    continued_line = "parentheses"
                                    unclosed_parens = num_parens
                                if num_attrib > 0:
                                    continued_line = "attributes"
                if gFixIndent:
                    fixed_lines.append(fixed_line)
                    if fixed_line != orig_line:
                        did_fix = True
            # end of check_style

        line_to_parse += line
        if line.endswith("\\\n"):
            line_to_parse = line_to_parse[:-2] + "\n"
            lines_to_parse.append("")  # preserve line numbering
        else:
            lines_to_parse.append(line_to_parse)
            line_to_parse = ""

    if continued_line == "parentheses":  # pragma: no cover
        warning("Unclosed parenthesis at end of file %s (opened on line %d)"
                % (filename, (gLineNo[-1]-paren_continued_lines)))

    # write file, if there are any changes
    if gFixIndent:
        fixfile = filename + ".fixed"
        if did_fix:
            with open(fixfile, "w") as out_file:
                out_file.writelines(fixed_lines)
            print("Wrote %s with fixed indentation" % fixfile)
        else:
            if os.path.exists(fixfile) and os.path.isfile(fixfile):
                print("%s did not need fixes to indentation; removing %s" % (filename, fixfile))
                os.remove(fixfile)
            elif gVerbose:
                print("%s did not need fixes to indentation" % filename)

    # empty the arrays
    fixed_lines = []
    lines_of_code = []

    # process the lines
    gLineNo[-1] = 0
    preProcessNaturesAndContribs(lines_to_parse)
    gLineNo[-1] = 0
    if gPreProcess:
        lines_of_code = preProcessLines(lines_to_parse)
    else:
        parseLines(lines_to_parse)

    if len(gIfDefStatus) > initial_ifdef:
        error("Missing `endif")
        while len(gIfDefStatus) > initial_ifdef:
            gIfDefStatus.pop()

    # done with this file
    gFileName.pop()
    gLineNo.pop()
    if len(gFileName) == 0 and len(gScopeList) > 0:
        error("Missing 'endmodule'")

    return lines_of_code
# end of parseFile


################################################################################
# main routine


class InfoAction(argparse.Action):
    """ Class to handle --info """
    def __init__(self, option_strings, dest, nargs=None, help=None):
        argparse.Action.__init__(self, option_strings=option_strings, dest=dest, nargs=nargs, help=help)
    def __call__(self, parser, namespace, values, option_string=None):
        print("""
This program runs checks on Verilog-A models to detect common problems, including:

- hidden state (variables used before assigned)
- unused parameters or variables
- division by zero for parameters (1/parm where parm's range allows 0)
- integer division (1/2 = 0) and domain errors for ln() and sqrt()
- incorrect ddx() usage
- ports without direction and/or discipline
- incorrect access functions for discipline
- unnamed noise sources
- poor coding style
""")
        sys.exit()


class VersionAction(argparse.Action):
    """ Class to handle --version """
    def __init__(self, option_strings, dest, nargs=None, help=None):
        argparse.Action.__init__(self, option_strings=option_strings, dest=dest, nargs=nargs, help=help)
    def __call__(self, parser, namespace, values, option_string=None):
        print("\nVAMPyRE version %s\n" % gVersionNumber)
        sys.exit()


#
#   Process command line arguments
#

argpar = argparse.ArgumentParser(description='Parse and run basic checks on Verilog-A models')
argpar.add_argument('main_file',              help='name of top-level Verilog-A file')
argpar.add_argument('-a', '--all',            help='equivalent to --max_num 0', action='store_true')
argpar.add_argument('-b', '--binning',        help='analyze binning equations', action='store_true')
argpar.add_argument('-d', '--debug',          help='turn on debug mode', action='store_true')
argpar.add_argument('-D', '--define',         help='define token', action='append', default=[])
argpar.add_argument('-f', '--fix_indent',     help='fix indentation problems', action='store_true')
argpar.add_argument('-i', '--info',           help='print detailed usage information', action=InfoAction, nargs=0)
argpar.add_argument('-I', '--inc_dir',        help='add search path for include', action='append', default=[])
argpar.add_argument(      '--max_num',        help='maximum number of each error/warning message type', type=int, \
                                              default=5)
argpar.add_argument('-n', '--no_style',       help='turn off style checking', action='store_false')
argpar.add_argument(      '--not_compact',    help='not a compact model', action='store_false')
argpar.add_argument(      '--preprocess',     help='run preprocessor only', action='store_true')
argpar.add_argument(      '--print_defaults', help='print default values for all parameters', action='store_true')
argpar.add_argument('-s', '--indent_spaces',  help='number of spaces to indent', type=int, default=4)
argpar.add_argument(      '--superfluous',    help='report superfluous assignments', action='store_true')
argpar.add_argument('-v', '--verbose',        help='turn on verbose mode', action='store_true')
argpar.add_argument(      '--version',        help='print version number', action=VersionAction, nargs=0)
cmdargs       = argpar.parse_args()
gIncDir       = cmdargs.inc_dir
gDebug        = cmdargs.debug
gDefines      = cmdargs.define
gBinning      = cmdargs.binning
gFixIndent    = cmdargs.fix_indent
gCompactModel = cmdargs.not_compact
gPrDefVals    = cmdargs.print_defaults
gPreProcess   = cmdargs.preprocess
gSuperfluous  = cmdargs.superfluous
if cmdargs.all:
    gMaxNum = 0
else:
    gMaxNum = cmdargs.max_num
gSpcPerInd = cmdargs.indent_spaces
gStyle     = cmdargs.no_style
gVerbose   = cmdargs.verbose
if gFixIndent and not gStyle:
    fatal("Cannot fix indentation without checking style")
if gVerbose:
    print("VAMPyRE version %s" % gVersionNumber)
initializeModule()
if gPreProcess:
    gFixIndent = False
    gStyle = False
    output = parseFile(cmdargs.main_file)
    for ln in output:
        print(ln)
else:
    parseFile(cmdargs.main_file)
    printModuleSummary()
